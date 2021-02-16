from .Constants import DNA, AA
from .Utils import cal_outlier_thres_by_iqr, read_seq_file_for_eval
from Bio.Align import PairwiseAligner, substitution_matrices
from multiprocessing import Pool
from pandas import DataFrame
import numpy as np
import os

class ClusterEval:
    _CLUSTER_ID_COL = 'Cluster Id'
    _CLUSTER_SIZE_COL = 'No. of sequences'
    _AVG_SEQ_IDENT_COL = 'Average sequence identity'
    _MIN_SEQ_IDENT_COL = 'Min. sequence identity'
    _CENTER_SEQ_COL = 'Center sequence'
    _SEQ_FOR_MIN_SEQ_IDENT_COL = 'Sequence for min. sequence identity'
    _EVAL_TABLE_COLS = [_CLUSTER_SIZE_COL, _AVG_SEQ_IDENT_COL, _MIN_SEQ_IDENT_COL, _CENTER_SEQ_COL,
                        _SEQ_FOR_MIN_SEQ_IDENT_COL]

    _is_init = False

    @classmethod
    def init(cls, config, num_of_threads):
        cls._diff_len_nt_aligner = PairwiseAligner(mode='global',
                                                   match_score=config.dna_match_score,
                                                   mismatch_score=config.dna_mismatch_penalty,
                                                   query_open_gap_score=config.dna_gap_open_penalty,
                                                   query_extend_gap_score=config.dna_gap_extend_penalty,
                                                   target_internal_open_gap_score=config.dna_gap_open_penalty,
                                                   target_internal_extend_gap_score=config.dna_gap_extend_penalty)

        cls._same_len_nt_aligner = PairwiseAligner(mode='global',
                                                   match_score=config.dna_match_score,
                                                   mismatch_score=config.dna_mismatch_penalty,
                                                   open_gap_score=config.dna_gap_open_penalty,
                                                   extend_gap_score=config.dna_gap_extend_penalty)

        protein_score_mtrx = substitution_matrices.load(config.protein_score_mtrx)

        cls._diff_len_aa_aligner = PairwiseAligner(mode='global',
                                                   substitution_matrix=protein_score_mtrx,
                                                   query_open_gap_score=config.protein_gap_open_penalty,
                                                   query_extend_gap_score=config.protein_gap_extend_penalty,
                                                   target_internal_open_gap_score=config.protein_gap_open_penalty,
                                                   target_internal_extend_gap_score=config.protein_gap_extend_penalty)

        cls._same_len_aa_aligner = PairwiseAligner(mode='global',
                                                   substitution_matrix=protein_score_mtrx,
                                                   open_gap_score=config.protein_gap_open_penalty,
                                                   extend_gap_score=config.protein_gap_extend_penalty)

        cls._num_of_threads = num_of_threads
        cls._is_init = True

    @staticmethod
    def _select_max_len_seq_rec(cluster_seq_recs):
        max_len_seq_rec = None
        max_seq_len = 0

        for seq_rec in cluster_seq_recs:
            seq_len = len(seq_rec.seq)

            if seq_len > max_seq_len:
                max_len_seq_rec = seq_rec
                max_seq_len = seq_len

        return max_len_seq_rec

    @staticmethod
    def _select_center_seq_rec(cluster_seq_recs, global_edge_weight_mtrx):
        cluster_seq_lens = np.empty(len(cluster_seq_recs))

        for i in range(len(cluster_seq_recs)):
            cluster_seq_lens[i] = len(cluster_seq_recs[i].seq)

        candidate_seq_len_range = cal_outlier_thres_by_iqr(cluster_seq_lens)
        candidate_seq_id_to_seq_rec_map = dict()
        cluster_seq_ids = list()

        for i in range(len(cluster_seq_recs)):
            seq_id = int(cluster_seq_recs[i].id)
            cluster_seq_ids.append(seq_id)

            if cluster_seq_lens[i] >= candidate_seq_len_range[0] and cluster_seq_lens[i] <= candidate_seq_len_range[1]:
                candidate_seq_id_to_seq_rec_map[seq_id] = cluster_seq_recs[i]

        candidate_seq_ids = list(candidate_seq_id_to_seq_rec_map.keys())
        candidate_seq_ids.sort()

        center_seq_index = np.argmax(np.nansum(global_edge_weight_mtrx[np.ix_(candidate_seq_ids, cluster_seq_ids)],
                                               axis=1))

        return candidate_seq_id_to_seq_rec_map[candidate_seq_ids[center_seq_index]]

    @classmethod
    def _generate_seq_pairs_for_align(cls, seq_cluster_ptrs, global_edge_weight_mtrx, seq_file_path):
        seq_id_to_non_singleton_cluster_id_map = dict()

        for cluster_id in range(np.max(seq_cluster_ptrs) + 1):
            cluster_seq_ids = np.argwhere(seq_cluster_ptrs == cluster_id).flatten()
            if cluster_seq_ids.size > 1:
                for seq_id in cluster_seq_ids:
                    seq_id_to_non_singleton_cluster_id_map[seq_id] = cluster_id

        if len(seq_id_to_non_singleton_cluster_id_map) == 0:
            return list()

        cluster_to_seq_recs_map = read_seq_file_for_eval(seq_file_path, seq_id_to_non_singleton_cluster_id_map)
        seq_rec_pairs_to_align = list()

        for cluster_id, cluster_seq_recs in cluster_to_seq_recs_map.items():
            if len(cluster_seq_recs) > 2:
                center_seq_rec = cls._select_center_seq_rec(cluster_seq_recs, global_edge_weight_mtrx)
            else:
                center_seq_rec = cls._select_max_len_seq_rec(cluster_seq_recs)

            for seq_rec in cluster_seq_recs:
                if seq_rec.id != center_seq_rec.id:
                    seq_rec_pairs_to_align.append((center_seq_rec, seq_rec, cluster_id))

        return seq_rec_pairs_to_align

    @classmethod
    def _cal_dna_seq_ident(cls, seq_rec_pair_to_align):
        seq1 = seq_rec_pair_to_align[0].seq
        seq2 = seq_rec_pair_to_align[1].seq
        seq1_len = len(seq1)
        seq2_len = len(seq2)

        if seq1_len > seq2_len:
            short_seq = seq2
            long_seq = seq1
            nt_aligner = cls._diff_len_nt_aligner
        elif seq2_len > seq1_len:
            short_seq = seq1
            long_seq = seq2
            nt_aligner = cls._diff_len_nt_aligner
        else:
            short_seq = seq1
            long_seq = seq2
            nt_aligner = cls._same_len_nt_aligner

        rev_short_seq = short_seq.reverse_complement()
        fwd_align_score = nt_aligner.score(short_seq, long_seq)
        rev_align_score = nt_aligner.score(rev_short_seq, long_seq)

        if rev_align_score > fwd_align_score:
            short_seq = rev_short_seq

        max_seq_ident = 0
        check_count = 0

        for seq_align in nt_aligner.align(short_seq, long_seq):
            seq_align_pattern = format(seq_align, 'fasta').split(os.linesep)[1]
            seq_ident = seq_align_pattern.count('|') / len(seq_align_pattern.strip('-'))

            if seq_ident > max_seq_ident:
                max_seq_ident = seq_ident
                check_count = 0
                continue

            check_count += 1
            if check_count > 4:
                break

        return str(seq_rec_pair_to_align[0].id), str(seq_rec_pair_to_align[1].id), seq_rec_pair_to_align[2], \
            max_seq_ident

    @classmethod
    def _cal_protein_seq_ident(cls, seq_rec_pair_to_align):
        seq1 = seq_rec_pair_to_align[0].seq
        seq2 = seq_rec_pair_to_align[1].seq
        seq1_len = len(seq1)
        seq2_len = len(seq2)

        if seq1_len > seq2_len:
            short_seq = seq2
            long_seq = seq1
            aa_aligner = cls._diff_len_aa_aligner
        elif seq2_len > seq1_len:
            short_seq = seq1
            long_seq = seq2
            aa_aligner = cls._diff_len_aa_aligner
        else:
            short_seq = seq1
            long_seq = seq2
            aa_aligner = cls._same_len_aa_aligner

        max_seq_ident = 0
        check_count = 0

        for seq_align in aa_aligner.align(short_seq, long_seq):
            seq_align_pattern = format(seq_align, 'fasta').split(os.linesep)[1]
            seq_ident = seq_align_pattern.count('|') / len(seq_align_pattern.strip('-'))

            if seq_ident > max_seq_ident:
                max_seq_ident = seq_ident
                check_count = 0
                continue

            check_count += 1
            if check_count > 4:
                break

        return str(seq_rec_pair_to_align[0].id), str(seq_rec_pair_to_align[1].id), seq_rec_pair_to_align[2], \
            max_seq_ident

    @classmethod
    def _set_cal_seq_ident_func(cls, seq_type):
        if seq_type == AA:
            return cls._cal_protein_seq_ident
        elif seq_type == DNA:
            return cls._cal_dna_seq_ident

    @classmethod
    def _convert_to_dataframe(cls, cluster_id_to_eval_output_map, seq_file_info):
        cluster_eval_output_df = DataFrame.from_dict(cluster_id_to_eval_output_map, orient='index',
                                                     columns=cls._EVAL_TABLE_COLS)
        cluster_eval_output_df.index.name = cls._CLUSTER_ID_COL
        cluster_eval_output_df.astype({cls._CLUSTER_SIZE_COL:int, cls._AVG_SEQ_IDENT_COL:float,
                                       cls._MIN_SEQ_IDENT_COL:float})
        cluster_eval_output_df[cls._AVG_SEQ_IDENT_COL] /= (cluster_eval_output_df[cls._CLUSTER_SIZE_COL] - 1)
        cluster_eval_output_df = cluster_eval_output_df.round({cls._AVG_SEQ_IDENT_COL:3, cls._MIN_SEQ_IDENT_COL:3})
        cluster_eval_output_df.replace({cls._CENTER_SEQ_COL: seq_file_info.seq_id_to_seq_name_map}, inplace=True)
        cluster_eval_output_df.replace({cls._SEQ_FOR_MIN_SEQ_IDENT_COL: seq_file_info.seq_id_to_seq_name_map},
                                       inplace=True)
        cluster_eval_output_df.index += 1

        return cluster_eval_output_df.sort_index()

    @classmethod
    def eval_clusters(cls, seq_cluster_ptrs, global_edge_weight_mtrx, seq_file_info):
        if not cls._is_init:
            return None

        seq_rec_pairs_to_align = cls._generate_seq_pairs_for_align(seq_cluster_ptrs, global_edge_weight_mtrx,
                                                                   seq_file_info.seq_file_path)

        cluster_id_to_eval_output_map = dict()
        seq_ident_func = cls._set_cal_seq_ident_func(seq_file_info.seq_type)

        with Pool(processes=cls._num_of_threads, maxtasksperchild=10) as pool:
            for seq_ident_tuple in pool.imap_unordered(seq_ident_func, seq_rec_pairs_to_align, 2000):
                center_seq_id, cluster_seq_id, cluster_id, seq_ident = seq_ident_tuple
                if cluster_id in cluster_id_to_eval_output_map:
                    cluster_eval_output = cluster_id_to_eval_output_map[cluster_id]
                    cluster_eval_output[0] += 1
                    cluster_eval_output[1] += seq_ident

                    if seq_ident < cluster_eval_output[2]:
                        cluster_eval_output[2] = seq_ident
                        cluster_eval_output[4] = cluster_seq_id
                else:
                    cluster_id_to_eval_output_map[cluster_id] = [2, seq_ident, seq_ident, center_seq_id, cluster_seq_id]

        return cls._convert_to_dataframe(cluster_id_to_eval_output_map, seq_file_info)

    @classmethod
    def eval_clusters_single_thread(cls, seq_cluster_ptrs, global_edge_weight_mtrx, seq_file_info):
        if not cls._is_init:
            return None

        seq_rec_pairs_to_align = cls._generate_seq_pairs_for_align(seq_cluster_ptrs, global_edge_weight_mtrx,
                                                                   seq_file_info.seq_file_path)

        cluster_id_to_eval_output_map = dict()
        seq_ident_func = cls._set_cal_seq_ident_func(seq_file_info.seq_type)

        for seq_ident_tuple in map(seq_ident_func, seq_rec_pairs_to_align):
            center_seq_id, cluster_seq_id, cluster_id, seq_ident = seq_ident_tuple
            if cluster_id in cluster_id_to_eval_output_map:
                cluster_eval_output = cluster_id_to_eval_output_map[cluster_id]
                cluster_eval_output[0] += 1
                cluster_eval_output[1] += seq_ident

                if seq_ident < cluster_eval_output[2]:
                    cluster_eval_output[2] = seq_ident
                    cluster_eval_output[4] = cluster_seq_id
            else:
                cluster_id_to_eval_output_map[cluster_id] = [2, seq_ident, seq_ident, center_seq_id, cluster_seq_id]

        return cls._convert_to_dataframe(cluster_id_to_eval_output_map, seq_file_info)
