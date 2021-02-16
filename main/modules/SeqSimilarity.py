from .Constants import AA
import numpy as np
import os
import re
import shlex
import subprocess

class SeqSimilarity:
    _dna_kmer_size = None
    _protein_kmer_size = None
    _dna_sketch_size = None
    _protein_sketch_size = None
    _seed = None
    _min_shared_hash_ratio = None
    _noise_filter_thres = None
    _num_of_threads = None
    _p_value = None
    _is_init = False

    @classmethod
    def init(cls, user_params, p_value=0.0001):
        if user_params.kmer_size is None:
            cls._dna_kmer_size = user_params.default_dna_kmer_size
            cls._protein_kmer_size = user_params.default_protein_kmer_size
        else:
            cls._dna_kmer_size = user_params.kmer_size
            cls._protein_kmer_size = user_params.kmer_size

        if user_params.sketch_size is None:
            cls._dna_sketch_size = user_params.default_dna_sketch_size
            cls._protein_sketch_size = user_params.default_protein_sketch_size
        else:
            cls._dna_sketch_size = user_params.sketch_size
            cls._protein_sketch_size = user_params.sketch_size

        cls._seed = user_params.seed
        cls._min_shared_hash_ratio = user_params.min_shared_hash_ratio
        cls._noise_filter_thres = user_params.noise_filter_thres
        cls._num_of_threads = user_params.num_of_threads
        cls._p_value = p_value
        cls._is_init = True

    @classmethod
    def set_to_run_in_single_thread(cls):
        cls._num_of_threads = 1

    @classmethod
    def _parse_mash_output(cls, fid, mash_seq_name_to_seq_id_map, seq_count):
        max_seq_id = seq_count - 1
        global_edge_weight_mtrx = np.zeros((seq_count, seq_count), dtype=np.float32)

        with os.fdopen(fid) as f:
            while True:
                line = f.readline()
                line_fields = line.rstrip().split('\t')
                seq_name1 = line_fields[0]
                seq_name2 = line_fields[1]
                seq_id1 = mash_seq_name_to_seq_id_map[seq_name1]
                seq_id2 = mash_seq_name_to_seq_id_map[seq_name2]

                if seq_id1 == max_seq_id and seq_id2 == max_seq_id:
                    break

                if seq_id1 == seq_id2:
                    continue

                if cls._min_shared_hash_ratio is not None:
                    m = re.match(r'(\d+)/(\d+)', line_fields[4])
                    if m:
                        if int(m.group(1)) / int(m.group(2)) < cls._min_shared_hash_ratio:
                            continue
                    else:
                        continue

                global_edge_weight_mtrx[seq_id1, seq_id2] = 1 - float(line_fields[2])

        return global_edge_weight_mtrx

    @classmethod
    def get_pairwise_similarity(cls, seq_file_info):
        if not cls._is_init:
            return None

        mash_command = 'mash dist -C -i -v {} -p {} {} {}'.format(cls._p_value, cls._num_of_threads,
                                                                  seq_file_info.seq_file_path,
                                                                  seq_file_info.seq_file_path)

        if seq_file_info.seq_type == AA:
            mash_command = '{} -a -k {} -s {}'.format(mash_command, cls._protein_kmer_size, cls._protein_sketch_size)
        else:
            mash_command = '{} -k {} -s {}'.format(mash_command, cls._dna_kmer_size, cls._dna_sketch_size)

        if cls._seed is not None:
            mash_command = '{} -S {}'.format(mash_command, cls._seed)

        fr, fw = os.pipe()

        with subprocess.Popen(args=shlex.split(mash_command), stdout=fw, stderr=subprocess.DEVNULL) as p:
            global_edge_weight_mtrx = cls._parse_mash_output(fr, seq_file_info.mash_seq_name_to_seq_id_map,
                                                             seq_file_info.seq_count)

        os.close(fw)

        global_edge_weight_mtrx[global_edge_weight_mtrx < cls._noise_filter_thres] = 0

        return global_edge_weight_mtrx
