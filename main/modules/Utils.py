import sys
sys.path.append('modules')

from .Constants import *
from Bio import SeqIO
from collections import namedtuple
import numpy as np
import os
import re

SEQ_TYPE_MATCH = 100
SEQ_TYPE_MISMATCH = 101
SEQ_TYPE_LOW_QUAL = 102
SEQ_TYPE_NO_TARGET = 103

def _check_seq_type(seq_str):
    if re.match(IUPAC_DNA_STR_PATTERN, seq_str):
        unambig_dna_seq_str = re.sub(IUPAC_AMBIG_DNA_BASES, '', seq_str)
        if len(unambig_dna_seq_str) / len(seq_str) >= SEQ_QUAL_THRES:
            return DNA
        else:
            return AA
    elif re.match(IUPAC_AA_STR_PATTERN, seq_str):
        unambig_protein_seq_str = re.sub(IUPAC_AMBIG_AA_BASES, '', seq_str)
        if len(unambig_protein_seq_str) / len(seq_str) >= SEQ_QUAL_THRES:
            return AA

    return None

def infer_seq_type(seq_file_path):
    dna_seq_count = 0
    protein_seq_count = 0
    unknown_seq_count = 0
    
    with open(seq_file_path, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_type = _check_seq_type(str(seq_record.seq))
            if seq_type != None:
                if seq_type == AA:
                    protein_seq_count += 1
                elif seq_type == DNA:
                    dna_seq_count += 1
            else:
                unknown_seq_count += 1
                    
    if unknown_seq_count > dna_seq_count and unknown_seq_count > protein_seq_count:
        return None, 'The majority of sequences are of unknown sequence type'
    
    if dna_seq_count == protein_seq_count:
        error_msg = 'Cannot determine the sequence type for the majority. \
            Consider specifying the sequence type via the \'-b\' option.'
        return None, error_msg
    
    if protein_seq_count > dna_seq_count:
        return AA, None
    
    if dna_seq_count > protein_seq_count:
        return DNA, None
    
    return None, 'Failed to infer sequence type.'

def match_seq_type(seq_str, target_seq_type=None):
    if target_seq_type == AA:
        if re.match(IUPAC_AA_STR_PATTERN, seq_str):
            unambig_protein_seq_str = re.sub(IUPAC_AMBIG_AA_BASES, '', seq_str)
            if len(unambig_protein_seq_str) / len(seq_str) >= SEQ_QUAL_THRES:
                return SEQ_TYPE_MATCH
            else:
                return SEQ_TYPE_LOW_QUAL
        else:
            return SEQ_TYPE_MISMATCH
    elif target_seq_type == DNA:        
        if re.match(IUPAC_DNA_STR_PATTERN, seq_str):
            unambig_dna_seq_str = re.sub(IUPAC_AMBIG_DNA_BASES, '', seq_str)
            if len(unambig_dna_seq_str) / len(seq_str) >= SEQ_QUAL_THRES:
                return SEQ_TYPE_MATCH
            else:
                return SEQ_TYPE_LOW_QUAL
        else:
            return SEQ_TYPE_MISMATCH
    
    return SEQ_TYPE_NO_TARGET

def _convert_to_mash_seq_name(seq_name):    
    m = re.match(FASTA_SEQ_HDR_PATTERN, seq_name)
    if m:
        if m.group(3) is None:
            return f'{m.group(1)}{MASH_COMMENT_FIELD_SEP}'
        else:
            return f'{m.group(1)}{MASH_COMMENT_FIELD_SEP}{m.group(3)}'

    return None

def check_seq_len(seq_type, seq_len, kmer_size, default_dna_kmer_size, default_protein_kmer_size):
    if kmer_size is not None:
        return seq_len < kmer_size, kmer_size
    elif seq_type == DNA:
        return seq_len < default_dna_kmer_size, default_dna_kmer_size
    else:
        return seq_len < default_protein_kmer_size, default_protein_kmer_size

def read_seq_file(seq_file_path, user_params=None, precluster_seq_type=None):
    SeqFileInfo = namedtuple('SeqFileInfo', ['mash_seq_name_to_seq_id_map', 'seq_id_to_seq_name_map',
                                             'error_log', 'seq_file_path', 'seq_type', 'mash_last_seq_name',
                                             'seq_count', 'max_seq_len'])
    mash_seq_name_to_seq_id_map = dict()
    seq_id_to_seq_name_map = dict()
    seq_error_log = list()
    target_seq_type = None
    last_seq_name = None
    mash_last_seq_name = None
    seq_count = 0
    max_seq_len = 0
    
    if user_params is None:
        is_infer_seq_type = False
        is_check_seq_len = False
    else:
        is_infer_seq_type = user_params.target_seq_type == 'auto'
        kmer_size = user_params.kmer_size
        default_dna_kmer_size = user_params.default_dna_kmer_size
        default_protein_kmer_size = user_params.default_protein_kmer_size
        is_check_seq_len = True
        
    if is_infer_seq_type:
        target_seq_type, error_msg = infer_seq_type(seq_file_path)
        if error_msg is not None:
            seq_error_log.append(error_msg)
            
            return SeqFileInfo(dict(), dict(), seq_error_log, seq_file_path, target_seq_type, '', 1, 0)
    else:
        if user_params is None:
            if precluster_seq_type is None:
                seq_error_log.append('Sequence type for sequence pre-clusters is missing')
                return SeqFileInfo(dict(), dict(), seq_error_log, seq_file_path, '', '', 1, 0)
            else:
                target_seq_type = precluster_seq_type
        else:            
            if user_params.target_seq_type == 'aa':
                target_seq_type = AA
            elif user_params.target_seq_type == 'dna':
                target_seq_type = DNA

    with open(seq_file_path, 'r') as f:
        mismatch_seq_type_msg = '\'{}\': Unidentifiable {} character(s) found'
        low_qual_seq_type_msg = '\'{}\': Less than ' + str(round(SEQ_QUAL_THRES * 100)) + \
            '% of sequence characters represent ordinary {}'
        short_seq_len_msg = '\'{}\': Sequence length shorter than the required k-mer size [{}]'
        invalid_seq_hdr_msg = '\'{}\': Sequence header consisting of more than one whitespace or flanking whitespace'

        for seq_record in SeqIO.parse(f, 'fasta'):
            last_seq_name = seq_record.description
            
            seq_type_match_status = match_seq_type(str(seq_record.seq), target_seq_type)
            if seq_type_match_status == SEQ_TYPE_MISMATCH:
                if target_seq_type == AA:
                    seq_error_log.append(mismatch_seq_type_msg.format(last_seq_name, 'amino acid'))
                elif target_seq_type == DNA:
                    seq_error_log.append(mismatch_seq_type_msg.format(last_seq_name, 'DNA base'))
            elif seq_type_match_status == SEQ_TYPE_LOW_QUAL:
                if target_seq_type == AA:
                    seq_error_log.append(low_qual_seq_type_msg.format(last_seq_name, 'amino acids'))
                elif target_seq_type == DNA:
                    seq_error_log.append(low_qual_seq_type_msg.format(last_seq_name, 'DNA bases'))

            seq_len = len(seq_record.seq)

            if is_check_seq_len:
                is_invalid_seq_len, true_kmer_size = check_seq_len(target_seq_type, seq_len, kmer_size,
                                                                   default_dna_kmer_size, default_protein_kmer_size)
                if is_invalid_seq_len:
                    seq_error_log.append(short_seq_len_msg.format(seq_record.description, true_kmer_size))

            if seq_len > max_seq_len:
                max_seq_len = seq_len

            mash_last_seq_name = _convert_to_mash_seq_name(last_seq_name)
            if mash_last_seq_name is None:
                seq_error_log.append(invalid_seq_hdr_msg.format(seq_record.description))

            mash_seq_name_to_seq_id_map[mash_last_seq_name] = seq_count
            seq_id_to_seq_name_map[str(seq_count)] = last_seq_name
            seq_count += 1

    return SeqFileInfo(mash_seq_name_to_seq_id_map, seq_id_to_seq_name_map, seq_error_log, seq_file_path, \
                       target_seq_type, mash_last_seq_name, seq_count, max_seq_len)

def read_seq_file_for_preclusters(seq_file_path, seq_name_to_precluster_map):
    precluster_to_seq_recs_map = dict()

    with open(seq_file_path, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            precluster_id = seq_name_to_precluster_map[seq_record.id]

            if precluster_id in precluster_to_seq_recs_map:
                precluster_to_seq_recs_map[precluster_id].append(seq_record)
            else:
                precluster_to_seq_recs_map[precluster_id] = [seq_record]

    return precluster_to_seq_recs_map

def read_seq_file_for_eval(seq_file_path, seq_id_to_non_singleton_cluster_id_map):
    cluster_to_seq_recs_map = dict()
    seq_id = 0

    with open(seq_file_path, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            if seq_id in seq_id_to_non_singleton_cluster_id_map:
                seq_record.description = ''
                seq_record.id = seq_id
                cluster_id = seq_id_to_non_singleton_cluster_id_map[seq_id]

                if cluster_id in cluster_to_seq_recs_map:
                    cluster_to_seq_recs_map[cluster_id].append(seq_record)
                else:
                    cluster_to_seq_recs_map[cluster_id] = [seq_record]

            seq_id += 1

    return cluster_to_seq_recs_map

def cal_outlier_thres_by_iqr(data_vals):
    if np.all(np.isnan(data_vals)):
        return np.array([0, 0])

    quartile_vals = np.nanpercentile(data_vals, [25, 75])
    inter_quartile_range = quartile_vals[1] - quartile_vals[0]

    return np.array([quartile_vals[0] - 1.5 * inter_quartile_range, quartile_vals[1] + 1.5 * inter_quartile_range])

def convert_to_seq_clusters(seq_cluster_ptrs, seq_id_to_seq_name_map):
    output_seq_clusters = list()

    for cluster_id in range(np.max(seq_cluster_ptrs) + 1):
        seq_cluster = list()
        for seq_id in np.argwhere(seq_cluster_ptrs == cluster_id).flatten():
            seq_cluster.append('{}{}'.format(seq_id_to_seq_name_map[str(seq_id)], os.linesep))

        output_seq_clusters.append(seq_cluster)

    return output_seq_clusters

def get_max_precision(*vals):
    max_precision = 2

    for val in vals:
        m = re.match(r'\d+\.(\d*[1-9])0*', str(val))
        if m:
            precision = len(m.group(1))
            if precision > max_precision:
                max_precision = precision

    return max_precision
