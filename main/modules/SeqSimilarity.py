from .Constants import AA, MASH_SKETCH_MSG_PATTERN, MASH_OUTPUT_PATTERN
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
    _max_dist = None
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
        cls._max_dist = 1 - user_params.noise_filter_thres
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
        mash_error_msg = None

        with os.fdopen(fid) as f:
            while True:
                mash_output = f.readline().rstrip()
                m = re.match(MASH_OUTPUT_PATTERN, mash_output)
                if not m:
                    if re.match(MASH_SKETCH_MSG_PATTERN, mash_output):
                        continue

                    mash_error_msg = 'Error occurred in Mash as follows:{}{}'.format(os.linesep, mash_output)
                    break

                seq_name1 = m.group(1)
                seq_name2 = m.group(2)

                if seq_name1 in mash_seq_name_to_seq_id_map:
                    seq_id1 = mash_seq_name_to_seq_id_map[seq_name1]
                else:
                    mash_error_msg = 'Failed to match \'{}\' from Mash'.format(seq_name1)
                    break

                if seq_name2 in mash_seq_name_to_seq_id_map:
                    seq_id2 = mash_seq_name_to_seq_id_map[seq_name2]
                else:
                    mash_error_msg = 'Failed to match \'{}\' from Mash'.format(seq_name2)
                    break

                if seq_id1 == max_seq_id and seq_id2 == max_seq_id:
                    break

                if seq_id1 == seq_id2:
                    continue

                if cls._min_shared_hash_ratio is not None:
                    if int(m.group(9)) / int(m.group(10)) < cls._min_shared_hash_ratio:
                        continue

                global_edge_weight_mtrx[seq_id1, seq_id2] = 1 - float(m.group(3))

        return global_edge_weight_mtrx, mash_error_msg

    @classmethod
    def get_pairwise_similarity(cls, seq_file_info):
        if not cls._is_init:
            return None

        mash_command = 'mash dist -C -i -v {} -d {} -p {} {} {}'
        mash_command = mash_command.format(cls._p_value, cls._max_dist, cls._num_of_threads,
                                           seq_file_info.seq_file_path, seq_file_info.seq_file_path)

        if seq_file_info.seq_type == AA:
            mash_command = '{} -a -k {} -s {}'.format(mash_command, cls._protein_kmer_size, cls._protein_sketch_size)
        else:
            mash_command = '{} -k {} -s {}'.format(mash_command, cls._dna_kmer_size, cls._dna_sketch_size)

        if cls._seed is not None:
            mash_command = '{} -S {}'.format(mash_command, cls._seed)

        fr, fw = os.pipe()

        with subprocess.Popen(args=shlex.split(mash_command), stdout=fw, stderr=fw) as p:
            global_edge_weight_mtrx, mash_error_msg = \
                cls._parse_mash_output(fr, seq_file_info.mash_seq_name_to_seq_id_map, seq_file_info.seq_count)

        os.close(fw)

        return global_edge_weight_mtrx, mash_error_msg
