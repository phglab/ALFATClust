from .Utils import read_seq_file_for_preclusters
from Bio import SeqIO
from shutil import rmtree
from tempfile import mkdtemp, mkstemp
import numpy as np
import os
import shlex
import subprocess

class Precluster:
    _temp_dir_path = None

    @staticmethod
    def _parse_precluster_results(seq_precluster_file_path, seq_file_path):
        seq_name_to_precluster_map = dict()
        precluster_sizes = list()
        last_center_seq_name = ''
        precluster_count = -1

        with open(seq_precluster_file_path, 'r') as f:
            for line in f.readlines():
                center_seq_name, seq_name = line.rstrip().split('\t')
                if center_seq_name != last_center_seq_name:
                    precluster_count += 1
                    seq_name_to_precluster_map[center_seq_name] = precluster_count
                    precluster_sizes.append(1)
                    last_center_seq_name = center_seq_name
                    if seq_name == center_seq_name:
                        continue

                seq_name_to_precluster_map[seq_name] = precluster_count
                precluster_sizes[-1] += 1

        return read_seq_file_for_preclusters(seq_file_path, seq_name_to_precluster_map), np.max(precluster_sizes)

    @classmethod
    def precluster_seq_file(cls, user_params, seq_file_path, max_seq_len):
        cls._temp_dir_path = mkdtemp()

        seq_db_file_name_prefix = 'seq_db'
        seq_precluster_db_file_name_prefix = 'seq_precluster_db'
        seq_db_file_path_prefix = os.path.join(cls._temp_dir_path, seq_db_file_name_prefix)
        seq_precluster_db_file_path_prefix = os.path.join(cls._temp_dir_path, seq_precluster_db_file_name_prefix)
        temp_cluster_dir_path = os.path.join(cls._temp_dir_path, 'temp')
        seq_precluster_file_path = os.path.join(cls._temp_dir_path, 'seq_preclusters.tsv')

        mmseqs_build_db_cmd = 'mmseqs createdb {} {}'.format(seq_file_path, seq_db_file_path_prefix)
        subprocess.run(args=shlex.split(mmseqs_build_db_cmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       universal_newlines=True)

        mmseqs_cluster_cmd = \
            'mmseqs cluster {} {} {} --min-seq-id {} -c {} -v 0'.format(seq_db_file_path_prefix,
                                                                        seq_precluster_db_file_path_prefix,
                                                                        temp_cluster_dir_path,
                                                                        user_params.res_param_end, 0.001)

        if max_seq_len > 65535:
            mmseqs_cluster_cmd = '{} --max-seq-len {}'.format(mmseqs_cluster_cmd, max_seq_len)

        subprocess.run(args=shlex.split(mmseqs_cluster_cmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       universal_newlines=True)

        mmseqs_convert_cmd = 'mmseqs createtsv {} {} {} {}'.format(seq_db_file_path_prefix, seq_db_file_path_prefix,
                                                                   seq_precluster_db_file_path_prefix,
                                                                   seq_precluster_file_path)
        subprocess.run(args=shlex.split(mmseqs_convert_cmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       universal_newlines=True)

        rmtree(temp_cluster_dir_path)

        for file_name in os.listdir(cls._temp_dir_path):
            if file_name.startswith(seq_db_file_name_prefix) or \
                file_name.startswith(seq_precluster_db_file_name_prefix):
                os.remove(os.path.join(cls._temp_dir_path, file_name))

        precluster_to_seq_recs_map, max_precluster_size = cls._parse_precluster_results(seq_precluster_file_path,
                                                                                        seq_file_path)
        cls.clear_temp_data()

        return precluster_to_seq_recs_map, max_precluster_size

    @classmethod
    def write_precluster_seq_records(cls, seq_records):
        f, output_file_path = mkstemp(dir=cls._temp_dir_path)
        SeqIO.write(seq_records, f, 'fasta')

        return output_file_path

    @classmethod
    def create_temp_dir(cls):
        if cls._temp_dir_path is None or not os.path.isdir(cls._temp_dir_path):
            cls._temp_dir_path = mkdtemp()

    @classmethod
    def clear_temp_data(cls):
        if cls._temp_dir_path is not None and os.path.isdir(cls._temp_dir_path):
            rmtree(cls._temp_dir_path)
            cls._temp_dir_path = None
