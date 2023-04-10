#!/usr/bin/env python3

import sys
sys.path.append('../modules')

from Config import Config
from Constants import *
from Utils import SEQ_TYPE_MATCH, SEQ_TYPE_MISMATCH, SEQ_TYPE_LOW_QUAL, check_seq_len, infer_seq_type, match_seq_type
from Bio import SeqIO
import argparse
import os

def generate_seq_type_err_msg(seq_hdr, seq_type_match_status, target_seq_type, seq_qual_thres_percent):        
    if seq_type_match_status == SEQ_TYPE_LOW_QUAL:
        if target_seq_type == AA:
            return f'\'{seq_hdr}\' => Less than {seq_qual_thres_percent}% of sequence characters represent ordinary amino acids'
        elif target_seq_type == DNA:
            return f'\'{seq_hdr}\' => Less than {seq_qual_thres_percent}% of sequence characters represent ordinary DNA bases'
    elif seq_type_match_status == SEQ_TYPE_MISMATCH:
        if target_seq_type == AA:
            return f'\'{seq_hdr}\' => Unidentifiable amino acid character(s) found'
        elif target_seq_type == DNA:
            return f'\'{seq_hdr}\' => Unidentifiable DNA base character(s) found'

    return f'\'{seq_hdr}\''

if __name__ == '__main__':
    main_dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    config_file_path = os.path.join(main_dir_path, 'settings.cfg')
    try:
        config = Config(config_file_path)
    except:
        sys.exit('Configuration file \'{}\' is not set properly'.format(config_file_path))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='seq_file_path', required=True,
                        help='DNA/Protein sequence FASTA file path')
    parser.add_argument('-o', '--output', action='store', dest='out_seq_file_path', required=True,
                        help='Output sequence FASTA file path')
    parser.add_argument('-e', '--error', action='store', dest='error_report_file_path',
                        help='Filtered sequences report path')
    parser.add_argument('-b', '--target', action='store', dest='target_seq_type',
                        help='Specify the sequence type (dna|aa) for input dataset or let the system decides [auto]',
                        choices=['aa', 'dna', 'auto'], default='auto')
    parser.add_argument('-k', '--kmer', type=int, help='K-mer size for Mash')
    args = parser.parse_args()
    
    if args.kmer is not None and args.kmer <= 0:
        sys.exit('K-mer size must be a positive integer')
    
    if args.target_seq_type == 'auto':
        target_seq_type, error_msg = infer_seq_type(args.seq_file_path)
        if error_msg is not None:
            sys.exit(error_msg)
    elif args.target_seq_type == 'aa':
        target_seq_type = AA
    elif args.target_seq_type == 'dna':
        target_seq_type = DNA
    
    if args.kmer is None:
        if target_seq_type == AA:
            if config.default_protein_kmer_size > 0:
                kmer_size = config.default_protein_kmer_size
            else:
                sys.exit('Default k-mer size for protein sequences must be a positive integer')
        elif target_seq_type == DNA:
            if config.default_dna_kmer_size > 0:
                kmer_size = config.default_dna_kmer_size
            else:
                sys.exit('Default k-mer size for DNA sequences must be a positive integer')
    else:
        kmer_size = args.kmer
        
    seq_qual_thres_percent = str(round(SEQ_QUAL_THRES * 100))
    
    if args.error_report_file_path is None:
        f_err = None
    else:
        f_err = open(args.error_report_file_path, 'w')
    
    proc_seq_count = 0
    exported_seq_count = 0
    filtered_seq_count = 0
    
    with open(args.seq_file_path, 'r') as f, open(args.out_seq_file_path, 'w') as fw:
        for seq_record in SeqIO.parse(f, 'fasta'):
            is_seq_filtered = False
            seq_type_match_status = match_seq_type(str(seq_record.seq), target_seq_type)
            
            if seq_type_match_status != SEQ_TYPE_MATCH:
                seq_type_err_msg = generate_seq_type_err_msg(seq_record.description, seq_type_match_status,
                                                             target_seq_type, seq_qual_thres_percent)
                if f_err is None:
                    print('                                      ', end='\r')
                    print(seq_type_err_msg)
                else:
                    f_err.write(f'{seq_type_err_msg}{os.linesep}')

                is_seq_filtered = True
                
            is_invalid_seq_len, true_kmer_size = \
                check_seq_len(target_seq_type, len(seq_record.seq), args.kmer, config.default_dna_kmer_size,
                              config.default_protein_kmer_size)
            
            if is_invalid_seq_len:
                seq_len_err_msg = f'\'{seq_record.description}\' => Sequence length shorter than the required k-mer size [{true_kmer_size}]'
                if f_err is None:
                    print('                                      ', end='\r')
                    print(seq_len_err_msg)
                else:
                    f_err.write(f'{seq_len_err_msg}{os.linesep}')

                is_seq_filtered = True
                
            if is_seq_filtered:
                filtered_seq_count += 1
            else:
                SeqIO.write(seq_record, fw, 'fasta')
                exported_seq_count += 1
                
            proc_seq_count +=1
            print(f'{proc_seq_count} sequences processed', end='\r')

    if f_err is not None:
        f_err.close()

    print('')
    print(f'Process completed')
    print(f'{exported_seq_count} sequences exported')
    print(f'{filtered_seq_count} sequences filtered')
