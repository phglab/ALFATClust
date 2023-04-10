#!/usr/bin/env python3

from Bio import SeqIO
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='seq_file_path', required=True,
                        help='DNA/Protein sequence FASTA file path')
    parser.add_argument('-o', '--output', action='store', dest='out_seq_file_path', required=True,
                        help='Output sequence FASTA file path')
    args = parser.parse_args()
    
    with open(args.seq_file_path, 'r') as f, open(args.out_seq_file_path, 'w') as fw:
        proc_seq_count = 0
        
        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.description.replace(' ', '_')
            seq_record.description = ''
            SeqIO.write(seq_record, fw, 'fasta')
    
            proc_seq_count += 1
            print(f'{proc_seq_count} sequences processed', end='\r')

        print('')
        print('Process completed')
