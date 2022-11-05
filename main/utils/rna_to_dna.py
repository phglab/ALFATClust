#!/usr/bin/env python3

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', dest='rna_seq_file_path', required=True,
                    help='RNA sequence FASTA file path')
parser.add_argument('-o', '--output', action='store', dest='output_dna_seq_file_path', required=True,
                    help='Output DNA sequence FASTA file path')
parser.add_argument('-r', '--reverse', action='store_true',
                    help='Perform reverse complement after reverse transcription')      
args = parser.parse_args()

if args.reverse:
    print('Reverse complement: ON')
else:
    print('Reverse complement: OFF')

print('Converting RNA sequences to DNA sequences...')
process_count = 0

with open(args.rna_seq_file_path, 'r') as f, open(args.output_dna_seq_file_path, 'w') as fw:
    for seq_record in SeqIO.parse(f, 'fasta'):
        seq_record.seq = seq_record.seq.back_transcribe().upper()
        if args.reverse:
            seq_record.seq = seq_record.seq.reverse_complement()

        SeqIO.write(seq_record, fw, 'fasta')

        process_count += 1
        print(f'Processed {process_count} sequences', end='\r')

print('')
print('Conversion completed')
