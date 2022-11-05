### **RNA to DNA conversion**

Convert RNA sequences to DNA sequences through reverse transcription

**Command**

./rna_to_dna.py [optional arguments] -i *\<RNA sequence file path>* -o *\<DNA sequence file path>*

**Mandatory arguments**

| Argument name                           | Description                                         |
| --------------------------------------- | --------------------------------------------------- |
| -i/--input *\<RNA sequence file path>*  | (full/relative) input RNA sequence FASTA file path  |
| -o/--output *\<DNA sequence file path>* | (full/relative) output DNA sequence FASTA file path |

**Optional arguments**

| Argument name | Description                                            |
| ------------- | ------------------------------------------------------ |
| -r/--reverse  | Perform reverse complement after reverse transcription |
| -h/--help     | show help message and exit                             |
