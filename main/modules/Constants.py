DNA = 1
AA = 2

IUPAC_DNA_STR_PATTERN = r'^[ACGTRYSWKMBDHVN]+$'
IUPAC_AA_STR_PATTERN = r'^[ABCDEFGHIKLMNPQRSTVWXYZ]+$'
IUPAC_AMBIG_DNA_BASES = r'[RYSWKMBDHVN]'
IUPAC_AMBIG_AA_BASES = r'[BXZ]'

FASTA_SEQ_NAME_WITH_COMMENT_PATTERN = r'^(\S*)\s?(.*)$'
MASH_SKETCH_MSG_PATTERN = r'^Sketching .+ \(provide sketch file made with "mash sketch" to skip\).+$'
MASH_OUTPUT_PATTERN = r'^(.+)\s+(.+)\s+(\d+(\.\d+)?(e\-\d+)?)\s+(\d+(\.\d+)?(e\-\d+)?)\s+(\d+)\/(\d+)$'
MASH_COMMENT_FIELD_SEP = ':'
