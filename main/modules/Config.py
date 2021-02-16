from Bio.Align import substitution_matrices
from configparser import ConfigParser

RES_PARAM_SECTION = 'ResolutionParameter'
START_VAL = 'Start'
END_VAL = 'End'
STEP_SIZE = 'StepSize'
THRESHOLD_SECTION = 'Threshold'
NOISE_FILTER = 'NoiseFilter'
PRECLUSTER = 'Precluster'
DNA_MASH_SECTION = 'DNAMash'
KMER_SIZE = 'Kmer'
SKETCH_SIZE = 'Sketch'
PROTEIN_MASH_SECTION = 'ProteinMash'
DNA_EVAL_SECTION = 'DNAEvaluation'
MATCH_SCORE = 'MatchScore'
MISMATCH_PENALTY = 'MismatchPenalty'
GAP_OPEN_PENALTY = 'GapOpeningPenalty'
GAP_EXTEND_PENALTY = 'GapExtensionPenalty'
PROTEIN_EVAL_SECTION = 'ProteinEvaluation'
SCORE_MATRIX = 'ScoreMatrix'

'''Object class for configuration file'''
class Config:
    @staticmethod
    def _get_score_matrix(score_mtrx_name):
        score_mtrx_name = score_mtrx_name.upper()
        if score_mtrx_name in set(substitution_matrices.load()):
            return score_mtrx_name
        else:
            return score_mtrx_name

    def __init__(self, config_file_path):
        config_settings = ConfigParser()
        config_settings.read(config_file_path)

        self._res_param_start = float(config_settings[RES_PARAM_SECTION][START_VAL])
        self._res_param_end = float(config_settings[RES_PARAM_SECTION][END_VAL])
        self._res_param_step_size = float(config_settings[RES_PARAM_SECTION][STEP_SIZE])

        self._noise_filter_thres = float(config_settings[THRESHOLD_SECTION][NOISE_FILTER])
        self._precluster_thres = int(config_settings[THRESHOLD_SECTION][PRECLUSTER])

        self._default_dna_kmer_size = int(config_settings[DNA_MASH_SECTION][KMER_SIZE])
        self._default_dna_sketch_size = int(config_settings[DNA_MASH_SECTION][SKETCH_SIZE])

        self._default_protein_kmer_size = int(config_settings[PROTEIN_MASH_SECTION][KMER_SIZE])
        self._default_protein_sketch_size = int(config_settings[PROTEIN_MASH_SECTION][SKETCH_SIZE])

        self._dna_match_score = float(config_settings[DNA_EVAL_SECTION][MATCH_SCORE])
        if self._dna_match_score < 0:
            self._dna_match_score *= -1

        self._dna_mismatch_penalty = float(config_settings[DNA_EVAL_SECTION][MISMATCH_PENALTY])
        if self._dna_mismatch_penalty > 0:
            self._dna_mismatch_penalty *= -1

        self._dna_gap_open_penalty = float(config_settings[DNA_EVAL_SECTION][GAP_OPEN_PENALTY])
        if self._dna_gap_open_penalty > 0:
            self._dna_gap_open_penalty *= -1

        self._dna_gap_extend_penalty = float(config_settings[DNA_EVAL_SECTION][GAP_EXTEND_PENALTY])
        if self._dna_gap_extend_penalty > 0:
            self._dna_gap_extend_penalty *= -1

        self._protein_score_mtrx = self._get_score_matrix(config_settings[PROTEIN_EVAL_SECTION][SCORE_MATRIX])

        self._protein_gap_open_penalty = float(config_settings[PROTEIN_EVAL_SECTION][GAP_OPEN_PENALTY])
        if self._protein_gap_open_penalty > 0:
            self._protein_gap_open_penalty *= -1

        self._protein_gap_extend_penalty = float(config_settings[PROTEIN_EVAL_SECTION][GAP_EXTEND_PENALTY])
        if self._protein_gap_extend_penalty > 0:
            self._protein_gap_extend_penalty *= -1

    @property
    def res_param_start(self):
        return self._res_param_start

    @property
    def res_param_end(self):
        return self._res_param_end

    @property
    def res_param_step_size(self):
        return self._res_param_step_size

    @property
    def noise_filter_thres(self):
        return self._noise_filter_thres

    @property
    def precluster_thres(self):
        return self._precluster_thres

    @property
    def default_dna_kmer_size(self):
        return self._default_dna_kmer_size

    @property
    def default_dna_sketch_size(self):
        return self._default_dna_sketch_size

    @property
    def default_protein_kmer_size(self):
        return self._default_protein_kmer_size

    @property
    def default_protein_sketch_size(self):
        return self._default_protein_sketch_size

    @property
    def dna_match_score(self):
        return self._dna_match_score

    @property
    def dna_mismatch_penalty(self):
        return self._dna_mismatch_penalty

    @property
    def dna_gap_open_penalty(self):
        return self._dna_gap_open_penalty

    @property
    def dna_gap_extend_penalty(self):
        return self._dna_gap_extend_penalty

    @property
    def protein_score_mtrx(self):
        return self._protein_score_mtrx

    @property
    def protein_gap_open_penalty(self):
        return self._protein_gap_open_penalty

    @property
    def protein_gap_extend_penalty(self):
        return self._protein_gap_extend_penalty
