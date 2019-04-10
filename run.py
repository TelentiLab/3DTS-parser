from parser import read_scores_into_ram, process_loci_file

_score_file_path = '/Users/yin.li/Downloads/archive/scores.js.gz'
_loci_file_path = '/Users/yin.li/Downloads/archive/structuralFeatures.loci.js.gz'
_result_file_path = './3dts_scores.tsv'
SCORE_LINES = 2273865
LOCI_LINES = 413021682

score_map = read_scores_into_ram(_score_file_path, SCORE_LINES)
process_loci_file(_loci_file_path, _result_file_path, score_map, LOCI_LINES)
