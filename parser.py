import json
import time
import gzip
from typing import Dict

_score_file_path = '/Users/yin.li/Downloads/archive/scores.js'
_loci_file_path = '/Users/yin.li/Downloads/archive/structuralFeatures.loci.js.gz'
SCORE_LINES = 2273865
LOCI_LINES = 413021682
FEATURE_KEYS = ['pdbId', 'pdbChain', 'uniprotFeatureName', 'pdbResidueMin', 'pdbResidueMax']


def read_score_into_ram(filename: str, total_lines: int = 0) -> Dict:
    start_time = time.time()
    score_map = {}
    counter = 0
    skipped_lines = []
    with open(filename, 'r') as file:
        for line in file:
            try:
                data: Dict = json.loads(line.strip())
            except:  # error parsing json
                skipped_lines.append(line)
                continue
            """
            "featureKey": { # this whole object is the key in this file
                "pdbId": "10GS",
                "pdbChain": "A",
                "uniprotFeatureName": "HELIX",
                "pdbResidueMin": "187",
                "pdbResidueMax": "194"
            }
            """
            feature = data.get('featureKey')
            if not feature:
                skipped_lines.append(line)
                continue
            key = '.'.join([feature[field] for field in FEATURE_KEYS])
            try:
                post_global_synonymous_rate = data.get(
                    'nsPostHeptamerIndependentChromosomeSpecificIntergenicRate'
                ).get('mean')
            except AttributeError:  # failed to find score
                skipped_lines.append(line)
                continue
            score_map[key] = post_global_synonymous_rate
            if total_lines:
                counter += 1
                ratio = counter / total_lines
                time_left = (time.time() - start_time) * (1 - ratio) / ratio
                print(f'reading score file({100 * ratio:.2f}%), estimated time left: {time_left:0.0f}s')
    end_time = time.time()
    print(f'# of lines skipped due to parsing error: {len(skipped_lines)}')
    print(f'score map created, time spent: {end_time - start_time}')
    return score_map


def process_loci_file(filename: str, total_lines: int = 0):
    start_time = time.time()
    score_map = read_score_into_ram(_score_file_path, SCORE_LINES)
    counter = 0
    skipped_lines = []
    with open('./scores_mapping.tsv', 'w') as file_out:
        file_out.write('#chr\tstart\tend\tscore\tfeature')
        with gzip.open(filename, 'rt') as file:
            for line in file:
                counter += 1
                data = json.loads(line.strip())
                feature = data.get('feature')
                key = '.'.join([feature.get(k) for k in FEATURE_KEYS])
                score = score_map.get(key)
                loci = data.get('locus')
                file_out.write(f"{loci}\t{score}\t{feature.get('uniprotFeatureName')}\n")
                if total_lines:
                    counter += 1
                    print(f'writing result({100 * counter / total_lines:.2f}%)')
                if counter >= 30:
                    break

    end_time = time.time()
    print(f'# of lines skipped due to parsing error: {len(skipped_lines)}')
    print(f'processing complete, time spent: {end_time - start_time}')


process_loci_file(_loci_file_path, LOCI_LINES)
