import json
import time
import gzip
import datetime
import logging
from typing import Dict

FEATURE_KEYS = ['pdbId', 'pdbChain', 'uniprotFeatureName', 'pdbResidueMin', 'pdbResidueMax']
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', level=logging.INFO)
_logger = logging.getLogger('3DTS_LOGGER')


def read_scores_into_ram(filename: str, total_lines: int = 0) -> Dict:
    """
    Process the scores.gz file and put (feature, score) mapping into RAM. The feature is serialized to a string
    separated by dot (.) in the format of:
    pdbId.pdbChain.uniprotFeatureName.pdbResidueMin.pdbResidueMax

    The feature object from the file looks like this:
    "featureKey": { # this whole object is the key in score file
        "pdbId": "10GS",
        "pdbChain": "A",
        "uniprotFeatureName": "HELIX",
        "pdbResidueMin": "187",
        "pdbResidueMax": "194"
    }
    :param filename: the full path to the score file
    :param total_lines: number of lines of the file, optional. If provided, progress and estimated time would be shown.
    :return: A dictionary containing feature to score mapping.
    """
    start_time = time.time()
    _logger.info('starting to read scores into RAM...')
    score_map = {}
    counter = 0
    skipped_lines = []
    with gzip.open(filename, 'rt') as file:
        for line in file:
            try:
                data: Dict = json.loads(line.strip())
            except:  # error parsing json
                skipped_lines.append(line)
                _logger.error(f'Failed to parse line into json: {line}')
                continue
            feature = data.get('featureKey')
            if not feature:
                skipped_lines.append(line)
                _logger.error(f'Failed to get featureKey from line: {line}')
                continue
            key = '.'.join([feature[field] for field in FEATURE_KEYS])
            try:
                post_global_synonymous_rate = data.get(
                    'nsPostHeptamerIndependentChromosomeSpecificIntergenicRate'
                ).get('mean')
            except AttributeError:  # failed to find score
                skipped_lines.append(line)
                _logger.error(f'Failed to get score from line: {line}')
                continue
            score_map[key] = post_global_synonymous_rate
            if total_lines:
                counter += 1
                if counter % (total_lines // 1000) == 0:
                    ratio = counter / total_lines
                    time_left = datetime.timedelta(seconds=(time.time() - start_time) * (1 - ratio) / ratio)
                    _logger.info(f'reading score file({100 * ratio:.2f}%), estimated time left: {time_left}')
    end_time = time.time()
    _logger.info(f'# of lines skipped due to parsing error: {len(skipped_lines)}')
    _logger.info(f'score map created, time spent: {end_time - start_time}')
    return score_map


def process_loci_file(loci_filename: str, output_filename: str, score_map: dict, total_lines: int = 0):
    start_time = time.time()
    _logger.info('start processing loci file...')
    counter = 0
    skipped_lines = []
    with open(output_filename, 'w') as file_out:
        file_out.write('#chr\tstart\tend\tscore\tfeature\n')
        with gzip.open(loci_filename, 'rt') as file:
            for line in file:
                try:
                    data = json.loads(line.strip())
                except:  # error parsing json
                    skipped_lines.append(line)
                    _logger.error(f'Failed to parse line into json: {line}')
                    continue
                feature = data.get('feature')
                if not feature:
                    skipped_lines.append(line)
                    _logger.error(f'Failed to get feature from line: {line}')
                    continue
                key = '.'.join([feature.get(k) for k in FEATURE_KEYS])
                score = score_map.get(key)
                if score is None:
                    # _logger.info(f'current loci does not map to any score: {key}')
                    skipped_lines.append(line)
                else:
                    loci = data.get('locus')
                    file_out.write(f"{loci}\t{score}\t{key}\n")
                if total_lines:
                    counter += 1
                    if counter % (total_lines // 5000) == 0:
                        ratio = counter / total_lines
                        time_left = datetime.timedelta(seconds=(time.time() - start_time) * (1 - ratio) / ratio)
                        _logger.info(f'writing result({100 * ratio:.2f}%), estimated time left: {time_left}')

    end_time = time.time()
    _logger.info(f'# of lines skipped due to parsing error: {len(skipped_lines)}')
    _logger.info(f'processing complete, time spent: {datetime.timedelta(seconds=end_time - start_time)}')
