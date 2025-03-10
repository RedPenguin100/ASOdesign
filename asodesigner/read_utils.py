from fuzzysearch import find_near_matches

from asodesigner.fold import get_trigger_mfe_scores_by_risearch, get_mfe_scores


def process_hybridization(args):
    i, l, target_seq, locus_to_data, target_cache_filename = args

    parsing_type = '2'
    scores = get_trigger_mfe_scores_by_risearch(target_seq[i: i + l], locus_to_data,
                                                minimum_score=900, neighborhood=l, parsing_type=parsing_type, target_file_cache=target_cache_filename)
    energy_scores = get_mfe_scores(scores, parsing_type=parsing_type)
    total_candidates = 0
    energy_sum = 0
    max_sum = 0
    binary_sum = 0
    for locus_scores in energy_scores:
        total_candidates += len(locus_scores)
        energy_sum += sum(locus_scores)
        max_sum += min(locus_scores)
        binary_sum += 1 if min(locus_scores) < -20 else 0
    return (i, l, total_candidates, energy_sum, max_sum, binary_sum)


def process_watson_crick_differences(args):
    i, l, target_seq, locus_to_data = args
    sense = target_seq[i:i + l]
    matches_per_distance = [0, 0, 0, 0]

    for locus_tag, locus_info in locus_to_data.items():
        matches = find_near_matches(sense, locus_info, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1
            if match.dist == 0:
                print(locus_tag)

    # Return a tuple containing the starting index, current l, and match counts
    return (i, l, matches_per_distance[0],
            matches_per_distance[1], matches_per_distance[2], matches_per_distance[3])
