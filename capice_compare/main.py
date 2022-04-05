def match_variants(capice37, capice38, consequences37, consequences38, selected_chromosomes, log, pathogenic_threshold):
    counts = {'total': 0}
    write_to_logs('build37\tbuild38\tconsequence\tscore build 37\tscore build 38\tdiff score\twarning(s)\n', log)
    # Scores increased and decreased in build 38 in comparison with build 37
    score_inc = 0
    score_dec = 0
    # All scores in build 37 and 38 to determine average score for model
    scores37 = []
    scores38 = []
    # Number of pathogenic scores in data for build 37 and 38
    pscores_37 = 0
    pscores_38 = 0

    with open(capice38) as b38:
        b38_data = b38.readlines()
        variants38, variants_list38, variants_per_pos38 = process_variants(b38_data, selected_chromosomes)
    with open(capice37) as b37:
        b37_data = b37.readlines()
        variants37, variants_list37, variants_per_pos37 = process_variants(b37_data, selected_chromosomes)

    for (v37, v38) in zip(variants_per_pos37, variants_per_pos38):
        variants_pos37 = variants_per_pos37[v37]
        variants_pos38 = variants_per_pos38[v38]

        for variant2match in variants_pos37:
            variant_info = variant2match.split('_')
            ref = variant_info[2]
            alt = variant_info[3]
            transcript = variant_info[6]

            for possible_match in variants_pos38:
                # Variants are sorted based on position, so we can assume that for the autosomal chromosomes the order
                # is the same. We also assume the ref, alt and transcript are the same in order to make matches
                # (checked for lp test set in vip and worked, but I realise it's not waterproof)
                if '_{}_{}_'.format(ref, alt) in possible_match and transcript in possible_match:
                    score37 = float(variants37[variant2match])
                    score38 = float(variants38[possible_match])
                    consequence37 = consequences37[variant2match]
                    consequence38 = consequences38[possible_match]

                    # Compare scores between build 37 and build 38
                    counts = compare_scores(score37, score38, consequence37, consequence38, counts, variant2match,
                                            possible_match, log, pathogenic_threshold)

                    scores37.append(score37)
                    scores38.append(score38)

                    if score37 > score38:
                        score_dec += 1
                    else:
                        score_inc += 1

                    if score37 > pathogenic_threshold:
                        pscores_37 += 1
                    if score38 > pathogenic_threshold:
                        pscores_38 += 1

    print_counts(counts, len(variants37), score_inc, score_dec, pscores_37, pscores_38, scores37, scores38)


def get_percentage(count, total):
    return round(100 * count / total, 2)


def generate_consequence_message(consequences, consequence, count, total):
    percentage = get_percentage(count, total)
    consequences.append((consequence, count, percentage))
    return consequences


def print_counts(counts, total, score_inc, score_dec, pscores_37, pscores_38, scores37, scores38):
    consequences = []
    consequences_warn = []

    # Print overall comparisons
    print('Big difference in scores: {}/{}\n'.format(counts['total'], total))

    print('Scores increased in 38: {}'.format(score_inc))
    print('Scores decreased in 38: {}\n'.format(score_dec))

    print('Number of pathogenic classifications in build 37: {}'.format(pscores_37))
    print('Number of pathogenic classifications in build 38: {}\n'.format(pscores_38))

    print('Average score 37: {}'.format(sum(scores37) / total))
    print('Average score 38: {}\n'.format(sum(scores38) / total))

    # Get consequence comparisons
    for count in counts:
        if count != 'total':
            consequences = generate_consequence_message(consequences, count, counts[count]['total'], total)
            consequences_warn = generate_consequence_message(consequences_warn, count, counts[count]['diff'],
                                                             counts['total'])

    # Print consequence comparisons for overall set, sorted on percentage
    print('Overall consequences:')
    for consequence in sorted(consequences, key=lambda x: -x[2]):
        print('{}: {} ({}%)'.format(*consequence))

    # Print consequence for scores that differ
    print('\nConsequences for warnings:')
    for consequence in sorted(consequences_warn, key=lambda x: -x[2]):
        if consequence[2] != 0:
            print('{}: {} ({}%)'.format(*consequence))


def compare_scores(score37, score38, consequence37, consequence38, counts, variant37, variant38, log,
                   pathogenic_threshold):
    # Get absolute difference between score in build 37 and build 38
    diff37_38 = get_score_diff(score37, score38)

    # Calculate outcome based on threshold
    outcome37 = calculate_outcome(score37, pathogenic_threshold)
    outcome38 = calculate_outcome(score38, pathogenic_threshold)

    # Determine whether diff is > 0.1 or if there's a B/P mismatch or if the consequences mismatches between the models
    diff_37_38_01 = diff37_38 > 0.1
    outcome_not_equal = outcome37 != outcome38
    consequence_mismatch = consequence37 != consequence38
    message = []
    if diff_37_38_01 or outcome_not_equal or consequence_mismatch:
        counts['total'] += 1
        write_to_logs(
            '{}\t{}\t{}\t{}({})\t{}({})\t{}\t'.format(variant37, variant38, consequence37, score37, outcome37, score38,
                                                      outcome38, diff37_38), log)
        if diff_37_38_01:
            message.append('Abs diff 37 - 38 > 0.1')
        if outcome_not_equal:
            message.append('Diff outcomes (B/P)')
        if consequence_mismatch:
            message.append('Consequences don\'t match!!!')
        write_to_logs('|'.join(message) + '\n', log)
    if not consequence_mismatch:
        # Count consequences for all scored transcripts and those with diff > 0.1 / B/P
        if consequence37 not in counts:
            counts[consequence37] = {'total': 1, 'diff': 0}
            if diff_37_38_01 or outcome_not_equal:
                counts[consequence37]['diff'] += 1
        else:
            counts[consequence37]['total'] += 1
            if diff_37_38_01 or outcome_not_equal:
                counts[consequence37]['diff'] += 1
    return counts


def process_variants(data, selected_chromosomes):
    variants = {}
    variants_list = []
    variants_per_pos = {}
    for line in data:
        columns_per_index = {0: 'chrom', 1: 'pos', 2: 'ref', 3: 'alt', 4: 'gene', 7: 'transcript', 9: 'score'}
        variant_id, variant_info = get_variant_from_line(line.strip('\n'), columns_per_index)
        chromosome = variant_info['chrom']
        if chromosome in selected_chromosomes or chromosome.replace('chr', '') in selected_chromosomes:
            variants_list.append(variant_id)
            variants[variant_id] = variant_info['score']
            chr_pos = '{}_{}'.format(chromosome, variant_info['pos'])
            if chr_pos in variants_per_pos:
                variants_per_pos[chr_pos].append(variant_id)
            else:
                variants_per_pos[chr_pos] = [variant_id]
    return variants, variants_list, variants_per_pos


def get_splitted_values(string2split, columns_per_index, split_char):
    splitted = string2split.split(split_char)
    return {column: splitted[index] for index, column in columns_per_index.items()}


def get_score_diff(score1, score2):
    return abs(score2 - score1)


def calculate_outcome(score, pathogenic_threshold):
    if score < pathogenic_threshold:
        return 'B'
    else:
        return 'P'


def get_variant_from_line(line, columns_per_index):
    variant_info = get_splitted_values(line, columns_per_index, '\t')
    variant_id = '{}_{}_{}_{}_{}_{}'.format(variant_info['chrom'], variant_info['pos'], variant_info['ref'],
                                            variant_info['alt'], variant_info['gene'], variant_info['transcript'])
    return variant_id, variant_info


def get_consequences(input_file):
    consequences = {}
    with open(input_file) as opened_file:
        input_data = opened_file.readlines()
        for line in input_data:
            columns_per_index = {0: 'chrom', 1: 'pos', 2: 'ref', 3: 'alt', 4: 'consequence', 5: 'gene', 8: 'transcript'}
            variant_id, variant_info = get_variant_from_line(line.strip('\n'), columns_per_index)
            consequences[variant_id] = variant_info['consequence']
    return consequences


def write_to_logs(text, log):
    log.write(text)


def main():
    # be careful, not sure if it works for non-autosomal chromosomes
    selected_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
                            '18', '19', '20', '21', '22']
    print('OUTCOME FOR CHR: ' + str(selected_chromosomes))
    log = open('log_chr{}.tsv'.format('_'.join(selected_chromosomes)), 'w')

    # Capice output
    b37_scores = 'input/capice_p37_capice.tsv'
    b38_scores = 'input/capice_p38_capice.tsv'

    # Capice input (vcf > VEP > BCF tools)
    input37 = 'input/capice_input_p37.tsv'
    input38 = 'input/capice_input_p38.tsv'

    pathogenic_threshold = 0.2

    # Get consequences from the input files to match to the output
    consequences37 = get_consequences(input37)
    consequences38 = get_consequences(input38)

    # Match scored transcripts
    match_variants(b37_scores, b38_scores, consequences37, consequences38, selected_chromosomes, log,
                   pathogenic_threshold)
    log.close()


if __name__ == '__main__':
    main()
