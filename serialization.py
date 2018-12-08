def dump_tn_position_scores(tn_position_scores, out_file_name):
    with open(out_file_name, 'w') as dumpfile:
        for chrom in sorted(tn_position_scores):
            for pos in sorted(tn_position_scores[chrom]):
                dumpfile.write('{}\t{}\t{}\n'.format(chrom, pos, tn_position_scores[chrom][pos]))    


def read_tn_position_scores_from_dump(file_name):
    tn_position_scores = dict()

    with open(file_name, 'r') as dumpfile:
        for line in dumpfile:
            tmp = line.strip().split('\t')

            chrom = tmp[0]
            if chrom not in tn_position_scores:
                tn_position_scores[chrom] = dict()

            tn_position_scores[chrom][int(tmp[1])] = float(tmp[2])

    return tn_position_scores








##############
## TODO: Test
##############




def read_clustering_borders_from_dump(file_name):
    borders = dict()

    with open(file_name, 'r') as dumpfile:
        for line in dumpfile:
            tmp = line.strip().split('\t')

            chrom = tmp[0]
            if chrom not in borders:
                borders[chrom] = []

            borders[chrom].append((int(tmp[1]), int(tmp[2])))

    return borders


