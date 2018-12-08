import numpy as np
import pandas as pd

import os.path

from Bio import SeqIO


import functions

import serialization
#import tree

#import STree


chrom_names = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
               'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'Mito']

chrom_names_utrs_table = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
                          '10', '11', '12', '13', '14', '15', '16', 'mt']

chrom_names_ref = ['NC_001133', 'NC_001134', 'NC_001135', 'NC_001136',
                   'NC_001137', 'NC_001138', 'NC_001139', 'NC_001140',
                   'NC_001141', 'NC_001142', 'NC_001143', 'NC_001144',
                   'NC_001145', 'NC_001146', 'NC_001147', 'NC_001148',
                   'NC_001224']


def convert_chromosome_id(chrom):
    return chrom_names[chrom_names_utrs_table.index(chrom)]


def strip_id(gene_instance):
    gene_instance['gene'] = gene_instance['Gene'].strip()
    return gene_instance


def split_orf_coordinates(gene_instance):
    elem = gene_instance['ORF coordinations'].strip().split(':')
    chrom = convert_chromosome_id(elem[0][3:])

    tmp = elem[1].split('-')
    start = int(tmp[0])

    tmp = elem[1][len(tmp[0]) + 1:].split('(')
    end = int(tmp[0])
    strand = tmp[1][0]

    gene_instance['chrom'], gene_instance['strand'] = chrom, strand
    # Converting original 0-based coordinates to 1-based coordinates
    gene_instance['gene_start'], gene_instance['gene_end'] = start + 1, end + 1
    
    return gene_instance


def convert_int_to_list(value):
    if value < 1000:
        return [value]

    values_list = []
    while value >= 1000:
        values_list.append(value % 1000)
        value = value // 1000
    values_list.append(value)

    values_list.reverse()
    return values_list


def check_compatibility(counts, lengths):
    if 0 in counts:
        raise ValueError('Counts contain 0')

    if 0 in lengths:
        raise ValueError('Lengths contain 0')

    if len(counts) != len(lengths):
        raise ValueError('Counts and Lengths do not match in size')


def select_max_count(counts, lengths):
    max_index, max_count = max(enumerate(counts), key=lambda x: x[1])
    max_length = lengths[max_index]

    return max_count, max_length


def crutch(gene_instance):
    if gene_instance['gene'] == 'YIR035C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 85
    elif gene_instance['gene'] == 'YDR155C':
        gene_instance['max_count'], gene_instance['utr_length'] = 1007, 88
    elif gene_instance['gene'] == 'YDR326C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 151
    elif gene_instance['gene'] == 'YDR424C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 50
    elif gene_instance['gene'] == 'YGL195W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 117
    elif gene_instance['gene'] == 'YGR108W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 173
    elif gene_instance['gene'] == 'YHR029C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 83
    elif gene_instance['gene'] == 'YJR112W-A':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 108
    elif gene_instance['gene'] == 'YKR087C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 124
    elif gene_instance['gene'] == 'YLL050C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 127
    elif gene_instance['gene'] == 'YLL040C':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 71
    elif gene_instance['gene'] == 'YLR056W':
        gene_instance['max_count'], gene_instance['utr_length'] = 20, 197
    elif gene_instance['gene'] == 'YLR430W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 96
    elif gene_instance['gene'] == 'YMR014W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 123
    elif gene_instance['gene'] == 'YNL292W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 137
    elif gene_instance['gene'] == 'YNL208W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 109
    elif gene_instance['gene'] == 'YNL087W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 205
    elif gene_instance['gene'] == 'YOL058W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 142
    elif gene_instance['gene'] == 'YOR153W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 115
    elif gene_instance['gene'] == 'YPL181W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 123
    elif gene_instance['gene'] == 'YPR048W':
        gene_instance['max_count'], gene_instance['utr_length'] = 10, 79

    return gene_instance


def find_utr_peak(gene_instance):
    elem_1 = gene_instance["3'UTR read counts"]
    elem_2 = gene_instance["3'UTR lengths"]

    if type(elem_1) is str:
        if type(elem_2) is str:
            counts = [int(n) for n in elem_1.strip().split(',')]
            lengths = [int(n) for n in elem_2.strip().split(',')]
            
            try:
                check_compatibility(counts, lengths)
                gene_instance['max_count'], gene_instance['utr_length'] = select_max_count(counts, lengths)
            except Exception as inst:
                print('Exception:', inst)
                print('Occurred in block 1 at index:', gene_instance.name)

            return gene_instance

        if type(elem_2) is int:
            counts = [int(n) for n in elem_1.strip().split(',')]
            lengths = convert_int_to_list(elem_2)

            try:
                check_compatibility(counts, lengths)
                gene_instance['max_count'], gene_instance['utr_length'] = select_max_count(counts, lengths)
            except Exception as inst:
                max_count, length = select_max_count(counts, lengths)
                if max_count > 0 and length > 0:
                    gene_instance['max_count'], gene_instance['utr_length'] = select_max_count(counts, lengths)
                else:
                    print('Exception:', inst)
                    print('Occurred in block 2 at index:', gene_instance.name)

            return gene_instance

        raise TypeError('Lengths are neither int nor str when Counts are str')
    
    if type(elem_1) is int:
        if type(elem_2) is str:
            counts = convert_int_to_list(elem_1)
            lengths = [int(n) for n in elem_2.strip().split(',')]

            try:
                check_compatibility(counts, lengths)
                gene_instance['max_count'], gene_instance['utr_length'] = select_max_count(counts, lengths)
            except Exception as inst:
                print('Exception:', inst)
                print('Occurred in block 3 at index:', gene_instance.name)

                gene_instance = crutch(gene_instance)

            return gene_instance

        if type(elem_2) is int:
            counts = convert_int_to_list(elem_1)
            lengths = convert_int_to_list(elem_2)

            try:
                check_compatibility(counts, lengths)
                gene_instance['max_count'], gene_instance['utr_length'] = select_max_count(counts, lengths)
            except Exception as inst:
                print('Exception:', inst)
                print('Occurred in block 4 at index:', gene_instance.name)

                gene_instance = crutch(gene_instance)

            return gene_instance

        raise TypeError('Lengths are neither int nor str when Counts are int')

    raise TypeError('Counts are neither int nor str')


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    return ''.join(bases)


def get_utr_sequences(gene_instance, records):
    for record in records:
        chrom_ref = chrom_names_ref[chrom_names.index(gene_instance['chrom'])]

        if record.id == 'ref|' + chrom_ref + '|':
            break

    utr = record.seq[gene_instance['utr_start']: gene_instance['utr_end'] + 1]
    if gene_instance['strand'] == '-':
        utr = reverse_complement(utr)
    else:
        utr = ''.join(utr)

    gene_instance['utr_seq'] = utr

    return gene_instance





######################################################





def blah():
    if assembly == 2:
        for f in db.features_of_type('CDS'):
            id = f.attributes['Name'][0][:-4]

            if id in gene_names:
                OK = True

                if f.start - 1 != utrs_table[utrs_table['gene'] == id]['start'].item():
                    OK = False

                if f.end - 1 != utrs_table[utrs_table['gene'] == id]['end'].item():
                    OK = False

                if f.strand != utrs_table[utrs_table['gene'] == id]['strand'].item():
                    OK = False

                if OK:
                    gene_names.remove(id)
                else:
                    print(id)
                    break

    elif assembly == 1:
        for f in db.features_of_type('CDS'):
            id = f.attributes['Name'][0]

            if id in gene_names:
                gene_names.remove(id)
                
    else:
        raise ValueError('Assembly version does not exist')


def create_ucsc_track_hub_file_for_3UTRs(aligner, utrs_table):
    bed_dir_name = '/data/mayrc/zheng/Lab/Buki/Yeast_By_Andrey/MyResults/clusters/{}/bed'.format(aligner)
    hub_file_name = bed_dir_name + '/3UTRs_track.bed'

    utrs = {}
    for chrom in sorted(chrom_names):
        if chrom == 'Mito':
            continue
        utrs[chrom] = []

    for index, row in utrs_table.iterrows():
        utrs[row['chrom']].append((row['utr_start'], row['utr_end']))
        

    with open(hub_file_name, 'w') as outfile:
        for chrom in sorted(chrom_names):
            if chrom == 'Mito':
                continue

            for pair in sorted(utrs[chrom], key=lambda x: x[0]):
                outfile.write('{}\t{}\t{}\t{}\t{}\n'.format('chr' + chrom, pair[0] - 1, pair[1], '3UTRs', 1000))


def main():
    original_table = pd.read_excel('../Yeast_3UTR_Peaks.3Seq.xlsx', sheet_name='Yeast_3UTR_Peaks.3Seq')

    #original_table.loc[original_table['Unnamed: 5'] == 97, "3'UTR lengths"] = 97
    #original_table.drop('Unnamed: 5', axis=1, inplace=True)
    #[6692 rows x 5 columns]
    
    utrs_info_table = original_table[pd.notnull(original_table["3'UTR read counts"])]
    #[5507 rows x 5 columns]
    
    utrs_info_table = utrs_info_table.assign(gene = np.nan)
    utrs_info_table = utrs_info_table.apply(strip_id, axis=1)
    utrs_info_table.drop('Gene', axis=1, inplace=True)

    utrs_info_table = utrs_info_table.assign(chrom = np.nan,
                                             gene_start = np.nan,
                                             gene_end = np.nan,
                                             strand = np.nan)
    # 0-base original coordinates will be converted to 1-based coordinates
    utrs_info_table = utrs_info_table.apply(split_orf_coordinates, axis=1)
    utrs_info_table.drop('ORF coordinations', axis=1, inplace=True)
    #[5507 rows x 8 columns]

    utrs_length_table = utrs_info_table.assign(max_count = np.nan,
                                               utr_length = np.nan)
    utrs_length_table = utrs_length_table.apply(find_utr_peak, axis=1)
    utrs_length_table = utrs_length_table[pd.notnull(utrs_length_table['utr_length'])]
    utrs_length_table[['max_count', 'utr_length']] = utrs_length_table[['max_count', 'utr_length']].astype(int)
    #[5494 rows x 10 columns]
    
    utrs_length_table.drop('ORF length', axis=1, inplace=True)
    utrs_length_table.drop("3'UTR read counts", axis=1, inplace=True)
    utrs_length_table.drop("3'UTR lengths", axis=1, inplace=True)
    #[5494 rows x 7 columns]

    utrs_table = utrs_length_table.assign(utr_start = np.nan,
                                          utr_end = np.nan,
                                          utr_seq = np.nan)
    
    gb = utrs_table.groupby(utrs_table['strand'])
    utrs_table.loc[gb.groups['+'], 'utr_start'] = utrs_table.loc[gb.groups['+'], 'gene_end'] + 1
    utrs_table.loc[gb.groups['+'], 'utr_end'] = utrs_table.loc[gb.groups['+'], 'gene_end'] + utrs_table.loc[gb.groups['+'], 'utr_length']
    utrs_table.loc[gb.groups['-'], 'utr_start'] = utrs_table.loc[gb.groups['-'], 'gene_start'] - utrs_table.loc[gb.groups['-'], 'utr_length']
    utrs_table.loc[gb.groups['-'], 'utr_end'] = utrs_table.loc[gb.groups['-'], 'gene_start'] - 1
    utrs_table[['utr_start', 'utr_end']] = utrs_table[['utr_start', 'utr_end']].astype(int)

    assembly = 1
    
    if assembly == 1:
        records = list(SeqIO.parse("../data/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa", "fasta"))
    elif assembly == 2:
        records = list(SeqIO.parse("../data/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa", "fasta"))
    else:
        raise ValueError('Assembly version does not exist')

    utrs_table = utrs_table.apply(get_utr_sequences, axis=1, records=records)
    #[5494 rows x 10 columns]

    print ("First test")
    aligner = 'hisat2'
    create_ucsc_track_hub_file_for_3UTRs(aligner, utrs_table)
    print ("First test")
   # return
    print ("First test")
    clusters_dir_name = '/data/mayrc/zheng/Lab/Buki/Yeast_By_Andrey/analysis/clusters/{}'.format(aligner)
    print ("First test")
    chromosomes_length = functions.get_chromosomes_length(assembly=1)
    print("Second test")
    print(chromosomes_length)
    
    for method in ['linear']:
       # for eps in range(15, 121, 5):
        for eps in range(15,52,2):
            for min_samples in range(6, 15, 2):
                dump_file_name = clusters_dir_name + '/dump_{}_tn_clusters_{}_{}.txt'.format(method, eps, min_samples)
                clustering_borders = serialization.read_clustering_borders_from_dump(dump_file_name)
                inverted_clustering_borders = functions.invert_clustering_borders(clustering_borders, chromosomes_length)

                j = 0
                chrom_prev = utrs_table.iloc[0]['chrom']

                essential_gene_ids = []

                essential_gene_ids = {0.50: [],
                                      0.55: [],
                                      0.60: [],
                                      0.65: [],
                                      0.70: [],
                                      0.75: [],
                                      0.80: [],
                                      0.85: [],
                                      0.90: [],
                                      0.95: []}
                
                for i in range(len(utrs_table)):
                    chrom = utrs_table.iloc[i]['chrom']

                    if chrom != chrom_prev:
                        chrom_prev = chrom
                        j = 0
                    
                    while inverted_clustering_borders[chrom][j][1] < utrs_table.iloc[i]['gene_start']:
                        j += 1                    
                
                    if inverted_clustering_borders[chrom][j][0] > utrs_table.iloc[i]['gene_end']:
                        continue
                    else:
                        k = j

                        overlap = 0
                        
                        while inverted_clustering_borders[chrom][k][0] <= utrs_table.iloc[i]['gene_end']:
                            if inverted_clustering_borders[chrom][k][0] >= utrs_table.iloc[i]['gene_start']:
                                if inverted_clustering_borders[chrom][k][1] <= utrs_table.iloc[i]['gene_end']:
                                    overlap += inverted_clustering_borders[chrom][k][1] - inverted_clustering_borders[chrom][k][0] + 1
                                else:
                                    overlap += utrs_table.iloc[i]['gene_end'] - inverted_clustering_borders[chrom][k][0] + 1
                            else:
                                if inverted_clustering_borders[chrom][k][1] <= utrs_table.iloc[i]['gene_end']:
                                    overlap += inverted_clustering_borders[chrom][k][1] - utrs_table.iloc[i]['gene_start'] + 1
                                else:
                                    overlap += utrs_table.iloc[i]['gene_end'] - utrs_table.iloc[i]['gene_start'] + 1

                            k += 1

                        for cutoff in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
                            if overlap / (utrs_table.iloc[i]['gene_end'] - utrs_table.iloc[i]['gene_start'] + 1) > cutoff:
                                essential_gene_ids[cutoff].append(i)

                for cutoff in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
                    essential_genes = utrs_table.iloc[essential_gene_ids[cutoff]]

                    essential_dir_name = '/data/mayrc/zheng/Lab/Buki/Yeast_By_Andrey/MyResults/essential'
                    out_file_name = essential_dir_name + '/essential_genes_{}_{}_{}_{}.csv'.format(method, eps, min_samples, cutoff)
                    
                    essential_genes.to_csv(out_file_name)
    
    """
    for method in ['linear']:
        print(method)
        print("Test")
        for eps in range(15, 121, 5):
            print(eps)
            for min_samples in range(6, 11, 2):
                print(eps,min_samples)
                dump_file_name = clusters_dir_name + '/dump_{}_tn_clusters_{}_{}.txt'.format(method, eps, min_samples)
                clustering_borders = serialization.read_clustering_borders_from_dump(dump_file_name)
                inverted_clustering_borders = functions.invert_clustering_borders(clustering_borders, chromosomes_length)

                j = 0
                chrom_prev = utrs_table.iloc[0]['chrom']

                essential_gene_ids = []

                essential_gene_ids = {0.80: [],
                                      0.85: [],
                                      0.90: [],
                                      0.95: []}
                
                for i in range(len(utrs_table)):
                    chrom = utrs_table.iloc[i]['chrom']

                    if chrom != chrom_prev:
                        chrom_prev = chrom
                        j = 0

                    while inverted_clustering_borders[chrom][j][1] < utrs_table.iloc[i]['utr_start']:
                        j += 1

                    if inverted_clustering_borders[chrom][j][0] > utrs_table.iloc[i]['utr_end']:
                        continue
                    else:
                        k = j

       	                overlap = 0

       	                while inverted_clustering_borders[chrom][k][0] <= utrs_table.iloc[i]['utr_end']:
                            if inverted_clustering_borders[chrom][k][0] >= utrs_table.iloc[i]['utr_start']:
                                if inverted_clustering_borders[chrom][k][1] <= utrs_table.iloc[i]['utr_end']:
                                    overlap += inverted_clustering_borders[chrom][k][1] - inverted_clustering_borders[chrom][k][0] + 1
                                else:
                                    overlap += utrs_table.iloc[i]['utr_end'] - inverted_clustering_borders[chrom][k][0] + 1
                            else:
                                if inverted_clustering_borders[chrom][k][1] <= utrs_table.iloc[i]['utr_end']:
                                    overlap += inverted_clustering_borders[chrom][k][1] - utrs_table.iloc[i]['utr_start'] + 1
                                else:
                                    overlap += utrs_table.iloc[i]['utr_end'] - utrs_table.iloc[i]['utr_start'] + 1

                            k += 1

                        for cutoff in [0.8, 0.85, 0.9, 0.95]:
                            if overlap / (utrs_table.iloc[i]['utr_end'] - utrs_table.iloc[i]['utr_start'] + 1) > cutoff:
                                essential_gene_ids[cutoff].append(i)

                for cutoff in [0.8, 0.85, 0.9, 0.95]:
                    essential_genes = utrs_table.iloc[essential_gene_ids[cutoff]]

                    essential_dir_name = '/data/mayrc/zheng/Lab/Buki/Yeast_By_Andrey/MyResults/essential'
                    out_file_name = essential_dir_name + '/essential_utrs_{}_{}_{}_{}.csv'.format(method, eps, min_samples, cutoff)

                    essential_genes.to_csv(out_file_name)  
    """




                
                
if __name__ == "__main__":
    main()
