import os.path

import gffutils

chrom_names = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
               'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'Mito']


def open_assembly_db(assembly=2):
    if assembly == 1:
        if os.path.exists('../analysis/S288C_R64-1-1.gff.db'):
            return gffutils.FeatureDB('../analysis/S288C_R64-1-1.gff.db',
                                      keep_order=True)
        else:
            gff_file_name = '../data/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff'
            return gffutils.create_db(gff_file_name,
                                      '../analysis/S288C_R64-1-1.gff.db')

    elif assembly == 2:
        if os.path.exists('../analysis/S288C_R64-2-1.gff.db'):
            return gffutils.FeatureDB('../analysis/S288C_R64-2-1.gff.db',
                                      keep_order=True)
        else:
            gff_file_name = '../data/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff'
            return gffutils.create_db(gff_file_name,
                                      '../analysis/S288C_R64-2-1.gff.db')

    else:
        raise ValueError('Assembly version does not exist')


def get_chromosomes_length(assembly=2):
    db = open_assembly_db(assembly)

    chromosomes = dict()

    for f in db.features_of_type('chromosome'):
        if f.id[3:] == 'mt':
            chromosomes['Mito'] = f.end
        else:
            chromosomes[f.id[3:]] = f.end

    return chromosomes


def invert_clustering_borders(clustering_borders, chromosomes_length):
    inverted_clustering_borders = dict()

    for chrom in sorted(chrom_names):
        if chrom == 'Mito':
            continue

        inverted_clustering_borders[chrom] = []

        sorted_borders = sorted(clustering_borders[chrom], key=lambda pair: pair[0])

        if sorted_borders[0][0] > 1:
            inverted_clustering_borders[chrom].append((1, sorted_borders[0][0] - 1))

        for i in range(len(sorted_borders) - 1):
            inverted_clustering_borders[chrom].append((sorted_borders[i][1] + 1, sorted_borders[i + 1][0] - 1))

        last_position = sorted_borders[len(sorted_borders) - 1][1]
        if last_position < chromosomes_length[chrom]:
            inverted_clustering_borders[chrom].append((last_position + 1, chromosomes_length[chrom]))
        
    return inverted_clustering_borders


def get_cds_coordinates(assembly=2):
    db = open_assembly_db(assembly)

    cds = dict()

    for f in db.features_of_type('CDS'):
        chrom = f.seqid[3:]
        if chrom not in cds:
            cds[chrom] = dict()

        if assembly == 2:
            gene_id = f.attributes['Name'][0][:-4]
        elif assembly == 1:
            gene_id = f.attributes['Name'][0]
        else:
            raise ValueError('Assembly version does not exist')

        # 1-based coordinates of regions
        if gene_id not in cds[chrom]:
            cds[chrom][gene_id] = [(f.start, f.end, f.strand)]
        else:
            cds[chrom][gene_id].append((f.start, f.end, f.strand))

    return cds


def get_genes_list(assembly=2):
    db = open_assembly_db(assembly)

    genes = {'known': {},
             'unknown': {}}
    
    for f in db.features_of_type('gene'):
        chrom = f.seqid[3:]
        if chrom not in genes['known']:
            genes['known'][chrom] = set()
            genes['unknown'][chrom] = set()
            
        if assembly == 2:
            gene_id = f.attributes['Name'][0][:-4]
        elif assembly == 1:
            if 'gene' not in f.attributes:
                genes['unknown'][chrom].add(f.attributes['ID'][0])
            else:
                genes['known'][chrom].add(f.attributes['gene'][0])
        else:
            raise ValueError('Assembly version does not exist')

    return genes


def get_genes_list_2(assembly=2):
    db = open_assembly_db(assembly)

    genes = {'known': set(),
             'unknown': set()}

    for f in db.features_of_type('gene'):
        #chrom = f.seqid[3:]
        #if chrom not in genes['known']:
        #    genes['known'][chrom] = set()
        #    genes['unknown'][chrom] = set()

        if assembly == 2:
            gene_id = f.attributes['Name'][0][:-4]
        elif assembly == 1:
            if 'gene' not in f.attributes:
                genes['unknown'].add(f.attributes['ID'][0])
            else:
                genes['known'].add(f.attributes['gene'][0])
        else:
            raise ValueError('Assembly version does not exist')

    return genes


def get_gene_name_id_converter(assembly=2):
    db = open_assembly_db(assembly)

    genes = {}

    for f in db.features_of_type('gene'):
        #chrom = f.seqid[3:]
        #if chrom not in genes:
        #    genes[chrom] = set()
        #    genes[chrom] = set()

        if assembly == 2:
            gene_id = f.attributes['Name'][0][:-4]
        elif assembly == 1:
            if 'gene' not in f.attributes:
                continue
            else:
                genes[f.attributes['gene'][0]] = f.attributes['ID'][0]
                genes[f.attributes['ID'][0]] = f.attributes['gene'][0]
        else:
            raise ValueError('Assembly version does not exist')

    return genes


def get_gene_name_id_converter_2(assembly=2):
    db = open_assembly_db(assembly)

    genes = {'name_to_id': {},
             'id_to_name': {}}

    for f in db.features_of_type('gene'):
        #chrom = f.seqid[3:]
        #if chrom not in genes:
        #    genes[chrom] = set()
        #    genes[chrom] = set()
        
        if assembly == 2:
            gene_id = f.attributes['Name'][0][:-4]
        elif assembly == 1:
            if 'gene' not in f.attributes:
                continue
            else:
                genes['name_to_id'][f.attributes['gene'][0]] = f.attributes['ID'][0]
                genes['id_to_name'][f.attributes['ID'][0]] = f.attributes['gene'][0]
        else:
            raise ValueError('Assembly version does not exist')

    return genes
