#!/usr/bin/env python
'''
Summarize MyDesign tsv file

Usage:
    SummarizeMyDesign.py -i myDesign.tsv -r ref.txt [-o myDesign_summary.txt]

Options:
    -i FILE, --mydesign=FILE  a tab-delimited myDesign file
    -r FILE, --refinfo=FILE   a white-space delimited file with two fields
                              per line (no header):
                              CHROMNAME LENGTH
    -o FILE, --out=FILE       a file that summarizes the myDesign file

'''
from collections import defaultdict

'''
headers:
Organism SNPId REF_STR SEQ CHR POS SNP_NEIGHBOR SNP_PRIORITY SNP_VAL CHR_TYPE
'''
def read_my_design (fi):

    mifi = open(fi, 'r')
    header = mifi.readline().rstrip('\n').split('\t')
    chrom_subtotals=defaultdict(lambda : 0)
    snp_dict={}
    for l in mifi:
        line = l.rstrip('\n').split('\t')

        snp_dict[line[1]] = [line[4], line[5], line[3]]
        chrom_subtotals[line[4]] += 1

    mifi.close()
    return snp_dict, chrom_subtotals

def read_ref_info(fi):

    mifi = open(fi, 'r')
    chrm={}
    chrm['TOTAL']=0
    for l in mifi:
        tmp = l.split()
        chrm[tmp[0]] = int(tmp[1])
        chrm['TOTAL'] += int(tmp[1])
    return chrm

def summarize_table(snps, subtots, chromlen, outf):

    total_len = 0.0
    covered_len = 0.0
    my_chroms = []
    out = open(outf, 'w')
    for chrom, length in sorted(chromlen.items(), key = lambda x: x[1], reverse = True):
        total_len+=length
        if chrom in subtots:
            covered_len += chromlen[chrom]
            my_chroms.append(chrom)

    print("Found at least one SNP probe in %i of %i Chromosomes, roughly %.2f%% of the genome." % 
          (len(my_chroms), len(chromlen.keys()), (covered_len/total_len*100) ) )
    out.write("CHROM\tNUM_PROBES\tCHOM_LENGTH\tSNPS_PER_kb\n")
    for c in my_chroms:
        kb_ratio = float(subtots[c]) / (float(chromlen[c])/1000.0)
        out.write("%s\t%i\t%i\t%.04f\n" % (c, subtots[c], chromlen[c], kb_ratio))

def main():
    from docopt import docopt
    args = docopt(__doc__)
    snps, subtotals = read_my_design(args['--mydesign'])
    chrom_lengths = read_ref_info(args['--refinfo'])
    if args['--out'] is None:
        outf = args['--mydesign']+'_summary.txt'
    else:
        outf = args['--out']

    summarize_table(snps, subtotals, chrom_lengths, outf)

if __name__ == '__main__':
    main()
