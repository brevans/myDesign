#!/usr/bin/env python
'''
Vcf2AxiomMyDesign: mining out SNPs from a vcf file and it's reference to a file ready for submission to the Axiom myDesign custom probe pipeline.

Usage:
    Vcf2AxiomMyDesign.py -s species -v file.vcf -r ref.fa [--minindiv 20] [-o outfile.tsv]

Options:
    -s STRING, --species=STRING  Your species name
    -v FILE, --vcf=FILE          A vcf file
    -r FILE, --fasta=FILE        The matching reference fasta file
    --minindiv=INT               Minimum number of individuals with genotypes
                                  for a position to be output [default: 1]
    -o FILE, --out=FILE          output file [default: myDesign.tsv]

'''

from docopt import docopt
import os
from Bio import SeqIO

class rec(object):
    def __init__ (self):
        '''
        container for previous, current, next objects
        '''
        self.up=0
        self.curr=0
        self.down=0

class vcf(object):
    '''
    fVcf: POSITION SORTED vcf file name
    vcf: a generator that walks over a vcf file

    sample VCF line:
    CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT sam1...samn
    '''
    def __init__(self, fVcf):
        self.fVcf = fVcf
        self.header_lines=[]
        self.vcf_handle=open(fVcf)
        self.rec = rec()

        l=''
        #initialize samples, strip off vcf header
        #stop at first line with real data
        while True:
            l = self.vcf_handle.readline()

            if l.startswith("#CHROM"):
                a=l.lstrip('#').rstrip().split()
                self.cols = a[:9]
                self.samples = a[9:]

            elif l.startswith("##"):
                self.header_lines.append(l.rstrip())

            else:
                self.rec.up = 0
                self.rec.curr = 0
                self.rec.down = self.parse_vcf_line(l)
                break
    
    def __iter__(self):
        return self

    def next(self):
        '''
        '''
        if self.rec.curr != None:
            self.rec.up = self.rec.curr
            self.rec.curr = self.rec.down
            self.rec.down = self.parse_vcf_line(self.vcf_handle.readline())
            return self.rec
        else:
            raise StopIteration()

    def parse_vcf_line(self, line):
        tmp = {}
        tmp['called_snps']=0
        vals = line.rstrip().split()
        if len(vals) != len(self.cols)+len(self.samples):
            return None
        for i,v in zip(self.cols,vals[:9]):
            tmp[i]=v
        for v in vals[9:]:
            if v != './.':
                tmp['called_snps']+=1
        return tmp

def check_snp(s, min_dist, min_sam):
    #end 
    if s.curr is None:
        print "End!"
        return False
    
    #minimum samples found
    if s.curr['called_snps'] < min_sam:
        print( "not enough SNPS called at " +s.curr['CHROM']+":"+s.curr['POS'] +" ("+str(s.curr['called_snps'])+")" )
        return False

    #only bi-allelic
    if len(s.curr['ALT']) !=1:
        print("too many alleles at "+ s.curr['CHROM']+":"+s.curr['POS'] + " ("+s.curr['ALT']+")")
        return False

    #check to make sure up and downstream snps are far enough away for probes to fit between
    for other in [s.up, s.down]:
        if other is None:
            continue
        if other == 0:
            continue
        if s.curr['CHROM'] == other['CHROM']:
            if abs( int(s.curr['POS']) - int(other['POS'])) < min_dist:
                print( "bad distance between " + s.curr['CHROM']+":"+s.curr['POS'] + " and " + other['CHROM']+":"+other['POS'] )
                return False
    return True


def generate_design(org, fVcf, fRef, min_sam, outfile):
    probe_size = 35
    min_dist = probe_size * 2 #35 bp probe for each position
    my_vcf = vcf(fVcf)
    ref = SeqIO.index(fRef, 'fasta')

    my_design=open(outfile, 'w')
    #write my design header
    my_design.write('Organism\tSNPId\tREF_STR\tSEQ\tCHR\tPOS\tSNP_NEIGHBOR\tSNP_PRIORITY\tSNP_VAL\tCHR_TYPE\n')
    ref_name = os.path.basename(fRef).rstrip('.fa')

    for snp in my_vcf:
        #if the snp is good
        snp_good = check_snp(snp, min_dist, min_sam)
        if snp_good:
            #write to my_design file
            p=int(snp.curr['POS'])-1 #(convert to 0-based position)
            c=snp.curr['CHROM']
            snpid = c + "_" + snp.curr['POS']
            ra = "[" + snp.curr['REF'] + "/" + snp.curr['ALT'] + ']'
            seq= str( ref[c][p-probe_size-1:p].seq + ra + ref[c][p+1:p+probe_size+1].seq )
            my_design.write( '\t'.join([org, snpid, ref_name, seq, c, str(p), '0', '1', 'Il', 'autosomal']) + '\n' )

def main():
    args = docopt(__doc__)
    generate_design(args['--species'], args['--vcf'], args['--fasta'], int(args['--minindiv']), args['--out'])
    
if __name__ == '__main__':
    main()
    
