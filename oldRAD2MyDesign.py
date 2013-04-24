#!/usr/bin/env python
import sys
import os
from Bio import SeqIO

inf = open(sys.argv[1])
reff = sys.argv[2]
ref_name = os.path.basename(reff).rstrip('.fa')
ref = SeqIO.index(reff, 'fasta')
out = open(sys.argv[3], 'w')
out.write('Organism\tSNPId\tREF_STR\tSEQ\tCHR\tPOS\tSNP_NEIGHBOR\tSNP_PRIORITY\tSNP_VAL\tCHR_TYPE\n')

probe_size = 35
org = 'Aedes_aegypti'

header = inf.next().split()
'''
stacks_id chr_num chr pos seq relpos ref alt
'''

for snp_line in inf:
    (stacks_id, chr_num, c, pos, seq, relpos, r, a) = snp_line.rstrip().split('\t')
    p = int(pos) + int(relpos)
    ra = '[%s/%s]' % (r, a)
    probe = str( ref[c][p-probe_size-1:p-1].seq + ra + ref[c][p+1:p+probe_size+1].seq )
    snpid = '%s_%d' % (c.replace('.', '_'), p)
    out.write( '\t'.join((org, snpid, ref_name, probe, c, str(p), '0', '1', 'Il', 'autosomal')) + '\n' )

out.close()
