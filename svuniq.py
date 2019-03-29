import os
import sys
import cyvcf2
import gzip
from pybedtools import BedTool

#def open_file(in_file):
#    if in_file[-3:] == '.gz':
#        opened  = gzip.open(vcf_file, 'r')
#    elif in_file[-3:] == 'vcf':
#        opened = open(vcf_file, 'r')
#    return opened

#def combine_beds(beds):
#    cmd = 'cat ' + ' '.join(beds) + ' | gsort /dev/stdin ~/genome_ref/genome.nochr.file > tmp.bed'
#    os.system(cmd)
    
def subtract_overlaps(key, values):
    tmpA = open('tmpA.bed', 'w')
    tmpB = open('tmpB.bed', 'w')
    sv = key.split(':')
    tmpA.write('\t'.join(sv) + '\n')
    tmpA.close()
    for overlap in values:
        tmpB.write('\t'.join(overlap) + '\n')
    tmpB.close()
    tmpAbed = BedTool('tmpA.bed')
    tmpBbed = BedTool('tmpB.bed')
    subtract = tmpAbed.subtract(tmpBbed)
    subtract.saveas('tmpsubtract.bed')
    chunks = open('tmpsubtract.bed')
    for chunk in chunks:
        fields = chunk.rstrip('\n').split('\t')
        cut = 'cut:' + '_'.join(sv[0:3])
        uniques.write('\t'.join(fields) + '\t' + cut + '\n')
    chunks.close()

vcf = cyvcf2.VCF(sys.argv[1], gts012=True)
bed = BedTool(sys.argv[2])#.split(',')
tmpbed = open('tmp.bed', 'w')

for v in vcf:
    chrom = v.CHROM
    pos = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'DEL' or svtype == 'DUP':
        if end > pos and chrom != 'hs37d5':
            out = [str(chrom), str(pos), str(end), svtype]
            tmpbed.write('\t'.join(out) + '\n')

tmpbed.close()

vcfbed = BedTool('tmp.bed')
intersect = vcfbed.intersect(bed, wao=True)
intersect.saveas('tmpintersect.bed')
intersects = open('tmpintersect.bed', 'r')

uniques = open('tmpuniqs.bed', 'w')
overlaps = {}
for line in intersects:
    fields = line.rstrip('\n').split('\t')
    chrom1,start1,end1,type1,chrom2,start2,end2,type2 = fields[0:8]
    sv = [str(chrom1),str(start1),str(end1),type1]
    if fields[4] == '.':
        uniques.write('\t'.join(sv) + '\t' + 'uncut' + '\n')
    else:
        key = ':'.join(sv)
        if key not in overlaps:
            overlaps[key] = []
        if type1 == type2:
            overlap = [str(chrom2),str(start2),str(end2),type2]
            overlaps[key].append(overlap)

intersects.close()

for key in overlaps:
    if len(overlaps[key]) == 0:
        sv = key.split(':')
        uniques.write('\t'.join(sv) + '\t' + 'uncut' + '\n')
    else:
        subtract_overlaps(key,overlaps[key])

uniques.close()

cleanup = open('tmpuniqs.bed', 'r')
svtypes = {}
for line in cleanup:
    fields = line.rstrip('\n').split('\t')
    chrom,start,end,svtype,cut = fields[0:5]
    if svtype not in svtypes:
        svtypes[svtype] = []
    svtypes[svtype].append(fields)

for svtype in svtypes:
    tmp = open('tmpclean.bed', 'w')
    for chunk in svtypes[svtype]:
        tmp.write('\t'.join(chunk) + '\n')
    tmp.close()
    clean = BedTool('tmpclean.bed')
    sort = clean.sort()
    merged = sort.merge()
    
