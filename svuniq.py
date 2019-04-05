import os
import sys
import cyvcf2
import gzip
from pybedtools import BedTool

def subtract_overlaps(key, values):
    tmpA = open('tmpA.bed', 'w')
    tmpB = open('tmpB.bed', 'w')
    sv = key.split(':')
    tmpA = [(sv[0],sv[1],sv[2],sv[3])]
    tmpAbed = BedTool(tmpA)
    tmpBbed = BedTool(values)
    AN = int(0)
    AC = int(0)
    for v in values:
        AN += int(v[4])
        AC += int(v[5])
    subtract = tmpAbed.subtract(tmpBbed)
    for interval in subtract:
        cut = 'cut:' + '_'.join(sv[0:3]) + ':' + str(AN) + ':' + str(AC)
        svlist = [interval[0],interval[1],interval[2],interval[3],cut]
        uniques[interval[3]].append(svlist)

vcf = cyvcf2.VCF(sys.argv[1], gts012=True)
popbed = BedTool(sys.argv[2])#.split(',')
tmpbed = []
for v in vcf:
    chrom = v.CHROM
    pos = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'DEL' or svtype == 'DUP' or svtype == 'INV':
        if end > pos:
            out = [str(chrom), str(pos), str(end), svtype]
            tmpbed.append(out)

vcfbed = BedTool(tmpbed)
intersect = vcfbed.intersect(popbed, wao=True)
uniques = {}
overlaps = {}
for interval in intersect:
    chrom1,start1,end1,type1,chrom2,start2,end2,gnomadid,type2,an,ac,homref,het,homalt,shared = interval
    sv = [str(chrom1),str(start1),str(end1),type1]
    if type1 not in uniques:
        uniques[type1] = []
    if chrom2 == '.':
        sv.append('uncut:' + '_'.join([str(chrom1),str(start1),str(end1)]))
        uniques[type1].append(sv)
    else:
        key = ':'.join(sv)
        if key not in overlaps:
            overlaps[key] = []
        if type1 == type2:
            overlap = [str(chrom2),str(start2),str(end2),type2,str(an),str(ac)]
            overlaps[key].append(overlap)

for key in overlaps:
    if len(overlaps[key]) == 0:
        sv = key.split(':')
        svlist = [str(sv[0]),str(sv[1]),str(sv[2]),sv[3],'uncut:' + '_'.join([str(sv[0]),str(sv[1]),str(sv[2])])]
        uniques[sv[3]].append(svlist)
    else:
        subtract_overlaps(key,overlaps[key])

merged = []
for svtype in uniques:
    svbed = BedTool(uniques[svtype])
    sortmerge = svbed.sort().merge(c="4,5", o="collapse")
    for interval in sortmerge:
        merged.append(interval)

mergedbed = BedTool(merged)
final = []
for interval in mergedbed.sort():
    chrom, start, end, svtypes, chunks = interval
    svtype = list(set(svtypes.split(',')))
    chunk = list(set(chunks.split(',')))
    out = [str(chrom), str(start), str(end), ','.join(svtype), ','.join(chunk)]
    final.append(out)

finalbed = BedTool(final)
finalbed.sort().saveas('uniques.bed')
