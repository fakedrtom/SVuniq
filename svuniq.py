import os
import sys
import cyvcf2
import gzip
from pybedtools import BedTool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--in',
                    metavar='INPUT VCF',
                    dest="i",
                    help='path to VCF')
#parser.add_argument('-f', '--minf',
#                    metavar='MINIMUM OVERLAP',
#                    dest="f",
#                    type=float,
#                    help='minimum reciprocal overlap required between SVs (default 0.01, must be between 0 and 1.0)')
parser.add_argument('-pop', '--pop',
                    metavar='PATH TO BED',
                    dest="pop",
                    help='path to SV bed file with population AFs')
parser.add_argument('-ci', '--ci',
                    metavar='USE CI BOUNDARIES',
                    dest="ci",
                    choices=['in','out'],
                    help='option to use out inner or outer confidence intervals (CIPOS95, CIEND95) for SV boundaries, must answer "in" or "out"')
parser.add_argument('-minf', '--minf',
                    metavar='MINIMUM ALLELE FREQUENCY',
                    dest="minf",
                    type=float,
                    help='option to only include SVs from the popoulation AFs bed that exceed a minimum AF threshold')

args = parser.parse_args()
if args.i is None:
    raise NameError('Must include path to VCF with option -i')
else:
    vcf = cyvcf2.VCF(args.i)
#if args.f is None:
#    minf = float(0.01)
#else:
#    minf = float(args.f)
#if minf > 1 or minf < 0:
#    raise NameError('minimum reciprocal overlap must be between 0 and 1.0')
if args.pop is not None:
    pop = gzip.open(args.pop, 'r')
    popbed = BedTool(pop)
else:
    raise NameError('please provide path to population AF bed file with -pop')
if args.minf is not None:
    minf = float(args.minf)
    if minf > 1 or minf < 0:
        raise NameError('minimum frequency must be between 0 and 1.0')

def subtract_overlaps(key, values):
#    tmpA = open('tmpA.bed', 'w')
#    tmpB = open('tmpB.bed', 'w')
    sv = key.split(':')
    tmpA = [(sv[0],sv[1],sv[2],sv[3])]
    tmpAbed = BedTool(tmpA)
    tmpBbed = BedTool(values)
    AN = int(0)
    AC = int(0)
    AFs = []
    for v in values:
        AN += int(v[4])
        AC += int(v[5])
        AFs.append(float(v[6]))
    maxaf = max(AFs)
    subtract = tmpAbed.subtract(tmpBbed)
    for interval in subtract:
        cut = 'cut:' + '_'.join(sv[0:3]) + ':' + str(maxaf) + ':' + str(AC) + ':' + str(AN)
        svlist = [interval[0],interval[1],interval[2],interval[3],cut]
        uniques[interval[3]].append(svlist)

#vcf = cyvcf2.VCF(sys.argv[1], gts012=True)
#popbed = BedTool(sys.argv[2])#.split(',')
tmp = []
for v in vcf:
    chrom = v.CHROM
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    start = int(v.POS)-1
    end = v.INFO.get('END')
    if v.INFO.get('SVTYPE') is not None:
        svtype=v.INFO.get('SVTYPE')
    if svtype == 'BND':
        end = int(v.POS)
    if args.ci == 'out':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[0])
        end = end + int(ciend[1])
    if args.ci == 'in':
        cipos = v.INFO.get('CIPOS95')
        ciend = v.INFO.get('CIEND95')
        start = start + int(cipos[1])
        end = end + int(ciend[0])
        if end < start:
            end = start + 1
    if svtype != 'BND':# or svtype == 'DUP' or svtype == 'INV':
#        if end > pos:
        out = [str(chrom), str(start), str(end), svtype]
        tmp.append(out)

tmpbed = BedTool(tmp)
intersect = tmpbed.intersect(popbed, wao = True) #, f = minf, r = True)
uniques = {}
overlaps = {}
uniqs = 0
overs = 0
for interval in intersect:
    chrom1,start1,end1,type1 = interval[0:4]
    chrom2,start2,end2,type2,af,ac,an,source = interval[4:12]
    sv = [str(chrom1),str(start1),str(end1),type1]
    if type1 not in uniques:
        uniques[type1] = []
    if chrom2 == '.':
        sv.append('uncut:' + '_'.join([str(chrom1),str(start1),str(end1)]))
        uniques[type1].append(sv)
        uniqs += 1
    else:
        key = ':'.join(sv)
        if key not in overlaps:
            overlaps[key] = []
        if type1 == type2 and args.minf is None:
            overlap = [str(chrom2),str(start2),str(end2),type2,str(an),str(ac),str(af)]
            overlaps[key].append(overlap)
            overs += 1
        elif type1 == type2 and args.minf is not None and float(af) > minf:
            overlap = [str(chrom2),str(start2),str(end2),type2,str(an),str(ac),str(af)]
            overlaps[key].append(overlap)
            overs += 1

for key in overlaps:
    if len(overlaps[key]) == 0:
        sv = key.split(':')
        svlist = [str(sv[0]),str(sv[1]),str(sv[2]),sv[3],'uncut:' + '_'.join([str(sv[0]),str(sv[1]),str(sv[2])])]
        uniques[sv[3]].append(svlist)
        uniqs += 1
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
counts = 0
for interval in mergedbed.sort():
    chrom, start, end, svtypes, chunks = interval
    svtype = list(set(svtypes.split(',')))
    chunk = list(set(chunks.split(',')))
    out = [str(chrom), str(start), str(end), ','.join(svtype), ','.join(chunk)]
    final.append(out)
    counts += 1

finalbed = BedTool(final)
finalbed.sort().saveas('uniques.bed')

print uniqs, 'completely unique SVs with no overlaps'
print overs, 'matching SV type overlaps found'
print counts, 'unique SV or SV chunks after merging listed in uniques.bed'
