#!/usr/bin/env python


# Making CSV-format variation DB for clinicians

import sys
import re


# INPUT -> VCF-file with dbNSFP annotation, fatHMM-MKL output file
# OUTPUT -> a CSV-formatted table with genotypes and predictions for all GATK-called variants
# Note that VCF'd better have SnpEff annotation, as well as dbNSFP predictions for SIFT, Polyphen2HVAR, MutationTaster, 
# Mutation Assessor and FATHMM. No multiallelic entries are allowed!


# Genotype SMART-formatter (0/0:45)
def gt(line):
    content = line.split('\t')
    genotypes = content[9:]
    smart_gt = []
    format = content[8].split(':')
    for i in genotypes:
        data = i.split(':')
        current_gt = data[0]
        for i, par in enumerate(format):
            if par == 'DP':
                current_dp = data[i]
                break
        this_genotype = current_gt + ':' + current_dp
        smart_gt.append(this_genotype)
    return ','.join(smart_gt)


# Pathogenicity annotation converters
sift_converter = {'D': 'Damaging', 'T': 'Tolerated', '-': '-'}
pph_converter = {'D': 'Damaging', 'P': 'Possibly damaging', 'B': 'Benign', '-': '-'}
mt_converter = {'A': 'Disease causing (auto)', 'D': 'Disease causing', 'N': 'Polymorphism', 'P' : 'Polymorphism (auto)', '-': '-'}
ma_converter = {'H': 'High', 'M': 'Medium', 'N': 'Neutral', 'L': 'Low', '-': '-'}
fhmm_converter = sift_converter
clnv_sig = {0: 'Uncertain significance', 1: 'Not provided', 2: 'Benign', 3: 'Likely benign', 4: 'Likely pathogenic', 5: 'Pathogenic', 6: 'Drug response', 7: 'MHC', 255: 'Other'} 

# Getting files
try:
    vcf, fathmm, output = sys.argv[1:]
except:
    print 'Invalid inputs specified\nUsage: csvKnit.py <vcf-file> <fathmm-mkl file> <output file>\n!'
    sys.exit()


# Checking input validity
try:
    with open(vcf, 'r') as a:
        pass
    with open(fathmm, 'r') as b:
        pass
except IOError:
    print 'Invalid inputs specified\nUsage: csvKnit.py <vcf-file> <fathmm-mkl file> <output file>\n!'
    sys.exit()


# Getting fathmm-MKL annotations for both coding and non-coding cases from its output
with open(fathmm, 'r') as fh:
    mkl_ann_coding = {}
    mkl_ann_noncoding = {}
    for line in fh:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        variation = ','.join(content[:4])
        if float(content[4]) > 0.5:
            mkl_ann_coding[variation] = 'Pathogenic'
        else:
            mkl_ann_coding[variation] = 'Benign'
        # Non-coding
        if float(content[6]) > 0.5:
            mkl_ann_noncoding[variation] = 'Pathogenic'
        else:
            mkl_ann_noncoding[variation] = 'Benign'


with open(vcf, 'r') as ifile:
    with open(output, 'w') as ofile:
        ofile.write('Gene,Chromosome,Position,rsID,REF,ALT,EffType,OMIM,ClinVar,SIFT,Polyphen2b,MutationTaster,MutationAssessor,FATHMM,fathmm-MKL,AF,')
        for line in ifile:
            if line.startswith('#'):
                if 'CHROM' in line:
                    content = line.strip().replace('#', '').split('\t')
                    samples = content[9:]
                    ofile.write(','.join(samples) + '\n')
                continue
            content = line.split('\t')
            # Checking once more for input validity
            if ',' in content[4]:
                print 'No multiallelic entries are allowed in the VCF-file!'
                sys.exit()
            # Taking out the dbNSFP annotations
            sift = re.findall('dbNSFP_SIFT_pred=[^:]*?(\w)', line)[0] if 'dbNSFP_SIFT' in line else '-'
            pph = re.findall('dbNSFP_Polyphen2_HVAR_pred=[^:]*?(\w)', line)[0] if 'dbNSFP_Polyp' in line else '-'
            mt = re.findall('dbNSFP_MutationTaster_pred=[^:]*?(\w)', line)[0] if 'dbNSFP_MutationTaster' in line else '-'
            ma = re.findall('dbNSFP_MutationAssessor_pred=[^:]*?(\w)', line)[0] if 'dbNSFP_MutationAssessor' in line else '-'
            fhmm = re.findall('dbNSFP_FATHMM_pred=[^:]*?(\w)', line)[0] if 'dbNSFP_FATHMM' in line else '-'
            # Getting FATHMM-MKL out of the pocket
            variation = ','.join(content[:2]) + ',' + content[3][0] + ',' + content[4][0]
            if variation in mkl_ann_coding:
                mkl = mkl_ann_coding[variation] if 'MODIFIER' not in line else mkl_ann_noncoding[variation]
            else:
                mkl = '-'
            # Grepping allele frequency for given variant and making final csv-string
            af = re.findall('AF=(\d.\d+)', line)[0]
            snpeff = re.findall('EFF=([A-Z0-9_]+)', line)[0]
            gene = content[7].split('|')[4]
            omim = '"=HYPERLINK(""http://www.omim.org/search?index=entry&sort=score+desc%2C+prefix_sort+desc&start=1&limit=10&search=' +  re.findall('rs\d+', content[2])[0] + '"",""OMIM record"")"' if (';OM' in line and 'rs' in content[2]) else '-'
            clinvar = int(re.findall('CLNSIG=(\d+)', line)[0]) if 'CLNSIG' in line else '-'
            clv_link = '"=HYPERLINK(""http://www.ncbi.nlm.nih.gov/clinvar/?term=' + re.findall('rs\d+', content[2])[0] + '"",""' + clnv_sig[clinvar] + '"")"' if ('CLNSIG' in line and 'rs' in content[2]) else '-'
            content[2] = content[2].replace(';', ':')
            content[2] = re.sub('rs(\d+)(.*)', 
                                '"=HYPERLINK(""http://ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=\g<1>"",""rs\g<1>\g<2>"")"', 
                                 content[2])
            output_csv_string = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (gene, ','.join(content[:5]), snpeff, omim, clv_link, sift_converter[sift], pph_converter[pph], mt_converter[mt], ma_converter[ma], fhmm_converter[fhmm], mkl, af, gt(line))
            ofile.write(output_csv_string)
