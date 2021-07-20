import os
import sys
import subprocess

binPath = f'{os.getcwd()}/cfDNA_analysis/bin'
#print("Hello World!")
#os.popen(f'python3 {binPath}/maf_based_fetal_fraction.py -r chr21 -bam graduate/JK-10/input/bam/TDP1807015083_cfDNA.recaled.bam --vcf graduate/JK-10/input/vcf/cfDNA/TDP1807015083_cfDNA.output-hc.vcf.filter.vcf.gz')
#cmd = f'python3 {binPath}/maf_based_fetal_fraction.py -r chr1 -bam graduate/JK-10/input/bam/TDP1807015083_cfDNA.recaled.bam --vcf small.vcf.gz'
#cmd = f'python3 {binPath}/nucleosome_score.py -bam graduate/JK-10/input/bam/TDP1807015083_cfDNA.recaled.bam -g /zfssz2/ST_MCHRI/BIGDATA/USER/guotingyou/youyou/hkgenes.all.txt'
cmd = f'python3 {binPath}/nucleosome_score.py -bam task/task5_transplantation/lung_3demo/bam_demo_30X/HZL_D10_18899.bam -g mhc.txt'
subprocess.run(cmd.split())
