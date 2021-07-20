#!/home/guotingyou/anaconda3/bin/
### Filename: True.step0_variant_filtration.sh ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

samtools='/zfssz2/ST_MCHRI/BIGDATA/USER/guotingyou/anaconda3/pkgs/samtools-1.9/samtools'
fasta='/hwfssz5/ST_MCHRI/REPRO/PROJECT/P18Z10200N0183/xujinjin/Direction/xujinjin/database/GRch37_index/GRch37.fasta'
out_name="cfDNA"
in_bam=$1
out_path=$2
position=$3
chrom=$4
out_file="${out_path}/${out_name}.${chrom}.pileup"

echo Start: `date`

time $samtools mpileup \
    -f $fasta \
    -r $chrom \
    -l $position \
    -C 50 \
    -q 30 \
    -Q 20 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP \
    --output ${out_file} \
    ${in_bam}

time gzip -f ${out_file}

echo End: `date`
echo "Practice_makes_perfect"
