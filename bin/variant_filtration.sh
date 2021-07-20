### Required Arguments ###

vcfFile=$1

### Optional Arguments ###

outdir=$2
if [ ! -n "$outdir" ]; then
    outdir=`pwd`;
elif [ ! -d $outdir ]; then
    echo "$outdir does not exist."
    exit
fi

prefix=$3
if [ ! -n "$prefix" ]; then
    perfix=$(echo `basename $vcfFile` | awk -F . '{print $1}')
fi

### Preparation ###

export PATH="/zfssz2/ST_MCHRI/BIGDATA/USER/qianxiaobo/southwest/DNA_Human_WGS_2018a/softs/java/jre1.8.0_101/bin:$PATH"
gatk=/zfssz2/ST_MCHRI/BIGDATA/USER/guotingyou/anaconda3/pkgs/gatk-4.1.4.0/gatk
tabix=/zfssz5/BC_PS/huangshujia/biosoft/extra/local/bin/tabix

# Compressed VCF should decompress
isCompressed=$(
    if ( file $vcfFile | grep -q "symbolic link" ); then
        if ( file $(readlink $1) | grep -q "gzip compressed" ); then
            echo "True"
        else
            echo "False"
        fi
    elif ( file $vcfFile | grep -q "gzip compressed" ); then
        echo "True"
    else
        echo "False"
    fi
)
if [ $isCompressed = "True" ]; then
    gzip -d --force $vcfFile && vcfFile=${vcfFile%.gz*}
fi

############
### Main ###
############

if [ -s "${outdir}/${prefix}.filtered.vcf.gz.finish" ]; then
    echo "[Exit] ${outdir}/${prefix}.filtered.vcf.gz.finish exists."
    exit
fi

# Step1-1: Select variants of type SNP
if [ -s "${outdir}/${prefix}.snp.vcf.gz" ]; then
    echo "[Skip] Step1-1: Select variants of type SNP"
else
    echo "[SNP Selection] Start `date`"
    time $gatk SelectVariants \
        -select-type SNP \
        -V $vcfFile \
        -O ${outdir}/${prefix}.snp.vcf.gz
    echo "[SNP Selection] End `date`"
fi

# Step1-2: Hard filtration for SNPs
if [ -s "${outdir}/${prefix}.snp.filter.vcf.gz" ]; then
    echo "[Skip] Step1-2: Hard filtration for SNPs"
else
    echo "[SNP Filtration] Start `date`"
    time $gatk VariantFiltration \
        -V ${outdir}/${prefix}.snp.vcf.gz \
        --filter-expression "QD < 2.0 || MQ < 50.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "SnpFilter" \
        -O ${outdir}/${prefix}.snp.filter.vcf.gz
    echo "[SNP Filtration] End `date`"
fi

# Step2-1: Select variants of type Indel
if [ -s "${outdir}/${prefix}.indel.vcf.gz" ]; then
    echo "[Skip] Step2-1: Select variants of type Indel"
else
    echo "[Indel Selection] Start `date`"
    time $gatk SelectVariants \
        -select-type INDEL \
        -V $vcfFile \
        -O ${outdir}/${prefix}.indel.vcf.gz
    echo "[Indel Selection] End `date`"
fi

# Step2-2: Hard filtration for Indels
if [ -s "${outdir}/${prefix}.indel.filter.vcf.gz" ]; then
    echo "[Skip] Step2-2: Hard filtration for Indels"
else
    echo "[Indel Filtration] Start `date`"
    time $gatk VariantFiltration \
        -V ${outdir}/${prefix}.indel.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "IndelFilter" \
        -O ${outdir}/${prefix}.indel.filter.vcf.gz
    echo "[Indel Filtration] End `date`"
fi 

# Step3: Merge
if [ -s "${outdir}/${prefix}.filtered.vcf.gz" ]; then
    echo "[Skip] Step3: Merge"
elif [ ! -s "${outdir}/${prefix}.snp.filter.vcf.gz" ]; then
    echo "[Error] ${outdir}/${prefix}.snp.filter.vcf.gz does not exist."
    exit
elif [ ! -s "${outdir}/${prefix}.indel.filter.vcf.gz" ]; then
    echo "[Error] ${outdir}/${prefix}.indel.filter.vcf.gz does not exist."
    exit
else
    echo "[Merge] Start `date`"
    time $gatk MergeVcfs \
        -I ${outdir}/${prefix}.snp.filter.vcf.gz \
        -I ${outdir}/${prefix}.indel.filter.vcf.gz \
        -O ${outdir}/${prefix}.filtered.vcf.gz
    echo "[Merge] End `date`"
fi

# Step4: Index and Final check
if [ ! -s "${outdir}/${prefix}.filtered.vcf.gz" ]; then
    echo "[Error] ${outdir}/${prefix}.filtered.vcf.gz does not exist."
    exit
else
    rm -f ${outdir}/${prefix}.snp.* ${outdir}/${prefix}.indel.*
    $tabix -f -p vcf ${outdir}/${prefix}.filtered.vcf.gz && \
    echo "Practice_makes_perfect" > ${outdir}/${prefix}.filtered.vcf.gz.finish
fi
