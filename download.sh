DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"


# Download the SILVA 138 pre-trained classifier for 16s sequences (99% OTUs from V4 region of sequences) 
# from https://docs.qiime2.org/2023.5/data-resources/ (accessed on 14 July 2023)
ls $DIR/16s/ | ( grep silva-138-99-515-806-nb-classifier.qza > /dev/null && echo "SILVA 138 pre-trained classifier already downloaded" ) \
|| ( echo "Download the SILVA 138 pre-trained classifier" && wget -P $DIR/16s/ https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza ) \
|| ( echo "Can't download the SILVA 138 pre-trained classifier." && exit )

# Download the UNITE pre-trained classifier for ITS sequences (99% OTUs from ITS1f/ITS2 region of sequences)
# from https://github.com/colinbrislawn/unite-train/releases (accessed on 14 July 2023)
ls $DIR/ITS/ | ( grep unite_ver9_99_all_29.11.2022-Q2-2023.5.qza > /dev/null && echo "UNITE pre-trained classifier already downloaded" ) \
|| ( echo "Download the UNITE pre-trained classifier" && wget -P $DIR/ITS/ https://github.com/colinbrislawn/unite-train/releases/download/9.0-qiime2-2023.5-demo/unite_ver9_99_all_29.11.2022-Q2-2023.5.qza ) \
|| ( echo "Can't download the UNITE pre-trained classifier." && exit )


# Download SRA files used in this article (16s data)
mkdir $DIR/16s/fastq

PRJNA980748_16s=(SRR24918031 SRR24918054 SRR24917972 SRR24918002 SRR24917969 SRR24918027 SRR24918024 SRR24917978 \
SRR24917958 SRR24917956 SRR24918053 SRR24918051 SRR24918049 SRR24918015 SRR24918013 SRR24918010 \
SRR24918042 SRR24918040 SRR24917976 SRR24917974 SRR24917955 SRR24917953 SRR24917951 SRR24917949 \
SRR24918047 SRR24918044 SRR24917997 SRR24917995 SRR24918006 SRR24918004 SRR24918001 SRR24917947 \
SRR24917945 SRR24918067 SRR24918065 SRR24917998 SRR24918039 SRR24918037 SRR24917993 SRR24917971 \
SRR24917968 SRR24917966 SRR24917964 SRR24918062 SRR24918060 SRR24918057 SRR24917989 SRR24918036 \
SRR24918034 SRR24918032)
PRJNA998468_16s=(SRR25433230 SRR25433229 SRR25433218 SRR25433207 SRR25433204 SRR25433203 SRR25433202 SRR25433201 \
SRR25433200 SRR25433199 SRR25433228 SRR25433227 SRR25433226 SRR25433225 SRR25433224 SRR25433223)
all_SRA=(${PRJNA980748_16s[*]} ${PRJNA998468_16s[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/16s/fastq/ && fasterq-dump $DIR/16s/fastq/$sra/$sra.sra -O $DIR/16s/fastq \
&& gzip $DIR/16s/fastq/$sra"_1.fastq" \
&& gzip $DIR/16s/fastq/$sra"_2.fastq" \
&& rm -R $DIR/16s/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Download SRA files used in this article (ITS data)

mkdir $DIR/ITS/fastq

PRJNA980748_ITS=(SRR24918026 SRR24918011 SRR24918045 SRR24917999 SRR24918058 SRR24918030 SRR24917979 SRR24917959 \
SRR24917957 SRR24918055 SRR24918052 SRR24918050 SRR24918048 SRR24918014 SRR24918012 SRR24918009 \
SRR24918041 SRR24917977 SRR24917975 SRR24917973 SRR24917954 SRR24917952 SRR24917950 SRR24917948 \
SRR24918046 SRR24918043 SRR24917996 SRR24918007 SRR24918005 SRR24918003 SRR24918000 SRR24917946 \
SRR24918068 SRR24918066 SRR24918064 SRR24918056 SRR24918038 SRR24917994 SRR24917992 SRR24917970 \
SRR24917967 SRR24917965 SRR24918063 SRR24918061 SRR24918059 SRR24917990 SRR24917988 SRR24918035 \
SRR24918033 SRR24917991)
PRJNA998468_ITS=(SRR25433222 SRR25433221 SRR25433220 SRR25433219 SRR25433217 SRR25433216 SRR25433215 SRR25433214 \
SRR25433213 SRR25433212 SRR25433211 SRR25433210 SRR25433209 SRR25433208 SRR25433206 SRR25433205)
all_SRA=(${PRJNA980748_ITS[*]} ${PRJNA998468_ITS[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/ITS/fastq/ && fasterq-dump $DIR/ITS/fastq/$sra/$sra.sra -O $DIR/ITS/fastq \
&& gzip $DIR/ITS/fastq/$sra"_1.fastq" \
&& gzip $DIR/ITS/fastq/$sra"_2.fastq" \
&& rm -R $DIR/ITS/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Rename SRA for input in qiime2 (16s data)
data_16s_1=($(ls $DIR/16s/fastq | grep _1.fastq.gz))
for oldname in ${data_16s_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

data_16s_2=($(ls $DIR/16s/fastq | grep _2.fastq.gz))
for oldname in ${data_16s_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

# Rename SRA for input in qiime2 (ITS data)
data_ITS_1=($(ls $DIR/ITS/fastq | grep _1.fastq.gz))
for oldname in ${data_ITS_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done

data_ITS_2=($(ls $DIR/ITS/fastq | grep _2.fastq.gz))
for oldname in ${data_ITS_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done