DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"


# Import SRA files in qiime2 (16s data)
echo "Import SRA files in qiime2 (16s data)..."
ls $DIR/16s/ | ( grep "16s.qza" > /dev/null && echo "Import already done (16s data)." ) || \
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $DIR/16s/fastq \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $DIR/16s/16s.qza

# Import SRA files in qiime2 (ITS data)
echo "Import SRA files in qiime2 (ITS data)..."
ls $DIR/ITS/ | ( grep "ITS.qza" > /dev/null && echo "Import already done (ITS data)." ) || \
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $DIR/ITS/fastq \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $DIR/ITS/ITS.qza

# DADA2 processing (16s data)
echo "DADA2 processing (16s data)..."
mkdir $DIR/16s/dada2
ls $DIR/16s/dada2 | ( grep -e "FeatureData\[Sequence\]_16s.qza" -e "FeatureTable\[Frequency\]_16s.qza" -e "16s-stats-dada2.qza" > /dev/null \
&& echo "DADA2 processing already done (16s data)." ) || \
qiime dada2 denoise-paired \
 --i-demultiplexed-seqs $DIR/16s/16s.qza \
  --p-trim-left-f 21 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
 --o-representative-sequences $DIR/16s/dada2/FeatureData[Sequence]_16s.qza \
 --o-table $DIR/16s/dada2/FeatureTable[Frequency]_16s.qza \
 --o-denoising-stats $DIR/16s/dada2/16s-stats-dada2.qza

# DADA2 processing (ITS data)
echo "DADA2 processing (ITS data)..."
mkdir $DIR/ITS/dada2
ls $DIR/ITS/dada2 | ( grep -e "FeatureData\[Sequence\]_ITS.qza" -e "FeatureTable\[Frequency\]_ITS.qza" -e "ITS-stats-dada2.qza" > /dev/null \
&& echo "DADA2 processing already done (ITS data)." ) || \
qiime dada2 denoise-paired \
 --i-demultiplexed-seqs $DIR/ITS/ITS.qza \
  --p-trim-left-f 22 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
 --o-representative-sequences $DIR/ITS/dada2/FeatureData[Sequence]_ITS.qza \
 --o-table $DIR/ITS/dada2/FeatureTable[Frequency]_ITS.qza \
 --o-denoising-stats $DIR/ITS/dada2/ITS-stats-dada2.qza

# Taxonomic identification of sequences using the QIIME 2 Scikit-learn algorithm (16s data)
echo "Taxonomic identification of sequences using the QIIME 2 Scikit-learn algorithm (16s data)..."
mkdir $DIR/16s/feature-classifier_classify-sklearn
ls $DIR/16s/feature-classifier_classify-sklearn | ( grep "FeatureData\[Taxonomy\]_16s.qza" \
&& echo "Taxonomic identification already done (16s data)." ) || \
qiime feature-classifier classify-sklearn \
  --i-classifier $DIR/16s/silva-138-99-515-806-nb-classifier.qza \
  --i-reads $DIR/16s/dada2/FeatureData[Sequence]_16s.qza \
  --o-classification $DIR/16s/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_16s.qza

# Taxonomic identification of sequences using the QIIME 2 Scikit-learn algorithm (ITS data)
echo "Taxonomic identification of sequences using the QIIME 2 Scikit-learn algorithm (ITS data)..."
mkdir $DIR/ITS/feature-classifier_classify-sklearn
ls $DIR/ITS/feature-classifier_classify-sklearn | ( grep "FeatureData\[Taxonomy\]_ITS.qza" > /dev/null \
&& echo "Taxonomic identification already done (ITS data)." ) || \
qiime feature-classifier classify-sklearn \
  --i-classifier $DIR/ITS/unite_ver9_99_all_29.11.2022-Q2-2023.5.qza \
  --i-reads $DIR/ITS/dada2/FeatureData[Sequence]_ITS.qza \
  --o-classification $DIR/ITS/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_ITS.qza


echo $DIR > dir.txt

# Processing 16s data using R
Rscript $DIR/16s/preprocess_16s.r
# Rscript $DIR/16s/barplot_16s.r
# Rscript $DIR/16s/upset_16s.r
# Rscript $DIR/16s/heatmap_16s.r
# Rscript $DIR/16s/alpha_16s.r
# Rscript $DIR/16s/beta_16s.r

# Processing ITS data using R
Rscript $DIR/ITS/preprocess_ITS.r
# Rscript $DIR/ITS/barplot_ITS.r
# Rscript $DIR/ITS/upset_ITS.r
# Rscript $DIR/ITS/heatmap_ITS.r
# Rscript $DIR/ITS/alpha_ITS.r
# Rscript $DIR/ITS/beta_ITS.r

# Microbiome network analysis
Rscript $DIR/network.r

rm dir.txt