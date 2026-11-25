#!/bin/bash
# Download dos cromossomos 13 e 19

wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz

# Conferindo se os dois arquivos foram salvos:

ls -lh chr13.fa.gz chr19.fa.gz

# Montar a referência hg19.fa

zcat chr13.fa.gz chr19.fa.gz | sed 's/^>chr/>/' > hg19.fa

# conferir os dois cabeçalhos de uma vez:

grep '^>' hg19.fa

# Indexar a referencia

samtools faidx hg19.fa

# Conferir o indice gerado

cat hg19.fa.fai

# Criar o dicionário da referencia

./gatk-4.2.2.0/gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict

# conferir as duas sequências no dicionário:

grep '^@SQ' hg19.dict

# Gerar a lista de intervalos

./gatk-4.2.2.0/gatk ScatterIntervalsByNs \
  -R hg19.fa \
  -O hg19.interval_list \
  -OT ACGT

# confirmar que os dois cromossomos têm intervalos

awk '!/^@/ {print $1}' hg19.interval_list | sort -n | uniq -c

# Conferir os arquivos das amostras extras

ls -lh *wp190* *wp191* af-only-gnomad-chr13-chr19.vcf.gz*

# Confirmar os IDs das amostras nos BAMs

samtools view -H tumor_wp190.bam | grep '^@RG'
samtools view -H tumor_wp191.bam | grep '^@RG'

# Rodar o Mutect2 para o par WP190/WP191

./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  -I tumor_wp191.bam \
  -normal WP191 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP190.somatic.vcf.gz

# Gerar as tabelas de contagem do tumor e do normal

./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp190.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP190.table

./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp191.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP191.table

# Conferindo os cabeçalhos e as primeiras contagens

head -n 5 WP190.table WP191.table

# Estimar a contaminação no tumor WP190

./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP190.table \
  -matched WP191.table \
  -O WP190.contamination.table

# Conferindo o resultado:

cat WP190.contamination.table

# Filtrar as chamadas do par WP190/WP191

./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP190.somatic.vcf.gz \
  --contamination-table WP190.contamination.table \
  -O WP190.filtered.vcf.gz

# contar as chamadas por valor da coluna FILTER:

zcat WP190.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c

# Conferir os BAMs do par WP017/WP018

ls -lh *wp017* *wp018*

# confira os IDs nos cabeçalhos:

samtools view -H tumor_wp017.bam | grep '^@RG'
samtools view -H tumor_wp018.bam | grep '^@RG'

# Rodar o Mutect2 para WP017/WP018

./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp017.bam \
  -I tumor_wp018.bam \
  -normal WP018 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP017.somatic.vcf.gz

# Gerar as tabelas de contagem do tumor e do normal

./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp017.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP017.table

./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp018.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP018.table

# Conferindo os IDs e as primeiras contagens:

head -n 5 WP017.table WP018.table

# Estimar a contaminação no tumor WP017

./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP017.table \
  -matched WP018.table \
  -O WP017.contamination.table

# Conferindo o resultado:

cat WP017.contamination.table

# Filtrar as chamadas do par WP017/WP018

./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP017.somatic.vcf.gz \
  --contamination-table WP017.contamination.table \
  -O WP017.filtered.vcf.gz

# contar as chamadas por valor de FILTER

zcat WP017.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c

# Identificar as três chamadas PASS

for vcf in WP190.filtered.vcf.gz WP017.filtered.vcf.gz; do
  printf '\n%s\n' "$vcf"
  zcat "$vcf" | awk '/^#CHROM/ || (!/^#/ && $7 == "PASS")'
done
