
# Atividade Amostras Extras: Pipeline de Variantes Somáticas

## Objetivo

Identificar variantes somáticas em amostras pareadas de tumor e normal, utilizando o GATK Mutect2, com posterior estimativa de contaminação, filtragem das chamadas e identificação das variantes elegíveis para anotação.

O pipeline contempla as seguintes etapas:

1. Preparação da referência genômica.
2. Chamada de variantes somáticas com Mutect2.
3. Estimativa de contaminação.
4. Filtragem das variantes.
5. Identificação das variantes classificadas como PASS.

## Contextualização

### Repositório de origem dos arquivos BAM

https://github.com/renatopuga/somatico

### Ferramentas utilizadas

- Bash
- SAMtools
- GATK v4.2.2.0
- Referência genômica hg19 (cromossomos 13 e 19)
- gnomAD como recurso de variantes germinativas

### Amostras analisadas

| Tumor | Normal pareado |
|-------|----------------|
| WP190 | WP191 |
| WP017 | WP018 |

---

## 1. Preparação da referência genômica

### Download dos cromossomos 13 e 19

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
```

### Conferência dos arquivos baixados

```bash
ls -lh chr13.fa.gz chr19.fa.gz
```

### Montagem da referência hg19.fa

```bash
zcat chr13.fa.gz chr19.fa.gz | sed 's/^>chr/>/' > hg19.fa
```

### Conferência dos cabeçalhos

```bash
grep '^>' hg19.fa
```

### Indexação da referência

```bash
samtools faidx hg19.fa
```

### Conferência do índice gerado

```bash
cat hg19.fa.fai
```

### Criação do dicionário da referência

```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
```

### Conferência das sequências no dicionário

```bash
grep '^@SQ' hg19.dict
```

### Geração da lista de intervalos

```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs \
  -R hg19.fa \
  -O hg19.interval_list \
  -OT ACGT
```

### Confirmação dos intervalos nos dois cromossomos

```bash
awk '!/^@/ {print $1}' hg19.interval_list | sort -n | uniq -c
```

---

## 2. Análise do par WP190/WP191

### Conferência dos arquivos das amostras extras

```bash
ls -lh *wp190* *wp191* af-only-gnomad-chr13-chr19.vcf.gz*
```

### Confirmação dos IDs das amostras nos BAMs

```bash
samtools view -H tumor_wp190.bam | grep '^@RG'
samtools view -H tumor_wp191.bam | grep '^@RG'
```

### Chamada de variantes somáticas com Mutect2

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  -I tumor_wp191.bam \
  -normal WP191 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP190.somatic.vcf.gz
```

### Geração das tabelas de contagem do tumor e do normal

```bash
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
```

### Conferência dos cabeçalhos e das primeiras contagens

```bash
head -n 5 WP190.table WP191.table
```

### Estimativa de contaminação no tumor WP190

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP190.table \
  -matched WP191.table \
  -O WP190.contamination.table
```

### Conferência do resultado

```bash
cat WP190.contamination.table
```

### Filtragem das chamadas

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP190.somatic.vcf.gz \
  --contamination-table WP190.contamination.table \
  -O WP190.filtered.vcf.gz
```

### Contagem das chamadas por valor da coluna FILTER

```bash
zcat WP190.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c
```

---

## 3. Análise do par WP017/WP018

### Conferência dos arquivos BAM

```bash
ls -lh *wp017* *wp018*
```

### Confirmação dos IDs nos cabeçalhos

```bash
samtools view -H tumor_wp017.bam | grep '^@RG'
samtools view -H tumor_wp018.bam | grep '^@RG'
```

### Chamada de variantes somáticas com Mutect2

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp017.bam \
  -I tumor_wp018.bam \
  -normal WP018 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19.interval_list \
  -O WP017.somatic.vcf.gz
```

### Geração das tabelas de contagem do tumor e do normal

```bash
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
```

### Conferência dos IDs e das primeiras contagens

```bash
head -n 5 WP017.table WP018.table
```

### Estimativa de contaminação no tumor WP017

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP017.table \
  -matched WP018.table \
  -O WP017.contamination.table
```

### Conferência do resultado

```bash
cat WP017.contamination.table
```

### Filtragem das chamadas

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP017.somatic.vcf.gz \
  --contamination-table WP017.contamination.table \
  -O WP017.filtered.vcf.gz
```

### Contagem das chamadas por valor da coluna FILTER

```bash
zcat WP017.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c
```

---

## 4. Identificação das variantes PASS

Extração das variantes que passaram pelos filtros aplicados pelo FilterMutectCalls nos dois pares de amostras.

```bash
for vcf in WP190.filtered.vcf.gz WP017.filtered.vcf.gz; do
  printf '\n%s\n' "$vcf"
  zcat "$vcf" | awk '/^#CHROM/ || (!/^#/ && $7 == "PASS")'
done
```

As variantes classificadas como PASS constituem as chamadas elegíveis para as etapas subsequentes de anotação e interpretação.

---

## Script completo

O pipeline também está disponível em formato Bash no arquivo [pipeline_somatico.sh](pipeline_somatico.sh).
