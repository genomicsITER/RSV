<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division -->
<!-- of the Institute of Technology an Renewable Energy (ITER) -->
<!-- Tenerife, Canary Islands, SPAIN -->
<!-- See the "Contact us" section to collaborate with us to growth -->
<!-- this repository. ;=) -->

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/RSV" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/logos_GH.png" width="auto" /> 
      </a>
</p>

# RSV
A public repository of resources related to RSV analysis maintained by ITER.

This is the result of an ongoing joint-effort of the following institutions and laboratories:
<ul>
 <li><a href="https://www3.gobiernodecanarias.org/sanidad/scs/organica.jsp?idCarpeta=10b3ea46-541b-11de-9665-998e1388f7ed">Servicio de Microbiología, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://fciisc.org/">Fundación Canaria Instituto de Investigación Sanitaria de Canarias</a> at the Research Unit, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://portalciencia.ull.es/grupos/6361/detalle">Laboratorio de Inmunología Celular y Viral</a>, Unidad de Farmacología, Facultad de Medicina, Universidad de La Laguna, 38200 San Cristóbal de La Laguna, Spain.</li>
 <li><a href="https://www.iter.es/areas/area-genomica/">Genomics Division, Instituto Tecnológico y de Energías Renovables</a>, 38600 Santa Cruz de Tenerife, Spain</li>
</ul>
<hr>
<!-- ------------------ SECTION ------------------ -->

# Code for Illumina short-reads processing #
## Table of contents ##
<ul>
  <li><a href="#1">1. Quality control of raw reads with FastQC</a></li>
  <li><a href="#2">2. Taxonomic classification of raw reads using Kraken2 with PlusPF database</a></li>
  <li><a href="#3">3. Adapter trimming with fastp</a></li>
  <li><a href="#4">4. Quality control of trimmed reads with FastQC</a></li>
  <li><a href="#5">5. Remove host reads using Kraken2 with humanDB database</a></li>
  <li><a href="#6">6. Multi-reference alignment with bbmap to select reference</a></li>
  <li><a href="#7">7. <i>De novo</i> assembly</a></li> 
  <li><a href="#8">8. Align reads</a></li>
  <li><a href="#9">9. Trim adapters with iVar before create consensus sequences</a></li>
  <li><a href="#10">10. Coverage analysis with Mosdepth</a></li>
  <li><a href="#11">11. Create consensus sequence with iVar</a></li>
  <li><a href="#12">12. Variant-calling with iVar</a></li>
  <li><a href="#13">13. Lineage classification with NextClade</a></li>
  <li><a href="#14">14. Multisample Alignment and Phylogenetic Analysis</a></li>
</ul>


<hr>

<a name="1"></a>
#### 1. Quality control of raw reads with FastQC:
```Bash
# In case you installed FastQC through conda:
conda activate fastqc

r1="sample_R1_001.fastq.gz"
r2="sample_R2_001.fastq.gz"
outdir="./fastqc_results"

# For each sample run following commands:
fastqc ${r1} --outdir ${outdir}
fastqc ${r2} --outdir ${outdir}

conda deactivate
```

<a name="2"></a>
#### 2. Taxonomic classification of raw reads using Kraken2 with PlusPF database:
```Bash
# In case you installed Kraken2 through conda:
conda activate kraken2

database="/path/to/kraken2/databases/k2_pluspf_20230605"

r1="sample_R1_001.fastq.gz"
r2="sample_R2_001.fastq.gz"

report="sample.kraken2-report.txt"

kraken2 --db ${database} \
  --report ${report} \
  --paired ${r1} ${r2}

# Reports could be loaded to Pavian server (https://fbreitwieser.shinyapps.io/pavian/)
# for an interactive analysis of the classification results.

conda deactivate
```

<a name="3"></a>
#### 3. Adapter trimming with fastp:
```Bash
# In case you installed fastp through conda:
conda activate fastp

threads=16

primers_fasta="RSV_A_B.primers.fasta"

r1="sample_R1_001.fastq.gz"
r2="sample_R2_001.fastq.gz"

out1="sample_R1_001.fastp.fastq.gz"
out2="sample_R2_001.fastp.fastq.gz"

json="sample.fastp.json"
html="sample.fastp.html"

fastp --thread ${threads} \
  --in1 ${r1} \
  --in2 ${r2} \
  --out1 ${out1} \
  --out2 ${out2} \
  --json ${json} \
  --html ${html} \
  --adapter_fasta ${primers_fasta}

conda deactivate
```

<a name="4"></a>
#### 4. Quality control of trimmed reads with FastQC:
```Bash
# In case you installed FastQC through conda:
conda activate fastqc

trimmed_r1="sample_R1_001.fastp.fastq.gz"
trimmed_r2="sample_R2_001.fastp.fastq.gz"
outdir="./fastqc_trimmed_results"

fastqc ${trimmed_r1} --outdir ${outdir}
fastqc ${trimmed_r2} --outdir ${outdir}

conda deactivate
```

<a name="5"></a>
#### 5. Remove host reads using Kraken2 with humanDB database:
```Bash
# In case you installed Kraken2 through conda:
conda activate kraken2

database="/path/to/kraken2/databases/kraken2-human-db"

trimmed_r1="sample_R1_001.fastp.fastq.gz"
trimmed_r2="sample_R2_001.fastp.fastq.gz"

classified="sample.classified#.fastq"
unclassified="sample.unclassified#.fastq"

report="sample.kraken2-report.txt"

kraken2 --db ${database} \
  --unclassified-out ${unclassified} \
  --classified-out ${classified} \
  --report ${report} \
  --paired ${r1} ${r2}

conda deactivate
```

<a name="6"></a>
#### 6. Multi-reference alignment with bbmap to select reference:
```Bash
# In case you installed bbmap through conda:
conda activate bbmap

# Reference with both RSV subtypes:
multifasta="hRSV_A_and_B.fasta"

r1="sample.unclassified_1.fastq"
r2="sample.unclassified_2.fastq"

out="sample.multifasta_alignment.sam"
cov_stats="sample.multifasta_alignment_bbmap_covstats.tsv"
all_stats="sample.multifasta_alignment_bbmap_allstats.txt"

threads=16

bbmap.sh \
  in=${r1} \
  in2=${r2} \
  outm=${out} \
  ref=${ref} \
  threads=${threads} \
  covstats=${cov_stats} \
  local=true interleaved=false maxindel=80 -Xmx32g > ${all_stats} 2>&1

# Select reference using this script from Graniger-Lab's Revica pipeline:
# https://github.com/greninger-lab/revica/blob/main/bin/select_reference.py
python3 select_reference.py \
  -bbmap_covstats ${covstats} \
  -b sample \
  -reads_num ${reads_used} \
  -mapped_reads ${mapped_reads} \
  -m 3 \
  -p 0

# This script will generate a file with the following name structure:
# sample_accession_reference_vid.txt

conda deactivate
```

<a name="7"></a>
#### 7. <i>De novo</i> assembly:
```Bash

# Unicycler uses SPAdes to assemble the Illumina reads into an assembly graph
# In case you installed Unicycler through conda:
conda activate spades

threads=16

#Define scrubbed FASTQ files (files without human-derived reads)
r1_scrubbed="sample.scrubbed_1.fastq"
r2_scrubbed="sample.scrubbed_1.fastq"

asmdir=${outdir}/assembly
mkdir -p ${asmdir}

spades.py -t ${threads} -1 ${r1_scrubbed} -2 ${r2_scrubbed} -o ${asmdir}

conda deactivate
```     

<a name="8"></a>
#### 8. Align reads:
```Bash
# Select reference file:
selected_reference=$( awk '{ print $3 }' sample_accession_reference_vid.txt )

if [ ${selected_reference} -eq "hRSV_A_England_397_2017" ]; then
  reference=/path/to/hRSV_A_England_397_2017.fasta
  type="RSV-A"
elif [ ${selected_reference} -eq "hRSV_B_Australia_VIC-RCH056_2019" ]; then
  reference=/path/to/hRSV_B_Australia_VIC-RCH056_2019.fasta
  type="RSV-B"
fi

r1=${sampledir}/${sample}.fastp.unclassified_1.fastq
r2=${sampledir}/${sample}.fastp.unclassified_2.fastq

# Align:
outfile=sample.aligned-to-${type}.sam
bwa mem -Y ${reference} ${r1} ${r2} > ${outfile}

# Convert to BAM, sort and index:
infile=sample.aligned-to-${type}.sam
outfile=sample.aligned-to-${type}.sorted.bam
samtools sort ${infile} > ${outfile}
samtools index ${outfile}

# QCs:
infile=sample.aligned-to-${type}.sorted.bam
outfile1=sample.aligned-to-${type}.sorted.flagstat
outfile2=sample.aligned-to-${type}.sorted.idxstats
samtools flagstat ${infile} > ${outfile1}
samtools idxstats ${infile} > ${outfile2}

# Discard unmapped reads
infile=sample.aligned-to-${type}.sorted.bam
outfile=sample.aligned-to-${type}.sorted.mapped.bam
samtools view -F 0x04 -b ${infile} > ${outfile}
samtools index ${outfile}
```

<a name="9"></a>
#### 9. Trim adapters with iVar before create consensus sequences:
```Bash
# In case you installed iVar through conda:
conda activate ivar

primers_bed="RSV_A_B.primers.bed"

# Adapter trimming:
infile=sample.aligned-to-${type}.sorted.mapped.bam
prefix=sample.aligned-to-${type}.sorted.mapped.trimmed
ivar trim -e \
  -i ${infile} \
  -b ${primer_bed} \
  -p ${prefix}

# Sort and index:
infile=sample.aligned-to-${type}.sorted.mapped.trimmed.bam
outfile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.bam
samtools sort ${infile} > ${outfile}
samtools index ${outfile}

conda deactivate
```

<a name="10"></a>
#### 10. Coverage analysis with Mosdepth and Samtools:
```Bash
# In case you installed Mosdepth through conda:
conda activate mosdepth

infile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.bam
prefix=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted
mosdepth --threads 4 --thresholds 1,10,20,50,100,500,1000 ${prefix} ${infile}

conda deactivate

infile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.bam
outfile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.depth
samtools depth -a ${infile} > ${outfile}
```

<a name="11"></a>
#### 11. Create consensus sequence with iVar:
```Bash
# In case you installed iVar through conda:
conda activate ivar

# Minimum depth to call consensus:
min_depth=5

infile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.bam
prefix=sample.aligned-to-${type}.ivar_consensus

samtools mpileup -A -Q 0 ${infile} | ivar consensus -p ${prefix} -q 10 -t 0 -m ${min_depth}

conda deactivate
```

<a name="12"></a>
#### 12. Variant-calling with iVar:
```Bash
# In case you installed iVar through conda:
conda activate ivar

infile=sample.aligned-to-${type}.sorted.mapped.trimmed.sorted.bam
prefix=sample.aligned-to-${type}.ivar_calling

samtools mpileup --reference ${reference} ${infile} | ivar variants -r ${reference} -p ${prefix}

conda deactivate
```

<a name="13"></a>
#### 13. Lineage classification with NextClade:
```Bash
strain=$( cat sample.strain.txt )

# Create a multisample FASTA with more sequences separated by reference:
sequences_dir="/path/to/more/sequences"
outfile="multifasta.${type}.fa"

cat ${sequences_dir}/*.aligned-to-${type}.ivar_consensus.fa ${infile} > ${outfile}

# You can upload this multifasta file to Nextclade (https://clades.nextstrain.org/) and select the
# correct pathogen to run clade assignment, mutation calling, and sequence quality checks.
```

<a name="14"></a>
#### 14. Multisample Alignment and Phylogenetic Analysis:
```Bash
# Align multisample FASTA using MAFFT:
# You should replace the --thread parameter with a proper value that suits your execution environment.
infile="multifasta.${type}.fa"
outfile="multifasta.${type}.mafft-aligned.fa"
logfile="multifasta.${type}.mafft-aligned.log"
mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 48 ${infile} 1> ${outfile} 2> ${logfile}

# Generate phylogenetic tree using IQ-TREE with best-fit model:
# Alternatively, you could choose a model like JC by using the -m parameter.
iqtree -s ${outfile} -nt AUTO

# You will end up with a phylogenetic tree in Newick format (*.treefile extension) which can be represented 
# using your favorite tool such as NCBI Tree Viewer.
```
  
<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
