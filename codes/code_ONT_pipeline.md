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

# Code for Oxford Nanopore Technologies (ONT) long-reads processing #
## Table of contents ##
<ul>
  <li><a href="#1">1. Quality control of raw reads with NanoPlot</a></li>
  <li><a href="#2">2. Taxonomic classification of raw reads using Kraken2 with PlusPF database</a></li>
  <li><a href="#3">3. Remove host reads using Kraken2 with humanDB database</a></li>
  <li><a href="#4">4. Quality control of dehosted reads with NanoPlot</a></li>
  <li><a href="#5">5. Select reference using IRMA and BLAST</a></li>
  <li><a href="#6">6. Align reads</a></li>
  <li><a href="#7">7. Trim adapters with artic <i>align_trim</i> tool</a></li> 
  <li><a href="#8">8. Coverage analysis</a></li>
  <li><a href="#9">9. Calling with Artic and Medaka</a></li>
  <li><a href="#10">10. Create consensus with Artic and BCFtools</a></li>
  <li><a href="#11">11. Lineage classification with NextClade</a></li>
</ul>


<hr>

<a name="1"></a>
#### 1. Quality control of raw reads with NanoPlot:
```Bash
# In case you installed NanoPlot through conda:
conda activate nanoplot

infile="SAMPLE.fastq.gz"
outdir="./nanoplot_raw"

threads=16

NanoPlot --threads ${threads} --verbose --fastq ${infile} --outdir ${outdir}

conda deactivate
```

<a name="2"></a>
#### 2. Taxonomic classification of raw reads using Kraken2 with PlusPF database:
```Bash
# In case you installed Kraken2 through conda:
conda activate kraken2

database="/path/to/kraken2/databases/k2_pluspf_20230605"

infile="SAMPLE.fastq.gz"
report="SAMPLE.kraken2-report.txt"
log="SAMPLE.kraken2.log"

kraken2 --db ${database} \
    --report ${report} \
    ${infile} &> ${log}

# Reports could be loaded to Pavian server (https://fbreitwieser.shinyapps.io/pavian/)
# for an interactive analysis of the classification results.

conda deactivate
```

<a name="3"></a>
#### 3. Remove host reads using Kraken2 with humanDB database:
```Bash
# In case you installed kraken2 through conda:
conda activate kraken2

database="/path/to/kraken2/databases/kraken2-human-db"

infile="SAMPLE.fastq.gz"

classified="SAMPLE.classified#.fastq"
unclassified="SAMPLE.unclassified#.fastq"

report="SAMPLE.dehosted.kraken2-report.txt"
log="SAMPLE.dehosted.kraken2.log"

kraken2 --db ${database} \
  --unclassified-out ${unclassified} \
  --classified-out ${classified} \
  --report ${report} \
  --paired ${infile} &> ${log}

conda deactivate
```

<a name="4"></a>
#### 4. Quality control of dehosted reads with NanoPlot:
```Bash
# In case you installed NanoPlot through conda:
conda activate nanoplot

infile="SAMPLE.unclassified.fastq"
outdir="./nanoplot_dehosted"

NanoPlot --threads ${threads} --verbose --fastq ${infile} --outdir ${outdir}

conda deactivate
```

<a name="5"></a>
#### 5. Select reference using IRMA and BLAST:
```Bash
export PATH=/path/to/IRMA/flu-amd:$PATH

module=RSV0
external_config=/path/to/IRMA/flu-amd/IRMA_RES/modules/RSV0/config/RSV0-minion.sh

infile="SAMPLE.unclassified.fastq"
outdir="./reference_selection"

IRMA ${module} -c ${external_config} ${infile} ${outdir}/${sample}

# BLAST:
conda activate blast

indir="./reference_selection"
query=$( ls -1 ${indir}/*.fasta )

outfile=${indir}/SAMPLE.blast.txt

blastn -remote -db nt -num_alignments 1 -outfmt "6 stitle score bitscore qcovs" -query ${query} -out ${outfile}

conda deactivate

# Select reference by parsing BLASTn results:
selected_reference=$( sort -t$'\t' -k3 -k4 -nr ${basedir}/${sample}/05_Reference_Selection/${sample}.blast.txt | head -1 )

if [[ ${selected_reference} == *"Human respiratory syncytial virus A"* ]]; then
   echo "RSV-A" > blast_selected_ref
elif [[ ${selected_reference} == *"Human respiratory syncytial virus B"* ]]; then
   echo "RSV-B" > blast_selected_ref
fi

# Reference files:
ref_type=$( cat blast_selected_ref )

if [[ ${ref_type} == "RSV-A" ]]; then
   selected_ref=${refdir}/hRSV_A_England_397_2017.fasta
   primer_bed_collapsed=${refdir}/hRSV_A_England_397_2017.primer_scheme.collapsed.bed
   genes_gf_bed=${refdir}/RSV_A_genes_GF.bed
elif [[ ${ref_type} == "RSV-B" ]]; then
   selected_ref=${refdir}/hRSV_B_Australia_VIC-RCH056_2019.fasta
   primer_bed_collapsed=${refdir}/hRSV_B_Australia_VIC-RCH056_2019.primer_scheme.collapsed.bed
   genes_gf_bed=${refdir}/RSV_B_genes_GF.bed
fi

```

<a name="6"></a>
#### 6. Align reads:
```Bash
# In case you installed minimap2 through conda:
conda activate minimap2

infile="SAMPLE.unclassified.fastq"

# Align:
outfile=SAMPLE.aligned-to-${ref_type}.sam
minimap2 -ax map-ont ${selected_ref} ${infile} > ${outfile}

# Convert to BAM, sort and index:
infile=SAMPLE.aligned-to-${ref_type}.sam
outfile=SAMPLE.aligned-to-${ref_type}.sorted.bam
samtools sort ${infile} > ${outfile}
samtools index ${outfile}

# Quality control:
infile=SAMPLE.aligned-to-${ref_type}.sorted.bam
outfile=SAMPLE.aligned-to-${ref_type}.sorted.flagstat
samtools flagstat ${infile} > ${outfile}

infile=SAMPLE.aligned-to-${ref_type}.sorted.bam
outfile=SAMPLE.aligned-to-${ref_type}.sorted.idxstats
samtools idxstats ${infile} > ${outfile}

# Discard unmapped reads:
infile=SAMPLE.aligned-to-${ref_type}.sorted.bam
outfile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.bam
samtools view -F 0x04 -b ${infile} > ${outfile}
samtools index ${outfile}

conda deactivate
```

<a name="7"></a>
#### 7. Trim adapters with Artic <i>align_trim</i> tool:
```Bash
# In case you installed Artic pipeline through conda (see here: https://artic.readthedocs.io/en/latest/installation/):
conda activate artic

infile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.bam
prefix=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed

align_trim ${primer_bed_collapsed} --remove-incorrect-pairs --report ${prefix}.alignreport.txt --verbose < ${infile} 2> ${prefix}.alignreport.er | \
samtools sort -T SAMPLE - -o ${prefix}.sorted.bam

# Convert to BAM, sort and index:
infile=${prefix}.bam
outfile=${prefix}.sorted.bam
samtools index ${outfile}

# Quality control:
infile=${prefix}.sorted.bam
outfile=${prefix}.sorted.flagstat
samtools flagstat ${infile} > ${outfile}

infile=${prefix}.sorted.bam
outfile=${prefix}.sorted.idxstats
samtools idxstats ${infile} > ${outfile}

conda deactivate
```     

<a name="8"></a>
#### 8. Coverage analysis:
```Bash
infile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed.sorted.bam
mosdepth_prefix=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed.sorted
samtools_outfile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed.sorted.depth

# In case you installed Mosdepth through conda:
conda activate mosdepth

# Coverage analysis with mosdepth:
mosdepth \
    --threads ${threads} \
    ${mosdepth_prefix} \
    ${infile}

# Coverage analysis with mosdepth per amplicon:
mosdepth \
    --threads ${threads} \
    --by ${primer_bed_collapsed} \
    --thresholds 1,10,20,50 \
    ${mosdepth_prefix}.all_amplicons \
    ${infile}

# Coverage in G and F genes:
mosdepth \
    --threads 4 \
    --by ${genes_gf_bed} \
    --thresholds 1,10,20,50 \
    ${mosdepth_prefix}.genes_gf \
    ${infile}

conda deactivate

# Coverage analysis with samtools:
samtools depth -a ${infile} > ${samtools_outfile}

# BAMdash plot:
# In case you installed bamdash through conda:
conda activate bamdash

ref_id=$( samtools view -H ${infile} | grep "@SQ" | awk '{ print $2 }' | cut -d":" -f2 )

bamdash -b ${infile} -r ${ref_id} -t ${primer_bed_collapsed} ${genes_gf_bed}

conda deactivate
```

<a name="9"></a>
#### 9. Calling with Artic and Medaka:
```Bash
# In case you installed Artic pipeline through conda:
conda activate artic

infile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed.sorted.bam

medaka_model="r1041_e82_400bps_sup_v4.3.0"

# 1) Variant calling on each read group (RG = pool):
for pool in $( awk '{ print $5 }' ${primer_bed_RSV_A_and_B} | sort | uniq ); do
    hdf=SAMPLE.aligned-to-${ref_type}.medaka_consensus.pool${pool}.hdf
    vcf=SAMPLE.aligned-to-${ref_type}.medaka_consensus.pool${pool}.vcf
    ann=SAMPLE.aligned-to-${ref_type}.medaka_consensus.pool${pool}.annotate.vcf

    medaka consensus --model ${medaka_model} --threads ${threads} --chunk_len 800 --chunk_ovlp 400 --RG ${pool} ${infile} ${hdf}

    medaka variant ${selected_ref} ${hdf} ${vcf}

    medaka tools annotate --pad 25 --RG ${pool} ${vcf} ${selected_ref} ${infile} ${ann}
done

# 2) Merge called variants for each read group:
prefix=SAMPLE.aligned-to-${ref_type}.medaka_consensus
vcf_pool1=SAMPLE.aligned-to-${ref_type}.medaka_consensus.pool1.annotate.vcf
vcf_pool2=SAMPLE.aligned-to-${ref_type}.medaka_consensus.pool2.annotate.vcf

artic_vcf_merge ${prefix} ${primer_bed} 1:${vcf_pool1} 2:${vcf_pool2} 2> ${prefix}.primersitereport.txt

# 3) Check and filter VCFs:
prefix=SAMPLE.aligned-to-${ref_type}.medaka_consensus.merged

bgzip -f ${prefix}.vcf
tabix -p vcf ${prefix}.vcf.gz
artic-tools check_vcf --dropPrimerVars --dropOverlapFails --vcfOut ${prefix}.filtered.vcf ${prefix}.vcf.gz ${primer_bed} 2> ${prefix}.vcfreport.txt

# 4) Filter the variants to produce PASS and FAIL lists, then index them:
artic_vcf_filter --medaka ${prefix}.filtered.vcf ${prefix}.filtered.pass.vcf ${prefix}.filtered.fail.vcf
bgzip -f ${prefix}.filtered.pass.vcf
tabix -p vcf ${prefix}.filtered.pass.vcf.gz

conda deactivate
```

<a name="10"></a>
#### 10. Create consensus with Artic and BCFtools:
```Bash
# In case you installed Artic pipeline through conda:
conda activate artic

infile=SAMPLE.aligned-to-${ref_type}.sorted.mapped.trimmed.sorted.bam
vcf_pass=SAMPLE.aligned-to-${ref_type}.medaka_consensus.merged.filtered.pass.vcf.gz
vcf_fail=SAMPLE.aligned-to-${ref_type}.medaka_consensus.merged.filtered.fail.vcf

# 1) Get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
artic_make_depth_mask --store-rg-depths ${selected_ref} ${infile} coverage_mask.txt
artic_mask ${selected_ref} ${coverage_mask} ${vcf_fail} preconsensus.fasta

# 2) Generate the consensus sequence:
bcftools consensus -f preconsensus.fasta ${vcf_pass} -m ${coverage_mask} -o consensus.fasta

# 3) Apply the header to the consensus sequence and run alignment agains the reference sequence:
fasta_header="./ARTIC/medaka"
muscle_infile=SAMPLE.muscle.in.fasta
muscle_outfile=SAMPLE.muscle.out.fasta

artic_fasta_header ${consensus} ${fasta_header}
cat ${consensus} ${selected_ref} > ${muscle_infile}
muscle -in ${muscle_infile} -out ${muscle_outfile}

conda deactivate
```

<a name="11"></a>
#### 11. Lineage classification with NextClade:
```Bash
strain=$( cat sample.strain.txt )

# Create a multisample FASTA with more sequences separated by reference:
sequences_dir="/path/to/more/sequences"
outfile="multifasta.${type}.fa"

cat ${sequences_dir}/*.aligned-to-${type}.ivar_consensus.fa ${infile} > ${outfile}

# You can upload this multifasta file to Nextclade (https://clades.nextstrain.org/) and select the
# correct pathogen to run clade assignment, mutation calling, and sequence quality checks.
```
  
<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
