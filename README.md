 <a name="RSV"></a> 
 
<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/RSV" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/logos_GH.png" width="auto" /> 
      </a>
</p>

# Respiratory Syncytial Virus genomic surveillance in the Canary Islands

The COVID-19 pandemic has shown the impact of genomic surveillance of emergent and re-emergent pathogens based on Next Generation Sequencing (NGS), as has been recognized by the World Health Organization [<a href="#References">1,2</a>]. Guiding the Public Health response has been accelerated thanks to the generalization of the NGS, allowing the identification and monitoring of emerging SARS-CoV-2 variants in a routine basis across the World.

Here we present a public repository of **Respiratory Syncytial Virus** related resources maintained by ITER-FIISC-HUNSC-ULL task force.

This is the result of a continuous collaborative effort of the following Institutions and Laboratories:
<ul>
 <li><a href="https://www3.gobiernodecanarias.org/sanidad/scs/organica.jsp?idCarpeta=10b3ea46-541b-11de-9665-998e1388f7ed">Servicio de Microbiología, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://fciisc.org/">Fundación Canaria Instituto de Investigación Sanitaria de Canarias</a> at the Research Unit, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://portalciencia.ull.es/grupos/6361/detalle">Laboratorio de Inmunología Celular y Viral</a>, Unidad de Farmacología, Facultad de Medicina, Universidad de La Laguna, 38200 San Cristóbal de La Laguna, Spain.</li>
 <li><a href="https://www.iter.es/areas/genomics/?lang=en">Genomics Division, Instituto Tecnológico y de Energías Renovables</a>, 38600 Santa Cruz de Tenerife, Spain.</li>
</ul>

<hr>
<!-- ------------------ SECTION ------------------ -->

# Table of contents #
<ul>
  <li><a href="#Virological-posts">A draft of the first Respiratory Syncytial Virus genomes from the Canary Islands, Spain, 2022-2023</a></li>
  <li><a href="#Protocols">Protocols for library preparation and sequencing of RSV viral genomes</a></li>
      <ul>
          <li><a href="#Illumina-protocol">Illumina-based protocol</a></li>
          <li><a href="#ONT-protocol">Oxford Nanopore Technologies-based protocol</a></li>
          <li><a href="#PCR-primers">PCR primers</a></li>
      </ul>
 <li><a href="#Bioinformatic pipelines">Bioinformatic pipelines</a></li>
      <ul>
           <li><a href="#Code-Illumina">Code for Illumina short-reads processing</a></li>
           <!-- <li><a href="#Code-ONT">Code for Nanopore long-reads processing and hybrid <i>de novo</i> assemby</a></li> -->
           <li><a href="#List-of-software">List of bioinformatic software used in our pipelines</a></li>
           <li><a href="#Useful-Files">Useful files for the pipelines (FASTA references, BED files, etc.)</a></li>
      </ul>
  <li><a href="#Sequences-and-Classification-Results">Sequences and Classification Results</a></li>
  <li><a href="#Other-repos">Other useful repositories with resources to study RSV</a></li>
  <li><a href="#References">References</a></li>
  <li><a href="#Acknowledgements">Acknowledgements</a></li>
  <li><a href="#License and Attribution">License and Attribution</a></li>
  <li><a href="#Participating">Participating</a></li>
  <li><a href="#How-to-cite">How to cite this work</a></li>
  <li><a href="#Update logs">Update logs</a></li>
</ul>

<hr>
<!-- ------------------ SECTION 1 ------------------ -->

<a name="Virological-posts"></a>
## A draft of the first Respiratory Syncytial Virus genomes from the Canary Islands, Spain, 2022-2023

The Respiratory Syncytial Virus (RSV) is the leading cause of acute lower respiratory infections in children [<a href="#References">3</a>]. The first genome sequences of RSV A and B subtypes described by us are phylogenetically related to the multiple virus genomes deposited in <a href="https://gisaid.org/">GISAID</a> that correspond to circulating RSV worldwide, as shown in Figures 1 and 2.

<p align="center">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/figures/Tree_hRSV-A_CanaryIslands_Spain_EN.png" width="75%" />
  </a>
</p>

**Figure 1**. A phylogenetic tree depicting the position of the genome draft of RSV-A sampled in February and August 2023, from patients in the Canary Islands along with NCBI GenBank publicly available sequences as computed by <a href="https://clades.nextstrain.org/">Nextstrain</a> using the the 'hRSV/A/England/397/2017' (EPI_ISL_412866) reference.

<p align="center">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/figures/Tree_hRSV-B_CanaryIslands_Spain_EN.png" width="75%" />
  </a>
</p>

**Figure 2**. A phylogenetic tree depicting the position of the genome draft of RSV-B sampled in the period November 2022-June 2023, from patients in the Canary Islands along with NCBI GenBank publicly available sequences as computed by <a href="https://clades.nextstrain.org/">Nextstrain</a> using the the 'hRSV/B/Australia/VIC-RCH056/2019' (EPI_ISL_1653999) reference.

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 2 ------------------ -->

<a name="Protocols"></a>
## Protocols for library preparation and sequencing of RSV genomes

<a name="Illumina-protocol"></a>
**Illumina-based protocol**

One of the sequencing strategies followed for SARS-CoV-2 surveillance is the use of amplicons derived from primer pools designed by the ARTIC community following a tiling approach [<a href="#References">4,5,6</a>] (see the <a href="#PCR-primers">PCR-primers</a> section).

Davina-Nunez et al. [<a href="#References">7</a>] have adapted the Illumina COVIDSeq™ Assay (RUO) kit to get the sequence of RSV A and B genomes. Their protocol uses a combination of two nested primer sets covering both RSV-A and RSV-B in the same reaction [<a href="#References">7,8,9,10</a>], followed by the Illumina COVIDSeq™ Assay protocol with minor modifications, and taking advantage of the same reagents included in the kit: <a href="https://www.protocols.io/view/whole-genome-amplification-of-respiratory-syncytia-cvpbw5in.html">Whole-Genome Amplification of Respiratory Syncytial Virus (RSV) using Illumina CovidSeq reagents for Next-Generation Sequencing V.1</a> at <a href="https://www.protocols.io/view/whole-genome-amplification-of-respiratory-syncytia-cvpbw5in.html">protocols.io</a>.

<hr>

<a name="ONT-protocol"></a>
**Oxford Nanopore Technologies-based protocol**

_Work in progress_. Come back by the end of December 2023 to find new stuff in this section.

<hr>
          
<a name="PCR-primers"></a>
**PCR primers**

PCR primers from Davina-Nunez et al. [<a href="#References">7</a>] (and adapted from Wang et al. [<a href="#References">8</a>]) have been used. These primers allows the splitting into two pools of non-consecutive amplicons (odd-numbered amplicon primers in pool-1; even-numbered amplicon primers in pool-2). <a href="https://github.com/genomicsITER/RSV/blob/main/primer_schemes/hRSV-A-B_Primers.tsv">Download this file</a> in tab-separated format.

<p align="center">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/figures/hRSV-A-B_Primers.png" width="500" />
  </a>
</p>

The coverage (x) for hRSV samples classified as RSV-A and RSV-B are shown in figure 3.

<p align="center">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/figures/Coverage_plots_one-sample.png" width="75%" />
  </a>
</p>

**Figure 3**. Coverage for hRSV A and B obtained with an Illumina MiSeq2 sequencer and MiSeq Reagent Kit v2 Nano (2x150 bp) collected from nasopharyngeal swabs. 

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 3 ------------------ -->

<a name="Bioinformatic pipelines"></a>
## Bioinformatic pipelines

The following diagram (Figure 3) represents a full pipeline used to derive the consensus FASTA sequence of RSV viruses using short-read Illumina sequencing.

In the upper part of the diagram, there is a typical pipeline to process short reads, from the basecalling to the final consensus FASTA sequence, and downstream analysis such as the phylogenetic inference.

In the lower part of the diagram, a typical pipeline to process long reads is shown. In addition, it shows how to perform a hybrid de novo genome assembly combining short and long reads.

Several consensus RSV A and B sequences derived from the pipeline based on mapping of Illumina short reads against a RSV (A or B) reference genome have been obtained so far. They have been deposited in GISAID EpiRSV (see <a href="#Sequences-and-Classification-Results">'Sequences'</a> section below).

<p align="center">
  <a href="https://github.com/genomicsITER/RSV" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/RSV/blob/main/figures/RSV_pipeline.png" width="auto" /> 
  </a>
</p>

**Figure 3**. Full bioinformatic pipeline to obtain the RSV sequences and to infer phylogenetic relationships with other RSV genomes available obtained from public repositories as provided by Nextstrain.

<hr>

<a name="Code-Illumina"></a>
**Code for Illumina short-reads processing**

<a href="https://github.com/genomicsITER/RSV/blob/main/codes/code_Illumina_pipeline.md"><img src="https://github.com/genomicsITER/RSV/blob/main/images/Code-Window-icon.png" width="32px" /></a>  See a detailed pipeline with examples of command usage for [Illumina short reads](https://github.com/genomicsITER/RSV/blob/main/codes/code_Illumina_pipeline.md).

<!-- 

<hr>

<a name="Code-ONT"></a>
**Code for Nanopore long-reads processing and hybrid *de novo* genome assemby**

<a href="https://github.com/genomicsITER/monkeypox/blob/main/codes/code_ONT-and-HybridAssembly_pipeline.md"><img src="https://github.com/genomicsITER/monkeypox/blob/main/images/Code-Window-icon.png" width="32px" /></a>  See a detailed pipeline with examples of command usage for [Oxford Nanopore Technology long-reads](https://github.com/genomicsITER/monkeypox/blob/main/codes/code_ONT-and-HybridAssembly_pipeline.md).

-->

<hr>

<a name="List-of-software"></a>
**List of bioinformatic software used in our pipelines**

<details>
<summary>Bioinformatic software (click to display):</summary>
<ul>
  <li>Conda manual for installation of numerous open-source tools used in these pipelines:<a href="https://docs.conda.io/en/latest/">Conda documentation</a></li>
  <li>Programming environment of general purpose: <a href="https://www.r-project.org/">R v.4.1.3</a></li>
  <li>Quality Control of Illumina reads: <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC v0.11.9</a></li>
  <li>Adapter trimming: <a href="https://github.com/OpenGene/fastp">fastp v0.23.2</a></li>
  <li>Remove Human mapping-reads from your FASTQ files: <a href="https://ccb.jhu.edu/software/kraken2/">Kraken2 v.2.1.2</a>. If you have issues when downloading the database indexes, try this <a href="https://benlangmead.github.io/aws-indexes/k2" >alternative site</a> from <a href="https://github.com/BenLangmead" >BenLangmead</a>.</li>
  <li>Visualization of Kraken2 reports: <a href="https://ccb.jhu.edu/software/pavian/">Pavian v.1.0</a></li>
  <li>Assemble of Illumina short-reads: <a href="https://github.com/rrwick/Unicycler">Unicycler v0.5.0</a></li>
  <li>Benchmarking and quality control of assemblies: <a href="http://quast.sourceforge.net/">QUAST v.5.0.2</a></li>
  <li>CLI tool to search in nucleotide databases using a nucleotide query: <a href="https://blast.ncbi.nlm.nih.gov/doc/blast-help/">BLAST+ v.2.12.0</a></li>
  <li>Mapping of short-reads: <a href="https://github.com/lh3/bwa/">BWA v.0.7.17-r1188</a></li>
  <li>Get mapping statistics, manipulate BAM files, and generate mpileups for FASTA consensus: <a href="https://github.com/samtools/samtools">SAMtools v.1.6</a></li>
  <li>Compute the depth of coverage and other statistics: <a href="https://github.com/brentp/mosdepth/">Mosdepth v.0.3.3</a></li>
  <li>Perform the variant calling and consensus: <a href="https://github.com/andersen-lab/ivar/">iVar v.1.3.1</a></li>
  <li>Multiple Sample Alignment: <a href="https://mafft.cbrc.jp/alignment/server/">MAFFT v.7.505</a></li>
  <li>Phylogenomic inference and tree computing: <a href="http://www.iqtree.org/">IQ-TREE v.2.2.0.3</a></li>
  <li>Framework for analyses and visualization of pathogen genome data (Nextstrain-Influenza in this case): <a href="https://nextstrain.org/influenza/">Nextstrain</a></li>
  <li>Visualization of phylogenetic trees: <a href="http://tree.bio.ed.ac.uk/software/figtree/">Figtree</a></li>
  <li>Visualization of phylogenetic trees: <a href="https://bioconductor.org/packages/release/bioc/html/ggtree.html/">ggtree 3.15</a></li>
  <li>Annotation of genomes: <a href="https://pcingola.github.io/SnpEff/">SnpEff v.5.1d</a></li>
</ul>
</details>

<hr>

<a name="Useful-Files"></a>
## Useful files for the pipelines

### Reference sequences

|RSV-A|RSV-B|
|:---:|:---:|
|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/12814/">hRSV/A/England/397/2017 EPI_ISL_412866 2017-01-01<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png"  style="width:32px;" title="Download FASTA" alt="Download FASTA" /></a>|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/12814/">hRSV/B/Australia/VIC-RCH056/2019 EPI_ISL_1653999 2019-03-04<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png"  style="width:32px;" title="Download FASTA" alt="Download FASTA" /></a>|

<br>

### BED and FASTA files

Primer schemes files (BED and FASTA) are required in the trimming step of PCR-primers.

Example of a <code>BED file</code> for hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01 using the <a href="#PCR-primers">primer-scheme</a>:

```
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	1	22	RSV_A_1_pool1_LEFT	60	-
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	1752	1779	RSV_A_1_pool1_RIGHT	60	+
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	1556	1576	RSV_A_2_pool2_LEFT	60	-
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	3377	3400	RSV_A_2_pool2_RIGHT	60	+
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	2897	2919	RSV_A_3_pool1_LEFT	60	-
hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01	4802	4826	RSV_A_3_pool1_RIGHT	60	+
...
```

Example of the <code>FASTA file</code> paired with the corresponding<code>BED file</code> for hRSV/A/England/397/2017|EPI_ISL_412866|2017-01-01 using the <a href="#PCR-primers">primer-scheme</a>:

```
>RSV_A_1_pool1_LEFT
ACGSGAAAAAATGCGTACAAC
>RSV_A_1_pool1_RIGHT
GAAGATTGTGCTATACCAAAATGAACA
>RSV_A_2_pool2_LEFT
ACAGGCATGACTCTCCTGAT
>RSV_A_2_pool2_RIGHT
TTGGGTGTGGATATTTGTTTCAC
>RSV_A_3_pool1_LEFT
GCYATGGCAAGACTYAGGAATG
>RSV_A_3_pool1_RIGHT
GTTTGCYGAGGCTATGAATATGAT
...
```

<br>

Please, download paired BED and FASTA files for each RSV subtype.

<div align="center">
 
|RSV Subtype|BED|FASTA|
|:---:|:---:|:---:|
|RSV-A|<a href="https://github.com/genomicsITER/RSV/blob/main/primer_schemes/A/RSV_A_primers.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file"/></a>|<a href="https://github.com/genomicsITER/RSV/blob/main/primer_schemes/A/RSV_A_primers.fasta"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download FASTA file" alt="Download FASTA file"/></a>|
|RSV-B|<a href="https://github.com/genomicsITER/RSV/blob/main/primer_schemes/B/RSV_B_primers.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file"/></a>|<a href="https://github.com/genomicsITER/RSV/blob/main/primer_schemes/B/RSV_B_primers.fasta"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download FASTA file" alt="Download FASTA file"/></a>|

</div>

<br>
 
<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 3 ------------------ -->

<a name="Sequences-and-Classification-Results"></a>
## Sequences and Classification Results ##

## Deposited sequences ##

Sequences are being deposited at <a href="https://gisaid.org/">GISAID</a>. You may search in GISAID by using the accession codes provided or proceed directly downloading our RSV sequences using the links provided below.

**Sequences of RSV-A**
  <ul>
    <li>Accesion 1: <a href="https://github.com/genomicsITER/RSV/blob/main/sequences/RSV-A/EPI_ISL_18321533" title="Download FASTA" >EPI_ISL_18321533</a></li>
    <li>Accesion 2: <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/RSV-A/EPI_ISL_18321558" title="Download FASTA" >EPI_ISL_18321558</a></li>
  </ul>

**Sequences of RSV-B**
  <ul>
    <li>Accesion 3: <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/RSV-B/EPI_ISL_18323795" title="Download FASTA" >EPI_ISL_18323795</li>
    <li>Accesion 4: <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/RSV-B/EPI_ISL_18323796" title="Download FASTA" >EPI_ISL_18323796</li>
    <li>Accesion 5: <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/RSV-B/EPI_ISL_18323797" title="Download FASTA" >EPI_ISL_18323797</li>
    <li>Accesion 6: <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/RSV-B/EPI_ISL_18323798" title="Download FASTA" >EPI_ISL_18323798</li>
  </ul>

  (*) NOTE: Some sequence/s may be incomplete.

  <br>

## Classification Results ##

|GISAID accession|Isolate name|RSV Subtype|Clade(G_clade)|Location|
|:---|:---:|:---:|:---:|:---:|
|EPI_ISL_18321533|A/Spain/CN-HUNSC_ITER|A|A23(G_clade:GA2.3.5)|Europe/Spain/Canary Islands|
|EPI_ISL_18321558|A/Spain/CN-HUNSC_ITER|A|A23(G_clade:GA2.3.5)|Europe/Spain/Canary Islands|
|EPI_ISL_18323795|A/Spain/CN-HUNSC_ITER|B|B6(G_clade:GB5.0.5a)|Europe/Spain/Canary Islands|
|EPI_ISL_18323796|A/Spain/CN-HUNSC_ITER|B|B6(G_clade:GB5.0.5a)|Europe/Spain/Canary Islands|
|EPI_ISL_18323797|A/Spain/CN-HUNSC_ITER|B|B6(G_clade:GB5.0.5a)|Europe/Spain/Canary Islands|
|EPI_ISL_18323798|A/Spain/CN-HUNSC_ITER|B|B6(G_clade:GB5.0.5a)|Europe/Spain/Canary Islands|

(*) NOTE: other metadata are available for these samples in <a href="https://gisaid.org/">GISAID</a> and from the authors upon a reasonable request.

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>


<!-- ------------------ SECTION 6 ------------------ -->

<a name="Other-repos"></a>
## Other useful repositories with resources to study RSV

<details>
<summary>Kudos to all research teams behind the scenes in all these repositories and web platforms (click to display):</summary>

<ul>
  <li><a href="">Repo</a></li>
  <li><a href="">Repo</a></li>
  <li><a href="">Repo</a></li>
  <li><a href="">Repo</a></li>
</ul>

</details>
 
<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 7 ------------------ -->

<a name="References"></a>
## References ##

<ol>
<li><a href="https://www.who.int/publications/i/item/9789240018440">Genomic sequencing of SARS-CoV-2. A guide to implementation for maximum impact on public health</a>, WHO, January 8, 2021.</li>
<li><a href="https://apps.who.int/iris/handle/10665/3"> Report “Global genomic surveillance strategy for pathogens with pandemic and epidemic potential, 2022-2032”</a>. Ginebra, WHO, 2022.</li>
<li>Heppe-Montero M., Walter S., Hernández-Barrera V. et al. Burden of respiratory syncytial virus-associated lower respiratory infections in children in Spain from 2012 to 2018. BMC Infect Dis 22, 315 (2022). <a href="https://doi.org/10.1186/s12879-022-07261-1">doi:10.1186/s12879-022-07261-1</a></li>
<li>Gohl DM, Garbe J, Grady P, et al. A rapid, cost-effective tailed amplicon method for sequencing SARS-CoV-2. BMC Genomics. 2020;21(1):863. Published 2020 Dec 4. <a href="https://doi.org/10.1186/s12864-020-07283-6">doi:10.1186/s12864-020-07283-6</a>.</li>
<li>Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M. Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. PLoS One. 2020;15(9):e0239403. Published 2020 Sep 18. <a href="https://doi.org/10.1371/journal.pone.0239403">doi:10.1371/journal.pone.0239403</a>.</li> 
<li>Koskela von Sydow A, Lindqvist CM, Asghar N, et al. Comparison of SARS-CoV-2 whole genome sequencing using tiled amplicon enrichment and bait hybridization. Sci Rep. 2023;13(1):6461. Published 2023 Apr 20. <a href="https://doi.org/10.1038%2Fs41598-023-33168-1">doi:10.1038/s41598-023-33168-1</a>.</li>
<li>Davina-Nunez C, Perez-Castro S, Godoy-Diz M, Regueiro-Garcia B 2023. Whole-Genome Amplification of Respiratory Syncytial Virus (RSV) using Illumina CovidSeq reagents for Next-Generation Sequencing. <a href="https://dx.doi.org/10.17504/protocols.io.eq2lyjzbrlx9/v1">protocols.io</a>.</li>
<li>Wang L, Ng TFF, Castro CJ, et al. Next-generation sequencing of human respiratory syncytial virus subgroups A and B genomes. J Virol Methods. 2022;299:114335. <a href="https://doi.org/10.1016/j.jviromet.2021.114335">doi:10.1016/j.jviromet.2021.114335</a>.</li>
<li>Goya S., Rojo G.L., Nabaes Jordar M.S., Valinotto L.E., Mistchenko A.S., Viegas M. Whole genome sequencing of respiratory syncytial (RSV) virus from clinical samples with low viral load. <a href="https://protocols.io/view/whole-genome-sequencing-of-respiratory-syncytial-r-bmhak32e">protocols.io</a>.</li>
<li>Lin Y., Koble J., Prashar P., Pottekat A., Middle C., Kuersten S., Oberholzer M., Brazas R., Whitlock D., Schlaberg R., Schroth G.P. A sequencing and subtyping protocol for Influenza A and B viruses using Illumina® COVIDSeq™ Assay Kit. <a href="https://protocols.io/view/a-sequencing-and-subtyping-protocol-for-influenza-crv3v68n">protocols.io</a>.</li>
</ol>
  
<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 8 ------------------ -->

<a name="Acknowledgements"></a>
## Acknowledgements ##

This study has been funded by Cabildo Insular de Tenerife (CGIEU0000219140 and "_Apuestas científicas del ITER para colaborar en la lucha contra la COVID-19_"); by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, epidemiological surveillance based on massive sequencing, Personalized Medicine and Biotechnology (OA17/008 and OA23/043); and by the agreement between Consejería de Educación, Universidades, Cultura y Deportes del Gobierno de Canarias y Cabildo Insular de Tenerife, 2022-2025 (AC0000014697).

This study is also an activity within the project Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network (101113109 - RELECOV 2.0) of the EU4Health Programme (EU4H) by the European Health and Digital Executive Agency (HaDEA), under the coordination of Instituto de Salud Carlos III (ISCIII).

We acknowledge the researchers and their institutions who released RSV sequences through NCBI GenBank, GISAID, and ENA that are being used in our studies. 

We also thank the authors, the laboratories that originated and submitted the genetic sequences and the metadata for sharing their work, as shown on Nextstrain, and:
<ul>
  <li>Hadfield <i>et al</i>, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics (2018).</li>
  <li>Sagulenko <i>et al</i>, TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017).</li>
</ul>

<!-- We would like to acknowledge the contributions of several researchers and laboratories who share their preliminary results through the [Virological](https://virological.org/) website. -->

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 9 ------------------ -->

<a name="License and Attribution"></a>
## License and Attribution ##

This repository and data exports are released under the CC BY 4.0 license. Please acknowledge the authors, the originating and submitting laboratories for the genetic sequences and metadata, and the open source software used in this work (third-party copyrights and licenses may apply).

Please cite this repository as: _"Influenza repository of the Reference Laboratory for Epidemiological Surveillance of Pathogens in the Canary Islands (accessed on YYYY-MM-DD)"_. And do not forget to <a href="#How-to-cite">cite the paper</a> (see the section "How to cite" below) when it becomes available. 

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 10 ------------------ -->

<a name="Participating"></a>
## Participating ##

> Want to share your relevant links? Place a Direct Message to @labcflores, @adrmunozb or @resocios on Twitter (see below).

 <p align="left">
  <a href="https://www.iter.es/areas/area-genomica/" title="Contact us at the Genomics Division of the Institute of Technology and Renewable Energy (ITER), Tenerife, Canary Islands, Spain">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/ITER_logo.png" width="30%" /> 
  </a>
</p>

By AMB <a href="https://twitter.com/adrmunozb" title="Follow to @resocios on Twitter" >@adrmunozb <img src="https://github.com/genomicsITER/influenza/blob/main/images/Twitter_Social_Icon_Circle_Color.png" width="32px" /></a> and JMLS <a href="https://twitter.com/resocios" title="Follow to @resocios on Twitter" >@resocios <img src="https://github.com/genomicsITER/influenza/blob/main/images/Twitter_Social_Icon_Circle_Color.png" width="32px" /></a>

Follow us on Twitter <a href="https://twitter.com/labcflores" title="Follow to @labcflores on Twitter" > @labcflores<img src="https://github.com/genomicsITER/influenza/blob/main/images/Twitter_Social_Icon_Circle_Color.png" width="32px" /></a>

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 11 ------------------ -->

<a name="How-to-cite"></a>
## How to cite this work ##

> This work has not been publised yet. See 'License and Attribution' section to cite this repository.

> To use the deposited sequences at GISAID, please, acknowledge this work as recommended by GISAID. Find the 'GISAID acknowledge tables' <a href="https://github.com/genomicsITER/RSV/tree/main/sequences/acknowledgements">here</a>.

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 12 ------------------ -->

<a name="Update logs"></a>
## Update logs ##

> September 27, 2023. Created the private version of this repository.

<p align="right">
  <a href="#RSV" title="Up">
    <img src="https://github.com/genomicsITER/RSV/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
