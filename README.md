 <a name="Influenza"></a> 
 
<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/influenza" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/logos_GH.png" width="auto" /> 
      </a>
</p>

# Influenza

The COVID-19 pandemic has shown the impact of genomic surveillance of emergent and re-emergent pathogens based on Next Generation Sequencing (NGS), as has been recognized by the World Health Organization [1,2]. Guiding the Public Health response has been accelerated thanks to the generalization of the NGS, allowing the identification and monitoring of emerging SARS-CoV-2 variants in a routine basis accross the World.

Here we present a public repository of Influenza (Inf) related resources maintained by ITER-FIISC-HUNSC- task force.

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
  <li><a href="#Virological-posts"><!-- Virological post: -->A draft of the first Influenza genomes from the Canary Islands, Spain, 2022-2023</a></li>
  <li><a href="#Protocols">Protocols for library preparation and sequencing of Influenza viral genomes</a></li>
      <ul>
          <li><a href="#Illumina-protocol">Illumina-based protocol</a></li>
          <li><a href="#ONT-protocol">ONT-based protocol</a></li>
          <li><a href="#PCR-universal-primers">PCR universal primers</a></li>
      </ul>
 <li><a href="#Bioinformatic pipelines">Bioinformatic pipelines</a></li>
      <ul>
           <li><a href="#Code-Illumina">Code for Illumina short-reads processing</a></li>
           <!-- <li><a href="#Code-ONT">Code for Nanopore long-reads processing and hybrid <i>de novo</i> assemby</a></li> -->
           <li><a href="#List-of-software">List of bioinformatic software used in our pipelines</a></li>
           <li><a href="#Useful-Files">Useful files for the pipelines (FASTA references, BED files, etc.)</a></li>
      </ul>
  <li><a href="#Sequences-and-Classification-Results">Sequences and Classification Results</a></li>
 
  <!-- <li><a href="#How-to-GenBank">How to download sequences and metadata from GenBank</a></li> -->
 
  <li><a href="#Other-repos">Other useful repositories with resources to study Influenza</a></li>
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
## <!-- Virological post: -->A draft of the first Influenza genomes from the Canary Islands, Spain, 2022-2023 ##

<!-- A technical post with the draft of the first Influenza genomes from the Canary Islands, Spain, will be shared through Virological shortly. Keep reading! [Here.](https://virological.org/) 
-->

The first genome sequences of Influenza virus A/H1N1, A/H3N2, and B (Victoria) described by us are phylogenetically related to the multiple viral genomes deposited in <a href="https://gisaid.org/">GISAID</a> that correspond to the past 2022-2023 seasonal flu wave in the Northern hemisphere, as shown in Figures 1 and 2.

<p align="center">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/figures/Tree_Influenza_A-H1N1_CanaryIslands_Spain_EN.png" width="auto" />
  </a>
</p>

**Figure 1**. A phylogenetic tree depicting the position of the genome draft of Influenza A/H1N1 isolated in the period November-December 2022, from patients from the Canary Islands along with NCBI GenBank publicly available sequences as computed by <a href="https://clades.nextstrain.org/">Nextstrain</a> using the HA gene and Influenza A H1N1pdm HA [<a href="https://www.ncbi.nlm.nih.gov/nuccore/CY121680">A/California/07/2009 (CY121680)</a>] as reference.


<p align="center">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/figures/Tree_Influenza_A-H3N2_CanaryIslands_Spain_EN.png" width="auto" />
  </a>
</p>

**Figure 2**. A phylogenetic tree depicting the position of the genome draft of Influenza A/H3N2 isolated in the period October-December 2022, from patients from the Canary Islands along with NCBI GenBank publicly available sequences as computed by <a href="https://clades.nextstrain.org/">Nextstrain</a> using the HA gene and Influenza A H3N2 HA [<a href="https://www.ncbi.nlm.nih.gov/nuccore/CY163680">A/Wisconsin/67/2005 (CY163680)</a>] as reference.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 2 ------------------ -->

<a name="Protocols"></a>
## Protocols for library preparation and sequencing of Influenza viral genomes ##

<a name="Illumina-protocol"></a>
**Illumina-based protocol**

One of the sequencing strategies followed for SARS-CoV-2 surveillance is the use of amplicons derived from primer pools designed by the ARTIC community following a tiling approach (reference). However, this approach is not suitable for influenza viruses because of their mutational burden and higher variability. Besides, it is possible to use the so called universal primers [Zhou references] taking advantage of the conserved promoter regions at the 5' and 3' ends of the influenza genome segments in order to amplify the entire genome using larger amplicons [references] (see the PCR-primers section).

Lin Y. et al. [5] have adapted the Illumina COVIDSeq™ Assay (RUO) kit to get the sequence of influenza A and B viral genomes. Their protocolo uses a combination of two primer sets (references), follows the Illumina COVIDSeq™ Assay protocol with minor modifications, and the take advantage of the same reagents included in the kit: <a href="https://www.protocols.io/view/a-sequencing-and-subtyping-protocol-for-influenza-n2bvj8mrxgk5/v1">A sequencing and subtyping protocol for Influenza A and B viruses using Illumina® COVIDSeq™ Assay Kit</a> at <a href="https://www.protocols.io/view/a-sequencing-and-subtyping-protocol-for-influenza-n2bvj8mrxgk5/v1">protocols.io</a>.

According to Lin Y. et al., this protocol provides accurate information for subtyping, lineage tracing, and antiviral resistance detection of influenza viruses.

<hr>

<a name="ONT-protocol"></a>
**ONT-based protocol**

Working in progress. Come back by the end of December 2023 to find new stuff in this section.

<hr>
          
<a name="PCR-universal-primers"></a>
**PCR universal primers**

PCR Universal Primers from Zhout et al. (<a href="https://link.springer.com/protocol/10.1007/978-1-61779-621-0_11">2012</a>, <a href="https://doi.org/10.1128/jcm.03265-13">2014</a>).


<p align="center">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/figures/Influenza_Zhou-2012-2014_Universal_Primers.png" width="75%" />
  </a>
</p>


<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 3 ------------------ -->

<a name="Bioinformatic pipelines"></a>
## Bioinformatic pipelines ##

The following diagram (Figure 3) represents a full pipeline used to derive the consensus FASTA sequence of influenza viruses using short-reads Illumina sequencing technology. A complementary pipeline is being designed and will be shown here by October 2023.

In the upper part of the diagram, there is a typical pipeline to process short-reads, from the basecalling to the final consensus FASTA sequence, and downstream analysis such as the phylogenetic inference.

In the lower part of the diagram, it is shown a typical pipeline to process long-reads. In addition, it shows how to perform a hybrid *de novo* assembly combining short- and long-reads.

Several consensus Influenza A/H1N1 and A/H3N2 sequences derived from the pipeline based on mapping of Illumina short-reads against an Influenza reference genome have been obtained so far. They have been deposited in GISAID EpiFlu (see <a href="#Sequences">'Sequences'</a> section below).

<p align="center">
  <a href="https://github.com/genomicsITER/influenza" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/influenza/blob/main/figures/Influenza_pipeline.png" width="auto" /> 
  </a>
</p>

**Figure 3**. Full bioinformatic pipeline to obtain the Influenza sequences and to infer phylogenetic relationships with other Influenza viral genomes available obtained from public repositories as provided by Nextstrain.

<hr>

<a name="Code-Illumina"></a>
**Code for Illumina short-reads processing**

<a href="https://github.com/genomicsITER/influenza/blob/main/codes/code_Illumina_pipeline.md"><img src="https://github.com/genomicsITER/influenza/blob/main/images/Code-Window-icon.png" width="32px" /></a>  See a detailed pipeline with examples of command usage for [Illumina short-reads](https://github.com/genomicsITER/influenza/blob/main/codes/code_Illumina_pipeline.md).

<!-- 

<hr>

<a name="Code-ONT"></a>
**Code for Nanopore long-reads processing and hybrid *de novo* assemby**

<a href="https://github.com/genomicsITER/monkeypox/blob/main/codes/code_ONT-and-HybridAssembly_pipeline.md"><img src="https://github.com/genomicsITER/monkeypox/blob/main/images/Code-Window-icon.png" width="32px" /></a>  See a detailed pipeline with examples of command usage for [Oxford Nanopore Technology long-reads](https://github.com/genomicsITER/monkeypox/blob/main/codes/code_ONT-and-HybridAssembly_pipeline.md).

-->

<hr>

<a name="List-of-software"></a>
**List of bioinformatic software used in our pipelines**
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

<hr>

<a name="Useful-Files"></a>
**Useful files for the pipelines**

**Reference sequences**

|A/H1N1|A/H3N2|B Victoria|B Yamagata|
|:---:|:---:|:---:|:---:|
|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/641809/">Influenza A virus (A/California/07/2009(H1N1))<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download FASTA" alt="Download FASTA" /></a>|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/393902/">Influenza A virus (A/Wisconsin/67/2005(H3N2))<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png"  style="width:32px;" title="Download FASTA" alt="Download FASTA" /></a>|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/604436/">Influenza B virus Victoria (B/Brisbane/60/2008)<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download FASTA" alt="Download FASTA"/></a>|<a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1089607/">Influenza B virus Yamagata (B/Wisconsin/01/2010)<br><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download FASTA" alt="Download FASTA" /></a>|

<br>

**BED files**

Primer schemes in BED format are required in the trimming step of PCR-primers.

Example of a BED file for segment 1 (FJ984387.1) of Influenza A virus (A/California/07/2009(H1N1)) using the <a href="#PCR-universal-primers">primer-scheme</a>:

```
FJ984387.1	1	12	Seg1_Uni12/Inf-1_LEFT	1	+	GGGGGGAGCAAAAGCAGG
FJ984387.1	1	12	Seg1_Uni12/Inf-3_LEFT	1	+	GGGGGGAGCGAAAGCAGG
FJ984387.1	2268	2280	Seg1_Uni13/Inf-1_RIGHT	1	-	CGGGTTATTAGTAGAAACAAGG
FJ984387.1	1	14	Seg1_B-Pbs-UniF_LEFT	1	+	GGGGGGAGCAGAAGCGGAGC
FJ984387.1	2266	2280	Seg1_B-Pbs-UniF_RIGHT	1	-	CCGGGTTATTAGTAGAAACACGAGC
```

Please, download the BED files **separately** (one file per Influenza segment).

|Virus strain|Seg-1|Seg-2|Seg-3|Seg-4|Seg-5|Seg-6|Seg-7|Seg-8|
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|Inf. A virus (A/California/07/2009(H1N1))|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg1.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png"  style="width:32px;" title="Download BED file" alt="Download BED file"/></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg2.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg3.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg4.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg5.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg6.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg7.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H1N1/Influenza_A_H1N1_A_California_07_2009.seg8.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a> | 
|Inf. A virus (A/Wisconsin/67/2005(H3N2))|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg1.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg2.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg3.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg4.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg5.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg6.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg7.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/A/H3N2/Influenza_A_H3N2_A_Wisconsin_67_2005.seg8.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|
|Inf. B virus Victoria (B/Brisbane/60/2008)|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg1.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg2.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg3.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg4.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg5.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg6.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg7.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Victoria/Influenza_B_Brisbane_60_2008.seg8.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|
|Inf. B virus Yamagata (B/Wisconsin/01/2010)|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg1.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg2.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg3.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg4.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg5.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg6.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg7.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|<a href="https://github.com/genomicsITER/influenza/blob/main/primer_schemes/B/Yamagata/Influenza_B_Wisconsin_01_2005.seg8.primer_scheme.bed"><br><img src="images/Pictogrammers-Material-Text-box-search-outline.64.png" style="width:32px;" title="Download BED file" alt="Download BED file" /></a>|



<br>
 
<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 3 ------------------ -->

<a name="Sequences-and-Classification-Results"></a>
## Sequences and Classification Results ##

## Deposited sequences ##

Sequences are being deposited at <a href="https://gisaid.org/">GISAID</a>. You may search in GISAID by using the accession codes provided or directly download our Influenza sequences using the links provided below.

**Sequences of Influenza A/H1N1**
  <ul>
    <li>Accesion 1: <a href="https://github.com/genomicsITER/influenza/blob/main/sequences/A/H1N1/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150368383_2022.fasta">EPI_ISL_18128205</a></li>
    <li>Accesion 2: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/A/H1N1/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150369675_2022.fasta">EPI_ISL_18308442</a></li>
    <li>Accesion 3: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/A/H1N1/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150372413_2022.fasta">EPI_ISL_18308501</a></li>
  </ul>

**Sequences of Influenza A/H3N2**
  <ul>
    <li>Accesion 4: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/A/H3N2/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150368334_2022.fasta">EPI_ISL_18313569</li>
    <li>Accesion 5: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/A/H3N2/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150369410_2022.fasta">EPI_ISL_18313571</li>
    <li>Accesion 6: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/A/H3N2/gisaid_epiflu_A_Spain_CN-HUNSC_ITER_150369722_2022.fasta">EPI_ISL_18313572</li>
  </ul>

**Sequences of Influenza B Victoria**
  <ul>
    <li>Accesion 7: <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/B/Victoria/gisaid_epiflu_B_Spain_CN-HUNSC_ITER_150387165_2023.fasta">EPI_ISL_18313574 </li>
  </ul>

  (*) NOTE: Some segment/s sequence/s may be incomplete.

## Classification Results ##

|GISAID accession|Isolate name|Subtype|Clade|Location|
|:---|:---:|:---:|:---:|:---:|
|EPI_ISL_18128205|A/Spain/CN-HUNSC_ITER|A/H1N1|6B.1A.5a.2a|Europe/Spain/Canary Islands|
|EPI_ISL_18308442|A/Spain/CN-HUNSC_ITER|A/H1N1|6B.1A.5a.2a.1|Europe/Spain/Canary Islands|
|EPI_ISL_18308501|A/Spain/CN-HUNSC_ITER|A/H1N1|6B.1A.5a.2a.1|Europe/Spain/Canary Islands|
|EPI_ISL_18313569|A/Spain/CN-HUNSC_ITER|A/H3N2|3C.2a1b.2a.2b|Europe/Spain/Canary Islands|
|EPI_ISL_18313571|A/Spain/CN-HUNSC_ITER|A/H3N2|3C.2a1b.2a.2b|Europe/Spain/Canary Islands|
|EPI_ISL_18313572|A/Spain/CN-HUNSC_ITER|A/H3N2|3C.2a1b.2a.2b|Europe/Spain/Canary Islands|
|EPI_ISL_18313574|B/Spain/CN-HUNSC_ITER|B|V1A.3a.2|Europe/Spain/Canary Islands|

(*) NOTE: other metadata are available for these samples in <a href="https://gisaid.org/">GISAID</a> and from the authors upon a reasonable request.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>


<!-- ------------------ SECTION 4 ------------------ -->
<!--
<a name="How-to-GenBank"></a>
## How to download sequences and metadata from GenBank ##

> Manual download

<ol>
<li>Browse to <a href="https://www.ncbi.nlm.nih.gov/genbank/">GenBank</a>.</li>
<li>Select '<i>Nucleotide</i>' from the combo box.</li>
<li>Fill in the accession code of the sequence you want to download (i.e. ON782054) or just write the name of the species (i.e. Monkeypox, and then clic on a certain accession code you are interested in).</li>
<li>Click on '<i>FASTA</i>' link</li>
<li>Click on '<i>Send to</i>' on the upper right part of the screen.</li>
<li>Select the option '<i>file</i>'.</li>
<li>Select '<i>FASTA</i>' as download format.</li>
<li>Click on '<i>Generate</i>' button.</li>
</ol>

> Programmatically download

<a href="https://github.com/genomicsITER/monkeypox/blob/main/codes/code_retrieve-seqs-from-GenBank.md"><img src="https://github.com/genomicsITER/monkeypox/blob/main/images/Code-Window-icon.png" width="32px" /></a>  We provide a full Python code to retrieve all sequences larger than 190,000 bases from GenBank as example. [See the code](https://github.com/genomicsITER/monkeypox/blob/main/codes/code_retrieve-seqs-from-GenBank.md).



<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr> -->


<!-- ------------------ SECTION 6 ------------------ -->

<a name="Other-repos"></a>
## Other useful repositories with resources to study Influenza ##

> Kudos to all research teams behind the scenes in all these repositories:

<ul>
  <li><a href="https://umasangumathi.github.io/FluLINE/">FluLINE</a></li>
  <li><a href="https://github.com/hkailee/FluSeq">FluSeq</a></li>
  <li><a href="https://insaflu.readthedocs.io/en/latest/index.html">INSaFLU</a></li>
  <li><a href="https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database">NCBI Influenza Virus Database</a></li>
  <li><a href="https://clades.nextstrain.org/">Nextclade</a></li>
  <li><a href="https://nextstrain.org/influenza/">Nextstrain</a></li>
  <li><a href="https://github.com/peterk87/nf-flu">nf-flu</a></li>
  <li><a href="https://github.com/greninger-lab/revica">REVICA</a></li>
</ul>

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 7 ------------------ -->

<a name="References"></a>
## References ##

<ol>
<li><a href="https://www.who.int/publications/i/item/9789240018440">Genomic sequencing of SARS-CoV-2. A guide to implementation for maximum impact on public health</a>, WHO, January 8, 2021.</li>
<li><a href="https://apps.who.int/iris/handle/10665/3"> Report “Global genomic surveillance strategy for pathogens with pandemic and epidemic potential, 2022-2032”</a>. Ginebra, WHO, 2022.</li>
<li>Zhou B, Wentworth DE. Influenza A virus molecular virology techniques. Methods Mol Biol. 2012;865:175-192. doi:<a href="https://link.springer.com/protocol/10.1007/978-1-61779-621-0_11">10.1007/978-1-61779-621-0_11</a>.</li>
<li>Zhou B, Lin X, Wang W, et al. Universal influenza B virus genomic amplification facilitates sequencing, diagnostics, and reverse genetics. J Clin Microbiol. 2014;52(5):1330-1337. doi:<a href="https://doi.org/10.1128/jcm.03265-13">10.1128/JCM.03265-13</a>.</li>
<li>Ying Lin, Jeffrey Koble, Priyanka Prashar, Anita Pottekat, Christina Middle, Scott Kuersten, Michael Oberholzer, Robert Brazas, Darcy Whitlock, Robert Schlaberg, Gary P. Schroth. A sequencing and subtyping protocol for Influenza A and B viruses using Illumina® COVIDSeq™ Assay Kit. Protocols.io. doi:<a href="https://dx.doi.org/10.17504/protocols.io.n2bvj8mrxgk5/v1">dx.doi.org/10.17504/protocols.io.n2bvj8mrxgk5/v1
</a></li>
</ol>
  
<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 8 ------------------ -->

<a name="Acknowledgements"></a>
## Acknowledgements ##

This study has been funded by Cabildo Insular de Tenerife (CGIEU0000219140 and "_Apuestas científicas del ITER para colaborar en la lucha contra la COVID-19_"), Instituto de Salud Carlos III (FI18/00230) cofunded by European Union (ERDF) "_A way of making Europe_", by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, Personalized Medicine and Biotechnology (OA17/008); and by the agreement between Consejería de Educación, Universidades, Cultura y Deportes del Gobierno de Canarias y Cabildo Insular de Tenerife, 2022-2025.

This study is also an activity within the project _Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network (RELECOV 2.0)_ of the EU4Health Programme (EU4H), under the coordination of Instituto de Salud Carlos III (ISCIII).

We acknowledge the researchers and their institutions who released Influenza sequences through NCBI GenBank, GISAID, and ENA that are being used in our studies. 

We also thank the authors, the laboratories that originated and submitted the genetic sequences and the metadata for sharing their work, as shown on Nextstrain, and:
<ul>
  <li>Hadfield <i>et al</i>, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics (2018).</li>
  <li>Sagulenko <i>et al</i>, TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017).</li>
</ul>

<!-- We would like to acknowledge the contributions of several researchers and laboratories who share their preliminary results through the [Virological](https://virological.org/) website. -->

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 9 ------------------ -->

<a name="License and Attribution"></a>
## License and Attribution ##

This repository and data exports are released under the CC BY 4.0 license. Please acknowledge the authors, the originating and submitting laboratories for the genetic sequences and metadata, and the open source software used in this work (third-party copyrights and licenses may apply).

Please cite this repository as: _"Influenza repository of the Reference Laboratory for Epidemiological Surveillance of Pathogens in the Canary Islands (accessed on YYYY-MM-DD)"_. And do not forget to <a href="#How-to-cite">cite the paper</a> (see the section "How to cite" below) when it becomes available. 

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
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
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 11 ------------------ -->

<a name="How-to-cite"></a>
## How to cite this work ##

> This work has not been publised yet. See 'License and Attribution' section to cite this repository.

> To use the deposited sequences at GISAID, please, acknowledge this work as recommended by GISAID. Find the 'GISAID acknowledge tables' <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/acknowledgements">here</a>.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 12 ------------------ -->

<a name="Update logs"></a>
## Update logs ##

> September 27, 2023. This repository became public. Enjoy the reading! ;=)

> September 26, 2023. Updated many sections: bioinformatic pipeline, primer-schemes (required BED files for the pipelines), deposited sequences, Influenza virus A and B reference sequences, and other useful repositories with resources to study Influenza.
 
> July 26, 2023. Created the private version of this repository.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
