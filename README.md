GeM-Pro and q-GeM
================

## IMPORTANT: Scripts and Example Data will be available after manuscript acceptance (under evaluation)

**Authors:** Mariano Torres Manno; María D. Pizarro; Marcos Prunello;
Christian Magni; Lucas D. Daurelo; Martín Espariz

## I. GeM-Pro

GeM-Pro (functional profiling by orthologous repertoire inferred by
Bayesian probability and Synteny analysis) is a tool for strain
profiling based on their repertoires of orthologous pathways of
interest.

### Requirements

UNIX like OS (Linux, MacOSX)

R, Blast, and optionally RStudio, OAT\_cmd and Prodigal programs should
be preinstalled. They can be downloaded from:

R. <https://cran.r-project.org/mirrors.html>

BLAST.
[ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](https://goo.gl/GuANzX)

RStudio. <https://www.rstudio.com/products/rstudio/download/>

OAT\_cmd. <https://www.ezbiocloud.net/tools/orthoani>

Prodigal <https://github.com/hyattpd/prodigal/releases/>

### GeM-Pro Quick Tutorial

Three cases are considered in this quick tutorial.

**Case 1:** The aim is to profile strains with known (and verified\!\!)
phylogenetic group. All strain genome sequences are well annotated
(their CDS are well predicted or known).

More information can be found in the [Example
Case](./Examples/Example_Case1.md)

*Step 1.* Create gene/pathways (Query\_pathways.csv) and phylogenetic
group (Phylogenetic\_groups.csv) files. See Section I.1. (Input files)
for reference.

*Step 2*. Convert feature files to RDS
files.

| Input files                                       | Script (section I.5) | Output Files               |
| ------------------------------------------------- | -------------------- | -------------------------- |
| strains.txt <br> Strain\_feature\_table.txt files | Subject\_feature.R   | Subject\_Feature\_list.RDS |

Example:

``` bash
Rscript Subject_feature.R Sample/strains.C1.txt Sample/Features/ _feature_table.txt
```

*Step 3*. Perform BLAST
searches.

| Input files                                | Script (section I.3) | Output Files       |
| ------------------------------------------ | -------------------- | ------------------ |
| strains.txt <br> strain\_protein.faa files | Blast.R              | output BLAST files |

Example:

``` bash
Rscript Blast.R FZB42 Sample/strains.C1.txt ALL Sample/Proteomes/ Sample/Blast_out/
```

*Step 4*. Run
GeM-Pro.

| Input files                                                                                                                                 | Script (section I.7) | Output Files                                                                                           |
| ------------------------------------------------------------------------------------------------------------------------------------------- | -------------------- | ------------------------------------------------------------------------------------------------------ |
| Phylogenetic\_groups.csv <br> output BLAST files <br> Query\_pathways.csv <br> Subject\_Feature\_list.RDS <br> REF\_STR\_feature\_table.txt | GeMPro.R             | IWNO.RDS <br> GeMPro\_DB.RDS <br> synteny.list.RDS <br> pathway.profiles.RDS <br> pathway.profiles.txt |

Example:

``` bash
Rscript GeMPro.R FZB42 Sample/Phylogenetic_groups.csv Sample/Blast_out/ Sample/Query_pathways.csv Subject_Feature_list.RDS Sample/Features/FZB42_feature_table.txt Sample/GeMPro_out/
```

*Step 5*. Perform
clustering.

| Input files          | Script (section I.8) | Output Files                       |
| -------------------- | -------------------- | ---------------------------------- |
| pathway.profiles.RDS | Clustering.R         | Clustering.RDS <br> Clustering.tre |

Example:

``` bash
Rscript Clustering.R pathway.profiles.RDS 100
```

**Case 2**: The aim is to profile strains with unknown phylogenetic
group. All strain genome sequences are well annotated (their CDS are
well predicted or known).

More information can be found in the [Example
Case](./Examples/Example_Case2.md)

*Step 1*. Create gene/pathways (Query\_pathways.csv) file. See Section
1.1 (Input files) for reference.

*Step 2*. Compute ANI values between strains and reference strains. (run
once for each reference
strain.)

| Input files                                       | Script (section I.2) | Output Files                   |
| ------------------------------------------------- | -------------------- | ------------------------------ |
| strains.txt <br> Strains\_name\_genomic.fna files | run\_ANI.R           | REF\_STR\_VS\_strain.ANI files |

Example:

``` bash
Rscript run_ANI.R 168 Sample/strains.C2.txt /usr/bin/OAT/ /usr/bin/ 4 Sample/Genomes/ Sample/ANI.out/
```

*Step 3*. Circumscribe strains in Phylogenetic
groups.

| Input files                                 | Script (section I.2) | Output Files                          |
| ------------------------------------------- | -------------------- | ------------------------------------- |
| REF\_STR\_VS\_strain.ANI files <br> ref.csv | ANI\_FUN.R           | Phylogenetic\_groups.csv <br> ANI.RDS |

Example:

``` bash
Rscript ANI_FUN.R Sample/ref.csv Sample/strains.C2.txt Sample/ANI.out/
```

After this point, this case is similar to Case 1. Go directly to Step 2
of **Case 1**\!

**Case 3**: The aim is to profile strains with known (and verified\!\!)
phylogenetic group. Some strain genome sequences are not annotated
(their CDS are not predicted).

More information can be found in the [Example
Case](./Examples/Example_Case3.md)

*Step 1*. Predict CDS and putative protein sequences.

| Input files         | Software (section I.4) | Output Files        |
| ------------------- | ---------------------- | ------------------- |
| Strain\_genomic.fna | prodigal               | Strain\_protein.faa |

Example:

``` bash
prodigal -i Sample/Genomes/37MA.fna -a Sample/Proteomes/37MA_protein.faa
```

*Step 2*. Convert feature files to RDS
files.

| Input files                                                                | Script (section I.5) | Output Files               |
| -------------------------------------------------------------------------- | -------------------- | -------------------------- |
| strains.txt <br> Strain\_feature\_table.txt files <br> Strain\_protein.faa | Subject\_feature.R   | Subject\_Feature\_list.RDS |

Example:

``` bash
Rscript Subject_feature.R Sample/strains.C3.txt Sample/Features/ _feature_table.txt Sample/Proteomes/ _protein.faa
```

After this point, this case is similar to **Case 1**. Go directly to
*Step 3* of **Case 1**\!

### GeM-Pro Tutorial

GeM-Pro classifies bacterial strains according to their repertoires of
genes/pathways of interest. The driving idea is that such genes/pathways
are correlated with a phenotype/behavior of interest (pathogenicity,
probiotic capability, plant growth promotion, etc.). Hence, the
genes/pathways of interest should be well characterized in a model
strain. GeM-Pro will predict the existence of isofunctional
genes/pathways in other strains based on the “orthology conjecture”. The
result will be the creation of a dendrogram describing the strains
gene/pathway content.

#### Section I.1. Input files.

In order to be run, GeM-Pro needs information about genes/pathways of
the model strain, as well as genomic information of each strain under
analysis. Also, a file indicating each strain phylogenetic group is
required.

Model strain genes/pathways information should be provided as a table in
comma-separated value format (CSV). It is suggested to name this table
Query\_pathways.csv. Gene accession version, name and pathway should be
indicated under the following column headers: Acc\_Ver, Gen and Pathway.

`Query_pathways.csv` example table:

| Acc\_Ver   | Gen   | Pathway      |
| ---------- | ----- | ------------ |
| ABS72764.1 | srfAA | surfactin    |
| ABS72765.1 | srfAB | surfactin    |
| ABS72767.1 | srfAC | surfactin    |
| ABS72768.1 | srfAD | surfactin    |
| ABS73214.1 | acoA  | Deg\_Acetoin |
| ABS73215.1 | acoB  | Deg\_Acetoin |
| ABS73216.1 | acoC  | Deg\_Acetoin |
| ABS73217.1 | acoL  | Deg\_Acetoin |
| ABS73218.1 | acoR  | Deg\_Acetoin |
| ABS75646.1 | alsD  | Syn\_Acetoin |
| ABS75647.1 | alsS  | Syn\_Acetoin |
| ABS75648.1 | alsR  | Syn\_Acetoin |

For each strain, their encoded protein sequences and feature
information, which indicates all CDS start and end coordinates in the
genome, are needed.

Hence, two input files are required for each strain to be analysed:

  - A fasta file with all protein sequences. It should be named
    name\_of\_strain\_protein.faa (i.e. `FZB42_protein.faa`)

  - A file describing the feature annotated in the genome to be
    analysed. It should be named name\_of\_strain\_feature\_table.txt
    (i.e. `FZB42_feature_table.txt`)

Bacteria should be assigned to a phylogenetic group (usually a species
or highly phylogenetically related strains). If the phylogenetic group
of each strain to be analysed is known, this information should be
provided in a table in comma-separated values format (CSV). It is
suggested to name this table Phylogenetic\_groups.csv. The name of the
strains to be analysed and their group name should be indicated under
the following column headers: Strains and Groups.

`Phylogenetic_groups.csv` file content example:

| Strains           | Groups        |
| ----------------- | ------------- |
| 41KF2b            | B.altitudinis |
| KACC\_16563       | B.altitudinis |
| BA06              | B.altitudinis |
| AH820             | B.anthracis   |
| A0248             | B.anthracis   |
| 7\_6\_55CFAA\_CT2 | B.cereus      |
| VD200             | B.cereus      |

If a strain under study has not been assigned to a phylogenetic group,
an extra file containing the genome sequence is required in order to
perform this classification. (See **Case 2**).

  - Genome sequence files in fasta format should be named
    name\_of\_strain\_genomic.fna (i.e. `FZB42_genomic.fna`)

As an example, for *B. velezensis* FZB42 files containing protein,
feature and genomic data could be downloaded from
[ftp.ncbi.nlm.nih.gov](https://goo.gl/3oVBvv) following the url links
indicated below.

  - [Protein](https://goo.gl/rxjba3)
  - [Feature](https://goo.gl/YkyJyW)
  - [Genomic](https://goo.gl/N2J3js)

Finally, a file indicating the list of strains to be analysed should be
provided. It is suggested to name this file strains.txt. Strains should
be listed one per line without header.

`strains.txt` file content example:

|                   |
| ----------------- |
| 41KF2b            |
| KACC\_16563       |
| BA06              |
| AH820             |
| A0248             |
| 7\_6\_55CFAA\_CT2 |
| VD200             |

#### Section I.2. Phylogenetic group circumsription (optional).

#### ANI computation

Phylogenetic groups are defined as those that have an ANI \> 94%
compared to the group reference strain. OAT\_cmd and BLAST programs
should be preinstalled in order to define phylogenetic groups based on
ANI values. The script run\_ANI.R calculates the ANI of a list of
strains against a reference strain. It should be run for each reference
strain individually, using the following command-line:

``` bash
Rscript run_ANI.R REF_STR strains.txt OAT/ /usr/bin/ 4 Genomic/ ANIs/
```

Where:

  - `run_ANI.R` is the script file.
  - `REF_STR` is reference strain name (i.e. FZB42) as it is recorded in
    the genome file (i.e. FZB42\_genomic.fna).
  - `strains.txt` is the file containing the list of strains to be
    analysed. It should contain one strain per line.
  - `OAT/` is the path of OAT\_cmd.jar.
  - `/usr/bin/` is the path where Blast is installed in your system.
  - `4` is the number of processors to be used. Only positive integers
    are accepted.
  - `Genomic/` is the path where the fasta files of the strains under
    study are located. They should be named as
    `Name_of_the_strain_genomic.fna` (i.e. `FZB42_genomic.fna`).
  - `ANIs/` is the output path.

Files with extension `.ANI` will be created in the output path specified
for each strain listed in `strains.txt`. The file name format will be
`REF_STR_VS_strain.ANI` (i.e. `41KF2b_VS_168.ANI` where `41KF2b` is the
reference strain and 168 the strain under analysis).

#### Definition of phylogenetic groups.

The script `ANI_FUN.R` should be run to create the phylogenetic groups.
All ANI files for all analysed reference strains should be located in
the same folder.

Command-line example:

``` bash
Rscript ANI_FUN.R ref.csv strains.txt ANIs/
```

Where:

  - `ANI_FUN.R` is the script file.
  - `ref.csv` is a CSV table containing reference strains and species
    names in the first and second columns, respectively, as it is shown
    in the example:

| Ref\_Strain          | Species             |
| -------------------- | ------------------- |
| 41KF2b               | B.altitudinis       |
| DSM\_7               | B.amyloliquefaciens |
| Ames                 | B.anthracis         |
| BSS                  | B.atrophaeus        |
| ATCC\_14579          | B.cereus            |
| ATCC\_14580\_DSM\_13 | B.licheniformis     |
| 168                  | B.subtilis          |
| NRRL\_B\_41580       | B.velezensis        |
| FSL\_W8\_0169        | B.wiedmannii        |

  - `strains.txt` is the file containing the list of strains to be
    analysed. It should contain one strain per line.
  - `ANIs/` is the folder where ANIs were saved in the previous step.

The result will be the creation of two files:

  - `Phylogenetic_groups.csv` is a table describing the strains present
    in each phylogenetic group.
  - `ANI_table.csv` is a table with all computed ANI values.

#### Section I.3. BLAST searches.

GeM-Pro needs BLAST output files of all protein sequences (not just
those of interest) of the model strain (i.e. *B. velezensis* FZB42)
against all protein sequences encoded in all the strains to be profiled.
Blast.R script runs BLASTP using fasta files containing the strains
protein sequences (proteomes in this context) and creates output files
in table format (BLASTp tables).

Command-line
example:

``` bash
Rscript Blast.R REF_STR strains.txt ALL Path/to/protein/fasta/ Blast_out/
```

Where:

  - `Blast.R` is the script file.
  - `REF_STR` is the name of the reference strain (i.e. FZB42).
  - `strains.txt` is the file containing the list of the strains to be
    analysed. It should contain one strain per line.
  - `ALL` is the argument that indicates that all model strains proteins
    will be used as query
  - `Path/to/protein/fasta/` is the path to the folder containing
    protein sequences in fasta format.
  - `Blast_out/` is the path where the output BLAST files (in table
    format) will be created.

The result will be the creation of BLAST output files in table
format.

#### Section I.4. Protein sequence prediction in non-annotated genomes (optional).

It case of non-annotated genomes, gene predictions could be performed
with Prodigal, using command line similar to the one written bellow.

``` bash
prodigal -i STRAIN.fna -a STRAIN_protein.faa
```

For more information please refer to
<https://github.com/hyattpd/prodigal/wiki>

#### Section I.5. Converting features files.

Synteny analysis performed by GeM-Pro and q-GeM, needs each CDS start
and end positions. The script Subject\_feature.R obtains this
information from the feature file of each genome and creates and R
object that is used by GeM-Pro and q-GeM.

Command-line example:

``` bash
Rscript Subject_feature.R strains.txt Features/ _feature_table.txt
```

Where:

  - `Subject_feature.R` is the script file.
  - `strains.txt` is the file containing the list of the strains to be
    analyzed. It should contain one strain per line.
  - `Features/` is the path to the fature files folder.
  - `_feature_table.txt` is the suffix of the feature files.

The result will be the creation of a RDS file named
`Subject_Feature_list.RDS` that will be used by GeMPro.R (see
below).

#### Section I.6. Feature information for non-annotated genomes (optional)

In case the Feature file is not available for a genome, the script
`Subject_feature.R` would be used to obtain CDS start and end sites
using the protein predictions performed by Prodigal (see above).
`Subject_feature.R` will automatically notice the absence of the feature
file and will extract the needed information from the Prodigal-predicted
protein files. However, the location of such files should be indicated.

Command-line
example:

``` bash
Rscript Subject_feature.R strains.txt Features/ _feature_table.txt Proteomes/ _protein.faa
```

Where:

  - `Subject_feature.R` is the script file.
  - `strains.txt` is the file containing the list of the strains to be
    analyzed. It should contain one strain per line.
  - `Features/` is the path to the feature files folder.
  - `_feature_table.txt` is the suffix of feature files.
  - `Proteomes/` is the path to the folder where the Prodigal-predicted
    protein fasta files are stored.
  - `_protein.faa` is the protein fasta file suffix.

The result will be the creation of a RDS file named `Subject_feature.R`
that will be used by `GeMPro.R` (see below).

#### Section I.7. GeM-Pro execution

The pipeline can be run for one phylogenetic group in particular or for
all together.

Command-line
example:

``` bash
Rscript GeMPro.R REF_STR Phylogenetic_groups.csv Blast_out/ Query_pathways.csv Subject_Feature_list.RDS REF_STR_feature_table.txt Out/
```

Where:

  - `GeMPro.R` is the script file.
  - `Phylogenetic_groups.csv` is the file that indicates the
    phylogenetic groups. See Section 1 (Input files).
  - `REF_STR` is the name of the model strain (i.e. FZB42).
  - `Blast_out/` is the path where BLAST outputs are stored.
  - `Query_pathways.csv` is the file where genes and pathways are listed
    in CSV format. See Section 1 (Input files).
  - `Subject_Feature_list.RDS` is the RDS file generated by
    Subject\_feature.R.
  - `REF_STR_feature_table.txt` is the feature file for the `REF_STR`
    (i.e. `FZB42_feature_table.txt`). This file could be downloaded form
    NCBI databases. See Section 1 (Input files).
  - `OUT/` is the folder where the output files will be stored.

The result will be the creation of the output files listed below. The
four RDS files could be read using R.

  - `IWNO.RDS`: It contains for each strain the applied cut-off (PID
    cut-offmean) as well as the list of othologues and non-orthologues
    proteins that results from the Internal witness of non-orthology
    (IWNO) verification. The strains are grouped by their phylogenetic
    groups.
  - `GeMPro_DB.RDS`: It contains a list with the following objects:
      - `orthoDensity` and `nonorthoDensity`: Othologues and
        non-orthologues density function for each phylogenetic group.
      - `orthoProbability` and `nonorthoProbability`: Probability values
        for each homologous pair.
      - `BayesTable`: It contains the information of the homologues
        proteins listed in `Query_pathways.csv` (accesion number of the
        query and subject protein, name of the subject organism, group
        name, %ID, % of coverage shared between query and subject
        proteins, probability of othology and non-orthology, Bayes
        Factor of Orthology, the GeM-Pro decision and the score define
        as \(P_{(Orthology)}/(P_{(Non-orthology)} × 100)\).
      - `BayesTable.all`: Similar to BayesTable but for all homologues
        proteins.
  - `synteny.list.RDS`: It contains the synteny analysis result for each
    strain.  
  - `pathway.profiles.RDS`: It contains the pathway identified in each
    strain under analysis.
  - `pathway.profiles.txt`: an extra ITOL file required for the
    representation of the pathways in the clustering. It contains the
    ITOL parameters and the pathway identified as well as their synteny
    score for each strain under analysis. It could be read by any text
    editor.

#### Section I.8. Hierarchical clustering.

`Clustering.R` script makes the hierarchical clustering of the strains
under study using the information stored in `pathway.profiles.RDS`.
`Clustering.R` requires the R packages pvclust and ape for hierarchical
clustering and tree construction. They must be downloaded if there were
not installed in your computer.

Command-line example:

``` bash
Rscript Clustering.R pathway.profiles.RDS nboot
```

Where:

  - `Clustering.R` is the script file.
  - `pathway.profiles.RDS` is a RDS file that contains the pathway
    identified in each strain under analysis. It is created by the
    script GeMPro.R.
  - `nboot` is the number of bootstrapping repetitions to make. See
    pvclust documentation for more details
    [here](https://cran.r-project.org/web/packages/pvclust/pvclust.pdf).

The result will be the creation of three files. `Clustering.RDS` could
be read using R whereas `Clustering.tre` and `Clustering.pdf` are
dendrograms that could be visualized with [iTOL](https://itol.embl.de/)
or any pdf viewer, respectively. To add shape plots to the dendrogram in
iTOL, the information in pathway.profiles.txt (created by GeMPro.R)
should be used. More information [here](https://itol.embl.de/help.cgi).

## II. q-GeM

q-GeM defines orthology classifiers based on `GeMPro_DB.RDS`
information. q-GeM could be used to profile new strains in a time and
computational efficient way and therefore would be applied in larger
datasets than GeM-Pro. However, as q-GeM uses GeM-Pro-DB the
phylogenetic group had to be analysed by GeM-Pro. Obviously, the strains
phylogenic group must be known (and verified).

### q-GeM Quick Tutorial

Two cases are considered in this quick tutorial.

**Case 4:** The aim is to profile strains with known (and verified\!\!)
phylogenetic group. The phylogenetic group was previously analyzed with
GeM-Pro.

More information can be found in the [Example
Case](./Examples/Example_Case3.md)

*Step 1*. Define Best
classifier.

| Input files    | Script (section II.2) | Output files                           |
| -------------- | --------------------- | -------------------------------------- |
| GeMPro\_DB.RDS | Best\_Classifier.R    | Name\_of\_group\_Best\_class.RDS files |

Example:

``` bash
Rscript Best_Classifier.R GeMPro_DB.RDS Sample/Best_Classifiers/
```

*Step 2*. Create gene/pathways (Query\_pathways.csv) and strains.txt.

*Step 3*. Perform BLAST
searches

| Input files                                                         | Script (section II.3) | Output files       |
| ------------------------------------------------------------------- | --------------------- | ------------------ |
| strains.txt <br> Query\_pathways.csv <br> Strain\_protein.faa files | Blast.R               | output BLAST files |

Example:

``` bash
Rscript Blast.R FZB42 Sample/strain.C4.txt Sample/Query_pathways.csv Sample/Proteomes/ Sample/Blast_out/
```

*Step 4*. Convert feature files to RDS
files.

| Input files                                       | Script (section I.5) | Output files               |
| ------------------------------------------------- | -------------------- | -------------------------- |
| strains.txt <br> Strain\_feature\_table.txt files | Subject\_feature.R   | Subject\_Feature\_list.RDS |

Example:

``` bash
Rscript Subject_feature.R Sample/strains.C4.txt Sample/Features/ _feature_table.txt
```

*Step 5*. Run
q-GeM.

| Input files                                                                                                                                                                             | Script (section II.5) | Output files                                                                         |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- | ------------------------------------------------------------------------------------ |
| Phylogenetic\_groups.csv <br> Query\_pathways.csv <br> output BLAST files <br> Name\_of\_group\_Best\_class.RDS files <br> Subject\_Feature\_list.RDS <br> REF\_STR\_feature\_table.txt | qGeM.R                | Orthologous.RDS <br> Synteny.RDS <br> pathway.profiles.RDS <br> pathway.profiles.txt |

Example:

``` bash
Rscript qGeM.R FZB42 Sample/Phylogenetic_groups.C4.csv Sample/Blast_out/ Sample/Query_pathways.csv Sample/Best_classifiers/ Sample/Subject_Feature_list.RDS Sample/Features/FZB42_feature_table.txt Sample/qGeM_out/
```

*Step 6*. Perform the
clustering.

| Input files          | Script (section I.8) | Output files                       |
| -------------------- | -------------------- | ---------------------------------- |
| pathway.profiles.RDS | Clustering.R         | Clustering.RDS <br> Clustering.tre |

Example:

``` bash
Rscript Clustering.R pathway.profiles.RDS 100
```

## q-GeM Tutorial

q-GeM pursues the same driving idea than GeM-Pro. q-GeM will predict the
existence of isofunctional genes/pathways regarding a model strain based
on the “orthology conjecture”. The result will be the creation of a
dendrogram describing strain gene/pathway content.

#### Section II.1. Input files.

Similar to GeM-Pro (See Section I.1. (Input files) for details), q-GeM
will need:

  - Information of the model strain genes/pathways that sould be
    provided in a table named `Query_pathways.csv`.
  - The list of strains to be analysed in a file named `strains.txt`.
  - A fasta file with all protein sequences for each strain (i.e.
    `FZB42_protein.faa`).
  - Features file for each strain (i.e. `FZB42_feature_table.txt`).
  - The information of each strain phylogenetic group in a table named
    `Phylogenetic_groups.csv`.
  - Additionally, q-GeM will need the file `GeMPro_DB.RDS` that should
    include the information of all phylogenetic groups to be analysed.
    See Section I.7.

#### Section II.2. Best Classifier selection.

The script `Best_Classifier.R` will define the best classifier for each
phylogenetic group. A `GeMPro_DB.RDS` file that includes information of
each phylogenetic group must be provided.

Command-line example:

``` bash
Rscript Best_Classifier.R GeMPro_DB.RDS Best_Classifiers/
```

Where:

  - `Best_Classifier.R` is the name of the script.
  - `GeMPro_DB.RDS` is the RDS file created by `GeMPro.R`.
  - `Best_Classifiers/` is the path where the output RDS files will be
    created.

The result will be the creation of many output files as phylogenetic
groups were included in the `GeMPro_DB.RDS` file. Each one will be named
`Name_of_group_Best_class.RDS` (i.e. `B.cereus_Best_class.RDS`). These
files will be used to predict orthology of genes of interest in new
strains.

#### Section II.3. BLAST searches.

**IMPORTANT**: q-GeM doesn’t need to compute all blasts of all protein
sequences against all of model strain proteins. It performs it decision
of orthology or non-orthology based on GeMPro-DB information. Hence, to
perform the classification only the Blast output of the genes listed in
`Query_pathways.csv` are needed. `Blast.R` script could be used to
obtain that.

Command-line
example:

``` bash
Rscript Blast.R REF_STR strain.txt Query_pathways.csv Path/to/protein/fasta/ Blast_out/
```

Where:

  - `Blast.R` is the script file
  - `REF_STR` is the name of the reference or model strain (i.e. FZB42).
  - `strains.txt` is the file containing the list of strains to be
    analysed. It should contain one strain per line.
  - `Query_pathways.csv` is the argument that indicates the model strain
    proteins (REF\_STR) that will be used as query.
  - `Path/to/protein/fasta/` is the path to the files containing the
    protein sequences in fasta format.
  - `Blast_out`/ is the path where the BLAST output files (in table
    format) will be created.

The result will be the creation of BLAST output files in table format.
In addition, a file named `sub_fasta.faa` will be created containing the
model strain protein sequence listed in `Query_pathways.csv`
(i.e. FZB42).

#### Section II.4. Converting features files.

Synteny analysis performed by q-GeM needs each CDS start and end
position. `Subject_feature.R` could be used as was described in Section
I.5. (Converting feature files).

#### Section II.5. q-GeM execution.

The q-GeM pipeline will use the best classifier for each strain
phylogenetic group.

Command-line
example:

``` bash
Rscript qGeM.R REF_STR Phylogenetic_groups.csv Blast_out/ Query_pathways.csv path/to/classifier/ Subject_Feature_list.RDS REF_STR_feature_table.txt
```

Where:

  - `qGeM.R` is the script file.
  - `REF_STR` is the name of the reference or model strain (i.e. FZB42).
  - `Phylogenetic_groups.csv` is the file that indicates the
    phylogenetic groups. See Section 1 (Input files).
  - `Blast_out/` is the path where the output BLAST are stored.
  - `Query_pathways.csv` is the file where genes and pathways are listed
    in CSV format. See Section I.1 (Input files).
  - `path/to/classifier/` is the path where the files the RDS file
    obtained with the `Best_Classifier.R` script
    (i.e.`B.cereus_Best_class.RDS`) are stored. See Section II.2. (Best
    Classifier selection).
  - `Subject_Feature_list.RDS` is the RDS file obtained with
    `Subject_feature.R` script. See Section I.5. (Converting features
    files).
  - `REF_STR_feature_table.txt` is the feature file for the `REF_STR`
    (i.e. `FZB42_feature_table.txt`). This file could be downloaded from
    NCBI databases. See Section 1 (Input files).

The result will be the creation of the files listed below. The three RDS
files could be read using R.

  - `Orthologous.RDS`: A list containing, as R dataframes, the blast
    output of the orthologous proteins of those listed in
    `Query_pathways.csv` for each strain.
  - `Synteny.RDS`: A list containing 2 objects:
      - `syn_all`:A list containing, as R dataframes, the results of the
        synteny analysis for each strain.
      - `syn_present`: A list containing, as R dataframes, the results
        of the synteny analysis for the pathways present for each
        strain.
  - `pathway.profiles.RDS`: It contains the identified pathway in each
    strain under analysis.
  - `pathway.profiles.csv`: an extra ITOL file required for the
    representation of the pathways in the clustering. It contains the
    ITOL parameters and the pathway identified as well as their synteny
    score for each strain under analysis. It could be read by any text
    editor.

#### Section II.6. Hierarchical clustering.

`Clustering.R` script cloud make the hierarchical clustering of the
strains under study using the information stored in
`pathway.profiles.RDS` as was described in Section I.8. (Hierarchical
clustering).
