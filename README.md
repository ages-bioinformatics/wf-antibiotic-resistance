# Antibiotics Resistence Gene Datenbase Workflow

This Nextflow-Workflow was developed for *in-silico* prediction of antimicrobial resistance genes.
The workflow is built arround [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder) 
and its main goal is to import all results into a dedicated database (currently either SQLite or MariaDB).  
  

## Requirements 
* [Nextflow](https://www.nextflow.io/)
* [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[Mamba](https://github.com/mamba-org/mamba)

Nextflow will prepare and install required software tools via conda. Approximately 5 GB of hdd-memory
 should be available for tool installation and 22 GB for sequence databases.  
   
Depending on the number of samples that should be processed, additional hdd-memory should be reserved for
 temporary and intermediate data (~50-100 GB for 3,000 samples in a single workflow). NextFlow
 integrates a large variety of parallelization-modes via [cloud or HPC-clusters](https://www.nextflow.io/docs/latest/executor.html),
  which can be configured in the `nextflow.config` file.
 A single sample needs about 0.2 core-hours on average. 

## Optional Database Server

It is recommended to use a database server (currently MariaDB/MySQL are supported). An empty database and database user
with SELECT/INSERT/UPDATE/DELETE/CREATE privileges should be prepared. The scheme (i.e. all tables) will be created
during the first `-entry upgrade` call, just after the amrfinder and resfinder sequence-databases have been downloaded.  
  
If no appropriate server is available, SQLite can be used,
 by adding `--use_sqlite` as parameter to the Workflow.

## Usage

### Set Nextflow secrets

To avoid putting sensitive information into code repositories, this workflow uses NextFlow-[Secrets](https://www.nextflow.io/docs/latest/secrets.html). To get started, you may set the secrets in your runtime as follows:

```
export NXF_ENABLE_SECRETS=true
nextflow secrets put -n amrdbpw -v PASSWORD
nextflow secrets put -n amrdbuser -v USERNAME
```
  
This only applies if a database is used as input or a database-server is used for the resistance-gene-database.

###  Install/Updata databases and install dependencies

The workflow is easy to set up. Before the first usage, run the upgrade-workflow to download all required
datbases and initialize the db-scheme for your selected database. Make sure to set up the `amrdb_name` and `amr_db`
in `nextflow.config` appropriately.  

```
nextflow run main.nf -entry upgrade [--db_dir path/to/target/dir] [--use_sqlite]
```

### Input Data

#### From Input-Directory
The main input type is genome assemblies (or complete genomes) in fasta-format. It is possible to provide fastq-files
(only used for resfinder) in addition.  

The default entry for the workflow uses `--in /path/to/data` parameter and recursively searches all (paired-)fastq and fasta-files,
 and matches both input types if possible.
The sample-name (basename without extensions of the fasta-file) will be used as the main sample identifier and imported
into the database under this name.  
  
#### Advanced: From other Database
Alternatively, a database connection and query can be specified in the `nextflow.config` 
file and the entry `-entry from_db` should be used. 

This entry expects 6 values per process which are:


| Field             | Type          | Description                       |
|-------------------|---------------|-----------------------------------|
| assembly_id       | int(11)       | unique identifier for input data  |
| process_id        | int(11)       | foreign key of reference table    |
| import_name       | varchar(200)  | name of the sample in database    |
| path_fasta        | varchar(1000) | path to fasta file                |
| path_fastq1       | varchar(1000) | path to forward reads             |
| path_fastq2       | varchar(1000) | path to reverse reads [nullable]  |

The `process_id` should point to a table `process` which has a `status` varchar column that can be set to "finished" in order
 to signal that this compute job has been finished, and will not return in the next SQL-query.

## Citations ##

This workflow relies heavily on work of others. This includes:  
- mlst: Seemann T, https://github.com/tseemann/mlst
- ResFinder: Bortolaia V et al., https://bitbucket.org/genomicepidemiology/resfinder
- AMRFinderPlus: Feldgarden M et al., https://github.com/ncbi/amr
- PhiSpy: Sajia A et al., https://github.com/linsalrob/PhiSpy
- PlasmidFinder: Carrattoli A et al., https://bitbucket.org/genomicepidemiology/plasmidfinder
- SpeciesFinder: Clausen PTLC et al., https://bitbucket.org/genomicepidemiology/speciesfinder
- Bakta: Schwengers O et al., https://github.com/oschwengers/bakta
- ISEScan: Xie Z and Tang H, https://github.com/xiezhq/ISEScan
- MOB-suite: Robertson J et al., https://github.com/phac-nml/mob-suite

### mlst ###
"This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust".

Seemann T, mlst Github https://github.com/tseemann/mlst


### ResFinder ###
ResFinder 4.0 for predictions of phenotypes from genotypes.
Bortolaia V, Kaas RS, Ruppe E, Roberts MC, Schwarz S, Cattoir V, Philippon A, Allesoe RL, Rebelo AR, Florensa AR, Fagelhauer L, Chakraborty T, Neumann B, Werner G, Bender JK, Stingl K, Nguyen M, Coppens J, Xavier BB, Malhotra-Kumar S, Westh H, Pinholt M, Anjum MF, Duggett NA, Kempf I, Nykasenoja S, Olkkola S, Wieczorek K, Amaro A, Clemente L, Mossong J, Losch S, Ragimbeau C, Lund O, Aarestrup FM. Journal of Antimicrobial Chemotherapy. 2020 Aug 11.
PMID: 32780112 doi: 10.1093/jac/dkaa345

### AMRFinderPlus ###
Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984.

### PhiSpy ###
Sajia Akhter, Katelyn McNair and Przemyslaw Decewicz. doi: https://doi.org/10.5281/zenodo.5945762

### PlasmidFinder ###
PlasmidFinder and pMLST: in silico detection and typing of plasmids. Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H. Antimicrob. Agents Chemother. 2014. April 28th. [Epub ahead of print]

### SpeciesFinder ###
Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, "Rapid and precise alignment of raw reads against redundant databases with KMA", BMC Bioinformatics, 2018;19:307.

### Bakta ###
Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685

### ISEScan ###
Zhiqun Xie, Haixu Tang. ISEScan: automated identification of Insertion Sequence Elements in prokaryotic genomes. Bioinformatics, 2017, 33(21): 3340-3347.

### MOB-Suite ###
Robertson, James et al. “Universal whole-sequence-based plasmid typing and its utility to prediction of host range and epidemiological surveillance.” Microbial genomics vol. 6,10 (2020): mgen000435. doi:10.1099/mgen.0.000435
