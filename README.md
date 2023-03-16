# Antibiotika Resistenzgen Datenbank Workflow

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
This entry expects 7 values per process which are:


| Field             | Type          | Description                       |
|-------------------|---------------|-----------------------------------|
| assembly_id       | int(11)       | unique identifier for input data  |
| process_id        | int(11)       | foreign key of reference table    |
| import_name       | varchar(200)  | name of the sample in database    |
| path_fasta        | varchar(1000) | path to fasta file                |
| path_fastq1       | varchar(1000) | path to forward reads             |
| path_fastq2       | varchar(1000) | path to reverse reads [nullable]  |
| species_name      | varchar(100)  | used for input in resfinder       |


The `process_id` should point to a table `process` which has a `status` varchar column that can be set to "finished" in order
 to signal that this compute job has been finished, and will not return in the next SQL-query.







<!--
### Update Databases
login on wsps1152 as user seqsphere

`cd /db/`

#### Amrfinder

`conda activate amrfinder`

`amrfinder -u`

copy files from '/proj/seqsphere/conda/mambaforge/envs/amrfinder/share/amrfinderplus/data/2023-04-17.1' to /db/amrfinder


#### Resfinder
rm -r resfinder_db/
rm -r pointfinder_db/
rm -r disinfinder_db/

git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/
-->

