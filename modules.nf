process run_resfinder_reads {
	label 'resfinder'
   input:
	tuple val(sample_id), path(reads), val(mlst_scheme)
   output:
        tuple val(sample_id), path("resfinder_reads_out")

   script:
   """
   organism=\$(translate_mlst_scheme_species.py --tool resfinder --name $mlst_scheme --translation_table ${projectDir}/data/mlst_schemes.txt)
   if [ ! -z "\$organism" ]; then
       python -m resfinder -ifq $reads -o resfinder_reads_out -s "\$organism" --acquired \
           --point -db_point $params.pointfinder_db -db_res $params.resfinder_db
   else
       python -m resfinder -ifq $reads -o resfinder_reads_out --acquired \
           -db_res $params.resfinder_db
   fi
   
   """
}

process run_resfinder_assembly {
   label 'resfinder'
   input:
        tuple val(sample_id), path(assembly), val(mlst_scheme)
   output:
        tuple val(sample_id), path(assembly), path("resfinder_assembly_out")

   script:
   """
   organism=\$(translate_mlst_scheme_species.py --tool resfinder --name $mlst_scheme --translation_table ${projectDir}/data/mlst_schemes.txt)
   if [ ! -z "\$organism" ]; then
       python -m resfinder -ifa $assembly -o resfinder_assembly_out -s "\$organism" --acquired \
            --point -db_point $params.pointfinder_db -db_res $params.resfinder_db 
   else
       python -m resfinder -ifa $assembly -o resfinder_assembly_out --acquired \
            -db_res $params.resfinder_db 
   fi
   """
}

process filter_assembly {
   errorStrategy "ignore" // fails in case no resistance genes found
   label 'amrdb'
   input:
	tuple val(sample_id), path(assembly), path(resfinder_output), path(armfinder_output), path(plasmidfinder_output)
   output:
        tuple val(sample_id), path("*.fa") optional true

   script:
   """
   new_assembly_name=${assembly.simpleName}.fa
   filter_contigs.py $assembly $resfinder_output/ResFinder_results_tab.txt $armfinder_output/amrfinder_results.txt $plasmidfinder_output > \$new_assembly_name

   if [ ! -s \$new_assembly_name ]; then
       rm \$new_assembly_name
   fi
   """
}

process run_mlst {
   label 'mlst'
   input:
        tuple val(sample_id), path(assembly)
   output:
        tuple val(sample_id), path("mlst_output.csv"), env(scheme)
   script:
   """
   mlst --quiet --csv $assembly > mlst_output.csv
   scheme=\$(cut -f2 -d, mlst_output.csv)
   """
}

process run_amrfinder {
   label 'amrfinder'
   cpus 4
   input:
        tuple val(sample_id), path(assembly), val(mlst_scheme)
   output:
        tuple val(sample_id), path("amrfinder_output")

   script:
   """
   organism=\$(translate_mlst_scheme_species.py --tool amrfinder --name $mlst_scheme --translation_table ${projectDir}/data/mlst_schemes.txt)
   mkdir amrfinder_output
   if [ ! -z "\$organism" ]; then
       amrfinder -n $assembly --database ${params.amrfinder_db}/latest --organism "\$organism" \
           -o amrfinder_output/amrfinder_results.txt  \
           --nucleotide_output amrfinder_output/amrfinder_nucleotides.fasta
   else
       amrfinder -n $assembly --database ${params.amrfinder_db}/latest \
           -o amrfinder_output/amrfinder_results.txt  \
           --nucleotide_output amrfinder_output/amrfinder_nucleotides.fasta
   fi
   """
}

process run_isescan {
   label 'isescan'
   cpus params.threads
   input:
        tuple val(sample_id), path(assembly)
   output:
        tuple val(sample_id), path("output/${assembly}.tsv") optional true

   script:
   """
   isescan.py --seqfile $assembly --output output --nthread ${params.threads}
   """
}

process run_bakta {
    label 'bakta'
    cpus params.threads
    input:
        tuple val(sample_id), path(assembly)
    output:
        tuple val(sample_id), path("bakta/${assembly}.tsv"), path("bakta/${assembly}.gbff")
    script:
    """
    bakta --db ${params.bakta_db} --keep-contig-headers --threads ${params.threads} --output bakta --prefix $assembly --force $assembly --skip-plot 
    """
}

process run_phispy {
    label 'phispy'
    cpus params.threads
    input:
        tuple val(sample_id), path(assembly), path(annotation)
    output:
        tuple val(sample_id), path("phispy/prophage_coordinates.tsv")
    script:
    """
	phispy $annotation -o phispy --threads ${params.threads}
    """
}

process run_speciesfinder {
    label 'speciesfinder'
    input:
        tuple val(sample_id), path(assembly) 
    output:
        tuple val(sample_id), path("data.txt")
    script:
    """
    mkdir -p tmp
    speciesfinder.py -i $assembly -o . -p ${params.speciesfinder_db}/ -t tmp/
    """
}

process run_plasmidfinder {
    label 'plasmidfinder'
    input:
        tuple val(sample_id), path(assembly) 
    output:
        tuple val(sample_id), path("results/results_tab.tsv")
    script:
    """
    mkdir -p results
	plasmidfinder.py -i $assembly  -o results -p ${params.plasmidfinder_db} -x 
    """
}
     
process add_result_to_database {
   label 'amrdb'
   secret 'amrdbpw'
   secret 'amrdbuser'

   input:
        tuple val(sample_id), val(method), path(input_path), val(db_version), val(tool_version)
   output:
        tuple val(sample_id), val("success")
   script:
   if (params.use_sqlite)
   """
   amrdb_add_results.py -d ${params.sqlite_db_path} \
       --method $method -i $input_path --external_id $sample_id --tool_version "$tool_version" \
       --db_version "$db_version"
   """
   else
   """
   amrdb_add_results.py -H ${params.amrdb_host} -d ${params.amrdb_name} -p "\$amrdbpw" -u "\$amrdbuser" \
       --method $method -i $input_path --external_id $sample_id --tool_version "$tool_version" \
       --db_version "$db_version"
   """
}


process add_assembly_to_database {
   label 'amrdb'
   secret 'amrdbpw'
   secret 'amrdbuser'

   input:
	tuple val(sample_id), val(sample_name), val(method), path(input_path), path(assembly), val(db_version), val(tool_version)
   output:
        val true
   script:
   if (params.use_sqlite)
   """
   amrdb_add_results.py -d ${params.sqlite_db_path} \
       --method $method -i $input_path --assembly $assembly --external_id $sample_id \
       --sample_name $sample_name --tool_version "$tool_version" --db_version "$db_version"
   """
   else
   """
   amrdb_add_results.py -H ${params.amrdb_host} -d ${params.amrdb_name} -p \$amrdbpw -u \$amrdbuser \
       --method $method -i $input_path --assembly $assembly --external_id $sample_id \
       --sample_name $sample_name --tool_version "$tool_version" --db_version "$db_version"
   """
}


process mob_suite {
   label 'mobsuite'
   input:
        tuple val(sample_id), path(assembly)
   output:
        tuple val(sample_id), path("mobtyper_results.txt")

   script:
   """
   mob_typer --multi --infile $assembly --out_file mobtyper_results.txt --database_directory ${params.mob_typer_db}
   """
}

process get_resfinder_version {
   label 'resfinder'
   output:
        tuple val("resfinder"), env(db_version), env(tool_version)

   script:
   """
   db_version=\$(cat ${params.resfinder_db}/VERSION)
   tool_version=\$(python -m resfinder --version)
   """

}

process get_speciesfinder_version {
   label 'speciesfinder'
   output:
        tuple val("speciesfinder"), env(db_version), env(tool_version)

   script:
   """
   db_version=\$(cat ${params.speciesfinder_db}/VERSION)
   tool_version="2022-05-30" #hardcoded, as script is currently distributed with workflow
   """

}

process get_amrfinder_version {
   label 'amrfinder'
   output:
        tuple val("amrfinder"), env(db_version), env(tool_version)

   script:
   """
   db_version=\$(cat ${params.amrfinder_db}/latest/version.txt)
   tool_version=\$(amrfinder --version)
   """
}

process get_isescan_version {
   label 'isescan'
   output:
        tuple val("isescan"), env(db_version), env(tool_version)

   script:
   """
   db_version="unknown"
   tool_version=\$(isescan.py --version)
   """
}

process get_bakta_version {
   label 'bakta'
   output:
        tuple val("bakta"), env(db_version), env(tool_version)
   script:
   """
   date=\$(cat ${params.bakta_db}"/version.json" | jq -r .date)
   v1=\$(cat ${params.bakta_db}"/version.json" | jq -r .major)
   v2=\$(cat ${params.bakta_db}"/version.json" | jq -r .minor)
   type=\$(cat ${params.bakta_db}"/version.json" | jq -r .type)
   db_version=(\$date"-"\$v1"."\$v2"-"\$type)
   tool_version=\$(bakta --version)
   """
}

process get_mlst_version {
   label 'mlst'
   output:
        tuple val("mlst"), env(db_version), env(tool_version)
   script:
   """
   # use install date of mlst as db_version
   db_version=\$(stat \$(which mlst) | grep Change | cut -f2 -d" ")
   tool_version=\$(mlst --version | cut -f2 -d" ")
   """
}

process get_mobtyper_version {
   label 'mobsuite'
   output:
        tuple val("mobtyper"), env(db_version), env(tool_version)

   script:
   """
   db_version=\$(cat $params.mob_typer_db"/status.txt" | grep -Eo 'Download date: [[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
   tool_version=\$(mob_typer --version)
   """
}

process get_plasmidfinder_version {
	label 'plasmidfinder'
	output:
		tuple val("plasmidfinder"), env(db_version), env(tool_version)
	script:
	"""
	db_version=\$(cat $params.plasmidfinder_db"/VERSION.txt")
	#tool_version=\$(mamba env export | grep '\\- plasmidfinder' | grep -E -o '[0-9]\\.[0-9]\\.[0-9]')
	tool_version="2.1.6"
	"""
}

process get_phispy_version {
	label 'phispy'
	output:
		tuple val("phispy"), env(db_version), env(tool_version)
	script:
	"""
	db_version="-"
	tool_version=\$(phispy -v)
	"""
}

process update_process_status {
    label 'amrdb'
    secret 'mysql_user'
    secret 'mysql_pw'
    input:
        tuple val(processId), val(status)
    script:
    """
    update_process_status.py -H ${params.db_host} -d ${params.db_name} -u \$mysql_user -p \$mysql_pw \
            --status $status --ids "${processId.collect { it.toString() }.join(",")}"
    """
}
