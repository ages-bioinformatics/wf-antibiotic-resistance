process create_resfinder_db {
	label "resfinder"
	output:
		val("success")
	script:
	"""
	rm -rf ${params.resfinder_db}
	mkdir -p ${params.resfinder_db}
	git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git ${params.resfinder_db}
	cd ${params.resfinder_db}
	python INSTALL.py
	"""
}


process create_pointfinder_db {
        label "resfinder"
	// input just defined to avoid simultaneous execution of all processes in parallel
        input:
                val(start)
        output:
                val("success")
	script:
	"""
	rm -rf ${params.pointfinder_db}
	mkdir -p ${params.pointfinder_db}
	git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git ${params.pointfinder_db}
	cd ${params.pointfinder_db}
	python INSTALL.py
	"""
}


process create_amrdb {
	label "amrdb"
	secret 'amrdbpw'
	secret 'amrdbuser'
	input:
		val(resfinder_db_created)
		val(amrfinder_db_created)
	script:
	if (params.use_sqlite == false)
	"""
	update_resfinder_database.py --database ${params.amrdb_name} --hostname ${params.amrdb_host} --user \$amrdbuser --password \$amrdbpw --resfinder_db ${params.resfinder_db} --amrfinder_db ${params.amrfinder_db}/latest
	"""
	else
	"""
	update_resfinder_database.py --database ${params.sqlite_db_path} --resfinder_db ${params.resfinder_db} --amrfinder_db ${params.amrfinder_db}/latest
	"""
}


process create_bakta_db {
	label "bakta"
	// input just defined to avoid simultaneous execution of all processes in parallel
        input:
                val(start)
        output:
                val("success")
	script:
	"""
	target_dir=\$(dirname ${params.bakta_db}
	mkdir -p \$target_dir
	bakta_db download --output \$target_dir --type light
	""" 
}


process create_mobtyper_db {
	label "mobsuite"
	cpus params.threads
	// input just defined to avoid simultaneous execution of all processes in parallel
        input:
                val(start)
        output:
                val("success")
	script:
	"""
	echo -e ">HelloWorld\nAACCGGTT" > test.fasta
	mob_typer -i test.fasta -o results -n $params.threads -d ${params.mob_typer_db}
	"""
}


process create_amrfinder_db {
	label "amrfinder"
	output:
		val("success")
	script:
	"""
	amrfinder_update -d ${params.amrfinder_db}
	"""
}

process create_speciesfinder_db {
	label "speciesfinder"
	// input just defined to avoid simultaneous execution of all processes in parallel
        input:
                val(start)
        output:
                val("success")
	script:
	"""
	rm -rf ${params.speciesfinder_db}
	mkdir -p ${params.speciesfinder_db}
	git clone https://bitbucket.org/genomicepidemiology/speciesfinder_db.git ${params.speciesfinder_db}
	wget -O - https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSUParc_tax_silva.fasta.gz | \
		gunzip -c > ${params.speciesfinder_db}/16srna_database.fasta
	pushd ${params.speciesfinder_db}
	python INSTALL.py
	echo "SILVA_138.1_SSUParc" > VERSION
	"""
}

