include {
    run_resfinder_reads;
    run_resfinder_assembly;
    run_amrfinder;
    run_isescan;
    run_bakta;
    run_phispy;
    run_mlst;
    run_speciesfinder;
    run_plasmidfinder;
    add_assembly_to_database;
    add_result_to_database;
    mob_suite;
    filter_assembly;
    get_resfinder_version;
    get_amrfinder_version;
    get_mobtyper_version;
    get_bakta_version;
    get_isescan_version;
    get_mlst_version;
    get_speciesfinder_version;
    get_plasmidfinder_version;
    get_phispy_version;
    update_process_status;
} from './modules.nf'

include {
    create_resfinder_db;
    create_amrdb;
    create_bakta_db;
    create_pointfinder_db;
    create_mobtyper_db;
    create_amrfinder_db;
    create_speciesfinder_db;
    create_plasmidfinder_db;
    create_phispy_db;
    create_mlst_db;
} from './setup.nf'


workflow from_db {
    // workflow entry is from database
    input_channel = Channel.sql.fromQuery(params.db_input_query, db: params.db_name)

    amr_detection(input_channel)
    amr_detection.out
	.groupTuple()
        .filter{it[1].every{ it == "success" }}
        .map { it -> [it[0], "finished"] }
	.set { successful_jobs }
    successful_jobs
        .join(input_channel.map{ it -> [it[0], it[2]] }, remainder: true)
	.filter { it[1] == null }
	.map { it -> [it[0], "failed"] }
        .set { failed_jobs }
    //assumes, that each assembly has one a single process running simulateously!
    successful_jobs
        .concat(failed_jobs)
        .join(input_channel.map { it -> [it[0], it[1]] })
        .map { it -> [it[2], it[1]] }
        .groupTuple(by: 1)
        .set { update_process_status_ch }
    update_process_status(update_process_status_ch)
}

workflow {
    // default create input channel from input path (tries to join readfiles on samplename)
    // fastq-input is optional
    Channel.fromPath(["${params.in}/**.fasta", "${params.in}/**.fa", "${params.in}/**.fna"])
        .map { it -> [it.getBaseName(), it] }
        .set { fasta_ch }
    fasta_ch.join(Channel.fromFilePairs(["${params.in}/**_R{1,2}*.f*q*", "${params.in}/**_{1,2}.f*q*"], flat: true), remainder:true)
        .map { it -> [0, 0, it[0], it[1], it[2], it[3], null] }
        .first()
        .view()
        .set { input_ch } 
    amr_detection(input_ch)
}

workflow amr_detection {
    take:
        input_channel
    main:
	// use mlst first to get an idea which species-parameter should be used
	run_mlst(input_channel.map{ it -> [it[0], it[3]] })
        // handle single end (currently not implemented for input-directory entry)
        input_channel.map { it -> [it[0], (it[5] == null ? [it[4]]:[it[4], it[5]])] }
            .filter{it[1].every{ it !== null }} //omit where no fastq in input
	    .join(run_mlst.out.map { it -> [it[0], it[2]] })
            .set { fastq_in }
        // main tools run
	input_channel
	    .map { it -> [it[0], it[3]] }
	    .join(run_mlst.out.map { it -> [it[0], it[2]] })
	    .set { assembly_species_tagged }
        run_resfinder_reads(fastq_in)
        run_resfinder_assembly(assembly_species_tagged)
        run_amrfinder(assembly_species_tagged)
        run_speciesfinder(input_channel.map{ it -> [it[0], it[3]] })
        run_plasmidfinder(input_channel.map{ it -> [it[0], it[3]] })

        // contigs with resfinder or amrfinderresults will be used for annotation
        // however, they will be cropped to +/- 50kbp around resistance genes
        filter_assembly(run_resfinder_assembly.out.join(run_amrfinder.out).join(run_plasmidfinder.out))
        run_isescan(filter_assembly.out)
        run_bakta(filter_assembly.out)
        mob_suite(filter_assembly.out)
		run_phispy(run_bakta.out)
		

        // check tool and database versions 
        versions = get_versions()

        // prepare input for amrfinder and resfinder (assembly file is added to import cmd
        // to be able to import contig size)
        input_channel.map { it -> [it[0], it[2]] }
            .join(run_resfinder_assembly.out)
            .map { it -> ["resfinder", it[0], it[1], it[3], it[2]] }
            .set { resfinder_import_ch }

        input_channel.map { it -> [it[0], it[2], it[3]] }
            .join(run_plasmidfinder.out)
            .map { it -> ["plasmidfinder", it[0], it[1], it[3], it[2]] }
            .set { plasmidfinder_import_ch }

        input_channel.map { it -> [it[0], it[2], it[3]] }
            .join(run_amrfinder.out)
            .map { it -> ["amrfinder", it[0], it[1], it[3], it[2]] }
            .concat(resfinder_import_ch)
            .concat(plasmidfinder_import_ch)
            .combine(versions, by: 0)
            .map { it -> [it[1], it[2], it[0], it[3], it[4], it[5], it[6]] }
            .set { assembly_import_ch }
            
            //[sample_id, sample_name, method, input_path, assembly, db_version, tool_version]
        // special form of input, which allows to import contigs before so that other
        // processes which might not create contig-db-entries can add results
        add_assembly_to_database(assembly_import_ch)


        // combine other outputs with version and prepare for db import
        run_bakta.out.map { it -> ["bakta", it[0], it[1]] }
           .concat(run_isescan.out.map { it -> ["isescan", it[0], it[1]] })
           .concat(run_resfinder_reads.out.map { it -> ["resfinder", it[0], it[1]] })
           .concat(mob_suite.out.map { it -> ["mobtyper", it[0], it[1]] })
           .concat(run_speciesfinder.out.map { it -> ["speciesfinder", it[0], it[1]] })
           .concat(run_mlst.out.map { it -> ["mlst", it[0], it[1]] })
           .concat(run_phispy.out.map { it -> ["phispy", it[0], it[1]] })
           .combine(versions, by:0).map { it -> [it[1], it[0], it[2], it[3], it[4],] }
           .set { add_db }
		add_db.view()
        add_result_to_database(add_db)
    // returns Channel consisting of [sample_id, "finished"], only needed for from_db entry
    emit:
        add_result_to_database.out
    }
    
workflow get_versions {
    get_resfinder_version()
    get_amrfinder_version()
    get_mobtyper_version()
    get_bakta_version()
    get_isescan_version()
    get_mlst_version()
    get_speciesfinder_version()
    get_plasmidfinder_version()
    get_phispy_version()

    versions_channel = get_resfinder_version.out
        .concat(get_amrfinder_version.out)
        .concat(get_mobtyper_version.out)
        .concat(get_bakta_version.out)
        .concat(get_isescan_version.out)
        .concat(get_mlst_version.out)
        .concat(get_speciesfinder_version.out)
        .concat(get_plasmidfinder_version.out)
        .concat(get_phispy_version.out)
    versions_channel.view()
    
    emit:
        versions_channel

}


workflow upgrade {
    // initializes all databases in specified params.db_dir and sets up/updates
    // the antibiotic resistance gene database
    create_resfinder_db()
    create_amrfinder_db()
    create_bakta_db(create_resfinder_db.out)
    create_pointfinder_db(create_resfinder_db.out)
    create_mobtyper_db(create_amrfinder_db.out)
    create_speciesfinder_db(create_resfinder_db.out)
    create_mlst_db(create_speciesfinder_db.out)
    create_amrdb(create_resfinder_db.out, create_amrfinder_db.out)
    create_plasmidfinder_db(create_mlst_db.out)
    create_phispy_db(create_plasmidfinder_db.out)
}
