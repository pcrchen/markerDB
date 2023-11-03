
configfile: "config.yaml"

# set api key 
os.environ["ENTREZ_KEY"] = config["ncbi_api"]

outdir = config["out_directory"]

rule all:
    input:
        "{outdir}/db/seqs.fasta".format(outdir = outdir),
        "{outdir}/db/seqs_nr.fasta".format(outdir = outdir),
        "{outdir}/db/seqs_nr.aln".format(outdir = outdir),
        "{outdir}/db/taxonomy.txt".format(outdir = outdir),
        "{outdir}/db/taxonomy_nr.txt".format(outdir = outdir)


rule sequence_search:
    output:
        taxa = "{outdir}/raw/raw_taxa.tsv".format(outdir = outdir),
        seqs = "{outdir}/raw/raw_seqs.fasta".format(outdir = outdir)
    params:
        taxdb_cache = "taxdb"
    conda: 
        "envs/sequence_search.yaml"
    script:
        "scripts/sequence_search.R"

rule cmscan:
    input: 
        rules.sequence_search.output.seqs
    output:
        "{outdir}/raw/cmscan.out".format(outdir = outdir)
    params: 
        models = "cm_model/marker_cm_database.cm"
    conda:
        "envs/cmscan.yaml"
    threads: 
        config["threads"]
    shell:
        """
        cmscan --rfam --cut_ga --nohmmonly --tblout {output} --fmt 2 --cpu {threads} {params.models} {input} > /dev/null 2>&1
        """
        
rule trim_seqs:
    input:
        hits = rules.cmscan.output,
        seqs = rules.sequence_search.output.seqs,
        taxa = rules.sequence_search.output.taxa
    output:
        seqs_final = "{outdir}/db/seqs.fasta".format(outdir = outdir),
        seqs_final_nr = "{outdir}/db/seqs_nr.fasta".format(outdir = outdir),
        taxa_final = "{outdir}/db/taxonomy.txt".format(outdir = outdir),
        taxa_final_nr = "{outdir}/db/taxonomy_nr.txt".format(outdir = outdir)
    conda: 
        "envs/trim_seqs.yaml"
    script:
        "scripts/trim_seqs.R"

rule align:
    input:
        rules.trim_seqs.output.seqs_final_nr
    output:
        "{outdir}/db/seqs_nr.aln".format(outdir = outdir)
    conda:
        "envs/align.yaml"
    threads:
        config["threads"]
    shell:
        """
        mafft --auto --reorder --thread {threads} {input} > {output}
        """

rule write_seqs:
    input:
        seqs = rules.trim_seqs.output.seqs_final,
        taxa = rules.trim_seqs.output.taxa_final,
    output:
        directory("{outdir}/formats/full".format(outdir = outdir))
    conda:
        "envs/write_seqs.yaml"
    script:
        "scripts/write_seqs.R"
        
rule write_seqs_nr:
    input:
        seqs = rules.trim_seqs.output.seqs_final_nr,
        taxa = rules.trim_seqs.output.taxa_final_nr,
        aln = rules.align.output
    output:
        directory("{outdir}/formats/nr".format(outdir = outdir))
    conda:
        "envs/write_seqs.yaml"
    script:
        "scripts/write_seqs.R"
        
rule run_app:
    conda:
        "envs/app.yaml"
    params:
        db_dir = outdir
    script:
        "scripts/app.R"
        



