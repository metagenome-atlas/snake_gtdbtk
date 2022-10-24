import os



GTDB_DATA_URL = "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz"
DBDIR = os.path.realpath(config["database_dir"])
GTDBTK_DATA_PATH = os.path.join(DBDIR, "GTDB_V07")






localrules:
    download_gtdb,


rule download_gtdb:
    output:
        temp(f"{GTDBTK_DATA_PATH}/gtdb_data.tar.gz"),
    conda:
        "../envs/gtdbtk.yaml"
    threads: 1
    resources:
        time=int(config.get("runtime", {"long": 10})["long"]),
    log:
        "logs/download/gtdbtk.log",
    shell:
        " wget {GTDB_DATA_URL} -O {output} &> {log} "


rule extract_gtdb:
    input:
        rules.download_gtdb.output,
    output:
        touch(os.path.join(GTDBTK_DATA_PATH, "downloaded_success")),
    conda:
        "../envs/gtdbtk.yaml"
    threads: 1
    resources:
        time=int(config.get("runtime", {"long": 10})["long"]),
    log:
        "logs/download/gtdbtk_untar.log",
    shell:
        'tar -xzvf {input} -C "{GTDBTK_DATA_PATH}" --strip 1 2> {log}; '
        'echo "Set the GTDBTK_DATA_PATH environment variable to {GTDBTK_DATA_PATH} " >> {log}; '
        "conda env config vars set GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} "


rule identify:
    input:
        flag=rules.extract_gtdb.output,
        dir= config["genome_dir"]
    output:
        directory(f"{gtdb_dir}/identify"),
    threads: config["threads"]
    resources:
        mem_mb = 1000* config["large_mem"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/identify.txt",
        f"{gtdb_dir}/gtdbtk.log",
        f"{gtdb_dir}/gtdbtk.warnings.log",
    params:
        extension="fasta",
        outdir = lambda wc, output: Path(output[0]).parent
    shell:
        "gtdbtk identify "
        "--genome_dir {input.dir} "
        " --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


checkpoint align:
    input:
        f"{gtdb_dir}/identify",
    output:
        directory(f"{gtdb_dir}/align"),
    threads: config["threads"]
    resources:
        mem_mb = 1000* config["large_mem"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/align.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
    shell:
        "gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=config["genome_dir"],
    output:
        directory(f"{gtdb_dir}/classify"),
    threads: config["threads"]  #pplacer needs much memory for not many threads
    resources:
        mem_mb = 1000* config["large_mem"],
        #time=config["runtime"]["long"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fasta",
    shell:
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        " --tmpdir {resources.tmpdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


rule combine_taxonomy:
    input:
        folder=f"{gtdb_dir}/classify",
    output:
        combined=f"{gtdb_dir}/gtdbtk.combined.summary.tsv",
        taxonomy="taxonomy/gtdb_taxonomy.tsv",
    log:
        "logs/taxonomy/gtdbtk/combine.txt",
    script:
        "../scripts/combine_taxonomy.py"


rule build_tree:
    input:
        f"{gtdb_dir}/align/{{msa}}.user_msa.fasta.gz",
    output:
        temp("tree/{msa}.unrooted.tree"),
    log:
        "tree/{msa}.log",
        "tree/{msa}.warnings.log"
    threads: max(config["threads"], 3)
    params:
        outdir = lambda wc, output: Path(output[0]).parent
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "gtdbtk infer --msa_file {input} "
        " --out_dir {params.outdir} "
        " --prefix {wildcards.msa} "
        " --cpus {threads} "
        "--tmpdir {resources.tmpdir} "
        

localrules:
    root_tree,


rule root_tree:
    input:
        tree=rules.build_tree.output,
    wildcard_constraints:
        msa="((?!unrooted).)*",
    output:
        tree="tree/{msa}.nwk",
    conda:
        "../envs/tree.yaml"
    threads: 1
    resources:
        mem_mb = 1000* config["simplejob_mem"],
        #time=config["runtime"]["simplejob"],
    log:
        "logs/tree/root_tree_{msa}.log",
    script:
        "../scripts/root_tree.py"


