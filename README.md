# GTDB-tk Snakemake Workflow


Run as follows:
```
dbdir="databases"

genome_dir="genomes"

snakemake --use-conda --conda-prefix "$dbdir/conda_envs" \
 --config database_dir="$dbdir" genome_dir="$genome_dir"

```

where `dbdir` is the path to a (shared) directory to place the GTDB database and conda envs.
`genome_dir` should be the folder containing all genome fastas.


