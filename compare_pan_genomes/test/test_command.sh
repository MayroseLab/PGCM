#!/bin/bash
#PBS -S /bin/bash
#PBS -N compare_pan_genomes
#PBS -r y
#PBS -q itaymaa
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/compare_pan_genomes/test/compare_pan_genomes.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/compare_pan_genomes/test/compare_pan_genomes.OU

source ~/.bashrc
hostname
conda activate snakemake
PATH="/groups/itay_mayrose/liorglic/miniconda3/envs/snakemake/bin:$PATH"

cd /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/compare_pan_genomes/test
snakefile="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/compare_pan_genomes/compare_pan_genomes.snakefile"
qsub_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/pbs_qsub_snakemake_wrapper.py"
job_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/jobscript.sh"
snakemake -s $snakefile --configfile conf.yml --cluster "python $qsub_script" --latency-wait 60 --use-conda -p -j 10 --jobscript "$job_script" >out 2>err
