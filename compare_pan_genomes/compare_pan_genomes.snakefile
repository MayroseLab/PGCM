from Bio import SeqIO
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *
import pandas as pd

#load_info_file
pan_genomes = pd.read_table(config['pan_genomes_info']).set_index("pan_genome_name", drop=False)
assert pan_genomes.shape[0] == 2, "Exactly two pan genomes should be provided"

LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = pipeline_dir + "/conda_env"
annotation_pipeline_dir = os.path.dirname(pipeline_dir) + '/annotation_pipeline'

onstart:
    write_config_file(config)
    if not os.path.isdir(LOGS_DIR):
        os.mkdir(LOGS_DIR)

onsuccess:
    print("%s pipeline finished, no error" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

onerror:
    print("%s pipeline failed" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

#------------------------------------
#                RULES              |
#------------------------------------

localrules: all

pg1 = pan_genomes.index[0]
pg2 = pan_genomes.index[1]

rule all:
    input:
        os.path.join(config["out_dir"], 'report.html'),
        os.path.join(config["out_dir"], 'discrepancies.tsv'),
        os.path.join(config["out_dir"], '{}_trans_vs_{}_assemblies/transcript_mapping.tsv'.format(pg1,pg2)),
        os.path.join(config["out_dir"], '{}_trans_vs_{}_assemblies/transcript_mapping.tsv'.format(pg2,pg1))

def get(wildcards, what, which=None):
    if not which or which == 'PG':
        pg = wildcards.PG
    elif which == 'PG1':
        pg = wildcards.PG1
    elif which == 'PG2':
        pg = wildcards.PG2
    sample_dir = pan_genomes.loc[pg, 'path']
    if what == 'prot':
        return sample_dir + '/all_samples/pan_genome/pan_proteome.fasta'
    elif what == 'trans':
        return sample_dir + '/all_samples/pan_genome/pan_transcripts.fasta'
    elif what == 'pav':
        return sample_dir + '/all_samples/pan_genome/pan_PAV.tsv'
    elif what == 'per_sample':
        return sample_dir + '/per_sample'

rule extract_non_ref_PG:
    """
    Extract non-ref proteins into a new fasta
    """
    input:
        lambda wc: get(wc, what='prot')
    output:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        nonref_records = []
        for rec in SeqIO.parse(input[0], 'fasta'):
            if rec.id.startswith('PanGene'):
                nonref_records.append(rec)
        SeqIO.write(nonref_records, output[0], 'fasta')

rule make_blast_db:
    """
    Create blast DB for non-ref proteins of each pan genome
    """
    input:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta")
    output:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta.phr")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blast.yml'
    shell:
        """
        makeblastdb -dbtype prot -in {input} -input_type fasta
        """

rule blast_non_ref:
    """
    blast non ref proteins from one PG against the other
    """
    input:
        pg1_fasta=os.path.join(config["out_dir"], "{PG1}_nonref.fasta"),
        pg2_fasta=os.path.join(config["out_dir"], "{PG2}_nonref.fasta"),
        pg2_db=os.path.join(config["out_dir"], "{PG2}_nonref.fasta.phr")
    output:
        os.path.join(config["out_dir"], "{PG1}_nonref_vs_{PG2}_nonref.blast6")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/blast.yml'
    shell:
        """
        blastp -query {input.pg1_fasta} -db {input.pg2_fasta} -out {output} -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        """

rule find_matches:
    """
    Find the maximum weight matching between pan genome proteins
    """
    input:
        pg1_fasta=os.path.join(config["out_dir"], "{PG1}_nonref.fasta"),
        pg2_fasta=os.path.join(config["out_dir"], "{PG2}_nonref.fasta"),
        fw=os.path.join(config["out_dir"], "{PG1}_nonref_vs_{PG2}_nonref.blast6"),
        rev=os.path.join(config["out_dir"], "{PG2}_nonref_vs_{PG1}_nonref.blast6"),
    output:
        os.path.join(config["out_dir"], '{PG1}_vs_{PG2}_max_weight_matches.tsv')
    params:
        find_matches_script=os.path.join(pipeline_dir, 'match_non_ref.py'),
        pg1_name=lambda w: w.PG1,
        pg2_name=lambda w: w.PG2,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/match_non_ref.yml'
    shell:
        """
        python {params.find_matches_script} {input.pg1_fasta} {input.pg2_fasta} {input.fw} {input.rev} {output} --set1_name {params.pg1_name} --set2_name {params.pg2_name} --normalize_weight --filter "pident>90 & qcov>0.9 & scov>0.9"
        """

rule make_nonref_lists:
    """
    Create lists of matched/unmatched
    nonref genes
    """
    input:
        matches_tsv=os.path.join(config["out_dir"], '{}_vs_{}_max_weight_matches.tsv'.format(pg1,pg2)),
        pg1_fasta=os.path.join(config["out_dir"], "{}_nonref.fasta".format(pg1)),
        pg2_fasta=os.path.join(config["out_dir"], "{}_nonref.fasta".format(pg2)),
    output:
        pg1_matched=os.path.join(config["out_dir"], '{}_vs_{}_nonref_matched'.format(pg1,pg2)),
        pg1_unmatched=os.path.join(config["out_dir"], '{}_vs_{}_nonref_unmatched'.format(pg1,pg2)),
        pg2_matched=os.path.join(config["out_dir"], '{}_vs_{}_nonref_matched'.format(pg2,pg1)),
        pg2_unmatched=os.path.join(config["out_dir"], '{}_vs_{}_nonref_unmatched'.format(pg2,pg1)),
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        grep -w -E $(cut -f1 {input.matches_tsv} | tr '\\n' '|' | sed '$ s/.$//') {input.pg1_fasta} | tr -d '>' > {output.pg1_matched}
        grep '>' {input.pg1_fasta} | grep -v -w -E $(cut -f1 {input.matches_tsv} | tr '\\n' '|' | sed '$ s/.$//') | tr -d '>' > {output.pg1_unmatched}
        grep -w -E $(cut -f2 {input.matches_tsv} | tr '\\n' '|' | sed '$ s/.$//') {input.pg2_fasta} | tr -d '>' > {output.pg2_matched}
        grep '>' {input.pg2_fasta} | grep -v -w -E $(cut -f2 {input.matches_tsv} | tr '\\n' '|' | sed '$ s/.$//') | tr -d '>' > {output.pg2_unmatched}
        """

rule map_transcripts:
    """
    Map transcripts of one PG to
    all assemblies of the other PG
    """
    input:
        lambda wc: get(wc, what='trans', which='PG1')
    output:
        os.path.join(config["out_dir"], '{PG1}_trans_vs_{PG2}_assemblies/done')
    params:
        pg2_per_sample_dir=lambda wc: get(wc, what='per_sample',which='PG2'),
        out_dir=os.path.join(config["out_dir"], '{PG1}_trans_vs_{PG2}_assemblies'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/minimap2.yml'
    shell:
        """
        for g in `ls -1 {params.pg2_per_sample_dir}/*/RG_assembly_*/ragtag_output/ragtag.scaffolds.fasta`
        do
            sampleName=`echo $g | sed 's@.*/per_sample/\([^/]*\)/.*@\\1@'`
            minimap2 -x splice:hq -uf $g {input} -t {params.ppn} --paf-no-hit > {params.out_dir}/$sampleName.paf
        done
        touch {output}
        """

rule find_transcript_mapping:
    """
    Parse paf files to find
    mappings of each transcript
    """
    input:
        os.path.join(config["out_dir"], '{PG1}_trans_vs_{PG2}_assemblies/done')
    output:
        os.path.join(config["out_dir"], '{PG1}_trans_vs_{PG2}_assemblies/transcript_mapping.tsv')
    params:
        paf_dir=os.path.join(config["out_dir"], '{PG1}_trans_vs_{PG2}_assemblies'),
        find_trans_mapping_script=os.path.join(pipeline_dir, "find_trans_mapping.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/match_non_ref.yml'
    shell:
        """
        python {params.find_trans_mapping_script} {params.paf_dir} {output}
        """

rule create_report_nb:
    """
    Create jupyter notebook of comparison report
    """
    input:
        pg1_pav=pan_genomes.loc[pg1]['path'] + '/all_samples/pan_genome/pan_PAV.tsv',
        pg2_pav=pan_genomes.loc[pg2]['path'] + '/all_samples/pan_genome/pan_PAV.tsv',
        pg1_vs_pg2_matches=os.path.join(config["out_dir"], '{}_vs_{}_max_weight_matches.tsv'.format(pg1,pg2)),
    output:
        os.path.join(config["out_dir"], 'report.ipynb')
    params:
        nb_template=os.path.join(pipeline_dir, 'report_template.ipynb'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed -e 's@<PG1_PAV>@{input.pg1_pav}@' -e 's@<PG2_PAV>@{input.pg2_pav}@' -e 's@<PG1_VS_PG2_NON_REF_MATCHES>@{input.pg1_vs_pg2_matches}@' -e 's@<PG1_NAME>@%s@' -e 's@<PG2_NAME>@%s@' {params.nb_template} > {output} 
        """ %(pg1, pg2)

rule create_report_html:
    """
    Convert notebook to HTML report
    """
    input:
        os.path.join(config["out_dir"], 'report.ipynb')
    output:
        report=os.path.join(config["out_dir"], 'report.html'),
        discrep=os.path.join(config["out_dir"], 'discrepancies.tsv')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        jupyter nbconvert {input} --output {output.report} --no-prompt --no-input --execute --NotebookClient.timeout=-1 --ExecutePreprocessor.timeout=-1 --to html
        """
