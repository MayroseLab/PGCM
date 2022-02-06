"""

"""

from Bio import SeqIO
import pandas as pd
import os

pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = pipeline_dir + "/conda_env"

onstart:
    if not os.path.isdir(config['out_dir']):
        os.mkdir(config['out_dir'])
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

rule all:
    input:
        report = os.path.join(config['out_dir'],'report.txt'),
        analysis_tsv = os.path.join(config['out_dir'],'nonref_analysis.tsv')

rule separate_ref_and_nonref:
    input:
        config['pan_proteome']
    output:
        ref = os.path.join(config['out_dir'], 'ref_proteins.fasta'),
        nonref = os.path.join(config['out_dir'], 'nonref_proteins.fasta')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        with open(output[0],'w') as ref, open(output[1],'w') as nonref:
            for rec in SeqIO.parse(input[0], 'fasta'):
                if rec.seq.endswith('*'):
                    rec.seq = rec.seq[:-1]
                if rec.id.startswith('PanGene'):
                    print(rec.format('fasta'), file=nonref)
                else:
                    print(rec.format('fasta'), file=ref)

rule blat_nonref_vs_ref:
    input:
        ref = os.path.join(config['out_dir'], 'ref_proteins.fasta'),
        nonref = os.path.join(config['out_dir'], 'nonref_proteins.fasta')
    output:
        os.path.join(config['out_dir'], 'nonref_vs_ref.psl')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blat.yml'
    shell:
        """
        blat {input.ref} {input.nonref} {output} -prot
        """

rule blat_nonref_vs_homology_db:
    input:
        nonref = os.path.join(config['out_dir'], 'nonref_proteins.fasta'),
        hom_db = config['homology_proteins_db']
    output:
        os.path.join(config['out_dir'], 'nonref_vs_DB.psl')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blat.yml'
    shell:
        """
        blat {input.hom_db} {input.nonref} {output} -prot
        """

rule blat_top_hit:
    """
    For each BLAT result,
    find the top hit per
    query
    """
    input:
        os.path.join(config['out_dir'], 'nonref_vs_{x}.psl')
    output:
        os.path.join(config['out_dir'], 'nonref_vs_{x}.psl.top')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        colnames = ["match","mismatch","rep._match","Ns","Q_gap_count","Q_gap_bases","T_gap_count","T_gap_bases","strand","Q_name","Q_size","Q_start","Q_end","T_name","T_size","T_start","T_end","block_count","blocK_sizes","q_starts","t_starts"]
        psl_df = pd.read_csv(input[0], sep='\t', skiprows=5, names=colnames)
        best_hits = psl_df.sort_values('match', ascending=False).drop_duplicates('Q_name')
        best_hits.to_csv(output[0], sep='\t', index=False)

rule create_nonref_list:
    """
    Create a plain list of nonref proteins
    """
    input:
        os.path.join(config['out_dir'], 'nonref_proteins.fasta')
    output:
        os.path.join(config['out_dir'], 'nonref_proteins.list')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        grep '>' {input} | tr -d '>' > {output}
        """

rule find_domains:
    input:
        os.path.join(config['out_dir'], 'nonref_proteins.fasta')
    output:
        os.path.join(config['out_dir'], 'nonref_domains.tsv')
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    shell:
        """
        module load java/jdk1.8.25
        module load interproscan/5.32-71
        interproscan.sh -i {input} -t p -dp -pa -appl Pfam,ProDom,SUPERFAMILY,PIRSF --goterms --iprlookup -o {output} -f TSV -cpu {params.ppn}
        """

rule bwa_index:
    input:
        config['pan_transcriptome']
    output:
        config['pan_transcriptome'] + '.bwt'
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa index {input}
        """

rule map_RNA_seq:
    input:
        r1 = config['RNA-seq_R1'],
        r2 = config['RNA-seq_R2'],
        trans = config['pan_transcriptome'],
        ind = config['pan_transcriptome'] + '.bwt'
    output:
        os.path.join(config['out_dir'], 'RNA-seq.sam')
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa mem {input.trans} {input.r1} {input.r2} -t {params.ppn} -o {output}
        """

rule sam_to_filtered_bam:
    input:
        os.path.join(config['out_dir'], 'RNA-seq.sam')
    output:
        os.path.join(config['out_dir'], 'RNA-seq.Q20.sort.bam')
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools view -bh -F 4 -q 20 {input} -@ {params.ppn} | samtools sort - -@ {params.ppn} -o {output}
        """

rule count_reads:
    input:
        os.path.join(config['out_dir'], 'RNA-seq.Q20.sort.bam')
    output:
        os.path.join(config['out_dir'], 'RNA-seq.read_counts.tsv')
    params:
        count_script=os.path.join(pipeline_dir, 'uniq_count.py'),
        queue=config['queue'],
        ppn=2,
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools view {input} | cut -f3 | python {params.count_script} > {output}
        """

rule transcript_lengths:
    input:
        config['pan_transcriptome']
    output:
        os.path.join(config['out_dir'], 'pan_transcriptome.fa.fai')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools faidx {input} --fai-idx {output}
        """

rule analyze_nonref:
    """
    Integrate results into a report
    """
    input:
        nonref_list = os.path.join(config['out_dir'], 'nonref_proteins.list'),
        nonref_vs_ref = os.path.join(config['out_dir'], 'nonref_vs_ref.psl.top'),
        nonref_vs_db = os.path.join(config['out_dir'], 'nonref_vs_DB.psl.top'),
        domains = os.path.join(config['out_dir'], 'nonref_domains.tsv'),
        read_counts = os.path.join(config['out_dir'], 'RNA-seq.read_counts.tsv'),
        trans_len = os.path.join(config['out_dir'], 'pan_transcriptome.fa.fai')
    output:
        report = os.path.join(config['out_dir'],'report.txt'),
        analysis_tsv = os.path.join(config['out_dir'],'nonref_analysis.tsv')
    params:
        identical_to_ref_min_identity = config['identical_to_ref_min_identity'],
        identical_to_ref_length_diff_ratio = config['identical_to_ref_length_diff_ratio'],
        truncated_ref_min_identity = config['truncated_ref_min_identity'],
        truncated_ref_max_subject_cov = config['truncated_ref_max_subject_cov'],
        homology_min_identity = config['homology_min_identity'],
        min_rpkm = config['min_RPKM'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        nonref_list = pd.read_csv(input.nonref_list, names=['gene'], squeeze=True)
        nonref_list.index = nonref_list.values

        nonref_vs_ref_df = pd.read_csv(input.nonref_vs_ref, sep='\t')
        identical_to_ref = nonref_vs_ref_df.query('match/Q_size >= {} & Q_size/T_size >= {} & Q_size/T_size <= {}'.format(params.identical_to_ref_min_identity, 1-params.identical_to_ref_length_diff_ratio,1+params.identical_to_ref_length_diff_ratio))['Q_name']
        identical_to_ref.index = identical_to_ref.values
        truncated = nonref_vs_ref_df.query('match/Q_size >= {} & Q_size/T_size < {}'.format(params.truncated_ref_min_identity,params.truncated_ref_max_subject_cov))['Q_name']
        truncated.index = truncated.values

        nonref_vs_db_df = pd.read_csv(input.nonref_vs_db, sep='\t')
        has_homologs = nonref_vs_db_df.query('match/Q_size >= {} & match/T_size >= {}'.format(params.homology_min_identity,params.homology_min_identity))['Q_name']
        has_homologs.index = has_homologs.values

        colnames = ['Protein_accession','Sequence_MD5','Sequence_length','Analysis','Signature_accession','Signature_description','Start_location','Stop_location','Score','Status','Date','InterPro_annotations_accession','InterPro_annotations_description','GO_annotations','Pathways_annotations']
        ips_domains_df = pd.read_csv(input.domains, sep='\t', names=colnames)
        has_domain = ips_domains_df.query('Score < 0.0001')['Protein_accession'].drop_duplicates()
        has_domain.index = has_domain.values

        trans_len_df = pd.read_csv(input.trans_len, sep='\t', usecols=[0,1], index_col=0, names=['gene','length'], squeeze=True)
        read_count_df = pd.read_csv(input.read_counts, sep='\t', index_col=0, names=['gene','count'], squeeze=True)
        gene_rpkm = pd.concat([trans_len_df, read_count_df], axis=1).fillna(0)
        total_reads = sum(gene_rpkm['count'])
        gene_rpkm['RPKM'] = 10**9 * gene_rpkm['count'] / (gene_rpkm['length'] * total_reads)
        min_rpkm = params.min_rpkm
        nonref_expressed = gene_rpkm.loc[nonref_list].query('RPKM >= @min_rpkm')['RPKM']

        nonref_analysis_df = pd.concat([nonref_list, identical_to_ref, truncated, has_homologs, has_domain, nonref_expressed], axis=1, sort=True)
        nonref_analysis_df.columns = ['nonref_gene','identical_to_ref','truncated_ref','has_homologs','has_domain', 'expressed']
        nonref_analysis_df = nonref_analysis_df.applymap(lambda x: 0 if pd.isna(x) else 1)
        nonref_analysis_df['nonref_gene'] = nonref_analysis_df.index
        nonref_analysis_df.to_csv(output.analysis_tsv, sep='\t', index=False)

        total_nonref = nonref_analysis_df.shape[0]
        with open(output.report,'w') as fo:
            print("Total nonreference: {}".format(total_nonref), file=fo)
            for s in ['identical_to_ref', 'truncated_ref', 'has_homologs', 'has_domain', 'expressed']:
                n = nonref_analysis_df.query('{} == 1'.format(s)).shape[0]
                p = n/total_nonref*100
                print("{}: {} ({}%)".format(s,n,p), file=fo)

            reliable = nonref_analysis_df.query('identical_to_ref == 0 & truncated_ref == 0 & (has_homologs == 1 | expressed == 1)').shape[0]
            reliable_perc = reliable/total_nonref*100
            print("Reliable nonreference: {} ({}%)".format(reliable, reliable_perc), file=fo)
