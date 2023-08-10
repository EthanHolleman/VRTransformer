
rule permute_sequences:
    conda:
        '../envs/py.yml'
    input:
        fasta=lambda w: seq_paths[w.plasmid]
    output:
        expand(
            'output/permutedSeqs/{plasmid}/perm-{perm_num}.fa',
            perm_num=range(TOTAL_PERMS), allow_missing=True
            )
    script:'../scripts/permute_sequence.py'
    

rule run_dr_t:
    conda:
        '../envs/DrT.yml'
    input:
        'output/permutedSeqs/{plasmid}/perm-{perm_num}.fa'
    output:
        directory('output/DrTransformer/{plasmid}/perm-{perm_num}')
    params:
        name=lambda w: f'{w.plasmid}-perm-{w.perm_num}'
    shell:'''
    cat {input} | DrTransformer --name {params.name} --outdir {output} --logfile
    '''


rule run_dr_t_og_seqs:
    conda:
        '../envs/DrT.yml'
    input:
        lambda w: seq_paths[w.plasmid]
    output:
        directory('output/DrTransformerOGSeqs/{plasmid}/')
    shell:'''
    cat {input} | DrTransformer --name {plasmid} --outdir {output} --logfile
    '''




