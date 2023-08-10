

rule collect_energy_values:
    conda:
        "../envs/py.yml"
    input:
        "output/DrTransformer/{plasmid}/perm-{perm_num}",
    output:
        "output/DrTransformer/{plasmid}/perm-{perm_num}/highestOccupancy.energy",
    script:
        "../scripts/primary_struct_extractor.py"


rule collect_structures:
    conda:
        "../envs/py.yml"
    input:
        energy_files=expand(
            "output/DrTransformer/{plasmid}/perm-{perm_num}/highestOccupancy.energy",
            perm_num=range(TOTAL_PERMS),
            allow_missing=True,
        ),
        plasmids=expand(
            "output/permutedSeqs/{plasmid}/perm-{perm_num}.fa",
            perm_num=range(TOTAL_PERMS),
            allow_missing=True,
        ),
    output:
        "output/collectedStructures/{plasmid}/{plasmid}StructureEnergy.tsv",
    script:
        "../scripts/struct_collector.py"


rule collect_energy_og_seqs:
    conda:
        "../envs/py.yml"
    input:
        "output/DrTransformerOGSeqs/{plasmid}/",
    output:
        "output/DrTransformerOGSeqs/{plasmid}/highestOccupancy.energy",
    script:
        "../scripts/primary_struct_extractor.py"



