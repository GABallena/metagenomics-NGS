
rule multiple_alignment:
    input:
        args=f"{DIRS['clustering']}/arg_clusters.txt",  # Required output from preprocessing
        card_db=f"{DIRS['blast']}/translated_card.faa"  # Required output from preprocessing
    output:
        aligned=f"{DIRS['phylogeny']}/args_aligned.fasta"
    conda: envs["muscle"]
    shell:
        """
        mkdir -p "{DIRS[phylogeny]}"
        muscle -in "{input.card_db}" -out "{output.aligned}"
        """

rule phylogenetic_analysis:
    input:
        alignment=f"{DIRS['phylogeny']}/args_aligned.fasta"
    output:
        tree=f"{DIRS['phylogeny']}/arg_tree.nwk"
    conda: envs["fasttree"]
    shell:
        """
        mkdir -p "{DIRS[phylogeny]}"
        megacc -a infer_ML_amino_acid.mao \
            -d "{input.alignment}" \
            -o "{output.tree}" \
            -n "{output.tree}" \
            --model JTT \
            --partial-deletion 95 \
            --bootstrap 50
        """

rule visualize_phylogeny:
    input:
        tree=f"{DIRS['phylogeny']}/arg_tree.nwk",
        clusters=f"{DIRS['clustering']}/arg_clusters.txt"
    output:
        plot=f"{DIRS['phylogeny']}/arg_phylogeny_heatmap.pdf"
    script:
        "scripts/plot_phylogeny.R"