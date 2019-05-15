configfile:"config.json"

rule get_genome:
    output:
        "output/reference_genomes/{organism}.fasta.gz"
    run:
        org = wildcards["organism"]
        if "genome_URL" in config["reference_genomes"][org].keys():
            url = config["reference_genomes"][org]["genome_URL"]
            cmd = "wget " + url + " -O {output}"
        if "genome_path" in config["reference_genomes"][org].keys():
            path = config["reference_genomes"][org]["genome_path"]
            cmd = "cp " + path + " {output}"
        shell(cmd)

rule decompress_reference_genome:
    input:
        "output/reference_genomes/{organism}.fasta.gz"
    output:
        "output/reference_genomes_decompress/{organism}.fasta"
    shell:
        "gzcat {input} > {output}"

def derive_genomes_input(wildcards):
    sample = wildcards["sample"]
    genome = wildcards["genome"]
    print("Sample=",sample)
    organism = config["metagenomes"][sample][genome]["reference_genome"]
    return ["output/reference_genomes_decompress/" + organism + ".fasta"]

##In this step, we can later add noise (mutation)
rule derive_genomes:
    input:
        unpack(derive_genomes_input)
    output:
        "output/derived_genomes/{sample}/{genome}.fasta"
    shell:
        "cp {input} {output}"

rule simulate_reads:
    input:
        ref="output/derived_genomes/{sample}/{genome}.fasta"
    output:
        forward="output/simulated_reads_separate/{sample}/{genome}_R1.fastq",
        reverse="output/simulated_reads_separate/{sample}/{genome}_R2.fastq"
    run:
        sample = wildcards.sample
        genome = wildcards.genome
        reads = config["metagenomes"][sample][genome]["reads"]
        cmd="randomreads.sh ref={input.ref} paired out1={output.forward} out2={output.reverse} reads=" + reads + " len=150" 
        shell(cmd)

def combined_reads_forward_input(wildcards):
    sample=wildcards.sample
    return ["output/simulated_reads_separate/" + sample + "/" + genome + "_R1.fastq" for genome in config["metagenomes"][sample].keys()]

rule combined_reads_forward:
    input:
        unpack(combined_reads_forward_input)
    output:
        forward="output/simulated_reads/{sample}_R1.fastq"
    shell:
        "cat {input} > {output}"

def combined_reads_reverse_input(wildcards):
    sample=wildcards.sample
    return ["output/simulated_reads_separate/" + sample + "/" + genome + "_R2.fastq" for genome in config["metagenomes"][sample].keys()]

rule combined_reads_reverse:
    input:
        unpack(combined_reads_reverse_input)
    output:
        reverse="output/simulated_reads/{sample}_R2.fastq"
    shell:
        "cat {input} > {output}"


