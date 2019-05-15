configfile:"config.json"

rule get_genome:
    output:
        "output/genomes/{organism}.fasta.gz"
    run:
        org = wildcards["organism"]
        if "genome_URL" in config["organisms"][org].keys():
            url = config["organisms"][org]["genome_URL"]
            cmd = "wget " + url + " -O {output}"
        if "genome_path" in config["organisms"][org].keys():
            path = config["organisms"][org]["genome_path"]
            cmd = "cp " + path + " {output}"
        shell(cmd)
