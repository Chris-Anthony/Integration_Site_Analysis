IDS, = glob_wildcards("input/{id}.fastq.gz")
BUILD = "bowtie2/build_mm38/GRCm38"

rule all:
  input:  
    expand("fastqc/{sample}_fastqc.html", sample=IDS),
    expand("homer/{sample}.txt", sample=IDS)

rule FastQC:
  input:
    "input/{sample}.fastq.gz"
  output:
    "fastqc/{sample}_fastqc.html"
  threads: 6
  shell:
    "fastqc -o ./fastqc "
    "--noextract -f fastq "
    "-t {threads} "
    "{input}"

rule cutadapt_1:
  input:
    "input/{sample}.fastq.gz"
  output:
    "cutadapt/{sample}_trim1.fastq.gz"
  threads: 4
  shell:
    "cutadapt --cores={threads} "
    "-g \"GTAAACTTCCGACTTCAACTG;min_overlap=8\" "
    "--minimum-length 20 "
    "-o {output} "
    "{input}"

rule bowtie_1:
  input:
    "cutadapt/{sample}_trim1.fastq.gz"
  output:
    "bowtie2/{sample}_trim1.sam"
  params: BUILD
  threads: 8
  shell:
     "bowtie2 -p {threads} "
     "-x {params} "
     "-U {input} "
     "-S {output}"

rule unmapped_1:
  input:
    "bowtie2/{sample}_trim1.sam"
  output:
    "cutadapt/{sample}_trim1_unmapped.fastq"
  threads: 1
  shell:
     "cat {input} | grep -v ^@ | awk \'$5==0 || length($10)>75\' | "
     "awk \'{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}\' > {output}"

rule cutadapt_2:
  input:
    "cutadapt/{sample}_trim1_unmapped.fastq"
  output:
    "cutadapt/{sample}_trim2.fastq.gz"
  threads: 4
  shell:
     "cutadapt --cores={threads} "
     "-a \"TTGAAGTCGGAAGTTTACA$;max_errors=4\" "
     "-g \"^TGTAAACTTCCGACTTCAACTG;max_errors=5\" "
     "--minimum-length 20 "
     "-o {output} "
     "{input}"

rule bowtie_2:
  input:
    "cutadapt/{sample}_trim2.fastq.gz"
  output:
    "bowtie2/{sample}_trim2.sam"
  params: BUILD
  threads: 8
  shell:
     "bowtie2 -p {threads} "
     "-x {params} "
     "-U {input} "
     "-S {output}"

rule unmapped_2:
  input:
    "bowtie2/{sample}_trim2.sam"
  output:
    "cutadapt/{sample}_trim2_unmapped.fastq"
  threads: 1
  shell:
     "cat {input} | grep -v ^@ | awk \'$5==0 || length($10)>75\' | "
     "awk \'{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}\' > {output}"

rule cutadapt_3:
  input:
    "cutadapt/{sample}_trim2_unmapped.fastq"
  output:
    "cutadapt/{sample}_trim3.fastq.gz"
  threads: 4
  shell:
     "cutadapt --cores={threads} "
     "-a \"TTGAAGTCGGAAGTTTACA$;max_errors=6\" "
     "-g \"^TGTAAACTTCCGACTTCAA;max_errors=6\" "
     "--minimum-length 20 "
     "-o {output} "
     "{input}"

rule bowtie_3:
  input:
    "cutadapt/{sample}_trim3.fastq.gz"
  output:
    "bowtie2/{sample}_trim3.sam"
  params: BUILD
  threads: 8
  shell:
     "bowtie2 -p {threads} "
     "-x {params} "
     "-U {input} "
     "-S {output}"

rule combine_mapped:
  input:
    "bowtie2/{sample}_trim1.sam",
    "bowtie2/{sample}_trim2.sam",
    "bowtie2/{sample}_trim3.sam"
  output:
    "bowtie2/{sample}.sam"
  params:
    gid="{sample}"
  threads: 1
  shell:
    """
    cat bowtie2/{params.gid}_trim1.sam | grep ^@ > bowtie2/{params.gid}.sam
    cat bowtie2/{params.gid}_trim1.sam | grep -v ^@ | awk 'length($10)<=75 && ($5 == 3 || $5 == 8 || $5 == 23 || $5 == 24 || $5 == 40 || $5 == 42)' >> bowtie2/{params.gid}.sam
    cat bowtie2/{params.gid}_trim2.sam | grep -v ^@ | awk 'length($10)<=75 && ($5 == 3 || $5 == 8 || $5 == 23 || $5 == 24 || $5 == 40 || $5 == 42)' >> bowtie2/{params.gid}.sam
    cat bowtie2/{params.gid}_trim3.sam | grep -v ^@ | awk '$5 == 3 || $5 == 8 || $5 == 23 || $5 == 24 || $5 == 40 || $5 == 42' >> bowtie2/{params.gid}.sam
    """

rule samtools:
  input:
    "bowtie2/{sample}.sam"
  output:
    "samtools/{sample}.sam",
    "samtools/{sample}.bam"
  params:
    gid="{sample}"
  threads: 2
  shell:
    """
    samtools sort -@ {threads} -o samtools/{params.gid}.sam {input}
    samtools sort -@ {threads} -o samtools/{params.gid}.bam {input}
    """

rule find_start_site:
  input:
    "samtools/{sample}.sam"
  output:
    "bedtools/{sample}_fwd.bed",
    "bedtools/{sample}_rev.bed"
  params:
    gid="{sample}"
  threads: 1
  shell:
    """
    ./scripts/find_start_site_v3.py -i {input} -c 5 -l 5 -o ./bedtools
    """

rule consensus_fasta:
    input:
      "samtools/{sample}.bam"
    output:
      "bcftools/{sample}_consensus.fa"
    params:
      build="mm38.fa",
      gid="{sample}"
    shell:
      """
      samtools mpileup -Ou -f bowtie2/{params.build} {input} | bcftools call -mv -Oz -o bcftools/{params.gid}_calls.vcf.gz
      bcftools norm -f bowtie2/{params.build} bcftools/{params.gid}_calls.vcf.gz -Ob -o bcftools/{params.gid}_calls.norm.bcf
      bcftools filter --IndelGap 5 bcftools/{params.gid}_calls.norm.bcf -Ob -o bcftools/{params.gid}_calls.norm.flt-indels.bcf
      bcftools index bcftools/{params.gid}_calls.norm.flt-indels.bcf
      bcftools consensus -f bowtie2/{params.build} bcftools/{params.gid}_calls.norm.flt-indels.bcf -o {output}
      """

rule recover_start_site:
    input:
      fwd="bedtools/{sample}_fwd.bed",
      rev="bedtools/{sample}_rev.bed",
      consensus="bcftools/{sample}_consensus.fa"
    output:
      "start_site_filter/{sample}_fwd_recovered.bed",
      "start_site_filter/{sample}_fwd_unrecovered.bed",
      "start_site_filter/{sample}_rev_recovered.bed",
      "start_site_filter/{sample}_rev_unrecovered.bed"
    params:
      seq="TA",
      chrom_size="mm10.chrom.sizes",
      gid="{sample}"
    threads: 1
    shell:
      """
      ./scripts/start_site_filter_v4.py -i {input.fwd} -x {params.seq} -s bowtie2/{params.chrom_size} -b {input.consensus} -t 5 -l 67 -o ./start_site_filter
      ./scripts/start_site_filter_v4.py -i {input.rev} -x {params.seq} -s bowtie2/{params.chrom_size} -b {input.consensus} -t 5 -l 67 -o ./start_site_filter
      """

rule homer_annotation:
    input:
      fwd="start_site_filter/{sample}_fwd_recovered.bed",
      rev="start_site_filter/{sample}_fwd_unrecovered.bed"
    output:
      "homer/{sample}.txt"
    params:
      gid="{sample}",
      build="../../homer/data/genomes/mm10"
    threads: 1
    shell:
      """
      awk \'{{if ($3-$2 >= 65 && $3-$2 <= 74) print}}\' {input.fwd} > temp_{params.gid}.bed
      #awk \'{{if ($3-$2 >= 65 && $3-$2 <= 74) print}}\' {input.rev} >> temp_{params.gid}.bed
      cat temp_{params.gid}.bed | bedtools sort > start_site_filter/{params.gid}_filtered.bed
      rm temp_{params.gid}.bed
      annotatePeaks.pl  start_site_filter/{params.gid}_filtered.bed {params.build} > {output}
      """