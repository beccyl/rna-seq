/*********************************************************************************
 ** This file defines all the RNASeq pipeline stages.
 **
 ** Author: Rebecca Evans <rebecca.evans@mcri.edu.au>
 ** Last Update: 25/02/2015
 ********************************************************************************/
//VERSION=0.1

PICARD_HOME="/home/rebeccae/dev/downloads/picard-tools-1.125"
GATK_HOME="/home/rebeccae/dev/downloads/GenomeAnalysisTK-3.3-0"

//codeBase = file(bpipe.Config.config.script).parentFile.absolutePath

/**********  Parameters that must be check by the user: ***********************/
/*** Modify them below, or set them when you run bpipe (with the -p option)   */

//readLayout="paired" //change to "single" or single-end reads

//threads=12 //Threads to use when running the pipeline on a single sample. ie. the total threads will be samples*threads
genomeDir="/mnt/storage/shared/genomes/hg38/vep/77_GRCh38/"
genomeFasta="/mnt/storage/shared/genomes/hg38/vep/77_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
commonParams="--runThreadN 12 --outSAMattributes All"
//fusionParams="--chimSegmentMin 15 --chimJunctionOverhangMin 15"
fusionParams=""
//filterParams="--outFilterIntronMotifs RemoveNoncanonicalUnannotated"
filterParams=""
gtfFile="/mnt/storage/shared/genomes/hg38/vep/77_GRCh38/Homo_sapiens.GRCh38.77.gtf"

// Input pattern (see bpipe documentation for how files are grouped and split )
// group on start, split on end. eg. on ReadsA_1.fastq.gz, ReadA_2.fastq.gz
// this would group the read files into pairs like we want.
//fastqInputFormat="%_*.fastq.gz"

//fastaSuffix="fasta"  
//fastaInputFormat="%."+fastaSuffix

/***************** Other configurables *************************************/

//Default output name
//outputName="star-analysis"


/** Stages **/

first_pass = {
    output.dir = "1pass"

    exec """
       mkdir -p $output.dir;
            
       STAR --genomeDir $genomeDir --readFilesIn \
       /home/rebeccae/dev/Ians_MiSeq_RNA/data/20140908/GMTR-14-2_S3_L001_R1_001.fastq.gz \
       /home/rebeccae/dev/Ians_MiSeq_RNA/data/20140908/GMTR-14-2_S3_L001_R2_001.fastq.gz \
       --outFileNamePrefix $output.dir/ \
       --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate $commonParams $fusionParams $filterParams
    """

}

genome_second_pass = {

    output.dir = "genome_2pass"

    exec """
      mkdir -p $output.dir
            
      STAR --runMode genomeGenerate --genomeDir $output.dir --genomeFastaFiles $genomeFasta --runThreadN 12 \
      --outFileNamePrefix $output.dir/ \
      --sjdbFileChrStartEnd ./1pass/SJ.out.tab --sjdbOverhang 100 --sjdbGTFfile $gtfFile
    """

}


second_pass = {

    output.dir = "2pass"

    exec """
      mkdir -p $output.dir
      STAR --genomeDir ./genome_2pass/ --readFilesIn
      /home/rebeccae/dev/Ians_MiSeq_RNA/data/20140908/GMTR-14-2_S3_L001_R1_001.fastq.gz
      /home/rebeccae/dev/Ians_MiSeq_RNA/data/20140908/GMTR-14-2_S3_L001_R2_001.fastq.gz
      --outFileNamePrefix $output.dir/ \
      --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate $commonParams $fusionParams
    """
}


addrg = {

//      #add read groups, sort
    exec """
      java -Xmx6g -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
      I=./2pass/Aligned.sortedByCoord.out.bam O=$output.bam \
      SO=coordinate RGID=948490 RGLB=SRP028554 RGPL=ILLUMINA RGPU=lane RGSM=SRR948490 VALIDATION_STRINGENCY=LENIENT
    """

}

dedupe = {

//      #mark duplicates, and create index
    exec """
      java -Xmx6g -jar $PICARD_HOME/picard.jar MarkDuplicates I=$input.bam O=$output.bam \
      CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
    """
}

splitncigar = {

//      #split'n'trim and reassign mapping qualities
    exec """
      java -Xmx6g -jar $GATK_HOME/GenomeAnalysisTK.jar -T SplitNCigarReads \
      -R $genomeFasta \
      -I $input.bam -o $output.bam \
      -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    """
}

call_variants = {
//      #variant calling
    exec """
      java -Xmx6g -jar $GATK_HOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
      -R $genomeFasta \
      -I $input.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 \
      -o $output.vcf
    """
}

filter_variants = {
//      #variant filtering
    exec """
      java -Xmx6g -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration \
      -R $genomeFasta \
      -V $input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $output.vcf
    """
}

Bpipe.run {
    first_pass + genome_second_pass + second_pass + addrg + dedupe + splitncigar + call_variants + filter_variants
}
