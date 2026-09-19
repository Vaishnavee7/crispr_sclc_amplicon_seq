

def alignment(SAMPLE, PATH, REFERENCE, FASTQ1, FASTQ2, BWA, OUTPUT, SAMTOOLS, THREADS_SORT, fil, real=True, CENTER = 'c1', LIBRARY = 'l1', LIBRARY_PREPARATION = 'p1', platform = 'Illumina', interleaved= False, options = '', THREADS_BWA = 16):
    import gzip
    with gzip.open(FASTQ1, 'rt', encoding='utf-8') as f:
        for line in f:
            instrument, run_id, flowcell_id, lane, tile, x, y = line[1:].strip().split('\t')[0].split(' ')[0].split(':')
            break

    READGROUP = f'@RG\\tID:{instrument}.{flowcell_id}.{lane}\\tPL:{platform}\\tPU:{instrument}.{flowcell_id}\\tLB:{SAMPLE}.{LIBRARY}.{LIBRARY_PREPARATION}\\tDS:na\\tSM:{SAMPLE}\\tCN:{CENTER}'


    if real:
        if FASTQ2:

            command = f"{BWA} mem -M -v3 -t {THREADS_BWA} {options} -R '{READGROUP}' '{REFERENCE}' '{FASTQ1}' '{FASTQ2}' | {SAMTOOLS} sort -@ {THREADS_SORT} -m 1G -o '{OUTPUT}' -" 
        elif interleaved:
            pass
        else:
            command = f"{BWA} mem -M -v3 -t {THREADS_BWA} {options} -R '{READGROUP}' '{REFERENCE}' '{FASTQ1}' | {SAMTOOLS} sort -@ {THREADS_SORT} -m 1G -o '{OUTPUT}' -"

        fil.caller(command)
        command = "%s index '%s'" %(SAMTOOLS, OUTPUT)
        fil.caller(command)
    return OUTPUT

def AnnovarAnnotation(VCF, ANNOVAR, ANNOVAR_DB, ANNOVAR_REF, fil, THREADS = 1, simple = False):
    annotation_output = VCF.strip('.vcf') + '_annotated'

    command = '''%s --thread %s -buildver %s -out %s -remove -protocol refGene,cosmic81_coding,cosmic81_noncoding,avsnp147,dbscsnv11,exac03,kaviar_20150923,gnomad_genome,gnomad_exome,hrcr1,gme,esp6500siv2_all -operation g,f,f,f,f,f,f,f,f,f,f,f -arg '-splicing_threshold 25',,,,,,,,,,, -vcfinput %s %s''' %(ANNOVAR, THREADS, ANNOVAR_REF, annotation_output, VCF, ANNOVAR_DB)
#        command = '''%s --thread %s -buildver hg38 -out %s -remove -protocol refGene,cosmic81_coding,cosmic81_noncoding,avsnp147,dbscsnv11,exac03,kaviar_20150923,gnomad_genome,gnomad_exome,hrcr1 -operation g,f,f,f,f,f,f,f,f,f -arg '-splicing_threshold 25',,,,,,,,, -polish -vcfinput %s %s''' %(ANNOVAR, THREADS, annotation_output, VCF, ANNODB)
    
    if simple:
        command = '''%s --thread %s -buildver %s -out %s -remove -protocol refGene -operation g -arg '-splicing_threshold 25' -vcfinput %s %s''' %(ANNOVAR, THREADS, ANNOVAR_REF, annotation_output, VCF, ANNOVAR_DB)

    fil.caller(command)

def CallCopyRatioSegments(GATK, JAVA, CN_SEGMENTATION, SOMATIC_CNV, fil):
    command = "'{}' -jar -Xmx20g '{}' CallCopyRatioSegments --input '{}' --output '{}' ".format(JAVA, GATK, CN_SEGMENTATION, SOMATIC_CNV)
    fil.caller(command)
    
def cleanUp(NAME, PATH, fil):
    command = 'rm %s%s/%s_bwa.bam' %(PATH, NAME, NAME)
    fil.caller(command)
    command = 'rm %s%s/%s_bwa.bam.bai' %(PATH, NAME, NAME)
    fil.caller(command)
    command = 'rm %s%s/%s_bwa_dedupped.bam' %(PATH, NAME, NAME)
    fil.caller(command)
    command = 'rm %s%s/%s_bwa_dedupped.bam.bai' %(PATH, NAME, NAME)
    fil.caller(command)    
    command = 'rm %s%s/%s_bwa_dedupped_realigned.bam' %(PATH, NAME, NAME)
    fil.caller(command)    
    command = 'rm %s%s/%s_bwa_dedupped_realigned.bai' %(PATH, NAME, NAME)
    fil.caller(command)
    command = 'rm %s%s/%s_bwa_dedupped_realigned_leftaligned.bam' %(PATH, NAME, NAME)
    fil.caller(command)
    command = 'rm %s%s/%s_bwa_dedupped_realigned_leftaligned.bai' %(PATH, NAME, NAME)
    fil.caller(command)    

def CollectAllelicCounts(GATK, JAVA, REFERENCE, preprocessed_intervals, BAM, OUTPUT_TSV, fil):
    command = "'{}' -jar -Xmx20g '{}' CollectAllelicCounts -R '{}' -I '{}' -O '{}' -L '{}' ".format(JAVA, GATK, REFERENCE, BAM, OUTPUT_TSV, preprocessed_intervals)
    fil.caller(command)
    
def CollectReadCounts(GATK, JAVA, REFERENCE, preprocessed_intervals, BAM, OUTPUT_TSV, fil):
    command = "'{}' -jar -Xmx6g '{}' CollectReadCounts -R '{}' -I '{}' -O '{}' -L '{}' --interval-merging-rule OVERLAPPING_ONLY --format TSV".format(JAVA, GATK, REFERENCE, BAM, OUTPUT_TSV, preprocessed_intervals)
    fil.caller(command) 
    
def CreateSomaticPanelOfNormals(GATK, REFERENCE, DB_PATH, PON_VCF, DBSNP, JAVA, fil):
    command = "'{}' -jar -Xmx20g '{}' CreateSomaticPanelOfNormals  -V gendb://{} -O '{}' --germline-resource '{}' -R '{}' ".format(JAVA, GATK, DB_PATH, PON_VCF, DBSNP, REFERENCE)
    fil.caller(command)
    
def CreateReadCountPanelOfNormals(GATK, REFERENCE, PON_COVERAGE, TSVs_FOR_COVERAGE_PON, interval_filter_table, JAVA, fil): 
    import pandas as pd
    import numpy as np
    PON_inputs = ''
    before_filtering = pd.concat([pd.read_csv(tsv, comment = '@', sep = '\t', dtype={'CONTIG':str, 'START':int, 'END':int, 'COUNT':int}).set_index(['CONTIG', 'START', 'END']) for tsv in TSVs_FOR_COVERAGE_PON], axis=1)
    only_autosomes = before_filtering.loc[before_filtering.index.get_level_values(0).intersection(np.arange(30).astype(str)).unique()]
    median_50 = only_autosomes.median()>=50
    without_low_coverage_samples = only_autosomes.T[median_50].T
    without_low_coverage = without_low_coverage_samples[(without_low_coverage_samples.min(axis=1)>=10) & (without_low_coverage_samples.median(axis=1)>=50)]
    without_low_coverage.to_csv(interval_filter_table)    
    for TSV in TSVs_FOR_COVERAGE_PON:
        filtered_TSV = '.'.join(TSV.split('.')[:-1] + ['filtered', 'tsv'])
        filterReadCountTable(TSV, filtered_TSV, interval_filter_table)
        PON_inputs += "-I '{}' ".format(filtered_TSV)
    
    command = "'{}' -jar -Xmx20g '{}' CreateReadCountPanelOfNormals  --minimum-interval-median-percentile 5.0 -O '{}' {}".format(JAVA, GATK, PON_COVERAGE, PON_inputs)
    fil.caller(command)

def dedup(NAME, PATH, REFERENCE, BAM_LIST, OUTPUT, SAMTOOLS, PICARD, JAVA, fil, real=True):
    import os

    DUPMETRICS = '%s%s/metrics/%s.dedup_metrics' %(PATH, NAME, NAME)
    metrics_folder = '%s%s/metrics/' %(PATH, NAME)
    os.makedirs(metrics_folder, exist_ok = True)

    TEMP_DIR = '%s%s/dedup_temp/' %(PATH, NAME)
    os.makedirs(TEMP_DIR, exist_ok = True)
    INPUTS = ''
    for input_bam in BAM_LIST:
        INPUTS += " INPUT= '{}'".format(input_bam)
    command = r"'{}' -jar -Xmx22g '{}' MarkDuplicatesWithMateCigar TMP_DIR='{}' MINIMUM_DISTANCE=200 {} OUTPUT='{}' METRICS_FILE='{}' ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT".format(JAVA, PICARD, TEMP_DIR, INPUTS, OUTPUT, DUPMETRICS)
    if real:
        fil.caller(command)
        command = "%s index '%s'" %(SAMTOOLS, OUTPUT)
        fil.caller(command)
    return OUTPUT

def DenoiseReadCountsWithPON(GATK, REFERENCE, PON_COVERAGE, TSV, interval_filter_table, JAVA, fil):
    filtered_TSV = '.'.join(TSV.split('.')[:-1] + ['filtered', 'tsv'])
    standardized_TSV = '.'.join(TSV.split('.')[:-1] + ['filtered', 'standardized','tsv'])
    denoised_TSV = '.'.join(TSV.split('.')[:-1] + ['filtered', 'denoised','tsv'])
    filterReadCountTable(TSV, filtered_TSV, interval_filter_table)
    command = "'{}' -jar -Xmx20g '{}' DenoiseReadCounts --count-panel-of-normals '{}' -I '{}' --standardized-copy-ratios '{}' --denoised-copy-ratios '{}'".format(JAVA, GATK, PON_COVERAGE, filtered_TSV, standardized_TSV, denoised_TSV)
    fil.caller(command)    

def FilterMutectCalls(GATK, REFERENCE, VCF, FILTERED_VCF, JAVA, fil):
    command = "'{}' -jar -Xmx20g '{}' FilterMutectCalls -R '{}' -V '{}' -O '{}' --create-output-variant-md5".format(JAVA, GATK, REFERENCE, VCF, FILTERED_VCF)
    fil.caller(command)  


def filterReadCountTable(TSV, filtered_TSV, interval_filter_table):
    import pandas as pd
    keep = set(pd.read_csv(interval_filter_table, dtype=str).set_index(['CONTIG', 'START', 'END']).index)
    with open(TSV) as f, open(filtered_TSV, 'wt') as w:
        for line in f:
            if line[0]=='@':
                w.write(line)
            else:
                 ls = line.split('\t')
                 if ls[0] == 'CONTIG' or tuple(ls[:3]) in keep:
                     w.write(line)
            


def genomicsDBImport(GATK, GVCF_LIST, DB_PATH, TMP_DIR, JAVA, REFERENCE, fil, L=''):
    command = "rm -r '{}'".format(TMP_DIR)
    fil.caller(command)
    command = "mkdir '{}'".format(TMP_DIR)
    fil.caller(command)
    command = "rm -r '{}'".format(DB_PATH)
    fil.caller(command)    
    GVCFs = ''
    for GVCF in GVCF_LIST:
        GVCFs += " -V '{}'".format(GVCF)
    if not L:
        L = ""
        with open(REFERENCE + '.fai') as f:
            for line in f:
                chromosome = line.split()[0]
                L += ' -L {}'.format(chromosome)
    command = "'{}' -jar -Xmx16g '{}' GenomicsDBImport {}  --merge-contigs-into-num-partitions 30 --genomicsdb-workspace-path '{}' --tmp-dir '{}' {}".format(JAVA, GATK, GVCFs, DB_PATH, TMP_DIR, L)
    fil.caller(command)
    
def genotypeGVCFs(GATK, REFERENCE, PATH, DB_PATH, TMP_DIR, JAVA, fil, label = ''):
    COHORT_VCF = '{}cohort{}.vcf'.format(PATH, label)
    command = "'{}' -jar -Xmx16g '{}' GenotypeGVCFs -R '{}' -V gendb://{} --tmp-dir '{}' -O '{}' --create-output-bam-md5 --create-output-variant-md5".format(JAVA, GATK, REFERENCE, DB_PATH, TMP_DIR, COHORT_VCF)
    fil.caller(command)
    

def haplotypeCaller(NAME, PATH, REFERENCE, BAM, SAMTOOLS, GATK, JAVA, fil, L='', real=True):

    GVCF = '{}{}/{}.g.vcf'.format(PATH, NAME, NAME)
    command = f"'{JAVA}' -jar -Xmx16g '{GATK}' HaplotypeCaller -R '{REFERENCE}' -I '{BAM}' --emit-ref-confidence GVCF -O '{GVCF}' {L}--create-output-bam-md5 --create-output-variant-md5"
    if real:
        fil.caller(command)
    return GVCF

def indexReference(SAMTOOLS, PICARD, REFERENCE, fil):
    ###ADD BWA INDEX
    command = '{} faidx {}'.format(SAMTOOLS, REFERENCE)
    fil.caller(command)
    command = "/usr/bin/java -jar -Xmx16g '{}' CreateSequenceDictionary R='{}' O='{}'".format(PICARD, REFERENCE, REFERENCE.replace('.fa', '.dict'))
    fil.caller(command)   
   
def leftAlignVariants(VCF, REFERENCE, GATK, JAVA, fil):
    OUTPUT_VCF = VCF.replace('.vcf', '_leftAligned.vcf')
    command = r"'{}' -jar -Xmx8g '{}' LeftAlignAndTrimVariants -R '{}' -O '{}' -V '{}' --dont-trim-alleles --keep-original-ac --split-multi-allelics --create-output-bam-md5 --create-output-variant-md5".format(JAVA, GATK, REFERENCE, OUTPUT_VCF, VCF)
    fil.caller(command)
    return OUTPUT_VCF

        
def ModelSegments(GATK, JAVA, REFERENCE, ALLELIC_COUNTS, SAMPLE_COVERAGE_TSV_DENOISED, PATH, PREFIX,  fil):
    if ALLELIC_COUNTS:
        command = "'{}' -jar -Xmx20g '{}' ModelSegments -O '{}' --output-prefix {} --allelic-counts '{}' --denoised-copy-ratios '{}' ".format(JAVA, GATK, PATH, PREFIX, ALLELIC_COUNTS, SAMPLE_COVERAGE_TSV_DENOISED)
    else:
        command = "'{}' -jar -Xmx20g '{}' ModelSegments -O '{}' --output-prefix {} --denoised-copy-ratios '{}' ".format(JAVA, GATK, PATH, PREFIX, SAMPLE_COVERAGE_TSV_DENOISED)
    fil.caller(command)
    
def Mutect2Normal(GATK, NAME, REFERENCE, BAM, VCF, intervals, JAVA, fil, THREADS = 4):
    command = "'{}' -jar -Xmx20g '{}' Mutect2 -R '{}' -I '{}' -O '{}' --max-mnp-distance 0 --create-output-variant-md5 --native-pair-hmm-threads {}".format(JAVA, GATK, REFERENCE, BAM, VCF, THREADS)
    if intervals:
        command += ' -L %s' %intervals
    fil.caller(command)
    
def Mutect2TumorOnly(GATK, REFERENCE, BAM, VCF, PON_VCF, intervals, JAVA, fil, THREADS = 4):
    
    #### removed use of dbsnp in the mouse - check for human
    command = "'{}' -jar -Xmx20g '{}' Mutect2 -R '{}' -I '{}' -O '{}' --max-mnp-distance 0 --panel-of-normals '{}' --create-output-variant-md5 --native-pair-hmm-threads {}".format(JAVA, GATK, REFERENCE, BAM, VCF, PON_VCF, THREADS)
    if intervals:
        command += ' -L %s' %intervals
    fil.caller(command)  

def Pileup(GATK, JAVA, REFERENCE, BAM, OUTPUT_PILEUP, intervals, fil):
    command = "'{}' -jar -Xmx6g '{}' Pileup -R '{}' -O '{}' -I '{}' -L '{}' ".format(JAVA, GATK, REFERENCE, OUTPUT_PILEUP, BAM, intervals)
    fil.caller(command) 
    
def prepareAnnovarFiles(REFERENCE, GTF, ANNODB, GENETOPRED, BUILDVERSION, ANNOVARSEQFASTA, fil):
    temp = ANNODB + BUILDVERSION + '_ensGene0.txt'
    extendedPred = ANNODB + BUILDVERSION + '_ensGene.txt'
    MRNA = ANNODB + BUILDVERSION + '_ensGeneMrna.fa'
    command = f'{GENETOPRED} -genePredExt {GTF} {temp}' 
    fil.caller(command)
    command = f'nl {temp} > {extendedPred}'
    fil.caller(command)
    command = f'{ANNOVARSEQFASTA} -format ensGene -seqfile {REFERENCE} --outfile {MRNA} {extendedPred}' 
    fil.caller(command)

def PreprocessIntervals(GATK, JAVA, REFERENCE, intervals, fil):
    intervals_name_list = intervals.split('.')
    output_intervals = '.'.join(intervals_name_list[:-1] + ['preprocessed', 'interval_list'])
    command = "'{}' -jar -Xmx6g '{}' PreprocessIntervals -R '{}' -O '{}' -L '{}' --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY ".format(JAVA, GATK, REFERENCE, output_intervals, intervals)
    fil.caller(command)  
    

def recalibrate(NAME, PATH, REFERENCE, BAM, GATK, THREADS, JAVA, fil, known = [], context = 1, L=''):

    KNOWN_VCFS = ''
    if known:
        for element in known:
            KNOWN_VCFS += ' --known-sites %s ' %element
    else:
        # KNOWN_VCFS += ' --known-sites /media/filippo/HD_1/Databases_1/GATK-bundle/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz '
        # KNOWN_VCFS += ' --known-sites /media/filippo/HD_1/Databases_1/GATK-bundle/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz '    
        KNOWN_VCFS += ' --known-sites /media/filippo/HD_1/Databases_1/dbsnp/00-common_all.vcf.gz' 
        
    RECAL = '%s%s/metrics/%s_recalibration_data.table' %(PATH, NAME, NAME)
    
    if BAM.split('_')[-4:] == 'bwa_dedupped_realigned_leftaligned.bam'.split('_'): 
        RECALIBRATED = '%s%s/%s_recalibrated_context.bam' %(PATH, NAME, NAME)
    else:
        RECALIBRATED = BAM.replace('.bam', '_recalibrated.bam' )
    
    
    # command = r"/usr/bin/java -jar -Xmx20g %s BaseRecalibrator -I %s -R %s -O %s %s -ics %s -mcs %s %s" %(GATK, BAM, REFERENCE, RECAL, KNOWN_VCFS, context, context, L)
    command = r"'{}' -jar -Xmx20g '{}' BaseRecalibrator -I '{}' -R '{}' -O '{}' {} {}".format(JAVA, GATK, BAM, REFERENCE, RECAL, KNOWN_VCFS, L)
    fil.caller(command)
    
    command = r"'{}' -jar -Xmx20g '{}' ApplyBQSR -I '{}' -R '{}' -O '{}' --bqsr-recal-file '{}' {} --create-output-bam-md5".format(JAVA, GATK, BAM, REFERENCE, RECALIBRATED, RECAL, L)
    fil.caller(command)

    return RECALIBRATED
    
