
sys.path.append(r'/run/user/1000/gvfs/smb-share:server=134.95.194.53,share=onkosig/scripts/Filippo/')
import CHEOPS.scripts.filPack as fil
import CHEOPS.scripts.PipelinePack_Cheops as pipa
import time
import pysam
import seaborn as sns
import pyBigWig as pbw
import gffutils
import pybedtools
from skbio import io, Sequence
import editdistance

prepare_databases = False
split_fastq = False
parse_bam = False
write_vcf = False
annotate = False


phastcons_35 = '/media/filippo/HD_1/Databases_1/UCSC/mm39.phastCons35way.bw'

gtf_db = '/home/filippo/Databases/Reference_Genomes/mouse/ensembl/GRCm39/110/Mus_musculus.GRCm39.110.gtf.gz.db'
gtf = '/home/filippo/Databases/Reference_Genomes/mouse/ensembl/GRCm39/110/Mus_musculus.GRCm39.110.gtf.gz'

REFERENCE = '/home/filippo/Databases/Reference_Genomes/mouse/ensembl/GRCm39/110/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz'
smaller_reference_fasta = '/home/filippo/Databases/Reference_Genomes/mouse/ensembl/GRCm39/20240308_CRISPR_reference.fa'

BWA = '/home/filippo/bioinformatics_tools/bwa-0.7.17/bwa'
SAMTOOLS = '/home/filippo/bioinformatics_tools/samtools-1.13/samtools'
GATK = '/home/filippo/bioinformatics_tools/gatk-4.4.0.0/gatk'
PICARD = '/home/filippo/bioinformatics_tools/picard-2.26.0/picard.jar' 
JAVA= 'java'
BBMERGE = '/home/filippo/bioinformatics_tools/bbmap/bbmerge.sh'

ANNOVAR = '/home/filippo/bioinformatics_tools/annovar_2020Jun8/annotate_variation.pl'
ANNODB = '/home/filippo/bioinformatics_tools/annovar_2020Jun8/mousedb/'
ANNOVARSEQFASTA = '/home/filippo/bioinformatics_tools/annovar_2020Jun8/retrieve_seq_from_fasta.pl'
TABLEANNOVAR = '/home/filippo/bioinformatics_tools/annovar_2020Jun8/table_annovar.pl'

BUILDVERSION = 'GRCm39'
GENETOPRED = '/home/filippo/bioinformatics_tools/gtfToGenePred'
GTF = '/home/filippo/Databases/Reference_Genomes/mouse/ensembl/GRCm39/110/Mus_musculus.GRCm39.110.gtf.gz'

DATE = '20240308'

fastq_folder = '/run/user/1000/gvfs/smb-share:server=134.95.194.53,share=onkosig/New Server_October 2016/Team_Data/Naja_Eckert_Data/Experiments/CrisprGuides/NGS AmpliconSequencing/Data/bastet2.ccg.uni-koeln.de/downloads/NGS_NE01_neckert_A006850319/'

main_folder = '/media/filippo/HD_1/Projects/SCLC/GLUTAMATE/CRISPR_guides/20240308_guides_NGS/'
fastq_folder = '/media/filippo/HD_1/Projects/SCLC/GLUTAMATE/CRISPR_guides/20240308_guides_NGS/NGS_NE03_neckert_A006200382/'


#samples = pd.read_csv(f'{fastq_folder}Sample_Names.tab', sep='\t')

fastqs = {  'RPMC':{  'F': f'{fastq_folder}A006200382_219739_S1_L000_R1_001.fastq.gz',
                            'R': f'{fastq_folder}A006200382_219739_S1_L000_R2_001.fastq.gz'},
}

amplicon_info = pd.read_csv('/media/filippo/HD_1/Projects/SCLC/GLUTAMATE/CRISPR_guides/20240308_guides_NGS/20240308.CRISPR_amplicon.information.csv', index_col = [3,5,6])
amplicon_info.loc[:, 'amplicon_start amplicon_end amplicon_sequence'.split()] = 0

if prepare_databases:
    command = f"{ANNOVAR} -buildver mm39 -downdb -webfrom annovar ensGene {ANNODB}"
    #fil.caller(command)

    command = f'mv {ANNODB}mm39_ensGene.txt {ANNODB}mm39_ensGene0.txt'
    #fil.caller(command)

    with open(f'{ANNODB}mm39_ensGene0.txt') as f, open(f'{ANNODB}mm39_ensGene.txt', 'wt') as w:
        for line in f:
            w.write(line.replace('chr', ''))

    command = f"{ANNOVARSEQFASTA} '{ANNODB}mm39_ensGene.txt' -seqfile {REFERENCE[:-3]} -format ensGene -outfile {ANNODB}mm39_ensGeneMrna.fa"
    #fil.caller(command)

    #pipa.prepareAnnovarFiles(REFERENCE[:-3], GTF, ANNODB, GENETOPRED, BUILDVERSION, ANNOVARSEQFASTA, fil)

    with pysam.FastaFile(REFERENCE) as pyref, open(smaller_reference_fasta, 'wt') as out_f:
        for (chromosome, strand, center), guide_data in amplicon_info.groupby(level = [0,1,2]):
            print(chromosome, strand, center)
            larger_sequence = pyref.fetch(str(chromosome), center-500, center+500)
            positions = []
            for sequence in guide_data.sequence:
                positions.append(max(larger_sequence.find(sequence.upper()), larger_sequence.find(fil.reverseComplement(sequence[10:].upper()))))
            start, end = sorted(list(set(positions)))
            excluded = int(sequence != larger_sequence[start:start+len(sequence)])
            end += len(guide_data.sequence.iloc[excluded])
            smaller_sequence = larger_sequence[start:end]
            out_f.write(f">{guide_data.gene.iloc[0]}_chr{chromosome}_{center-500 + start}_{center-500 + end}\n{smaller_sequence}\n")
            amplicon_info.loc[(chromosome, strand, center), 'amplicon_start'] = int(center-500 + start)
            amplicon_info.loc[(chromosome, strand, center), 'amplicon_end'] = int(center-500 + end)
            amplicon_info.loc[(chromosome, strand, center), 'amplicon_sequence'] = smaller_sequence
    
    amplicon_info['sample_id'] = (amplicon_info['Mouse ID'] + '_' + amplicon_info['Sample'])
    amplicon_info.to_csv('/media/filippo/HD_1/Projects/SCLC/GLUTAMATE/CRISPR_guides/20240308_guides_NGS/20240308.CRISPR_amplicon.information_updated.csv')
    command = f'bwa index {smaller_reference_fasta}'
    fil.caller(command)

    #db = gffutils.create_db(gtf, dbfn=gtf_db, force=True, keep_order=True,merge_strategy='merge', disable_infer_genes=True, disable_infer_transcripts=True, sort_attribute_values=True)

def apply_cigar_and_mask_N(read, reference_sequence):
    modified_seq = []
    read_seq = read.query_sequence
    seq_index = 0
    ref_index = 0
    
    for cigartuple in read.cigartuples:
        operation, length = cigartuple
        if operation == 0:  # M
            match_sequence = read_seq[seq_index:seq_index+length]
            if 'N' in match_sequence:
                ref_match_sequence = reference_sequence[ref_index:ref_index+length]
                match_sequence = ''.join([base if base != 'N' else ref_match_sequence[i] for i, base in enumerate(match_sequence)])
            modified_seq.append(match_sequence)
            seq_index += length
            ref_index += length
        elif operation == 2:  # D
            #modified_seq.append('n' * length)
            ref_index += length
        elif operation == 1:  # I
            modified_seq.append(read_seq[seq_index:seq_index+length])
            seq_index += length
        elif operation == 4:  # S (soft-clipping)
            # Increment seq_index without adding to modified_seq to exclude soft-clipped bases
            seq_index += length
        # Handle other operations as needed

    return ''.join(modified_seq)

def trim(ref, alt):
    if ref == alt:
        return ref, alt, 0, 0
    for i, (ref_base, alt_base) in enumerate(zip(ref, alt)):
        #print (i, ref_base, alt_base, 'forward')
        if ref_base != alt_base:
            break

    i = max(i-1, 0)
    ref, alt = ref[i:], alt[i:]

    for i2, (ref_base, alt_base) in enumerate(zip(ref[::-1], alt[::-1])):
        #print (i2, ref_base, alt_base, 'reverse')
        if ref_base != alt_base:
            break
    
    if i2:
        ref, alt = ref[:-i2], alt[:-i2]

    return ref, alt, i, i2

amplicon_info = pd.read_csv('/media/filippo/HD_1/Projects/SCLC/GLUTAMATE/CRISPR_guides/20240308_guides_NGS/20240308.CRISPR_amplicon.information_updated.csv')

subsamples = amplicon_info[~pd.isna(amplicon_info.sample_id)].sort_values('sample_id').sample_id.unique()
for sample in ['RPMC']:
    if split_fastq:
        FQ1 = fastqs[sample]['F']
        FQ2 = fastqs[sample]['R']
        fq_merged =  f'{main_folder}{sample}_{DATE}_CRISPR_amplicon.merged.fq.gz'
        fq_unmerged_forward = f'{main_folder}{sample}_{DATE}_CRISPR_amplicon.unmerged.forward.fq.gz'
        fq_unmerged_reverse = f'{main_folder}{sample}_{DATE}_CRISPR_amplicon.unmerged.reverse.fq.gz'
        command = f'{BBMERGE} in1={FQ1} in2={FQ2} out={fq_merged} outu={fq_unmerged_forward} outu2={fq_unmerged_reverse}'
        #fil.caller(command)

        barcode_index = amplicon_info[~pd.isna(amplicon_info.sample_id)].set_index('barcode').sample_id.to_dict()
        
        barcode_index_reverse = {fil.reverseComplement(barcode): barcode_index[barcode] for barcode in barcode_index}

        subsample_fastqs = {subsample:open(f'{main_folder}{subsample}_{DATE}_CRISPR_amplicon.fq', 'w') for subsample in barcode_index.values()}


        def read_four_lines_at_a_time(file_path):
            with open(file_path, 'r') as file:
                lines_group = []
                for line in file:
                    lines_group.append(line.rstrip())  # Remove newline characters
                    if len(lines_group) == 4:
                        # Process the group of 4 lines here
                        yield lines_group
                        lines_group = []  # Reset the list for the next group of lines

        found, not_found = 0, 0
        for header, sequence, empty, qualities in read_four_lines_at_a_time(fq_merged.replace('.gz', '')):
            
            subsample = False
            barcode = sequence[:10]
            try:
                subsample = barcode_index[barcode]
                trimmed_quality_scores = qualities[10:]
                trimmed_sequence = sequence[10:]

            except KeyError:
                try:
                    barcode = sequence[-10:]
                    subsample = barcode_index_reverse[barcode]
                    trimmed_quality_scores = qualities[:-10]
                    trimmed_sequence = sequence[:-10]

                except KeyError:
                    not_found += 1
                    continue

            found += 1
            subsample_fastqs[subsample].write(f'{header}\n{trimmed_sequence}\n+\n{trimmed_quality_scores}\n')

            if not found % 1000:
                print (found, found / (found + not_found), end = '\r')


        for subsample in subsample_fastqs:
            subsample_fastqs[subsample].close()
    

for subsample in subsamples:

    print ('\n\nSTARTING', subsample, '\n\n')
    sample_folder = f'{main_folder}{subsample}/'
    os.makedirs(f'{sample_folder}', exist_ok=True)
    sample_fastq = f'{main_folder}{subsample}_{DATE}_CRISPR_amplicon.fq'
    sample_fastq_gz = f'{sample_fastq}.gz'
    bam = f'{sample_folder}{subsample}.{DATE}.CRISPR_amplicon.bam'
    vcf = f'{sample_folder}{subsample}.{DATE}.CRISPR_amplicon.vcf'

    if parse_bam:
        command = f'gzip -k -f "{sample_fastq}"'
        fil.caller(command)
        print('\n\naligment started for',subsample, '\n\n')
        BAM = pipa.alignment(subsample, fastq_folder, smaller_reference_fasta, sample_fastq_gz, False, BWA, bam, SAMTOOLS, 20, fil, real=True, CENTER = 'CCG', LIBRARY = '1', LIBRARY_PREPARATION = 'CRISPR_amplicons', platform = 'Illumina', options = '-w 300 -A 3 -B 12 -O 18 -E 1 -L 150')

        pybam = pysam.AlignmentFile(bam, 'rb')
        pyref = pysam.FastaFile(smaller_reference_fasta)

        variants = {}
        count = 0
        # Iterate over each alignment in BAM file
        for read in pybam.fetch(until_eof=True):
            # Skip if read is unmapped
            chromosome = read.reference_name
            if read.mapping_quality <10:
                continue
            chromosome, start, end = read.reference_name, read.reference_start, read.reference_end
            reference_sequence = pyref.fetch(chromosome, start, end)
            #aligned_sequence = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            aligned_sequence = apply_cigar_and_mask_N(read, reference_sequence)
            if 'N' in aligned_sequence:
                continue

            reference_position = f'{chromosome}:{start}-{end}'
            
            trim_ref, trim_alt, i, i2 = trim(reference_sequence, aligned_sequence)
            start += i
            end -= i2

            variant_key = (chromosome, start, end, trim_ref, trim_alt)
            variants[variant_key] = variants.get(variant_key, 0) + 1
            
            count += 1
            if not count % 100000:
                print(count, flush = True, end = ' ')

        #create a multiindex dataframe with the variants

        variants_df = pd.DataFrame.from_dict(variants, orient='index', columns=['count'])
        variants_df.index = pd.MultiIndex.from_tuples(variants_df.index, names='chromosome start end reference aligned'.split())
        variants_df.to_csv(f'{sample_folder}/{subsample}.{DATE}.CRISPR_amplicon.raw_variants.csv')

        pybam.close()
        pyref.close()

    if write_vcf:
        variants_df = pd.read_csv(f'{sample_folder}/{subsample}.{DATE}.CRISPR_amplicon.raw_variants.csv')
        variants_df['chromosome'] = variants_df.chromosome.str.replace('.0', '')
        variants_df = variants_df.set_index('chromosome start end reference aligned'.split())
        depth_dic =  variants_df.groupby(level=0).sum().to_dict()['count']

        # Open a file to write the output
        with open(vcf, 'wt') as vcf_file, pysam.FastaFile(REFERENCE) as pyref:
            # Write the header lines
            vcf_file.write("##fileformat=VCFv4.2\n")
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for key, count in variants_df.iterrows():
                chromosome, start, end, ref, alt = key
                real_chromosome = chromosome.split('_')[1].replace('chr', '')
                depth = depth_dic[chromosome]
                start = int(start) + int(chromosome.split('_')[2])
                end = int(end) + int(chromosome.split('_')[2])
                assert pyref.fetch(real_chromosome, start, end) == ref
                #write a line of the VCF file, using the 'counts' column of the variants_df as the depth of the variant allele and the variable "depth" for the total depth, using the FORMAT column
                vcf_file.write(f'{real_chromosome}\t{start}\t.\t{ref}\t{alt}\t.\t.\tDP={depth};AD={count.iloc[0]};target={chromosome}\n')


    if annotate:
        annotated = f'{sample_folder}{subsample}.{DATE}.CRISPR_amplicon.annotated.vcf'

        command = f"{TABLEANNOVAR} --thread 2 --buildver mm39 --out '{annotated}' --remove --protocol ensGene --operation g --vcfinput '{vcf}' {ANNODB}"
        fil.caller(command)

        multianno = f'{sample_folder}{subsample}.{DATE}.CRISPR_amplicon.annotated.vcf.mm39_multianno.txt'
        multianno_df = pd.read_csv(multianno, sep='\t')
        multianno_df['target'] = multianno_df.Otherinfo11.str.split('target=').str[1]
        multianno_df['depth'] = multianno_df.Otherinfo11.str.split('DP=').str[1].str.split(';').str[0].astype(int)
        multianno_df['alt_depth'] = multianno_df.Otherinfo11.str.split('AD=').str[1].str.split(';').str[0].astype(int)
        multianno_df['fraction'] = multianno_df.alt_depth / multianno_df.depth
        
        wts = multianno_df['Otherinfo7'] == multianno_df['Otherinfo8']
        multianno_df.loc[multianno_df[wts].index,'ExonicFunc.ensGene'] = 'wildtype'
        not_exonic = (multianno_df['Func.ensGene'] == 'intronic') | (multianno_df['Func.ensGene'] == 'intergenic') & (multianno_df['ExonicFunc.ensGene'] == '.')
        multianno_df.loc[multianno_df[not_exonic].index,'ExonicFunc.ensGene'] = 'wildtype'
        splicing = (multianno_df['Func.ensGene'] == 'splicing') & (multianno_df['ExonicFunc.ensGene'] == '.')
        multianno_df.loc[multianno_df[splicing].index,'ExonicFunc.ensGene'] = 'splicing'

        multianno_df['Gene.ensGene'] = multianno_df['Gene.ensGene'].str.split(';').str[0]
        multianno_df['conservation'] = 'na'

        db = gffutils.FeatureDB(gtf_db, )
        phastcons = pbw.open(phastcons_35)

        genes = {'ENSMUSG00000042453':'Reln',
        'ENSMUSG00000060534':'Dcc',
        'ENSMUSG00000024211':'Grm8'}

        results = []
        for gene in genes:
            gene_info = db[gene]
            for i in db.children(gene_info, featuretype='exon', order_by='start'):
                exon_conservation = phastcons.values('chr'+i.chrom, i.start-100, i.end+100)
                results.append(pd.Series(index = pd.MultiIndex.from_arrays([np.array([i.chrom]*len(exon_conservation)), np.arange(i.start-100, i.end+100)], names=['chr', 'pos']), data=exon_conservation))

        results_df = pd.concat(results, axis=0).sort_index()
        conservations_list = []
        for chromosome, start, end, consequence in zip(multianno_df.Chr, multianno_df.Start, multianno_df.End, multianno_df['ExonicFunc.ensGene']):
            if consequence[:13] == 'nonsynonymous' or consequence[:13] == 'nonframeshift':
                conservation = np.round(results_df.loc[str(chromosome), start:end].max(),2)
                conservations_list.append(conservation)                
            else:
                conservations_list.append('na')

        multianno_df['conservation'] = conservations_list
        multianno_df.loc[:,'fraction target Gene.ensGene ExonicFunc.ensGene conservation'.split()].groupby('target Gene.ensGene ExonicFunc.ensGene conservation'.split()).sum().to_csv(f'{sample_folder}/{subsample}.{DATE}.CRISPR_amplicon.summary.csv')


results = {}

def grouper(sample_results):
    grouped_for_image = []
    for element, conservation in zip(sample_results['ExonicFunc.ensGene'], sample_results['conservation']):
        if 'wildtype' in element:
            grouped_for_image.append( 'wildtype / synonymous')
        elif 'nonframeshift' in element or 'nonsynonymous' in element:
            if conservation != 'na' and float(conservation) ==1:
                grouped_for_image.append('nonsynonymous conserved')
            else:
                grouped_for_image.append( 'nonsynonymous not conserved')
        elif 'frameshift' in element:
            grouped_for_image.append( 'frameshift')
        elif  'splicing' in element or 'stopgain' in element:
            grouped_for_image.append( 'stop / splicing')
        elif 'synonymous' in element:
            grouped_for_image.append( 'wildtype / synonymous')
        else:
            grouped_for_image.append( element)
    return grouped_for_image

for subsample in subsamples:
    sample_folder = f'{main_folder}{subsample}/'
    sample_results = pd.read_csv(f'{sample_folder}/{subsample}.{DATE}.CRISPR_amplicon.summary.csv')
    sample_results['grouped_for_image'] = grouper(sample_results)

    sample_results.loc[:,'gene_symbol'] = sample_results['target'].str.split('_').str[0]
    sample_results.loc[:,'chromosome'] = sample_results['target'].str.split('_').str[1].str.replace('chr', '')
    sample_results.loc[:,'amplicon_start'] = sample_results['target'].str.split('_').str[2].astype(int)
    sample_results.loc[:,'amplicon_end'] = sample_results['target'].str.split('_').str[3].astype(int)
    results[subsample] = sample_results.groupby('chromosome amplicon_start amplicon_end gene_symbol grouped_for_image'.split()).sum().loc[:, 'fraction']
    
    
results = pd.DataFrame(results).reset_index().set_index('gene_symbol amplicon_start'.split()).sort_index()
plot_data = results.set_index('grouped_for_image').drop(['chromosome', 'amplicon_end'], axis=1).loc[:, subsamples].T
rcParams, cm, fontsize = fil.setupForFigure()

categories = 'frameshift, stop / splicing, nonsynonymous conserved, nonsynonymous not conserved, wildtype / synonymous'.split(', ')
colors = [sns.color_palette('rainbow',20)[element] for element in [19,18, 17, 14,3]]

fig, ax = fil.figureOfSize(9,9, 1, 1, 2.5,2)

current = pd.Series(index = subsamples, data = 0).sort_index()
#create stacked bar plot with matplotlib

for category, color in zip(categories, colors):
    category_data = plot_data[category].loc[subsamples]
    print (category, '\n', current, '\n', category_data, '\n\n')

    #x_position has an empty space every three elements
    ax.bar(np.arange(len(subsamples)), category_data.values, bottom=current, label=category, color=color)
    current += category_data.loc[current.index].values

ax.set_ylim(0,1.05)
ax.set_ylabel('fraction of reads', fontsize = fontsize)

#add legend for the bar plot colors
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=2, fontsize = fontsize)
fig.savefig(f'{main_folder}RPMC_Grm8.png', dpi = 600)
