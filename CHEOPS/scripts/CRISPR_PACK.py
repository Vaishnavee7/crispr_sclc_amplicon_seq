#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:32:44 2020

@author: filippo
"""
import os
import sys

sys.path.append(
    r'/run/user/1000/gvfs/smb-share:server=134.95.194.53,share=onkosig/scripts/Filippo/'
)
import RAMSES.scripts.filPack as fil
import gc
import bisect
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as stmu


def analyze_guide_distribution(values,
                               value_label,
                               results_label,
                               results_folder,
                               guide_numbers=[3, 4, 5]):

    columns_dic = {
        'sgrna': 'sgrna',
        'sgRNA Target Sequence': 'sgrna',
        'Target Gene Symbol': 'gene_symbol',
        'gene_symbol': 'gene_symbol'
    }

    values.index.name = columns_dic[values.index.name]

    median_by_gene = values.groupby(level=0).median()
    guides_by_gene = values.groupby(level=0).count()

    results_by_gene = pd.DataFrame({
        'n_guides': guides_by_gene,
        f'median_{value_label}': median_by_gene,
        'p_val': np.nan,
        'q_val': np.nan
    }).sort_values(by=f'median_{value_label}')

    newfig = f'{results_folder}/empirical_distributions_{value_label}_{results_label}.png'
    figure, ax = fil.figureOfSize(9, 6, 1, 1, 1, 1)
    ax.set_xlabel(f'median guide {value_label}')
    ax.set_ylabel('millions of random samplings')

    print('\nsimulating distributions', value_label, results_label, flush=True)
    print('\nn of guides per gene', guide_numbers, flush=True)
    distributions = {}
    gc.collect()

    for guide_number in guide_numbers:
        distributions[guide_number] = []
    for i in range(1500):
        if not i % 100:
            print(i, end=' ', flush=True)

        for guide_number in distributions:
            distributions[guide_number] += values.sample(
                frac=1).rolling(guide_number).median()[guide_number -
                                                       1:].tolist()

    print('\nsorting and printing simulated distributions', flush=True)

    for distribution in distributions:
        print(distribution, end=' ', flush=True)
        distributions[distribution] = np.sort(
            distributions[distribution][:100000000])

        #pd.Series(distributions[distribution]).to_hdf(f'/home/filippo/Downloads/temp/distribution_{value_label}_{results_label}_{distribution}.hdf', 'distribution', complevel=9, index = None)

        ax.hist(distributions[distribution],
                bins=250,
                label=f'{distribution} random guides',
                alpha=0.5)
    ax.legend()
    figure.savefig(newfig, dpi=300)
    plt.close('all')

    print('\ncomputing p-values', flush=True)
    for gi, (gene, (distribution, median_value, pval,
                    qval)) in enumerate(results_by_gene.iterrows()):
        if distribution in distributions:
            insertion_equal_less = bisect.bisect(distributions[distribution],
                                                 median_value)
            insertion_less = bisect.bisect_left(distributions[distribution],
                                                median_value)
            results_by_gene.loc[gene, 'p_val'] = (
                min(insertion_equal_less, 100000000 - insertion_less) +
                1) / 50000000
        if not gi % 100:
            print(f'{gi}     ', end='\r', flush=True)

    no_nan_p = results_by_gene.p_val.dropna()
    no_nan_q = stmu.multipletests(no_nan_p, alpha=0.05, method='fdr_bh')[1]
    results_by_gene.q_val = pd.Series(no_nan_q, index=no_nan_p.index).reindex(
        results_by_gene.index)
    results_by_gene.to_csv(
        f'{results_folder}results_{value_label}_{results_label}.csv')


def compareTwoSamples(sample1,
                      sample2,
                      labels,
                      guides,
                      PATH,
                      precomputed=False):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import statsmodels.stats.multitest as stmu
    import bisect
    import gc

    label_name = '_'.join(labels)
    label1, label2 = labels[:2]

    try:
        shift = (sample2 - sample1).to_pandas()
    except AttributeError:
        shift = sample2 - sample1

    guide_grouper = guides.loc[shift.index].reset_index().set_index(
        ['sgrna', 'gene_symbol']).index.unique()
    shift_guides_by_gene = shift.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).count()
    guide_numbers = sorted(set(shift_guides_by_gene.tolist()))

    shift_mean_by_gene = shift.loc[guide_grouper.get_level_values(0)].groupby(
        guide_grouper.get_level_values(1)).mean()
    sample1_mean_by_gene = sample1.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).mean()
    sample2_mean_by_gene = sample2.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).mean()
    results_by_gene = pd.DataFrame({
        'n_guides': shift_guides_by_gene,
        'mean_shift': shift_mean_by_gene,
        label1 + '_mean': sample1_mean_by_gene,
        label2 + '_mean': sample2_mean_by_gene,
        'p_val': 1,
        'q_val': 1
    }).sort_values(by='mean_shift')

    if not precomputed:
        #        newfig = '%scomparisons/empirical_distributions/%s.png' %(PATH, label_name)
        #        fig = plt.figure(figsize=(10,5))
        #        fig.set_tight_layout(False)
        #        ax = fig.add_axes([0.15,0.15,0.7,0.7])
        #        ax.set_xlabel('mean guide shift', size = 16)
        #        ax.set_ylabel('millions of random samplings', size = 16)

        print('\nsimulating distributions', labels, flush=True)
        print('\nn of guides per gene', guide_numbers, flush=True)
        distributions = {}
        gc.collect()
        for guide_number in guide_numbers:
            distributions[guide_number] = []
        for i in range(1500):
            if not i % 100:
                print(i, end=' ', flush=True)
            for guide_number in distributions:
                distributions[guide_number] += shift.sample(
                    frac=1).rolling(guide_number).mean()[guide_number -
                                                         1:].tolist()

        print('\nsorting and printing simulated distributions', flush=True)

        for distribution in distributions:
            print(distribution, end=' ', flush=True)
            distributions[distribution] = np.sort(
                np.array(distributions[distribution][:100000000]))
#                sns.distplot(distributions[distribution], ax=ax, bins = 250, label = '%s random guides' %kdistribution)

        for gene, (distribution, mean_shift, first, second, pval,
                   qval) in results_by_gene.iterrows():
            if distribution in distributions:
                insertion_equal_less = bisect.bisect(
                    distributions[distribution], mean_shift)
                insertion_less = bisect.bisect_left(
                    distributions[distribution], mean_shift)
                results_by_gene.loc[gene, 'p_val'] = (
                    min(insertion_equal_less, 100000000 - insertion_less) +
                    1) / 50000000


#        ax.legend(fontsize =16)
#        fig.savefig(newfig, dpi=300)
#        plt.close('all')

    results_by_gene.q_val = stmu.multipletests(results_by_gene.p_val,
                                               alpha=0.05,
                                               method='fdr_bh')[1]
    results_by_gene.to_csv('%scomparisons/results_two_%s.csv' %
                           (PATH, label_name))


def generateFileLabel(labels):
    name = ''
    for label in labels:
        name += '_' + label
    return name


def generateShuffledGuidesIndex(guide_by_sample, labels, PATH):
    import pandas as pd
    import numpy as np
    guides_ilocs = pd.Series(np.arange(len(guide_by_sample)))
    shuffled = []
    for i in range(1500):
        if not i % 100:
            print(i, end=' ', flush=True)
        shuffled += guides_ilocs.sample(frac=1).tolist()
    label_name = '_' + '_'.join(labels)
    pd.Series(shuffled, name='guide_index').to_hdf(
        PATH + 'guideShufflingIndex' + label_name + '.hdf',
        'shuffling_index',
        complevel=9,
        index=None)


def compareMultipleSamples(guide_by_sample,
                           sample_index_1,
                           sample_index_2,
                           guides,
                           labels,
                           PATH='/HD1/Projects/SCLC/CRISPR_SCREEN/',
                           bin_resolution=100000,
                           repetitions=100000000):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import statsmodels.stats.multitest as stmu
    import bisect
    import gc
    import itertools

    label_name = '_'.join(labels)
    label1, label2 = labels[:2]
    shuffled_guide_index = pd.read_hdf(
        PATH + 'guideShufflingIndex_' + label_name + '.hdf',
        'shuffling_index')[:repetitions + 100]
    sample1, sample2 = guide_by_sample.loc[:, sample_index_1].median(
        axis=1), guide_by_sample.loc[:, sample_index_2].median(axis=1)
    shift = sample2 - sample1
    max_bin = (guide_by_sample.max(axis=1) - guide_by_sample.min(axis=1)).max()
    binner = lambda x: (x * bin_resolution) / max_bin

    all_samples = list(sample_index_1) + list(sample_index_2)
    combination_of_samples = list(
        itertools.combinations(all_samples, len(sample_index_1)))
    shift_guide_combinations = {}
    shuffling_count = -1
    for element in combination_of_samples:
        random_index_2 = sorted([
            subelement for subelement in all_samples
            if subelement not in element
        ])
        random_index_1 = sorted(element)
        random_index_key_list = sorted(tuple([random_index_1, random_index_2]))
        random_index_key = tuple(random_index_key_list[0] +
                                 random_index_key_list[1])
        if random_index_key not in shift_guide_combinations:
            shuffling_count += 1
            shift_guide_combinations[
                random_index_key] = guide_by_sample.loc[:, random_index_1].median(
                    axis=1) - guide_by_sample.loc[:, random_index_2].median(
                        axis=1)

    guide_grouper = guides.loc[shift.index].reset_index().set_index(
        ['sgrna', 'gene_symbol']).index.unique()
    shift_guides_by_gene = shift.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).count()
    guide_numbers = sorted(set(shift_guides_by_gene.tolist()))

    shift_mean_by_gene = shift.loc[guide_grouper.get_level_values(0)].groupby(
        guide_grouper.get_level_values(1)).mean()
    sample1_mean_by_gene = sample1.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).mean()
    sample2_mean_by_gene = sample2.loc[guide_grouper.get_level_values(
        0)].groupby(guide_grouper.get_level_values(1)).mean()

    results_by_gene = pd.DataFrame({
        'n_guides': shift_guides_by_gene,
        'mean_shift': shift_mean_by_gene,
        label1 + '_mean': sample1_mean_by_gene,
        label2 + '_mean': sample2_mean_by_gene,
        'p_val': 1,
        'q_val': 1
    })

    distributions = {
        guide_number: pd.Series(index=np.arange(bin_resolution + 1), data=0)
        for guide_number in guide_numbers
    }

    newfig = '%scomparisons/empirical_distributions/distribution_%s.png' % (
        PATH, label_name)
    fig = plt.figure(figsize=(10, 5))
    fig.set_tight_layout(False)
    ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])
    ax.set_xlabel('mean guide shift', size=16)
    ax.set_ylabel('random samplings', size=16)
    print('\nsimulating distributions', label_name, flush=True)

    for guide_number in guide_numbers:
        print('\tn of guides per gene', guide_number, flush=True)
        for random_index_key in shift_guide_combinations:
            print('\t\t', random_index_key, flush=True)

            gc.collect()
            binned_counts = binner(
                abs(shift_guide_combinations[random_index_key].
                    iloc[shuffled_guide_index].rolling(guide_number).mean()
                    [guide_number - 1:repetitions + guide_number -
                     1])).astype(int).to_frame(0).groupby(0).size()
            distributions[guide_number].iloc[
                binned_counts.index] += binned_counts
    distribution_df = pd.DataFrame(distributions)
    distribution_df.to_csv('%scomparisons/distribution_%s.csv' %
                           (PATH, label_name))

    distribution_pvals = (1 + distribution_df[::-1].cumsum()[::-1]) / (
        len(shift_guide_combinations) * repetitions + 1)
    distribution_pvals.to_csv('%scomparisons/distribution_pvals_%s.csv' %
                              (PATH, label_name))

    print('\nprinting simulated distributions', flush=True)

    for distribution in distributions:
        print(distribution, end=' ', flush=True)
        ax.plot(distributions[distribution],
                label=str(distribution) + ' guides')
    for gene, (distribution, mean_shift, first, second, pval,
               qval) in results_by_gene.iterrows():
        gene_bin = int(binner(abs(mean_shift)))
        results_by_gene.loc[gene,
                            'p_val'] = distribution_pvals.loc[gene_bin,
                                                              distribution]

    ax.legend(fontsize=16)
    fig.savefig(newfig, dpi=300)
    plt.close('all')

    results_by_gene.q_val = stmu.multipletests(results_by_gene.p_val,
                                               alpha=0.05,
                                               method='fdr_bh')[1]
    results_by_gene.to_csv('%scomparisons/results_multiple_%s.csv' %
                           (PATH, label_name))


def getGeneGuides(gene, guides):
    return guides[guides.loc[:, 'gene_symbol'] == gene].index


def processFastq(fastq, lines=None):

    def process(lines=None):
        ks = ['name', 'sequence', 'optional', 'quality']
        return {k: v for k, v in zip(ks, lines)}

    n = 4

    with gz.open(fastq, 'rt') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                yield lines
                lines = []


def read_library_fastq_helper(arguments):
    FASTQ, folder = arguments
    print(FASTQ, flush=True)

    log = '%s_log.txt' % (FASTQ)
    results_by_guide = '%s_df.csv' % (FASTQ)

    sample_info_df = sample_information[sample_information.folder == folder]
    results_dic = {}
    sample_dic = {}
    for ind, row in sample_info_df.iterrows():
        sample = row.sample_ID
        results_dic[sample] = pd.Series(index=all_guides,
                                        dtype=float).fillna(0).to_dict()
        sample_dic[(row.folder, row.BarcodeF)] = sample

    with open(log, 'wt') as w:
        w.write(
            'fastq, total, no_barcode, barcode_f, barcode_r guide_f, guide_r, no_guide, timestamp\n'
            .replace(', ', '\t'))
        barcodes = sample_info_df.BarcodeF.unique()
        total, no_barcode, barcode_f, guide_f, no_guide = 0, 0, 0, 0, 0
        if FASTQ[-3:] != '.gz':
            print('zipping', FASTQ, flush=True)
            command = 'gzip %s' % FASTQ
            fil.caller(command)

        fq_reader = processFastq(FASTQ)
        fastq_name = FASTQ.split('/')[-1].split(
            '_')[0] + '_na_VectorBuilder_ut'
        for record in fq_reader:
            total += 1
            barcode_found = False
            sequence = record[1]

            if folder == 'VectorBuilder':
                guide = fil.reverseComplement(sequence[115:135])
                sample = fastq_name
                if guide in all_guides_set:
                    results_dic[sample][guide] += 1
                    guide_f += 1
                else:
                    no_guide += 1
            else:
                for barcode in barcodes:
                    match = re.search(barcode, sequence)
                    if match:
                        barcode_found = True
                        barcode_f += 1
                        sample = sample_dic[(folder, barcode)]
                        guide = sequence[match.start() + 32:match.start() + 52]
                        if guide in all_guides_set:
                            results_dic[sample][guide] += 1
                            guide_f += 1
                        else:
                            no_guide += 1
                        break

                if not barcode_found:
                    no_barcode += 1
            if not total % 2000000:
                print(FASTQ,
                      '\n\t',
                      total,
                      no_barcode,
                      barcode_f,
                      guide_f,
                      no_guide,
                      flush=True)
        newline = '\t'.join([
            str(element) for element in [
                FASTQ, total, no_barcode, barcode_f, guide_f, no_guide,
                int(time.time())
            ]
        ])
        print(newline, flush=True)
        w.write(newline + '\n')
    pd.DataFrame.from_dict(results_dic).to_csv(results_by_guide)
