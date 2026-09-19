#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:59:29 2017

@author: filippo
"""
import pandas as pd
import numpy as np


def annotateGene(annotations_df, ax, annotate_outside=False):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    def search_free_spot(point_x, point_y, height, width, annotation_places):
        position = np.array([
            point_x + 5 * ((point_x > 50) * 2 - 1),
            point_y + 5 * ((point_y > 50) * 2 - 1)
        ]).astype(int)
        test = True
        checked = set()

        direction = np.array([1, 0])
        while True:
            if len(checked) > 10000:
                position, ha, va, range_x, range_y = np.array([
                    point_x + 5, point_y + 5
                ]), 'center', 'center', (point_x + 5, point_x + 5 +
                                         width), (point_y + 5,
                                                  point_y + 5 + height)
                break
            turn = turn_right(direction)

            if tuple(position + turn) in checked:
                position = position + direction
            else:
                position = position + turn
                direction = turn
            position_check = tuple(position)
            if position_check in checked:
                continue
            checked.add(position_check)

            ha, va = 'center', 'center'

            range_x, range_y = (position[0] - width / 2 - 1, position[0] +
                                width / 2 + 1), (position[1] - height / 2 - 1,
                                                 position[1] + height / 2 + 1)
            if min(range_x[0], range_y[0]) < 0 or max(range_x[1],
                                                      range_y[1]) > 100:
                continue

            test = ((
                annotation_places.loc[range_x[0]:range_x[1],
                                      range_y[0]:range_y[1]]).sum().sum() == 0)

            # retest = annotation_places.copy()
            # retest.loc[range_x[0]:range_x[1], range_y[0]:range_y[1]] += 2
            # plt.imshow(retest.T, origin= 'lower')
            # plt.show()

            if test:
                break

        return position, ha, va, range_x, range_y

    turn_right = lambda x: np.array([x[1], -x[0]])

    ax_position = ax.get_position()
    xrange = ax_position.x1 - ax_position.x0
    yrange = ax_position.y1 - ax_position.y0
    start_x, end_x = ax.get_xlim()
    start_y, end_y = ax.get_ylim()
    ax_x_range = end_x - start_x
    ax_y_range = end_y - start_y

    anno_start_x = start_x - (ax_position.x0 * ax_x_range /
                              xrange) * annotate_outside
    anno_end_x = end_x + (
        (1 - ax_position.x1) * ax_x_range / xrange) * annotate_outside
    anno_start_y = start_y - (ax_position.y0 * ax_y_range /
                              yrange) * annotate_outside
    anno_end_y = end_y + (
        (1 - ax_position.y1) * ax_y_range / yrange) * annotate_outside

    anno_x_range = anno_end_x - anno_start_x
    anno_y_range = anno_end_y - anno_start_y

    to0100x = lambda x: int(round((x - anno_start_x) / anno_x_range * 100, 0))
    to0100y = lambda x: int(round((x - anno_start_y) / anno_y_range * 100, 0))
    toDataX = lambda x: (x / 100 * anno_x_range) + anno_start_x
    toDataY = lambda x: (x / 100 * anno_y_range) + anno_start_y

    annotation_places = pd.DataFrame(index=np.arange(100),
                                     columns=np.arange(100),
                                     data=0)
    for label, (x, y) in annotations_df.loc[:, ['x', 'y']].iterrows():
        to_100_x, to_100_y = to0100x(x), to0100y(y)
        annotation_places.loc[to_100_x - 2:to_100_x + 3,
                              to_100_y - 2:to_100_y + 3] = 1

    r = ax.figure.canvas.get_renderer()
    annotations_df.loc[:, 'priority_2'] = annotations_df.index.str.len()
    for label, (x, y, fontsize, color) in annotations_df.sort_values(
            'priority', kind='stable').sort_values(
                'priority',
                kind='stable').loc[:, 'x y fontsize color'.split()].iterrows():

        t = ax.text(0, 0, label, fontsize=fontsize)
        tbox = t.get_window_extent(r)
        dbox = tbox.transformed(ax.transData.inverted())
        del ax.texts[-1]
        point_x, point_y, width, height = to0100x(x), to0100y(
            y
        ), dbox.width / anno_x_range * 100, dbox.height / anno_y_range * 100

        (ann_x, ann_y), ha, va, range_x, range_y = search_free_spot(
            point_x, point_y, height, width, annotation_places)
        annotation_places.loc[range_x[0]:range_x[1], range_y[0]:range_y[1]] = 1

        bbox_props = dict(facecolor='white',
                          edgecolor=color,
                          boxstyle='round,pad=0.1',
                          alpha=0.8)
        arrow_props = dict(arrowstyle="-",
                           connectionstyle="arc3, rad=-0.0",
                           relpos=(0.5, 0.5),
                           linewidth=1,
                           edgecolor=color,
                           alpha=0.8,
                           shrinkA=0)

        ax.annotate(label,
                    xy=(x, y),
                    xycoords='data',
                    xytext=(toDataX(ann_x), toDataY(ann_y)),
                    textcoords='data',
                    arrowprops=arrow_props,
                    bbox=bbox_props,
                    size=fontsize,
                    fontstyle='italic',
                    color=color,
                    ha=ha,
                    va=va)


def annotateSignificance(ax,
                         x0,
                         y0,
                         x1,
                         y1,
                         y_offset,
                         x_offsets=[0, 0],
                         label='',
                         linewidth=1,
                         fontsize=7,
                         bbox_props='standard',
                         height=0):
    ###assuming x0 is lower
    if bbox_props == 'standard':
        bbox_props = dict(facecolor='white',
                          edgecolor='white',
                          boxstyle='round,pad=0',
                          alpha=1,
                          linewidth=0.5)
    if bbox_props == 'no_border':
        bbox_props = dict(facecolor='white',
                          edgecolor='white',
                          boxstyle='round,pad=0',
                          alpha=0,
                          linewidth=0.5)

    y_max = max(y0, y1)
    y_top = y_max + y_offset * (2 + height * 2)
    y_bottom_0, y_bottom_1 = y0 + y_offset, y1 + y_offset
    x_0_off, x_1_off = x0 + x_offsets[0], x1 + x_offsets[1]

    ax.plot([x_0_off, x_0_off, x_1_off, x_1_off],
            [y_bottom_0, y_top, y_top, y_bottom_1],
            clip_on=False,
            linewidth=1,
            color='k')
    ax.text((x_0_off + x_1_off) / 2,
            y_top + y_offset,
            label,
            ha='center',
            va='bottom',
            bbox=bbox_props,
            clip_on=False,
            fontsize=fontsize)


def biomartFlatten(biomart, pd, np):
    CHROMOSOMES = pd.Index(
        '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT'.
        split())

    df0 = pd.read_csv(biomart, index_col=2, sep='\t').fillna(-1000)
    df0 = df0.loc[CHROMOSOMES.intersection(
        df0.index.unique())].reset_index().set_index('Gene name')
    df = df0.loc[:, 'Genomic coding start,Genomic coding end'.split(',')]
    has_coding = (df.loc[:, 'Genomic coding start']
                  > 0) & (df.loc[:, 'Genomic coding end'] > 0)
    df = df[has_coding].astype(int).sort_values(
        by='Genomic coding end').sort_values(
            by='Genomic coding start',
            kind='mergesort').sort_index(kind='mergesort')

    dfTranscripts = df0.loc[:, 'Exon region start (bp),Exon region end (bp)'.
                            split(',')].astype(int).sort_values(
                                by='Exon region end (bp)').sort_values(
                                    by='Exon region start (bp)',
                                    kind='mergesort').sort_index(
                                        kind='mergesort')

    df.loc[:, 'Genomic coding end'] += 3
    df.loc[:, 'Genomic coding start'] -= 2
    dfTranscripts.loc[:, 'Exon region end (bp)'] += 1

    temp = '/home/filippo/Databases/Reference_Genomes/temp/temp'
    merged_temp = '/home/filippo/Databases/Reference_Genomes/temp/merged_temp'
    flattened = biomart.split('.')[0] + '_flattened.txt'

    df.to_csv(temp, sep='\t', header=None)
    command = 'bedtools merge -i %s > %s' % (temp, merged_temp)
    caller(command)
    merged = pd.read_csv(merged_temp,
                         index_col=0,
                         header=None,
                         names='start,end'.split(','),
                         sep='\t')
    counts = (merged.loc[:, 'end'] -
              merged.loc[:, 'start']).groupby(level=0).sum()

    dfTranscripts.to_csv(temp, sep='\t', header=None)
    command = 'bedtools merge -i %s > %s' % (temp, merged_temp)
    caller(command)
    mergedTranscripts = pd.read_csv(merged_temp,
                                    index_col=0,
                                    header=None,
                                    names='start,end'.split(','),
                                    sep='\t')
    countsTranscripts = (mergedTranscripts.loc[:, 'end'] -
                         mergedTranscripts.loc[:, 'start']).groupby(
                             level=0).sum()

    final = pd.DataFrame(
        index=counts.index,
        columns=
        'gene_ensembl_ID chromosome gene_start gene_end total_genomic_length_of_coding_sequence total_length_of_transcripts strand'
        .split()).fillna(0)
    final.loc[:, 'total_genomic_length_of_coding_sequence'] = counts
    final.loc[:, 'total_length_of_transcripts'] = countsTranscripts
    gene_positions = df0.groupby(level=0).first()
    final.loc[:, 'chromosome'] = gene_positions.loc[:,
                                                    'Chromosome/scaffold name']
    final.loc[:, 'gene_start'] = gene_positions.loc[:, 'Gene start (bp)']
    final.loc[:, 'gene_end'] = gene_positions.loc[:, 'Gene end (bp)']
    final.loc[:, 'gene_ensembl_ID'] = gene_positions.loc[:, 'Gene stable ID']
    final.loc[:, 'strand'] = gene_positions.loc[:, 'Strand']

    final.to_csv(flattened, index_label='gene_name', sep='\t')


#import pandas as pd
#import numpy as np
#biomart = '/home/filippo/Databases/Reference_Genomes/20190425_Homo_sapiens_ensembl_biomart_exons_filterCCDS.txt'
##biomart = '/home/filippo/Databases/Reference_Genomes/20190425_Mus_musculus_ensembl_biomart_exons.txt'
#biomartFlatten(biomart, pd,np)


def calculate_non_overlapping_significant_genesets(plot_path,
                                                   all_significant_GO,
                                                   datasets_label,
                                                   datasets,
                                                   children,
                                                   parents,
                                                   terms_metadata_dic,
                                                   q_max=0.1,
                                                   min_overlap_index=0,
                                                   min_ratio=1):
    import numpy as np
    is_relative_of_hit = {}
    dataset_filtered = all_significant_GO.loc[datasets]
    significance_filter = (dataset_filtered.observed_over_expected
                           >= min_ratio) & (dataset_filtered.overlap_index
                                            >= min_overlap_index) & (
                                                dataset_filtered.q <= q_max)
    filtered_overlaps = dataset_filtered[significance_filter]
    summary_df = filtered_overlaps.groupby('geneset_2').count().loc[:, ['q']]
    summary_df.loc[:, 'total_q'] = -np.log10(
        filtered_overlaps.groupby('geneset_2').q.prod())
    summary_df.loc[:, 'term_id_1'] = filtered_overlaps.groupby(
        'geneset_2').first()['term_id_1']
    summary_df.loc[:, 'term_id_2'] = filtered_overlaps.groupby(
        'geneset_2').first()['term_id_2']
    summary_df.loc[:, 'score'] = summary_df.q + summary_df.total_q / (
        summary_df.total_q.max() + 1)
    for dataset in datasets:
        summary_df[dataset + '_q'] = filtered_overlaps.loc[
            dataset, 'geneset_2 q'.split()].set_index('geneset_2')
    results_df = summary_df.sort_values('score', ascending=False).copy()
    count = 0
    ranks = []
    for geneset_2, geneset_data in results_df.iterrows():
        total_q, score = geneset_data.total_q, geneset_data.score
        term_id = terms_metadata_dic[geneset_2]['term']
        if term_id in is_relative_of_hit:
            # print (geneset_2, is_relative_of_hit[term_id])
            rank = 'SubsetOverset'
        else:
            rank = count
            print(geneset_2, rank, round(total_q, 2), round(score, 2))
            count += 1
            for child in children[term_id]:
                is_relative_of_hit[child] = is_relative_of_hit.get(
                    child, []) + [term_id]
            for parent in parents[term_id]:
                is_relative_of_hit[parent] = is_relative_of_hit.get(
                    parent, []) + [term_id]

        ranks.append(rank)
    results_df.loc[:, 'rank'] = ranks
    results_df.to_csv('{}/{}_q{}_oi{}_ratio{}.csv'.format(
        plot_path, datasets_label, q_max, min_overlap_index, min_ratio))


def caller(command, verbose=True):
    import subprocess
    import sys
    if verbose:
        print()
        print(command)
    try:
        retcode = subprocess.call(command, shell=True)
        if retcode < 0:
            print(sys.stderr, "Child was terminated by signal", -retcode)
        else:
            if verbose:
                print(sys.stderr, "Child returned", retcode)
    except OSError as e:
        print(sys.stderr, "Execution failed:", e)
    try:
        LOG.write(command + '\n')
    except NameError:
        a = 1
    if verbose:
        print
        print


def callerReturn(command):
    import subprocess
    try:
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True)
        output, retcode = process.communicate()
    except OSError as e:
        print("Execution failed:", e)
    return output


def complement(s):
    dic = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g',
        'n': 'n',
        '+': '+',
        '-': '-',
        '/': '/'
    }
    newS = ''
    for i in s:
        newS += dic[i]
    return newS


def fadeColor(
        c1,
        c2,
        np,
        mix=0
):  #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    assert len(c1) == len(c2)
    assert mix >= 0 and mix <= 1, 'mix=' + str(mix)
    rgb1 = np.array([int(c1[ii:ii + 2], 16) for ii in range(1, len(c1), 2)])
    rgb2 = np.array([int(c2[ii:ii + 2], 16) for ii in range(1, len(c2), 2)])
    rgb = ((1 - mix) * rgb1 + mix * rgb2).astype(int)
    c = '#' + ''.join([hex(a)[2:] for a in rgb])
    return c


class FigureWithAnnotations():

    def __init__(self,
                 figsize,
                 annotation_color='k',
                 ha='center',
                 va='center',
                 fontsize=7,
                 linewidth=1,
                 DPI=600):
        import matplotlib.pyplot as plt
        import matplotlib
        import datetime
        import numpy as np
        print('''
interactive figure started
HOVER on a point to see the label and position
CLICK on a point to annotate it
CLICK AND HOLD an annotation to drag it around
PRESS x to activate or deactivate label eraser, when activated
    labels will deleted instead of dragged
PRESS y to save figure and figure data
               ''')
        self.timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        self.figsize = figsize
        self.pltFigure = plt.figure(figsize=figsize)
        self.pltFigure.set_tight_layout(False)
        self.open = True
        self.annotations = {}
        self.latest_hover = ''
        self.deleter = False
        self.cm = 1 / 2.54
        self.linewidth = linewidth
        self.arrow_props = dict(arrowstyle="-",
                                connectionstyle="arc3, rad=-0.0",
                                relpos=(0.5, 0.5),
                                linewidth=self.linewidth,
                                edgecolor='k',
                                alpha=0.8,
                                shrinkA=0,
                                shrinkB=0)
        self.bbox_props = dict(facecolor='white',
                               edgecolor='k',
                               boxstyle='round,pad=0.3',
                               alpha=0.8,
                               linewidth=self.linewidth)

        self.canvas = self.pltFigure.canvas
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('key_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.checkPathCollection = matplotlib.collections.PathCollection
        self.checkAnnotation = matplotlib.text.Annotation
        self.annotation_color = annotation_color
        self.ha = ha
        self.va = va
        self.fontsize = fontsize
        self.array = np.array
        self.DPI = DPI

    def save_figure_and_data(self):
        print('\nPRESSED y: saving figure and figure data')
        print(self.figure_path)
        print(self.figure_info_path)
        print(self.dump_path, '\n')

        self.pltFigure.set_size_inches(self.figsize)
        self.pltFigure.savefig(fname=self.figure_path,
                               dpi=self.DPI,
                               pad_inches=0)

        with open(self.figure_info_path,
                  'wt') as w1, open(self.dump_path, 'wt') as w2:
            header = 'label,point_x,point_y,label_x,label_y,annotation_color,ha,va,fontsize\n'
            w1.write(header)
            w2.write(header)

            for annotated_point_key in self.annotations:
                annotation = self.annotations[annotated_point_key][
                    'annotation'].ref_artist
                point_x, point_y = annotation.tracker
                label_x, label_y = annotation.get_position()
                label = annotation.get_text()
                label_2, annotation_color, ha, va, fontsize = self.annotations[
                    annotated_point_key]['writer']
                assert label == label_2
                newline = ','.join([
                    str(element) for element in [
                        label, point_x, point_y, label_x, label_y,
                        annotation_color, ha, va, fontsize
                    ]
                ]) + '\n'
                w1.write(newline)
                w2.write(newline)
        self.open = False

    def annotate_point(self,
                       ind,
                       point_x,
                       point_y,
                       label_x,
                       label_y,
                       annotation_color,
                       ha,
                       va,
                       fontsize,
                       saved_label=None):

        if (point_x, point_y) in self.annotations:
            self.annotations[(point_x, point_y)]['overlay_point'].remove()
            self.annotations[(point_x,
                              point_y)]['annotation'].ref_artist.remove()
            del self.annotations[(point_x, point_y)]

        try:
            size = self.sizes[ind]
        except IndexError:
            size = self.sizes[0]
        color = self.scatter.get_facecolors()[ind]
        label = saved_label if saved_label is not None else self.labels[ind]

        self.annotations[(point_x, point_y)] = {
            'overlay_point':
            self.aax.scatter(point_x,
                             point_y,
                             zorder=self.current_zorder,
                             s=size,
                             color=color,
                             linewidth=self.linewidth,
                             edgecolors=annotation_color),
            'annotation':
            self.aax.annotate(label,
                              xy=(point_x, point_y),
                              xycoords='data',
                              xytext=(label_x, label_y),
                              textcoords='data',
                              arrowprops=self.arrow_props,
                              bbox=self.bbox_props,
                              fontsize=fontsize,
                              fontstyle='italic',
                              color=annotation_color,
                              ha=ha,
                              va=va,
                              zorder=self.current_zorder - 0.1).draggable(),
            'writer': (label, annotation_color, ha, va, fontsize)
        }
        self.annotations[(
            point_x, point_y)]['annotation'].ref_artist.tracker = (point_x,
                                                                   point_y)
        self.canvas.draw_idle()
        self.current_zorder += 1

    def set_arrow_props(self, arrow_props):
        self.arrow_props = arrow_props

    def set_bbox_props(self, bbox_props):
        self.bbox_props = bbox_props

    def annotation_ax(self, ax):
        from matplotlib.widgets import Button
        self.aax = ax
        self.save_button = Button(self.aax, 'save')

    def on_pick(self, event):
        artist = event.artist
        self.current_event = event
        if type(artist) == self.checkPathCollection:
            self.point_picked(event, artist)
        elif type(artist) == self.checkAnnotation:
            self.annotation_picked(event, artist)

    def on_hover(self, event):
        self.latest_hover
        if event.inaxes == self.aax:
            cont, ind = self.scatter.contains(event)
            if cont:
                index = ind['ind'][0]
                label = self.labels[index]
                if label != self.latest_hover:
                    print(label, toPrecision(self.x[index], 2),
                          toPrecision(self.y[index], 2))
                    self.latest_hover = label

    def on_press(self, event):
        if event.key == 'x':
            self.deleter = not self.deleter
            print('pressed x:', ['not deleting', 'deleting'][self.deleter])
        if event.key == 'y':
            print('pressed y: saving')
            self.save_figure_and_data()

    def annotation_picked(self, event, artist):
        if self.deleter:
            artist_tracker = artist.tracker
            artist.remove()
            self.annotations[artist_tracker]['overlay_point'].remove()
            del self.annotations[artist_tracker]
            self.canvas.draw_idle()
        else:
            self.latest_annotation = artist

    def closest_point_to_click(self, event):
        ind = event.ind
        if len(ind) == 1:
            return ind[0]
        else:
            return self.closest_point_to_xy(event.mouseevent.xdata,
                                            event.mouseevent.ydata)

    def closest_point_to_xy(self, x, y):

        return np.nanargmin(((abs(self.offsets.data -
                                  self.array([[x, y]]))**2).sum(axis=1)**0.5))

    def point_picked(self, event, artist):
        ind = self.closest_point_to_click(event)
        point_x, point_y = self.x[ind], self.y[ind]
        label_x, label_y = point_x + self.offset_x, point_y + self.offset_y
        self.annotate_point(ind, point_x, point_y, label_x, label_y,
                            self.annotation_color, self.ha, self.va,
                            self.fontsize)

    def add_axes(self, *args, **kwargs):
        self.aax = self.pltFigure.add_axes(*args, **kwargs)
        return self.aax

    def add_scatter_data(self, scatter, labels, predetermined_annotations,
                         figure_path, dump_path):
        import os
        import pandas as pd
        import numpy as np
        os.makedirs(dump_path, exist_ok=True)
        self.scatter = scatter
        self.offsets = self.scatter.get_offsets()
        self.x, self.y = self.offsets.data.T
        self.labels = np.array(labels)
        self.figure_path = figure_path
        self.figure_info_path = '.'.join(
            figure_path.split('.')[:-1]) + '_label_positions.csv'
        self.dump_path = f'{dump_path}{self.timestamp}_{self.figure_info_path.split("/")[-1]}'
        self.annotations = {}
        self.current_zorder = self.scatter.get_zorder() + 1
        self.sizes = self.scatter.get_sizes()
        self.previously_plotted_labels = set()
        ax_size_inches = self.aax.get_window_extent().transformed(
            self.pltFigure.dpi_scale_trans.inverted())
        ax_size_data = self.aax.get_window_extent().transformed(
            self.aax.transData.inverted())
        self.offset_x = 0.1 * (
            1 / (ax_size_inches.width / self.cm)) * ax_size_data.width
        self.offset_y = 0.1 * (
            1 / (ax_size_inches.height / self.cm)) * ax_size_data.height

        if os.path.exists(self.figure_info_path):
            previous_annotations = pd.read_csv(self.figure_info_path,
                                               index_col=0)
            for label, (point_x, point_y, label_x, label_y, annotation_color,
                        ha, va, fontsize) in previous_annotations.iterrows():
                ind = self.closest_point_to_xy(point_x, point_y)

                # FIX: Snap to exact array coordinates to prevent floating-point key mismatch
                exact_x, exact_y = self.x[ind], self.y[ind]

                self.annotate_point(ind,
                                    exact_x,
                                    exact_y,
                                    label_x,
                                    label_y,
                                    annotation_color,
                                    ha,
                                    va,
                                    fontsize,
                                    saved_label=label)

        for label in predetermined_annotations:
            matched_indices = np.where(self.labels == label)[0]

            if len(matched_indices) > 0:
                ind = matched_indices[0]

                point_x, point_y = self.x[ind], self.y[ind]
                print(point_x, point_y)

                if (point_x, point_y) not in self.annotations:
                    label_x, label_y = point_x + self.offset_x, point_y + self.offset_y

                    self.annotate_point(ind,
                                        point_x,
                                        point_y,
                                        label_x,
                                        label_y,
                                        self.annotation_color,
                                        self.ha,
                                        self.va,
                                        self.fontsize,
                                        saved_label=label)

                    self.previously_plotted_labels.add(label)


def figureOfSize(fig_x=9,
                 fig_y=9,
                 left=1,
                 right=1,
                 down=1,
                 up=1,
                 figure_with_annotations=False,
                 bbox_props=None,
                 arrow_props=None,
                 DPI=600,
                 no_ax=False,
                 fontsize=7):
    cm = 1 / 2.54
    fig_x, fig_y = fig_x * cm, fig_y * cm
    left, right, down, up = left * cm / fig_x, right * cm / fig_x, down * cm / fig_y, up * cm / fig_y
    ax_x_start, ax_y_start = left, down
    ax_x_length = 1 - ax_x_start - right
    ax_y_length = 1 - ax_y_start - up
    if figure_with_annotations:
        figure = FigureWithAnnotations(figsize=[fig_x, fig_y],
                                       DPI=DPI,
                                       fontsize=fontsize)
        if bbox_props:
            figure.set_bbox_props(bbox_props)
        if arrow_props:
            figure.set_bbox_props(arrow_props)
    else:
        import matplotlib.pyplot as plt
        figure = plt.figure(figsize=[fig_x, fig_y])
        figure.set_tight_layout(False)
    if no_ax:
        return figure
    else:
        ax = figure.add_axes(
            [ax_x_start, ax_y_start, ax_x_length, ax_y_length])
        return figure, ax


def fisherHelper(contingency_key):
    import scipy.stats as st

    a, b, c, d = eval(contingency_key)
    contingency_table = ((a, b), (c, d))
    p = st.fisher_exact(contingency_table)[1]
    # print(contingency, p, flush = True)
    return contingency_key, p


def flatten_sample_id(string):
    #check if string is an iterable
    if isinstance(string, pd.Series) or isinstance(string, list) or isinstance(
            string, tuple) or isinstance(string, np.ndarray):
        return [flatten_sample_id(element) for element in string]
    for non_common_sign in '.,:;!@#$%^&*()_+-=[]{}|\\/':
        string = string.replace(non_common_sign, '_')
    string = string.lower()
    string = string.replace(' ', '_')
    return string


def gencodeReferenceByGene(gencode_gtf):
    print('OUTDATED, use parseGTF instead')


def geneticCode(with_rna=False):
    dic = {
        'CTT': 'L',
        'ATG': 'M',
        'AAG': 'K',
        'AAA': 'K',
        'ATC': 'I',
        'AAC': 'N',
        'ATA': 'I',
        'AGG': 'R',
        'CCT': 'P',
        'ACT': 'T',
        'AGC': 'S',
        'ACA': 'T',
        'AGA': 'R',
        'CAT': 'H',
        'AAT': 'N',
        'ATT': 'I',
        'CTG': 'L',
        'CTA': 'L',
        'CTC': 'L',
        'CAC': 'H',
        'ACG': 'T',
        'CAA': 'Q',
        'AGT': 'S',
        'CAG': 'Q',
        'CCG': 'P',
        'CCC': 'P',
        'TAG': '*',
        'TAT': 'Y',
        'GGT': 'G',
        'TGT': 'C',
        'CGA': 'R',
        'CCA': 'P',
        'CGC': 'R',
        'GAT': 'D',
        'CGG': 'R',
        'TTT': 'F',
        'TGC': 'C',
        'GGG': 'G',
        'TGA': '*',
        'GGA': 'G',
        'TGG': 'W',
        'GGC': 'G',
        'TAC': 'Y',
        'GAG': 'E',
        'TCG': 'S',
        'TTA': 'L',
        'GAC': 'D',
        'CGT': 'R',
        'GAA': 'E',
        'TCA': 'S',
        'GCA': 'A',
        'GTA': 'V',
        'GCC': 'A',
        'GTC': 'V',
        'GCG': 'A',
        'GTG': 'V',
        'TTC': 'F',
        'GTT': 'V',
        'GCT': 'A',
        'ACC': 'T',
        'TTG': 'L',
        'TCC': 'S',
        'TAA': '*',
        'TCT': 'S'
    }
    if with_rna:
        for codon in sorted(list(dic.keys())):
            rna_codon = ''.join([[base, 'U'][base == 'T'] for base in codon])
            dic[rna_codon] = dic[codon]
    return dic


def getGOParentsAndChildren(GO_parents, GO_children):
    parents = dict()
    with open(GO_parents) as f:
        for line in f:
            ls = line.rstrip('\n').split('\t')
            parents[ls[0]] = eval(ls[1])
    children = dict()
    with open(GO_children) as f:
        for line in f:
            ls = line.rstrip('\n').split('\t')
            children[ls[0]] = eval(ls[1])
    return parents, children


def getMouseHumanOrthologous(gencode_mouse,
                             gencode_human,
                             hcop_orthology_15,
                             timestamp=None):
    import pandas as pd
    import numpy as np

    def splitter(group):
        total = set()
        for element in group:
            total.update(element.split(','))
        return len(total)

    def chooseBest(grouped):
        index, group = grouped
        best = group.database_support.max()
        only_one_best = (group.database_support == best).sum() == 1
        if only_one_best:
            return [index] + list(
                group[group.database_support == best].values[0])
        else:
            return [np.nan, np.nan, np.nan]

    gencodes_files = {'human': gencode_human, 'mouse': gencode_mouse}
    if not timestamp:
        timestamp = str(int(hcop_orthology_15.split('/')[-1].split('_')[0]))
    mouse_version = gencode_mouse.split('gencode.')[1].split('.')[0]
    human_version = gencode_human.split('gencode.')[1].split('.')[0]
    orthology_label = '{}_{}_{}'.format(timestamp, human_version,
                                        mouse_version)

    genecodes, unraveled, geneset, protein_coding, ensemble_to_symbol = {}, {}, {}, {}, {}
    for species in 'human mouse'.split():
        print(species, flush=True)
        genecodes[species] = pd.read_csv(gencodes_files[species],
                                         sep='\t',
                                         index_col=2,
                                         comment='#',
                                         header=None).loc['exon']
        unraveled[species] = genecodes[species][
            genecodes[species].loc[:,
                                   0] != 'chrM'].loc[:,
                                                     [0, 3, 4, 6, 8]].rename(
                                                         columns={
                                                             0: 'chromosome',
                                                             3: 'start',
                                                             4: 'end',
                                                             6: 'strand',
                                                             8: 'information'
                                                         })
        gencode_strings = unraveled[species].information.str.replace(
            ' ', '').str.split(';')
        print(species, 'step 1', flush=True)
        info_index = dict((element.split('"')[0], dici)
                          for dici, element in enumerate(gencode_strings[0]))
        unraveled[
            species].loc[:, 'protein_coding_transcript'] = gencode_strings.str[
                info_index['transcript_type']].str[:-1].str.split(
                    '"').str[1] == 'protein_coding'
        #        unraveled[species].loc[:,'protein_coding_transcript'] = genecodes[species].iloc[:,7].str.contains('tag "CCDS"')
        unraveled[species].loc[:, 'transcript_id'] = gencode_strings.str[
            info_index['transcript_id']].str.split('"').str[1].str.split(
                '.').str[0]
        unraveled[species].loc[:, 'gene_id'] = gencode_strings.str[
            info_index['gene_id']].str.split('"').str[1].str.split('.').str[0]
        unraveled[species].loc[:, 'gene_symbol'] = gencode_strings.str[
            info_index['gene_name']].str[:-1].str.split('"').str[1]
        print(species, 'step 2', flush=True)

        geneset[species] = pd.Index(unraveled[species].loc[:,
                                                           'gene_id'].unique())
        protein_coding[species] = pd.Index(
            unraveled[species][unraveled[species].protein_coding_transcript].
            loc[:, 'gene_id'].unique())
        ensemble_to_symbol[species] = unraveled[
            species].loc[:, ['gene_id', 'gene_symbol']].drop_duplicates(
            ).set_index('gene_id')
        unraveled[species].to_csv(gencodes_files[species][:-4] + '_exons.csv',
                                  sep='\t',
                                  index=None)

    raw_orthology = pd.read_csv(
        hcop_orthology_15,
        sep='\t').loc[:,
                      'human_ensembl_gene mouse_ensembl_gene support'.split()]
    raw_orthology.columns = raw_orthology.columns[:-1].tolist() + [
        'database_support'
    ]
    raw_orthology = raw_orthology[(raw_orthology.human_ensembl_gene != '-')
                                  & (raw_orthology.mouse_ensembl_gene != '-')]
    raw_orthology = raw_orthology.set_index('human_ensembl_gene', drop=False)

    for only_protein_coding in [True, False]:
        protein_coding_flag = ['allGenes',
                               'proteinCoding'][only_protein_coding]
        print(protein_coding_flag, flush=True)
        intersector = [geneset, protein_coding][only_protein_coding]
        refined_orthology = '/media/filippo/HD_1/Databases/HCOP_orthology/{}_hcop_refined_{}.txt'.format(
            orthology_label, protein_coding_flag)

        orthology = raw_orthology.loc[
            raw_orthology.index.unique().intersection(
                intersector['human'])].set_index('mouse_ensembl_gene')
        orthology = orthology.loc[orthology.index.unique().intersection(
            intersector['mouse'])].reset_index().groupby(
                ['human_ensembl_gene',
                 'mouse_ensembl_gene']).agg(splitter).reset_index()
        grouped = orthology.set_index('human_ensembl_gene').groupby(level=0)
        best_for_human = pd.DataFrame(
            [chooseBest(element) for element in grouped],
            columns=[
                'human_ensembl_gene', 'mouse_ensembl_gene', 'database_support'
            ]).dropna()

        grouped_2 = best_for_human.set_index('mouse_ensembl_gene').groupby(
            level=0)
        best_for_mouse = pd.DataFrame(
            [chooseBest(element) for element in grouped_2],
            columns=[
                'mouse_ensembl_gene', 'human_ensembl_gene', 'database_support'
            ]).dropna()

        human_mouse_ortho = best_for_mouse.copy()
        human_mouse_ortho.loc[:, 'mouse_symbol'] = ensemble_to_symbol[
            'mouse'].loc[human_mouse_ortho.loc[:, 'mouse_ensembl_gene'],
                         'gene_symbol'].values
        human_mouse_ortho.loc[:, 'human_symbol'] = ensemble_to_symbol[
            'human'].loc[human_mouse_ortho.loc[:, 'human_ensembl_gene'],
                         'gene_symbol'].values
        human_mouse_ortho.loc[:,
                              'human_ensembl_gene mouse_ensembl_gene human_symbol mouse_symbol database_support'
                              .split()].to_csv(refined_orthology, index=None)


class IntervalIndexer:

    def __init__(self):
        self.interval_matrix = {}

    def loadIntervals(self, gene_information, bisect, np, pd):
        for chromosome in sorted(list(gene_information.keys())):
            print(chromosome, end=' ')
            all_positions = gene_information[chromosome][:, 1:].reshape(
                -1).astype(int)
            names = gene_information[chromosome][:, 0]
            sorter = np.argsort(all_positions)
            sorted_positions = all_positions[sorter]
            resorter = np.argsort(sorter)
            reshaped_resorter = resorter.reshape((-1, 2))
            self.interval_matrix[chromosome] = pd.DataFrame(
                index=names,
                columns=sorted(list(set(sorted_positions)))).fillna(0)
            for name, (start_index, end_index) in zip(names,
                                                      reshaped_resorter):
                self.interval_matrix[chromosome].loc[
                    name, sorted_positions[start_index:end_index]] = 1

            self.interval_matrix[chromosome] = self.interval_matrix[
                chromosome].astype(bool)

    def readRefSeq(self,
                   refseq,
                   bisect,
                   np,
                   pd,
                   fil,
                   os,
                   from_file=False,
                   filter_coding=True):
        store_folder = '/media/filippo/HD_1/Databases/GENE_MATRIX_HD5/%s/' % refseq.split(
            '/')[-1]
        command = 'mkdir %s' % store_folder
        fil.caller(command)
        if from_file:
            for chromosome_h5 in os.listdir(store_folder):
                chromosome = chromosome_h5.split('.')[0]
                self.interval_matrix[chromosome] = pd.read_hdf(
                    store_folder + chromosome_h5, 'chr' + chromosome)
        else:
            gene_information = {}
            header_dic = {}
            with open(refseq) as f:
                header = next(f).lstrip('#').rstrip('\n').split('\t')
                for hi, head in enumerate(header):
                    header_dic[head] = hi
                for head in header:
                    if 'chr' in head.lower():
                        print('chromosome header:', head)
                        chromosome_index = header_dic[head]
                    if head.lower() == 'name2' or head.lower() == 'gene':
                        print('gene header:', head)
                        gene_index = header_dic[head]
                    if head.lower() == 'name' or head.lower() == 'transcript':
                        print('transcript header:', head)
                        transcript_index = header_dic[head]
                    if 'start' in head.lower():
                        print('interval start:', head)
                        start_index = header_dic[head]
                    if 'end' in head.lower():
                        print('interval end:', head)
                        end_index = header_dic[head]
                for line in f:
                    ls = line.rstrip('\n').split('\t')
                    chromosome = ls[chromosome_index].lstrip('chr')
                    gene = ls[gene_index].upper()
                    transcript = ls[transcript_index].upper()
                    start = int(ls[start_index])
                    end = int(ls[end_index])
                    if filter_coding and not 'M' in transcript:
                        continue
                    if chromosome not in gene_information:
                        gene_information[chromosome] = {}
                    if gene not in gene_information[chromosome]:
                        gene_information[chromosome][gene] = [start, end]
                    gene_information[chromosome][gene][0] = min(
                        gene_information[chromosome][gene][0], start)
                    gene_information[chromosome][gene][1] = max(
                        gene_information[chromosome][gene][1], end)
            for chromosome in gene_information:
                gene_information[chromosome] = np.array([[
                    gene, gene_information[chromosome][gene][0],
                    gene_information[chromosome][gene][1]
                ] for gene in gene_information[chromosome]])
            self.loadIntervals(gene_information, bisect, np, pd)

            for chromosome in self.interval_matrix:
                self.interval_matrix[chromosome].to_hdf(
                    store_folder + chromosome + '.h5', 'chr' + chromosome)


def overlap_among_sets(
    sets_1_df,
    sets_2_df,
    n_all,
    terms_metadata_df,
    label,
    threads=24,
    results_folder='/media/filippo/HD_2/Human_studies/metaanalysis_results/pathway_analysis/'
):
    import pandas as pd
    import numpy as np
    from multiprocessing import Pool
    import statsmodels.stats.multitest as stmu
    import scipy.stats as st

    def formatContingencies(contingencies_subset_df):
        invert_max = contingencies_subset_df.iloc[:, [1, 2]].max(axis=1)
        invert_min = contingencies_subset_df.iloc[:, [1, 2]].min(axis=1)
        contingencies_subset_df.iloc[:, 1] = invert_max
        contingencies_subset_df.iloc[:, 2] = invert_min
        contingency_tuples = contingencies_subset_df.to_records(index=False)
        return [str(element) for element in contingency_tuples]

    print('calculating overlap bewteen ', label, flush=True)
    print('\tmultiple testing is performed for each set in sets_1_dic',
          flush=True)
    seen_contingencies = pd.read_csv(
        '/media/filippo/HD_1/Databases/seen_contingency_tables.tsv',
        sep='\t',
        index_col=0)['0'].to_dict()
    filtered_sets_1_df, filtered_sets_2_df = sets_1_df[sets_1_df.sum(
        axis=1).astype(bool)].fillna(0).astype(bool), sets_2_df[sets_2_df.sum(
            axis=1).astype(bool)].fillna(0).astype(bool)

    possible_overlapping_genes = filtered_sets_1_df.index.intersection(
        filtered_sets_2_df.index)
    nested_overlap_dic = {
        element: {
            element2: 0
            for element2 in filtered_sets_2_df.columns
        }
        for element in filtered_sets_1_df.columns
    }
    dic_1, dic_2 = filtered_sets_1_df.T.to_dict(
        orient='series'), filtered_sets_2_df.T.to_dict(orient='series')

    for gene in possible_overlapping_genes:
        gene_index_1 = dic_1[gene][dic_1[gene]].index
        gene_index_2 = dic_2[gene][dic_2[gene]].index

        for dic_key_1 in gene_index_1:
            for dic_key_2 in gene_index_2:
                nested_overlap_dic[dic_key_1][dic_key_2] += 1
    contingencies_df = pd.DataFrame({
        'size_overlap':
        pd.DataFrame.from_dict(nested_overlap_dic, orient='index').stack()
    })
    non_self_geneset = contingencies_df.index.get_level_values(
        0) != contingencies_df.index.get_level_values(1)
    contingencies_df = contingencies_df[non_self_geneset]
    contingencies_df.index.names = 'geneset_1 geneset_2'.split()

    # contingencies_df = pd.DataFrame(index = combinations_index, columns = 'n_remaining n1 n2 size_overlap'.split(), data = 0).astype(int)
    #
    # for gene in possible_overlapping_genes:
    #     contingencies_df.loc[(dic_1[gene][dic_1[gene]].index, dic_2[gene][dic_2[gene]].index), 'size_overlap'] += 1

    contingencies_df.loc[:, 'size_1'] = filtered_sets_1_df.sum().loc[
        contingencies_df.index.get_level_values(0)].values
    contingencies_df.loc[:, 'size_2'] = filtered_sets_2_df.sum().loc[
        contingencies_df.index.get_level_values(1)].values
    contingencies_df.loc[:,
                         'n1'] = contingencies_df.size_1 - contingencies_df.size_overlap.values
    contingencies_df.loc[:,
                         'n2'] = contingencies_df.size_2 - contingencies_df.size_overlap.values
    contingencies_df.loc[:,
                         'n_remaining'] = n_all - contingencies_df.n1 - contingencies_df.n2 - contingencies_df.size_overlap
    contingencies_df.loc[:, 'contingency_key'] = formatContingencies(
        contingencies_df.loc[:, 'n_remaining n1 n2 size_overlap'.split()])
    new_contingencies = [
        element for element in contingencies_df.contingency_key.unique()
        if not element in seen_contingencies
    ]

    print(len(new_contingencies),
          'contingency tables never seen before',
          flush=True)
    if len(new_contingencies):
        p = Pool(threads)
        for ti in range(int(len(new_contingencies) / 100) + 1):
            updater = dict(
                p.map(fisherHelper,
                      new_contingencies[ti * 110:(ti + 1) * 110],
                      chunksize=1))
            seen_contingencies.update(updater)
            print(ti * 100, flush=True, end=' ')
        p.close()
        p.join()
        pd.Series(seen_contingencies).to_csv(
            '/media/filippo/HD_1/Databases/seen_contingency_tables.tsv',
            sep='\t')

    contingencies_df.loc[:, 'q'] = np.nan

    contingencies_df.loc[:, 'fisher_p'] = contingencies_df.contingency_key.map(
        seen_contingencies)

    for label_1 in contingencies_df.index.get_level_values(0).unique():
        qs = stmu.multipletests(contingencies_df.loc[label_1].fisher_p,
                                alpha=0.05,
                                method='fdr_bh')[1]
        contingencies_df.loc[label_1, 'q'] = qs

    contingencies_df.loc[:,
                         'expected'] = contingencies_df.size_1 * contingencies_df.size_2 / n_all
    contingencies_df.loc[:,
                         'overlap_1'] = contingencies_df.size_overlap / contingencies_df.size_1
    contingencies_df.loc[:,
                         'overlap_2'] = contingencies_df.size_overlap / contingencies_df.size_2
    contingencies_df.loc[:,
                         'observed_over_expected'] = contingencies_df.size_overlap / contingencies_df.expected

    contingencies_df.loc[:, 'overlap_index'] = contingencies_df.loc[:, [
        'overlap_1', 'overlap_2'
    ]].max(axis=1)

    contingencies_df.loc[:, 'term_id_1'] = terms_metadata_df.reindex(
        contingencies_df.index.get_level_values(0)).term.values
    contingencies_df.loc[:, 'term_id_2'] = terms_metadata_df.reindex(
        contingencies_df.index.get_level_values(1)).term.values

    contingencies_df.to_csv('{}genesets_overlaps_{}.csv'.format(
        results_folder, label))

    contingencies_df_unique_index = pd.MultiIndex.from_arrays(
        np.array([
            tuple(sorted(list(element))) for element in contingencies_df.index
        ]).T)
    unique_contingencies = contingencies_df[~contingencies_df_unique_index.
                                            duplicated(keep='first')].copy()
    print('shapes', contingencies_df.shape, unique_contingencies.shape)
    unique_contingencies.loc[:, 'q'] = stmu.multipletests(
        unique_contingencies.fisher_p, alpha=0.05, method='fdr_bh')[1]
    unique_contingencies.to_csv('{}genesets_overlaps_{}_unique.csv'.format(
        results_folder, label))

    filtered_sets_1_df.astype(int).to_csv(
        '{}genesets_overlaps_{}_1.csv'.format(results_folder, label))
    filtered_sets_2_df.astype(int).to_csv(
        '{}genesets_overlaps_{}_2.csv'.format(results_folder, label))


def parseGTF(gtf,
             species='human',
             transcript_information=True,
             gene_name='gene_name',
             gene_id='gene_id',
             transcript_id='transcript_id',
             biotype='gene_biotype',
             gene_information=True):
    import pandas as pd
    import pandas as pd
    import re
    import numpy as np

    re_gene = re.compile(gene_id + ' "(.*?)";')
    re_transcript = re.compile(transcript_id + ' "(.*?)";')
    re_symbol = re.compile(gene_name + ' "(.*?)";')
    re_type = re.compile(biotype + ' "(.*?)";')

    gtf_df = pd.read_csv(gtf, comment='#', sep='\t', index_col=2, header=None)
    exons = gtf_df.loc['exon'].copy()
    exons.loc[:, 'transcript'] = exons.iloc[:,
                                            7].str.extract(re_transcript,
                                                           expand=False).values
    exons.loc[:, 'length'] = exons.iloc[:, 3] - exons.iloc[:, 2]
    exons.loc[:, f'ensembl_gene_id_{species}'] = exons.iloc[:, 7].str.extract(
        re_gene, expand=False).values
    grouped = exons.loc[:, f'transcript length ensembl_gene_id_{species}'.
                        split()].groupby(
                            [f'ensembl_gene_id_{species}',
                             'transcript']).sum()
    if gene_information:
        genes = gtf_df.loc['gene'].copy()
        genes.index = genes.iloc[:, 7].str.extract(re_gene, expand=False)
        genes.index.name = f'ensembl_gene_id_{species}'
        results = pd.DataFrame({
            f'gene_symbol_{species}':
            genes.iloc[:, 7].str.extract(re_symbol, expand=False),
            'biotype':
            genes.iloc[:, 7].str.extract(re_type, expand=False),
            'chromosome':
            genes.iloc[:, 0],
            'start':
            genes.iloc[:, 2],
            'end':
            genes.iloc[:, 3],
            'mean_transcript_length':
            grouped.groupby(level=0).mean().length.loc[genes.index],
            'median_transcript_length':
            grouped.groupby(level=0).median().length.loc[genes.index],
        })

        results.index.name = f'ensembl_gene_id_{species}'
        results.to_csv(
            gtf.replace('.gz', '').replace('.gtf', '_geneReferenceData.csv'))

    if transcript_information:
        transcripts = gtf_df.loc['transcript'].copy()
        transcripts.index = transcripts.iloc[:, 7].str.extract(re_transcript,
                                                               expand=False)
        transcripts.index.name = f'ensembl_transcript_id_{species}'
        transcripts_results = pd.DataFrame({
            f'gene_symbol_{species}':
            transcripts.iloc[:, 7].str.extract(re_symbol, expand=False),
            f'ensembl_gene_id_{species}':
            transcripts.iloc[:, 7].str.extract(re_gene, expand=False),
            'biotype':
            transcripts.iloc[:, 7].str.extract(re_type, expand=False),
            'chromosome':
            transcripts.iloc[:, 0],
            'start':
            transcripts.iloc[:, 2],
            'end':
            transcripts.iloc[:, 3],
            'genomic_length':
            transcripts.iloc[:, 3] - transcripts.iloc[:, 2],
            'transcript_length':
            grouped.groupby(level=1).sum().length.loc[transcripts.index],
        })
        transcripts_results.index.name = f'ensembl_transcript_id_{species}'
        transcripts_results.to_csv(
            gtf.replace('.gz', '').replace('.gtf',
                                           '_transcriptReferenceData.csv'))


def permutation_test(a, b, permutations, np):
    n, k = len(a), 0
    observed_difference = np.abs(np.mean(a) - np.mean(b))
    pooled_sample = np.concatenate([a, b])
    for j in range(permutations):
        np.random.shuffle(pooled_sample)
        k += (np.abs(np.mean(pooled_sample[:n]) - np.mean(pooled_sample[n:]))
              >= observed_difference)
    return k / permutations


def project3Don2D(ypos, xpos, xneg, np):
    center_to_upper = 1 - (1 / (np.sqrt(3) * 2))
    center_to_lower = 1 / (np.sqrt(3) * 2)
    lower_center_to_lower_apex = 0.5
    x = np.zeros_like(ypos)
    y = np.zeros_like(ypos)
    x += np.array((xpos * lower_center_to_lower_apex).tolist())
    x -= np.array((xneg * lower_center_to_lower_apex).tolist())
    y -= np.array((xpos * center_to_lower).tolist())
    y -= np.array((xneg * center_to_lower).tolist())
    y += np.array((ypos * center_to_upper).tolist())
    ax_x = np.array([0., 0., 0., 0.5, 0., -0.5, 0.])
    ax_y = np.array([0., 0.71132487, 0., -0.28867513, 0., -0.28867513, 0.])
    return x, y, ax_x, ax_y


def pseudonymize(
    real_ids,
    previous_table='/home/filippo/wfiler2_onkosig/New Server_October 2016/Team_Data/Filippo_Beleggia_Data/Databases/pseudonymization/'
):
    import time
    import pandas as pd
    import numpy as np
    import random
    import string
    import os
    import hashlib
    np.random.seed(0)
    os.makedirs(previous_table, exist_ok=True)
    existing_tables = []
    timestamp = time.strftime('%Y%m%d_%H%M%S')

    print('\npseudonymization')

    flattened_strings = [
        flatten_sample_id(string) for string in real_ids.index
    ]
    for string, flattened_string in zip(real_ids.index, flattened_strings):
        if string != flattened_string:
            print('FLATTENED', string, flattened_string)
    real_ids.index = flattened_strings

    for table in os.listdir(previous_table):
        if 'pseudonymization_table.csv' in table:
            existing_tables.append(table)

    if existing_tables:
        latest_table = max(existing_tables)
        print('pseudonymization continues from', latest_table)
        previous_table_df = pd.read_csv(
            os.path.join(previous_table, latest_table))
    else:
        previous_table_df = pd.DataFrame(columns=['original', 'pseudonymized'])
    new_columns = previous_table_df.columns.difference(real_ids.columns)

    for column in new_columns:
        previous_table_df[column] = ''

    for string in real_ids.index:
        if string in previous_table_df.index:
            old_data = previous_table_df.loc[string]
            for column in old_data.columns.intersection(real_ids.columns):
                if old_data[column] != real_ids.loc[string, column]:
                    raise ValueError(
                        f'{string} has different {column} in the old and new table'
                    )
            previous_table_df.loc[string,
                                  new_columns] = real_ids.loc[string,
                                                              new_columns]
        else:
            previous_table_df.loc[string] = real_ids.loc[string]
            previous_table_df.loc[string, 'pseudonymized'] = hashlib.md5(
                (string + timestamp).encode()).hexdigest()

    previous_table_df.to_csv(
        os.path.join(previous_table,
                     f'{timestamp}_pseudonymization_table.csv'))
    return previous_table_df.loc[flattened_strings]


def recurRegion(size, BASES, start=''):
    if len(start) == size:
        yield start
    else:
        for base in BASES:
            for rec in recurRegion(size, BASES, (start + base)):
                yield rec


def recurRegionDic(size, BASES, start=''):
    gen = recurRegion(size, BASES, start)
    dic, dic2 = {}, {}
    for i, sequence in enumerate(gen):
        dic[i] = sequence
        dic2[sequence] = i
    return dic, dic2


def reverse(s):
    newS = ''
    for i in range(len(s) - 1, -1, -1):
        newS += s[i]
    return newS


def reverseComplement(s):
    newS = ''
    dic = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g',
        'n': 'n',
        '+': '+',
        '-': '-',
        '/': '/'
    }
    for i in range(len(s) - 1, -1, -1):
        newS += dic[s[i]]
    return newS


def searchForScript(keywords,
                    folder='/home/filippo/wfiler2_onkosig/scripts/Filippo'
                    ):  #'/home/filippo/Scripts/'):
    import os

    results = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if not name.endswith('.py'):
                continue
            remake_path = '/'.join([root, name])
            with open(remake_path, encoding="ISO-8859-1") as f:
                data = f.read()
                title, content = 0, 0
                for keyword in keywords:
                    if keyword in name:
                        title += 1
                    if keyword in data:
                        content += 1
                if title or content:
                    results.append((title, content, remake_path))

    for element in sorted(results):
        print(element)


def setupForFigure(fontsize=7):
    import matplotlib
    from matplotlib import font_manager
    msft_fonts = '/usr/share/fonts/truetype/msttcorefonts/*.ttf'
    matplotlib_font = matplotlib.__file__.replace('__init__.py',
                                                  'mpl-data/fonts/ttf/')

    matplotlib.rcdefaults()
    command = f'cp {msft_fonts} {matplotlib_font}'
    caller(command)

    command = 'rm /home/filippo/.cache/matplotlib/* -r'
    caller(command)

    arial_path = matplotlib.font_manager.findfont('Arial')
    arial_object = matplotlib.font_manager.get_font(arial_path)

    matplotlib.rcParams.update(
        {  #"text.usetex": True,
            # "font.family": "sans-serif",
            "font.family": "Arial",
            "font.sans-serif": 'Arial',
            "font.serif": 'Arial',
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": 'none',
            'text.usetex': False,
            'figure.autolayout': False,
            'font.size': fontsize,
            "mathtext.fontset": "custom",
            "mathtext.it": "Arial:italic",
            "mathtext.rm": "Arial"
        }, )
    cm = 1 / 2.54
    return matplotlib.rcParams, cm, fontsize


def toPrecision(x, n=3, with_zero=True, return_string=True):
    import pandas as pd
    if type(x) == pd.core.frame.DataFrame:
        return x.applymap(toPrecision,
                          n=n,
                          with_zero=with_zero,
                          return_string=return_string)
    if type(x) == pd.core.frame.Series:
        if x.dtype == 'str':
            return x
        return x.apply(toPrecision,
                       n=n,
                       with_zero=with_zero,
                       return_string=return_string)
    string_format = '{:.%sg}' % n

    if type(x) == str:
        if return_string:
            return x
        x = float(x)

    try:
        iter(x)
        return [
            toPrecision(element,
                        n,
                        with_zero=with_zero,
                        return_string=return_string) for element in x
        ]
    except TypeError:
        formatted = string_format.format(x)
        if return_string:
            if with_zero:
                return formatted
            else:
                return formatted.replace('0.', '.')
        else:
            return float(formatted)


class ZarrArray:

    def __init__(self, zarr_path_columns):
        import dask.array as da
        import pandas as pd
        self.zarr_path_columns = zarr_path_columns
        self.zarr_path_rows = zarr_path_columns.replace(
            '_columns.zarr', '_rows.zarr')
        self.Xr = da.from_zarr(self.zarr_path_rows)
        self.Xc = da.from_zarr(self.zarr_path_columns)
        obs_path = self.zarr_path_columns.replace('_columns.zarr',
                                                  '.zarr_obs.csv')
        self.obs = pd.read_csv(obs_path, index_col=0)
        var_path = self.zarr_path_columns.replace('_columns.zarr',
                                                  '.zarr_var.csv')
        self.var = pd.read_csv(var_path, index_col=0)
        assert self.obs.shape[0] == self.Xr.shape[0] == self.Xc.shape[0]
        assert self.var.shape[0] == self.Xr.shape[1] == self.Xc.shape[1]
