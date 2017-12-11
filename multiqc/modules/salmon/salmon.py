#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
from scipy import spatial

from multiqc import config
from multiqc.plots import linegraph
from multiqc.plots import heatmap
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.GCModel import GCModel
from multiqc.modules.SEQModel import SEQModel

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
        href='http://combine-lab.github.io/salmon/',
        info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        self.gc_bias = False
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename( os.path.dirname(f['root']) )
            s_name = self.clean_s_name(s_name, f['root'])
            self.salmon_meta[s_name] = json.loads(f['f'])
        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        self.gc_bias_path_list = []
        self.seq_bias_path_list = []
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename( os.path.dirname(f['root']) )
                s_name = self.clean_s_name(s_name, f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed
                '''
                Check the meta_info.json file to check whether the salmon output is ran with gc bias or not.
                If ran with gc_bias then add its absolute path to the list of paths where the salmon output is ran with gc_bias
                '''
                meta_json_file_path = os.path.join(os.path.dirname(f['root']),'aux_info','meta_info.json')
                gc_bias_base_dir = os.path.dirname(f['root'])
                with open(meta_json_file_path,'r') as meta_data_file:
                    meta_info_data = json.load(meta_data_file)
                self.gc_bias = meta_info_data['gc_bias_correct']
                self.seq_bias = meta_info_data['seq_bias_correct']
                if self.gc_bias:
                    self.gc_bias_path_list.append(os.path.abspath(gc_bias_base_dir))
                if self.seq_bias:
                    self.seq_bias_path_list.append(os.path.abspath(gc_bias_base_dir))


        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning

        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
            'title': '% Aligned',
            'description': '% Mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['num_mapped'] = {
            'title': 'M Aligned',
            'description': 'Mapped reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Fragment Length Distribution',
            'ylab': 'Fraction',
            'xlab': 'Fragment Length (bp)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        '''
            Iterate over the list of paths where each path has salmon output ran with gc_bias.
            Using the GCModel class's utlity functions compute observed array, expected array and the weights.
            Multiply the observed array and expected array with the corresponding weights and create a Ordered Dictionary,
            containing the ratio of observed by expected array. Plot that Ordered Dict using matplotlib.
        '''
        self.gc_first_model_ratio = dict()
        self.gc_second_model_ratio = dict()
        self.gc_third_model_ratio = dict()
        self.gc_avg_ratio = dict()
        self.seq_three_prime = dict()
        self.seq_five_prime = dict()
        self.gc_average_data = []
        self.gc_heatmap_labels = []
        self.gc_heatmap_data = []

        for path_var in self.gc_bias_path_list:
            gc_model = GCModel()
            gc_model.from_file(path_var)
            obs_array = gc_model.obs_.tolist()
            exp_array = gc_model.exp_.tolist()
            obs_weights = list(gc_model.obs_weights_)
            exp_weights = list(gc_model.exp_weights_)
            self.path_var = path_var.split('/')[-2]
            ratio_dict = dict()
            avg_ratio_dict = OrderedDict()
            for i in range(len(obs_array)):
                obs = obs_array[i]
                exp = exp_array[i]
                obs_weight = obs_weights[i]
                exp_weight = exp_weights[i]
                ratio_value = OrderedDict()
                j = 1
                for o,e in zip(obs,exp):
                    ratio = (o*obs_weight)/(e*exp_weight)
                    ratio_value[j] = ratio
                    try:
                        avg_ratio_dict[j] += ratio
                    except:
                        avg_ratio_dict[j] = ratio
                    j += 1
                ratio_dict[i] = ratio_value

            for k in list(avg_ratio_dict.keys()):
                avg_ratio_dict[k] /= len(obs_array)

            self.gc_first_model_ratio[self.path_var] = ratio_dict[0]
            self.gc_second_model_ratio[self.path_var] = ratio_dict[1]
            self.gc_third_model_ratio[self.path_var] = ratio_dict[2]
            self.gc_avg_ratio[self.path_var] = avg_ratio_dict
            self.gc_average_data.append(list(avg_ratio_dict.values()))
            self.gc_heatmap_labels.append(self.path_var)

        for avg_data1 in self.gc_average_data:
            cosine_distance_vector = []
            for avg_data2 in self.gc_average_data:
                cosine_distance_vector.append(spatial.distance.cosine(avg_data1,avg_data2))
            self.gc_heatmap_data.append(cosine_distance_vector)

        files = list(self.gc_first_model_ratio.keys())
        self.model_ratios = dict()

        firstModelAvg = OrderedDict()
        secondModelAvg = OrderedDict()
        thirdModelAvg = OrderedDict()

        for k in files:
            firstModel = self.gc_first_model_ratio[k]
            secondModel = self.gc_second_model_ratio[k]
            thirdModel = self.gc_third_model_ratio[k]

            for key in list(firstModel.keys()):
                try:
                    firstModelAvg[key] += firstModel[key]
                    secondModelAvg[key] += secondModel[key]
                    thirdModelAvg[key] += thirdModel[key]
                except:
                    firstModelAvg[key] = firstModel[key]
                    secondModelAvg[key] = secondModel[key]
                    thirdModelAvg[key] = thirdModel[key]


        for k in list(firstModelAvg.keys()):
            firstModelAvg[k] = float(firstModelAvg[k]/len(files))
            secondModelAvg[k] = float(secondModelAvg[k]/len(files))
            thirdModelAvg[k] = float(thirdModelAvg[k]/len(files))

        modelAvg = {"First Model": firstModelAvg, "Second Model": secondModelAvg, "Third Model": thirdModelAvg}


        self.seq_3prime_ratio = dict()
        self.seq_5prime_ratio = dict()
        self.nucleotides = ['A','C','G','T']
        self.seq_3prime_avg_data = []
        self.seq_5prime_avg_data = []
        for path_var in self.seq_bias_path_list:
            seq_model = SEQModel()
            seq_model.from_file(path_var)
            obs3_array = seq_model.obs3_prime.tolist()
            exp3_array = seq_model.exp3_prime.tolist()
            obs5_array = seq_model.obs5_prime.tolist()
            exp5_array = seq_model.exp5_prime.tolist()
            self.path_var = path_var.split('/')[-2]
            ratio_dict_3prime = dict()
            ratio_dict_5prime = dict()
            avg_3prime_array = [0]*len(obs3_array[0])
            avg_5prime_array = [0]*len(obs5_array[0])
            for i in range(len(self.nucleotides)):
                obs_3prime = obs3_array[i]
                exp_3prime = exp3_array[i]
                obs_5prime = obs5_array[i]
                exp_5prime = exp5_array[i]
                ratio_3prime_dict = OrderedDict()
                ratio_5prime_dict = OrderedDict()
                j = 1
                for o,e in zip(obs_3prime, exp_3prime):
                    ratio = o/e
                    ratio_3prime_dict[j] = ratio
                    avg_3prime_array[j-1] += ratio
                    j += 1
                ratio_dict_3prime[self.nucleotides[i]] = ratio_3prime_dict 

                j = 1
                for o,e in zip(obs_5prime, exp_5prime):
                    ratio = o/e
                    ratio_5prime_dict[j] = ratio
                    avg_5prime_array[j-1] = ratio
                    j += 1
                ratio_dict_5prime[self.nucleotides[i]] = ratio_5prime_dict

            self.seq_3prime_avg_data.append([x/len(self.nucleotides) for x in avg_3prime_array])
            self.seq_5prime_avg_data.append([x/len(self.nucleotides) for x in avg_5prime_array])
            self.seq_3prime_ratio[self.path_var] = ratio_dict_3prime
            self.seq_5prime_ratio[self.path_var] = ratio_dict_5prime
        self.seq_3prime_heatmap_data = []
        self.seq_5prime_heatmap_data = []
        for avg_data1 in self.seq_3prime_avg_data:
            cosine_distance_vector = []
            for avg_data2 in self.seq_3prime_avg_data:
                cosine_distance_vector.append(spatial.distance.cosine(avg_data1,avg_data2))
            self.seq_3prime_heatmap_data.append(cosine_distance_vector)

        for avg_data1 in self.seq_5prime_avg_data:
            cosine_distance_vector = []
            for avg_data2 in self.seq_5prime_avg_data:
                cosine_distance_vector.append(spatial.distance.cosine(avg_data1,avg_data2))
            self.seq_5prime_heatmap_data.append(cosine_distance_vector)

        seq_heat_map_labels = [x.split('/')[-2] for x in self.seq_bias_path_list]

        A3_dict = dict()
        C3_dict = dict()
        G3_dict = dict()
        T3_dict = dict()
        A5_dict = dict()
        C5_dict = dict()
        T5_dict = dict()
        G5_dict = dict()

        for k in list(self.seq_3prime_ratio.keys()):
            A3_dict[k] = self.seq_3prime_ratio[k]['A']
            C3_dict[k] = self.seq_3prime_ratio[k]['C']
            G3_dict[k] = self.seq_3prime_ratio[k]['G']
            T3_dict[k] = self.seq_3prime_ratio[k]['T']
        
        for k in list(self.seq_5prime_ratio.keys()):
            A5_dict[k] = self.seq_5prime_ratio[k]['A']
            C5_dict[k] = self.seq_5prime_ratio[k]['C']
            G5_dict[k] = self.seq_5prime_ratio[k]['G']
            T5_dict[k] = self.seq_5prime_ratio[k]['T']

        A3_avg = dict()
        C3_avg = dict()
        G3_avg = dict()
        T3_avg = dict()
        A5_avg = dict()
        C5_avg = dict()
        T5_avg = dict()
        G5_avg = dict()
        files_count = len(self.seq_3prime_ratio.keys())

        for k in list(self.seq_3prime_ratio.keys()):
            A3 = A3_dict[k]
            C3 = C3_dict[k]
            G3 = G3_dict[k]
            T3 = T3_dict[k]
            A5 = A5_dict[k]
            C5 = C5_dict[k]
            G5 = G5_dict[k]
            T5 = T5_dict[k]
            for key in list(A3.keys()):
                try:
                    A3_avg[key] += A3[key]
                    C3_avg[key] += C3[key]
                    G3_avg[key] += G3[key]
                    T3_avg[key] += T3[key]
                    A5_avg[key] += A5[key]
                    C5_avg[key] += C5[key]
                    G5_avg[key] += G5[key]
                    T5_avg[key] += T5[key]
                except:
                    A3_avg[key] = A3[key]
                    C3_avg[key] = C3[key]
                    G3_avg[key] = G3[key]
                    T3_avg[key] = T3[key]
                    A5_avg[key] = A5[key]
                    C5_avg[key] = C5[key]
                    G5_avg[key] = G5[key]
                    T5_avg[key] = T5[key]

        for k in list(A3_avg.keys()):
            A3_avg[key] = A3_avg[key]/files_count
            C3_avg[key] = C3_avg[key]/files_count
            G3_avg[key] = G3_avg[key]/files_count
            T3_avg[key] = T3_avg[key]/files_count
            A5_avg[key] = A5_avg[key]/files_count
            C5_avg[key] = C5_avg[key]/files_count
            G5_avg[key] = G5_avg[key]/files_count
            T5_avg[key] = T5_avg[key]/files_count

        self.seq_bias_avg = {"A3":A3_avg, "C3":C3_avg, "G3":G3_avg, "T3":T3_avg, "A5":A5_avg, "C5":C5_avg, "G5":G5_avg, "T5":T5_avg}


        fconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: GC Bias Distribution in first model for different samples',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.gc_first_model_ratio, fconfig) )

        sconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: GC Bias Distribution in second model for different samples',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.gc_second_model_ratio, sconfig) )

        tconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: GC Bias Distribution in third model for different samples',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.gc_third_model_ratio, tconfig) )

        avgconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Avg GC Bias Distribution for across all samples',
        'ylab': 'Average Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(modelAvg, avgconfig) )

        taprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 3\' prime end for nucleotide A',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
       'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(A3_dict, taprimeconfig) )        

        tcprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 3\' prime end for nucleotide C',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
       'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(C3_dict, tcprimeconfig) )

        tgprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 3\' prime end for nucleotide G',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
       'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(G3_dict, tgprimeconfig) )

        ttprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 3\' prime end for nucleotide T',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
       'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(T3_dict, ttprimeconfig) )

        faprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 5\' end for nucleotide A',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(A5_dict, faprimeconfig) )

        fcprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 5\' end for nucleotide C',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(C5_dict, fcprimeconfig) )

        fgprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 5\' end for nucleotide G',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(G5_dict, fgprimeconfig) )

        ftprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 5\' end for nucleotide T',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(T5_dict, ftprimeconfig) )

        seqavgconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Avg Sequential Bias for each base across all samples for both 3\' and 5\' ends',
        'ylab': 'Average Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        gcheatmapconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Average GC Bias similarity',
        'ylab': 'Average Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        seq3primeheatmappconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Average Sequential Bias (3 Prime) similarity',
        'ylab': 'Average Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        seq5sprimeheatmappconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Average Sequential Bias (5 Prime) similarity',
        'ylab': 'Average Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        self.add_section( plot = linegraph.plot(self.seq_bias_avg, seqavgconfig) )
        self.add_section( plot = heatmap.plot(self.gc_heatmap_data, self.gc_heatmap_labels,self.gc_heatmap_labels,gcheatmapconfig))
        self.add_section( plot = heatmap.plot(self.seq_3prime_heatmap_data,seq_heat_map_labels,seq_heat_map_labels,seq3primeheatmappconfig))
        self.add_section( plot = heatmap.plot(self.seq_5prime_heatmap_data,seq_heat_map_labels,seq_heat_map_labels,seq5sprimeheatmappconfig))
        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig))


        '''
        Cosine Similarity for the Average GC Bias accross the samples
        Cosine Similarity for the Average Sequential Bias (3 Prime) accross the samples
        Cosine Similarity for the Average Sequential Bias (5 Prime) accross the samples

        '''