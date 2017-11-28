#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os

from multiqc import config
from multiqc.plots import linegraph
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
        self.seq_three_prime = dict()
        self.seq_five_prime = dict()

        for path_var in self.gc_bias_path_list:
            gc_model = GCModel()
            gc_model.from_file(path_var)
            obs_array = gc_model.obs_.tolist()
            exp_array = gc_model.exp_.tolist()
            obs_weights = list(gc_model.obs_weights_)
            exp_weights = list(gc_model.exp_weights_)
            self.path_var = path_var.split('/')[-2]
            ratio_dict = dict()
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
                    j += 1
                ratio_dict[i] = ratio_value
            self.gc_first_model_ratio[self.path_var] = ratio_dict[0]
            self.gc_second_model_ratio[self.path_var] = ratio_dict[1]
            self.gc_third_model_ratio[self.path_var] = ratio_dict[2]

        self.seq_3prime_ratio = dict()
        self.seq_5prime_ratio = dict()

        for path_var in self.seq_bias_path_list:
            seq_model = SEQModel()
            seq_model.from_file(path_var)
            obs3_array = seq_model.obs3_prime.tolist()
            exp3_array = seq_model.exp3_prime.tolist()
            obs5_array = seq_model.obs5_prime.tolist()
            exp5_array = seq_model.exp5_prime.tolist()
            self.path_var = path_var.split('/')[-2]
            ratio_dict_3prime = OrderedDict()
            i = 1
            for o,e in zip(obs3_array, exp3_array):
                ratio = o/e
                ratio_dict_3prime[i] = ratio
                i += 1

            ratio_dict_5prime = OrderedDict()
            j = 1
            for o,e in zip(obs5_array, exp5_array):
                ratio = o/e
                ratio_dict_5prime[j] = ratio
                j += 1

            self.seq_3prime_ratio[self.path_var] = ratio_dict_3prime
            self.seq_5prime_ratio[self.path_var] = ratio_dict_5prime


        fconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: GC Bias Distribution in first model for different experiments',
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
        'title': 'Salmon: GC Bias Distribution in second model for different experiments',
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
        'title': 'Salmon: GC Bias Distribution in third model for different experiments',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.gc_third_model_ratio, tconfig) )

        tprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 3\' end',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.seq_3prime_ratio, tprimeconfig) )        
        
        fprimeconfig = {
        'smooth_points': 500,
        'id': 'salmon_plot',
        'title': 'Salmon: Sequence Bias Distribution for different experiments measured from 5\' end',
        'ylab': 'Ratio (Observed/Expected)',
        'xlab': 'Read count',
        'ymin': 0,
        'xmin': 0,
        'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.seq_5prime_ratio, fprimeconfig) )


        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig) )
