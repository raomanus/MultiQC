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
        self.bias_path_list = []
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
                meta_json_file_path = os.path.join(os.path.dirname(f['root']),'aux_info','meta_info.json')
                gc_bias_base_dir = os.path.dirname(f['root'])
                with open(meta_json_file_path,'r') as meta_data_file:
                    meta_info_data = json.load(meta_data_file)
                self.gc_bias = meta_info_data['gc_bias_correct']
                if self.gc_bias:
                    self.bias_path_list.append(os.path.abspath(gc_bias_base_dir))


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
        for k,v in self.__dict__.items():
            print(k,':',v)
        '''

        for path_var in self.bias_path_list:
            gc_model = GCModel()
            gc_model.from_file(path_var)
            obs_array = gc_model.obs_.tolist()
            exp_array = gc_model.exp_.tolist()
            obs_weights = list(gc_model.obs_weights_)
            exp_weights = list(gc_model.exp_weights_)
            self.path_var = path_var.split('/')[-2]
            self.ratio_dict = dict()
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
                self.ratio_dict[i] = ratio_value
            rconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: GC Bias Distribution for '+ self.path_var,
            'ylab': 'Ratio (Observed/Expected)',
            'xlab': 'Read count',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
            }
            self.add_section( plot = linegraph.plot(self.ratio_dict, rconfig) )

        
        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig) )
