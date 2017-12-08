import pandas as pd
import numpy as np

class SEQModel:
    def __init__(self):
        self.obs3_prime = None
        self.exp3_prime = None
        self.obs5_prime = None
        self.exp5_prime = None
        self.dims_ = None
        self.valid_ = False

    def populate_model_(self, data_):
        import struct
        from numpy.linalg import norm

        model = None
        offset = 0
        int_struct = struct.Struct('@i')
        long_struct = struct.Struct('@q')

        context_length = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        left_context = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        right_context = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        offset += 3*(context_length*int_struct.size)

        vlmm_rows = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        vlmm_cols = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        vlmm_struct = struct.Struct('@' + vlmm_rows * vlmm_cols * 'd')
        offset += vlmm_struct.size

        nrows = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        ncols = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        model_struct = struct.Struct('@' + nrows * ncols * 'd')
        model = model_struct.unpack_from(data_[offset:])
        model = np.array(model)
        model = model.reshape(ncols,nrows).T
        return model
        
    '''
        dname is the root directory of salmon output
        Search for observed and expected gc files inside aux_info directory under dname.
        Read those files and populate the two matrices of order 3 X 25 and two weighted matrices of order 1 X 3
    '''
    def from_file(self, dname):
        import os
        import gzip
        obs3_name = os.path.sep.join([dname, 'aux_info', 'obs3_seq.gz'])
        exp3_name = os.path.sep.join([dname, 'aux_info', 'exp3_seq.gz'])
        obs5_name = os.path.sep.join([dname, 'aux_info', 'obs5_seq.gz'])
        exp5_name = os.path.sep.join([dname, 'aux_info', 'exp5_seq.gz'])
        obs_dat, exp_dat = None, None
        try:
            with gzip.open(obs3_name) as obs_file:
                obs_dat = obs_file.read()
            self.obs3_prime = self.populate_model_(obs_dat)
        except IOError:
            print("Could not open file {}".format(obs3_name))
            return False

        try:
            with gzip.open(exp3_name) as exp_file:
                exp_dat = exp_file.read()
            self.exp3_prime = self.populate_model_(exp_dat)
        except IOError:
            print("Could not open file {}".format(exp3_name))
            return False
            
        obs_dat, exp_dat = None, None
        try:
            with gzip.open(obs5_name) as obs_file:
                obs_dat = obs_file.read()
            self.obs5_prime = self.populate_model_(obs_dat)
        except IOError:
            print("Could not open file {}".format(obs5_name))
            return False

        try:
            with gzip.open(exp5_name) as exp_file:
                exp_dat = exp_file.read()
            self.exp5_prime = self.populate_model_(exp_dat)
        except IOError:
            print("Could not open file {}".format(exp5_name))
            return False

        self.valid_ = True
        return True