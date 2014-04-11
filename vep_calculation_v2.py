# -*- coding: utf-8 -*-
from __future__ import print_function, division
# from posturograph import *
from obci.analysis.obci_signal_processing import read_manager 
from obci.analysis.obci_signal_processing.tags import tags_file_reader as TFR

# import scipy.stats as stats
import matplotlib.pyplot as py
import numpy as np
from rys_10_20 import rys
import os.path, ast
from scipy.signal import filtfilt, cheby2, butter, resample
import random
from collections import defaultdict

class VEPAnalysis(object):
	def __init__(self, file_path, f_name_ind, f_name_group):
		"""
		format zwracanych danych frags['O1']['compatible-go'][[1024],[1024],[1024],...]
		warunki ['compatible-go','compatible-no-go','incompatible-go','incompatible-no-go']

		input:
		file_path
		f_name_ind   -- str   -- individual data file's name (a common part of the name for info, 
								tags and raw files)   
		f_name_group -- str   -- group data file's name (a common part of the name for info, 
								tags and raw files) 

		public interface:

		"""
		super(VEPAnalysis, self).__init__()
		self.channels = [u'Fp1',u'Fpz',u'Fp2',u'F7',u'F3',u'Fz',u'F4',u'F8',u'T3',u'C3',u'Cz',u'C4',u'T4',u'T5', 
		  u'P3', u'Pz', u'P4', u'T6', u'O1', u'Oz', u'O2']

		self.frags_group = self._cut_fragments(os.path.join(os.path.normpath(file_path), f_name_group))
		self.frags_ind = self._cut_fragments(os.path.join(os.path.normpath(file_path), f_name_ind))

	def _get_current_condition(self, file_name):
		f = os.path.basename(file_name).split('_')
		return 'blue' if f[2].lower() == 'r' or f[2].lower() == 'rb' else 'green'

	def _linked_ears_montage(self):
		sig = self.mgr.get_samples()
		ref_sig = []
		for i in xrange(23):
			new_sig = sig[i]-0.5*(self.mgr.get_channel_samples('A1') + self.mgr.get_channel_samples('A2'))
			ref_sig.append(new_sig - np.mean(new_sig))
		filtered_sig = self._filter_signal(ref_sig)
		return filtered_sig

	def _filter_signal(self, signal):
		sig_f = []
		fn = self.fs/2.
		# b,a = cheby2(1, 10, 0.5/fn, btype = "highpass")
		# b,a = butter(2, np.array([49,51])/fn, btype = "bandstop")
		# d,c = butter(2, 0.1/fn, btype = "highpass")
		d,c = butter(2, np.array([1,12])/fn, btype = "bandpass")
		for i in xrange(len(signal)):
			# sig_f.append(filtfilt(b,a,filtfilt(d, c, signal[i])))
			sig_f.append(filtfilt(d, c, signal[i]))
		return sig_f

	def _get_timestamps(self, tag_file):
		proper_key = self._get_proper_key()
		other_condition = 'green' if self.condition == 'blue' else 'blue'
		keys = ['green-left', 'green-right', 'blue-left', 'blue-right']
		times = {}
		for k in keys:
			times[k] = []
		reader = TFR.TagsFileReader(tag_file)
		tags = reader.get_tags()
		proper_cond = [key for key in keys if self.condition in key]
		wrong_cond = list(set(keys) - set(proper_cond))
		for tag in tags:
			for c in wrong_cond:
				if tag['name'] == c:
					if isinstance(tag['desc']['key'], basestring) and tag['desc']['key'] != proper_key:
						times[c].append(float(tag['start_timestamp']))
					elif isinstance(tag['desc']['key'], list):
						k = ast.literal_eval(tag['desc']['key'])
						if proper_key not in k:
							times[c].append(float(tag['start_timestamp']))
			for c in proper_cond:
				if tag['name'] == c:
					if isinstance(tag['desc']['key'], basestring) and tag['desc']['key'] == proper_key:
						times[c].append(float(tag['start_timestamp']))
					elif isinstance(tag['desc']['key'], list):
						k = ast.literal_eval(tag['desc']['key'])
						if proper_key in k:
							times[c].append(float(tag['start_timestamp']))				
		return times

	def _get_proper_key(self):
		if self.condition == 'blue':
			proper_key = 'right'  
		else: 
			proper_key = 'left'
		return proper_key

	def _normalize_keys(self, file_name, frags):
		proper_key = self._get_proper_key()
		other_key = 'left' if self.condition == 'blue' else 'right'
		other_condition = 'green' if self.condition == 'blue' else 'blue'
		frags_norm = defaultdict(dict)
		for channel in frags.keys():
			for key in frags[channel].keys():
				if key == self.condition+'-'+proper_key:
					frags_norm[channel]['compatible-go'] = frags[channel][key]
				elif key == self.condition+'-'+other_key:
					frags_norm[channel]['incompatible-go'] = frags[channel][key]
				elif key == other_condition+'-'+proper_key:
					frags_norm[channel]['compatible-no-go'] = frags[channel][key]
				elif key == other_condition+'-'+other_key:
					frags_norm[channel]['incompatible-no-go'] = frags[channel][key]
		return frags_norm

	def _cut_fragments(self, file_name):
		'''
		format zwracanych danych frags['O1']['blue-right'][[1024],[1024],[1024],...]
		'''
		self.condition = self._get_current_condition(file_name)
		self.mgr = read_manager.ReadManager(file_name+'.xml',
									   file_name+'.raw',
									   file_name+'.tag')
		channels = self.mgr.get_param('channels_names')
		pointsPerMikroV = self.mgr.get_param('channels_gains')
		self.fs = float(self.mgr.get_param('sampling_frequency'))
		# sig = self.mgr.get_samples()		### uncomment if signal is already filtered
		sig = self._linked_ears_montage()
		# sig = _one_ear_montage(self.mgr, self.fs)
		times = self._get_timestamps(file_name+'.tag')
		stim = {}
		for key in times.keys():
			stim[key] = [t*self.fs for t in times[key]]
		frags = {}
		for ch in channels[:23]:
			frags[ch] = {}
			for key in times.keys():
				frags[ch][key] = []
		for i in xrange(len(sig[:23])):
			sig[i] = sig[i]*float(pointsPerMikroV[i])
			for key in stim.keys():
				for j in stim[key]:
					frags[channels[i]][key].append(sig[i][j:j+int(self.fs)]-np.mean(sig[i][j-int(0.1*self.fs):j+int(self.fs)]))
		frags.pop('A1', None)
		frags.pop('A2', None)
		frags = self._normalize_keys(file_name,frags)
		return frags

if __name__ == '__main__':
	path = './badania_part3/'
	# f_list_group = ['21_K_RB_01']
	# f_list_ind = ['21_K_R_01']
	v = VEPAnalysis(path,'23_K_R_01','23_K_RB_01')