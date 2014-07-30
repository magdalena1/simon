# -*- coding: utf-8 -*-
from __future__ import print_function, division

import matplotlib.pyplot as py
import numpy as np
from collections import defaultdict
import os.path

class MPDataSerialization(object):
	def __init__(self, signal, fs, condition):

		super(MPDataSerialization, self).__init__()
		self.fs = fs
		self.condition = condition
		self.signal = signal
		length = len(self.signal['Pz']['compatible-go'][0])

		self.downsampled_signal = self._downsample_signal(length)

	def _downsample_signal(self, length, factor=4):
		''' 
		downsampling signal for binary serialization 
		'''
		self.length = int(length/factor)
		downsampled_signal = defaultdict(list)
		for channel in self.signal.keys():
			for trial in xrange(len(self.signal[channel][self.condition])):
				n = []
				for i in xrange(int(length)):
					if i%factor == 0:
						n.append(self.signal[channel][self.condition][trial][i])
				downsampled_signal[channel].append(n)
		return downsampled_signal

	def _downsample_averaged_signal(self, factor=8):
		''' 
		downsampling averaged signal
		'''
		self.length = int(self.fs/factor)
		downsampled_signal = defaultdict(dict)
		for channel in self.signal.keys():
			for condition in self.signal[channel].keys():
				for mode in self.signal[channel][condition].keys():
					n = []
					for i in xrange(int(self.fs)):
						if i%factor == 0:
							n.append(self.signal[channel][condition][mode][0][i])
					try:
						downsampled_signal[channel][condition][mode] = n
					except KeyError:
						downsampled_signal[channel][condition] = defaultdict(list)
						downsampled_signal[channel][condition][mode] = n
		return downsampled_signal

	def write_binary_MMP11(self, file_name, trials):
		'''
			MMP11 - multitrial, multichannel MP
		'''
		x = np.zeros([self.length,len(self.downsampled_signal),trials], dtype='<f')
		with open(file_name,'wb') as f:
			for trial in xrange(trials):
				for sample in xrange(128):
					for i,channel in enumerate(self.downsampled_signal.keys()):
						x[sample,i,trial] = self.downsampled_signal[channel][trial][sample]
						f.write(x[sample,i,trial])
		print('MMP11, downsampled: ',x.shape)
		return x

	def write_binary_MMP1_multitrial(self, file_name, channel, trials):
		'''
			MMP1 - multichannel MP (multitrial version) - sum of products maximized in all channels 
				   (only amplitudes varies across channels)
		'''
		x = np.zeros([self.length,trials], dtype='<f')
		for trial in xrange(trials):
			x[:,i] = self.downsampled_signal[channel][trial]
		with open(file_name,'wb') as f:
			x.tofile(f)
		print('MMP1 (multitrial), downsampled: ',x.shape)
		return x

	def write_binary_MMP1(self, file_name, trial):
		'''
			MMP1 - multichannel MP (for single trial) - sum of products maximized in all channels 
				   (only amplitudes varies across channels)
		'''
		x = np.zeros([self.length,len(self.downsampled_signal)], dtype='<f')
		for channel in self.downsampled_signal.keys():
			x[:,i] = self.downsampled_signal[channel][trial]
		with open(file_name,'wb') as f:
			x.tofile(f)
		return x

	def write_binary_MMP1_averaged_trials(self, file_name, trials, channel, averaged_nr):
		'''
			MMP1 - multichannel MP (for averaged trials) - sum of products maximized in all channels 
				   (only amplitudes varies across channels)
		'''
		x = np.zeros([self.length,trials/averaged_nr], dtype='<f')
		# x = np.zeros([int(0.45*self.length)-int(0.2*self.length),trials/averaged_nr], dtype='<f')
		mean_x = []
		i = 0
		for k in xrange(trials+1):
			mean_x.append(self.downsampled_signal[channel][k])
			if (averaged_nr != 1) and (k%averaged_nr == 0):
				if k != 0:
					x[:,i] = np.mean(mean_x, axis=0)
					# py.figure()
					# py.plot(x[:,i])
					i += 1
					mean_x = []
			elif averaged_nr == 1:
				try:
					x[:,k] = self.downsampled_signal[channel][k]#[int(0.2*self.length):int(0.45*self.length)]
				except IndexError:
					break
		with open(file_name,'wb') as f:
			x.tofile(f)
		# py.show()
		return x