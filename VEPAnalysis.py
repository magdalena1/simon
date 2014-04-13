# -*- coding: utf-8 -*-
import matplotlib.pyplot as py
import numpy as np
from rys_10_20 import rys
from collections import defaultdict

class VEPAnalysis(object):
	def __init__(self, signal_ind, signal_group, fs):
		"""
		komentarz...

		input:
		signal_ind   -- str   -- ind data   
		signal_group -- str   -- group data 

		public interface:
		draw_VEPs()

		"""
		super(VEPAnalysis, self).__init__()
		self.signal_ind = signal_ind
		self.signal_group = signal_group
		self.channels = self.signal_ind.keys()
		self.conditions = self.signal_ind[self.channels[0]].keys()
		self.fs = fs

		self.signal_ind_mean, self.signal_group_mean = self._estimate_mean_VEPs()
		# self.draw_VEPs(self.channels)

	def _estimate_mean_VEPs(self):
		signal_ind_mean = defaultdict(dict)
		signal_group_mean = defaultdict(dict)
		for channel in self.channels:
			for condition in self.signal_ind[channel].keys():
				signal_ind_mean[channel][condition] = [np.mean(self.signal_ind[channel][condition], axis=0),np.var(self.signal_ind[channel][condition], axis=0)/len(self.signal_ind[channel][condition])] 
				signal_group_mean[channel][condition] = [np.mean(self.signal_group[channel][condition], axis=0),np.var(self.signal_group[channel][condition], axis=0)/len(self.signal_group[channel][condition])] 
		return signal_ind_mean, signal_group_mean

	def find_p300(self, channel):
		p300 = dict()
		p300_ind = self._find_p300_amplitude(channel, self.signal_ind_mean)
		p300_group = self._find_p300_amplitude(channel, self.signal_group_mean)
		for condition in self.conditions:
			p300[condition] = [p300_ind[condition], p300_group[condition]]
		return p300

	def _find_p300_amplitude(self, channel, signal):
		P3_amplitudes = defaultdict(dict)
		for condition in self.conditions:
			deriv = np.gradient(np.array(signal[channel][condition][0][:0.5*self.fs]))
			ids = []
			amps = []
			for i in xrange(len(deriv)):
				if i != (0,1):
					if (deriv[i-1] > 0 and deriv[i] < 0):
						ids.append(('max',i-1))
					elif (deriv[i-1] < 0 and deriv[i] > 0):
						ids.append(('min',i-1))
			for i in xrange(len(ids)):
				if i != 0:
					if (ids[i-1][0] == 'min' and ids[i][0] == 'max'):
						amps.append([abs(signal[channel][condition][0][ids[i][1]]-signal[channel][condition][0][ids[i-1][1]]),
								    signal[channel][condition][1][ids[i][1]]+signal[channel][condition][1][ids[i-1][1]]])
			max_id = np.array(amps).argmax(axis=0)[0]
			P3_amplitudes[condition] = amps[max_id]
		return P3_amplitudes

	def draw_VEPs(self, channels):
		modes = ['go', 'no-go']
		t = np.linspace(0,1,self.fs)
		if isinstance(channels, list):
			for mode in modes:
				rys_go = rys(mode+' trials')
				for ch in channels:
					r1 = rys_go[ch].plot(t, self.signal_ind_mean[ch]['compatible-'+mode][0], 'r')
					rys_go[ch].fill_between(t, self.signal_ind_mean[ch]['compatible-'+mode][0]+self.signal_ind_mean[ch]['compatible-'+mode][1],
											   self.signal_ind_mean[ch]['compatible-'+mode][0]-self.signal_ind_mean[ch]['compatible-'+mode][1],
											   facecolor='r', alpha=0.5)
					r2 = rys_go[ch].plot(t, self.signal_group_mean[ch]['compatible-'+mode][0], 'g')
					rys_go[ch].fill_between(t, self.signal_group_mean[ch]['compatible-'+mode][0]+self.signal_group_mean[ch]['compatible-'+mode][1],
											   self.signal_group_mean[ch]['compatible-'+mode][0]-self.signal_group_mean[ch]['compatible-'+mode][1],
											   facecolor='g', alpha=0.5)
					r3 = rys_go[ch].plot(t, self.signal_ind_mean[ch]['incompatible-'+mode][0], 'c')
					rys_go[ch].fill_between(t, self.signal_ind_mean[ch]['incompatible-'+mode][0]+self.signal_ind_mean[ch]['incompatible-'+mode][1],
											   self.signal_ind_mean[ch]['incompatible-'+mode][0]-self.signal_ind_mean[ch]['incompatible-'+mode][1],
											   facecolor='c', alpha=0.5)
					r4 = rys_go[ch].plot(t, self.signal_group_mean[ch]['incompatible-'+mode][0], 'm')
					rys_go[ch].fill_between(t, self.signal_group_mean[ch]['incompatible-'+mode][0]+self.signal_group_mean[ch]['incompatible-'+mode][1],
											   self.signal_group_mean[ch]['incompatible-'+mode][0]-self.signal_group_mean[ch]['incompatible-'+mode][1],
											   facecolor='m', alpha=0.5)
					for r in [r1, r2]:
						py.setp(r,linewidth=2)
					rys_go[ch].set_xlim(0,0.6)
					# rys_go[ch].set_ylim(-6,10)
				rys_go['legend'].plot([],[],'r',label='Individual Compatible')
				rys_go['legend'].plot([],[],'g',label='Group Compatible')
				rys_go['legend'].plot([],[],'c',label='Individual Incompatible')
				rys_go['legend'].plot([],[],'m',label='Group Incompatible')
				rys_go['legend'].legend()
				# py.savefig('./rys/'+mode+'.png')#, bbox_inches='tight')
		else:
			for mode in modes:
				fig = py.figure()
				py.title(mode+' trials')
				r1 = py.plot(t, self.signal_ind_mean[channels]['compatible-'+mode][0], 'r', label='Individual Compatible')
				py.fill_between(t, self.signal_ind_mean[channels]['compatible-'+mode][0]+self.signal_ind_mean[channels]['compatible-'+mode][1],
								self.signal_ind_mean[channels]['compatible-'+mode][0]-self.signal_ind_mean[channels]['compatible-'+mode][1],
								facecolor='r', alpha=0.5)
				r2 = py.plot(t, self.signal_group_mean[channels]['compatible-'+mode][0], 'g', label='Group Compatible')
				py.fill_between(t, self.signal_group_mean[channels]['compatible-'+mode][0]+self.signal_group_mean[channels]['compatible-'+mode][1],
								self.signal_group_mean[channels]['compatible-'+mode][0]-self.signal_group_mean[channels]['compatible-'+mode][1],
								facecolor='g', alpha=0.5)
				r3 = py.plot(t, self.signal_ind_mean[channels]['incompatible-'+mode][0], 'c', label='Individual Incompatible')
				py.fill_between(t, self.signal_ind_mean[channels]['incompatible-'+mode][0]+self.signal_ind_mean[channels]['incompatible-'+mode][1],
								self.signal_ind_mean[channels]['incompatible-'+mode][0]-self.signal_ind_mean[channels]['incompatible-'+mode][1],
								facecolor='c', alpha=0.5)
				r4 = py.plot(t, self.signal_group_mean[channels]['incompatible-'+mode][0], 'm', label='Group Incompatible')
				py.fill_between(t, self.signal_group_mean[channels]['incompatible-'+mode][0]+self.signal_group_mean[channels]['incompatible-'+mode][1],
								self.signal_group_mean[channels]['incompatible-'+mode][0]-self.signal_group_mean[channels]['incompatible-'+mode][1],
								facecolor='m', alpha=0.5)
				for r in [r1, r2]:
					py.setp(r,linewidth=2)
				py.legend()
		py.show()