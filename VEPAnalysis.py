# -*- coding: utf-8 -*-
from __future__ import print_function, division

import matplotlib.pyplot as py
import matplotlib.mlab as mlab
import numpy as np
from rys_10_20 import rys
from collections import defaultdict
import csv
import os.path
import scipy.signal as ss

class VEPAnalysis(object):
	def __init__(self, signal_ind, signal_group, fs, mode='group'):
		"""
		komentarz...

		input:
		signal_ind   -- str   -- ind data   
		signal_group -- str   -- group data 

		public interface:
		draw_VEPs()
		find_p300(channel)

		"""
		super(VEPAnalysis, self).__init__()
		self.signal_ind = signal_ind
		self.signal_group = signal_group
		self.channels = self.signal_ind.keys()
		self.conditions = self.signal_ind[self.channels[0]].keys()
		self.fs = fs

		self.signal_ind_mean, self.signal_group_mean = self._estimate_mean_VEPs()
		self.signal = self._normalize_dict(mode)

		# self._find_veps_amplitude_means('Cz',self.signal_group_mean)

	def _normalize_dict(self, mode):
		signal = defaultdict(dict)
		if mode == 'group':
			signal_temp_ind = self.signal_ind_mean
			signal_temp_group = self.signal_group_mean
		else:
			signal_temp_ind = self.signal_ind
			signal_temp_group = self.signal_group
		for i in signal_temp_ind.keys():
			for j in signal_temp_ind[i].keys():
				signal[i][j] = {'individual': signal_temp_ind[i][j], 'group': signal_temp_group[i][j]}
		return signal

	def _estimate_mean_VEPs(self):
		signal_ind_mean = defaultdict(dict)
		signal_group_mean = defaultdict(dict)
		for channel in self.channels:
			for condition in self.signal_ind[channel].keys():
				signal_ind_mean[channel][condition] = [np.mean(self.signal_ind[channel][condition], axis=0),np.var(self.signal_ind[channel][condition], axis=0)/(len(self.signal_ind[channel][condition])-1)] 
				signal_group_mean[channel][condition] = [np.mean(self.signal_group[channel][condition], axis=0),np.var(self.signal_group[channel][condition], axis=0)/(len(self.signal_group[channel][condition])-1)] 
		return signal_ind_mean, signal_group_mean

	# def _write_to_file(self, erp_dict, name):
	# 	if not os.path.isfile(name):
	# 		with open(name,'wb') as f:
	# 			w = csv.writer(f)
	# 			for k,v in erp_dict.iteritems():
	# 				for k1,v1 in v.iteritems():
	# 					w.writerow([k,k1,v1[0],v1[1]])
	# 	else:
	# 		with open(name,'a') as f:
	# 			w = csv.writer(f)
	# 			for k,v in erp_dict.iteritems():
	# 				for k1,v1 in v.iteritems():
	# 					w.writerow([k,k1,v1[0],v1[1]])

	def _write_to_file_mean(self, erp_dict, name):
		if not os.path.isfile(name):
			with open(name,'wb') as f:
				w = csv.writer(f)
				for k,v in erp_dict.iteritems():
					for k1,v1 in v.iteritems():
						w.writerow([k,k1,v1[0]])
		else:
			with open(name,'a') as f:
				w = csv.writer(f)
				for k,v in erp_dict.iteritems():
					for k1,v1 in v.iteritems():
						w.writerow([k,k1,v1[0]])

	def find_VEPs_amplitudes(self,channel,write_to_file,vep):
		path = './csv/improved/'
		if vep in (['P3','P300']):
			p300 = defaultdict(dict)
			p300_ind = self._find_P3_amplitude(channel,self.signal_ind_mean)
			p300_group = self._find_P3_amplitude(channel,self.signal_group_mean)
			for condition in self.conditions:
				p300[condition] = {'individual':p300_ind[condition],'group':p300_group[condition]}
			if write_to_file:
				name = 'p300_'+channel+'_mean'+'.csv'
				self._write_to_file_mean(p300,os.path.join(path,name))
			return p300
		elif vep in (['N1','N100']):
			n100 = defaultdict(dict)
			n100_ind = self._find_N1_amplitude(channel,self.signal_ind_mean)
			n100_group = self._find_N1_amplitude(channel,self.signal_group_mean)
			for condition in self.conditions:
				n100[condition] = {'individual':n100_ind[condition],'group':n100_group[condition]}
			if write_to_file:
				name = 'n100_'+channel+'_mean'+'.csv'
				self._write_to_file_mean(n100,os.path.join(path,name))
			return n100
		elif vep in (['P1','P100']):
			p100 = defaultdict(dict)
			p100_ind = self._find_P1_amplitude(channel,self.signal_ind_mean)
			p100_group = self._find_P1_amplitude(channel,self.signal_group_mean)
			for condition in self.conditions:
				p100[condition] = {'individual':p100_ind[condition],'group':p100_group[condition]}
			if write_to_file:
				name = 'p100_'+channel+'_mean'+'.csv'
				self._write_to_file_mean(p100,os.path.join(path,name))
			return p100

	def _find_P1_amplitude(self,channel,signal):
		P1_amplitudes = defaultdict(dict)
		for condition in self.conditions:
			ids = self._find_extrema(channel,signal,condition)
			amps = []
			for i in xrange(len(ids)):
				if (ids[i-1][0] == 'max' and ids[i][0] == 'min'):
					if (ids[i-1][1] < 150 and ids[i-1][1] > 50):
						# amps.append([abs(signal[channel][condition][0][ids[i][1]]-signal[channel][condition][0][ids[i-1][1]]),
						# 	    	 signal[channel][condition][1][ids[i][1]]+signal[channel][condition][1][ids[i-1][1]]])
						amps.append(abs(np.mean(signal[channel][condition][0][ids[i-1][1]-0.025*self.fs:ids[i-1][1]+0.025*self.fs])))			
						break
			if len(amps):
				P1_amplitudes[condition] = amps
			else:
				P1_amplitudes[condition] = [0]
			# print ('***************',condition,P1_amplitudes[condition])
			# py.figure()
			# py.plot(signal[channel][condition][0])
			# py.vlines(ids[i-1][1]-0.025*self.fs,-6,10,colors='m')
			# py.vlines(ids[i-1][1]+0.025*self.fs,-6,10,colors='m')
			# py.xlim(0,0.5*self.fs)
			# py.show()
		return P1_amplitudes

	def _find_P3_amplitude(self,channel,signal):
		'''
			returns mean signal value 60 ms around P300 peak 
		'''
		P3_amplitudes = defaultdict(dict)
		for condition in self.conditions:
			ids = self._find_extrema(channel,signal,condition)
			amps = []
			for i in xrange(len(ids)):
				if i != 0:
					if (ids[i-1][0] == 'min' and ids[i][0] == 'max'):
						if (ids[i][1] < 380 and ids[i][1] > 250):
							amps.append(np.mean(signal[channel][condition][0][ids[i][1]-0.03*self.fs:ids[i][1]+0.03*self.fs]))
							break
			if len(amps):
				P3_amplitudes[condition] = amps
			else:
				P3_amplitudes[condition] = [0]
			# py.figure()
			# py.plot(signal[channel][condition][0])
			# py.vlines(ids[i][1]-0.03*self.fs,-6,10,colors='m')
			# py.vlines(ids[i][1]+0.03*self.fs,-6,10,colors='m')
			# py.xlim(0,0.5*self.fs)
			# print ('***************',condition,P3_amplitudes[condition])
			# py.show()
		return P3_amplitudes

	def _find_N1_amplitude(self,channel,signal):
		N1_amplitudes = defaultdict(dict)
		for condition in self.conditions:
			ids = self._find_extrema(channel,signal,condition)
			amps = []
			for i in xrange(len(ids)):
				if (ids[i-1][0] == 'max' and ids[i][0] == 'min'):
					if (ids[i-1][1] < 150 and ids[i-1][1] > 50):
						# amps.append([abs(signal[channel][condition][0][ids[i][1]]-signal[channel][condition][0][ids[i-1][1]]),
						# 	    	 signal[channel][condition][1][ids[i][1]]+signal[channel][condition][1][ids[i-1][1]]])
						amps.append(abs(np.mean(signal[channel][condition][0][ids[i][1]-0.025*self.fs:ids[i][1]+0.025*self.fs])))			
						break
			if len(amps):
				N1_amplitudes[condition] = amps
			else:
				N1_amplitudes[condition] = [0]
			# print ('***************',condition,N1_amplitudes[condition])
			# py.figure()
			# py.plot(signal[channel][condition][0])
			# py.vlines(ids[i][1]-0.025*self.fs,-6,10,colors='m')
			# py.vlines(ids[i][1]+0.025*self.fs,-6,10,colors='m')
			# py.xlim(0,0.5*self.fs)
			# py.show()
		return N1_amplitudes

	def _find_extrema(self,channel,signal,condition):
		deriv = np.gradient(np.array(signal[channel][condition][0][:0.5*self.fs]))
		ids = []
		for i in xrange(len(deriv)):
			if not (i == 0 or i == 1):
				if (deriv[i-1] > 0 and deriv[i] < 0):
					ids.append(('max',i-1))
				elif (deriv[i-1] < 0 and deriv[i] > 0):
					ids.append(('min',i-1))
		return ids

	def draw_VEPs(self, channels, signal):
		modes = ['go', 'no-go']
		t = np.linspace(-0.5,1,len(signal['Cz']['compatible-go']['individual'][0]))
		if isinstance(channels, list):
			for mode in modes:
				rys_go = rys(mode+' trials')
				for ch in channels:
					r1 = rys_go[ch].plot(t, signal[ch]['compatible-'+mode]['individual'][0], 'r')
					# rys_go[ch].fill_between(t, signal[ch]['compatible-'+mode]['individual'][0]+signal[ch]['compatible-'+mode]['individual'][1],
					# 						   signal[ch]['compatible-'+mode]['individual'][0]-signal[ch]['compatible-'+mode]['individual'][1],
					# 						   facecolor='r', alpha=0.5)
					r2 = rys_go[ch].plot(t, signal[ch]['compatible-'+mode]['group'][0], 'g')
					# rys_go[ch].fill_between(t, signal[ch]['compatible-'+mode]['group'][0]+signal[ch]['compatible-'+mode]['group'][1],
					# 						   signal[ch]['compatible-'+mode]['group'][0]-signal[ch]['compatible-'+mode]['group'][1],
					# 						   facecolor='g', alpha=0.5)
					r3 = rys_go[ch].plot(t, signal[ch]['incompatible-'+mode]['individual'][0], 'c')
					# rys_go[ch].fill_between(t, signal[ch]['incompatible-'+mode]['individual'][0]+signal[ch]['incompatible-'+mode]['individual'][1],
					# 						   signal[ch]['incompatible-'+mode]['individual'][0]-signal[ch]['incompatible-'+mode]['individual'][1],
					# 						   facecolor='c', alpha=0.5)
					r4 = rys_go[ch].plot(t, signal[ch]['incompatible-'+mode]['group'][0], 'm')
					# rys_go[ch].fill_between(t, signal[ch]['incompatible-'+mode]['group'][0]+signal[ch]['incompatible-'+mode]['group'][1],
					# 						   signal[ch]['incompatible-'+mode]['group'][0]-signal[ch]['incompatible-'+mode]['group'][1],
					# 						   facecolor='m', alpha=0.5)
					for r in [r1, r2]:
						py.setp(r,linewidth=2)
					rys_go[ch].set_xlim(-0.1,0.6)
					rys_go[ch].set_ylim(-6,16)
					rys_go[ch].axvline(x=0,color='r')
				rys_go['legend'].plot([],[],'r',label='Indywidualny zgodny')
				rys_go['legend'].plot([],[],'g',label='Grupowy zgodny')
				rys_go['legend'].plot([],[],'c',label='Indywidualny niezgodny')
				rys_go['legend'].plot([],[],'m',label='Grupowy niezgodny')
				rys_go['legend'].legend()
				# py.savefig('./rys/'+mode+'.png')#, bbox_inches='tight')
		else: 
			for mode in modes:
				fig = py.figure()
				py.title(mode+' trials')
				r1 = py.plot(t, signal[channels]['compatible-'+mode]['individual'][0], 'r', label='Indywidualny zgodny')
				# py.fill_between(t, signal[channels]['compatible-'+mode]['individual'][0]+signal[channels]['compatible-'+mode]['individual'][1],
				# 				signal[channels]['compatible-'+mode]['individual'][0]-signal[channels]['compatible-'+mode]['individual'][1],
				# 				facecolor='r', alpha=0.5)
				r2 = py.plot(t, signal[channels]['compatible-'+mode]['group'][0], 'g', label='Grupowy zgodny')
				# py.fill_between(t, signal[channels]['compatible-'+mode]['group'][0]+signal[channels]['compatible-'+mode]['group'][1],
				# 				signal[channels]['compatible-'+mode]['group'][0]-signal[channels]['compatible-'+mode]['group'][1],
				# 				facecolor='g', alpha=0.5)
				r3 = py.plot(t, signal[channels]['incompatible-'+mode]['individual'][0], 'c', label='Indywidualny niezgodny')
				# py.fill_between(t, signal[channels]['incompatible-'+mode]['individual'][0]+signal[channels]['incompatible-'+mode]['individual'][1],
				# 				signal[channels]['incompatible-'+mode]['individual'][0]-signal[channels]['incompatible-'+mode]['individual'][1],
				# 				facecolor='c', alpha=0.5)
				r4 = py.plot(t, signal[channels]['incompatible-'+mode]['group'][0], 'm', label='Grupowy niezgodny')
				# py.fill_between(t, signal[channels]['incompatible-'+mode]['group'][0]+signal[channels]['incompatible-'+mode]['group'][1],
				# 				signal[channels]['incompatible-'+mode]['group'][0]-signal[channels]['incompatible-'+mode]['group'][1],
				# 				facecolor='m', alpha=0.5)
				for r in [r1, r2]:
					py.setp(r,linewidth=2)
				py.xlim(0,0.6)
				py.legend()
		py.show()