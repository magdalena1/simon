#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import numpy as np
import collections
import matplotlib.pyplot as py
import matplotlib.gridspec as gridspec
import sys
from scipy.stats.stats import pearsonr
from scipy import stats
from operator import itemgetter

class BookImporter(object):
	def __init__(self, book_file, id_min=1, id_max=1, t_min=0, t_max=1, atom_min=0, atom_max=5):
		'''
		input arguments:
			book_file
			id_min 		- index of minimal trial to plot
			id_max		- index of maximal trial to plot
			t_min		- down edge of time period to reconstruct
			t_max		- upper edge of time period to reconstruct
			atom_min	- minimal atom to plot
			atom_max	- maximal atom to plot

		comment:
			if we want to reconstruct only one atom then atom_min should equal atom_max
		'''
		super(BookImporter, self).__init__()

		f = open(book_file,'rb')
		data, signals, atoms, epoch_s = self._read_book(f)
		self.epoch_s = epoch_s
		self.atoms = atoms
		self.signals = signals
		self.fs = data[5]['Fs']
		self.ptspmV = 0.0715
		self.id_min = id_min
		self.id_max = id_max
		self.t_min = t_min
		self.t_max = t_max
		self.atom_min = atom_min
		self.atom_max = atom_max

		amps = self._get_atoms_peak2peak_amplitudes()
		
		# self.draw_reconstructed_signal()
		# self.draw_mean_reconstructed_signal(atoms, signals)

		# amps = self._get_atoms_amplitudes()
		self.amps = amps
		# print(np.median(amps))
		# self.perform_linear_regression(amps)

		# t, f, E_a, sigt, signal_a, signal_reconstruction_a, alpha_mean = self.calculate_mean_map(df = 0.2, dt = 0.008, f_a = [2, 20])
		# self.alpha_mean = alpha_mean
		# print(pearsonr(amps,alpha_mean)) ### (Pearson’s correlation coefficient, 2-tailed p-value)
		# self.perform_linear_regression(alpha_mean)
		# print(np.median(alpha_mean))

		# t, f, E_a,sigt, signal_a, signal_reconstruction_a = self.calculate_map(self.atoms[2],self.signals[1][1][1], df = 0.2, dt = 0.008, f_a = [2, 20])

		# self.only_draw_map(t, f, E_a, sigt, signal_a, signal_reconstruction_a)

		py.show()

	def _get_type(self, ident, f):
		'''generacja dtype dla odczytania fragmentów, odpalana od razu,
		po odczytaniu identyfikatora'''
		# print(ident)
		if ident == 1:
			com_s = np.fromfile(f, '>u4', count=1)[0]
			if not com_s==0: ## comment
				return np.dtype([('comment', 'S'+str(com_s))])
			else:
				return None
				
		elif ident == 2: ## header
			head_s = np.fromfile(f, '>u4', count=1)[0]
			return None
			
		elif ident == 3: ## www adress
			www_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('www', 'S'+str(www_s))])
			
		elif ident == 4: ## data stworzenia pliku
			date_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('date', 'S'+str(date_s))])
			
		elif ident == 5: ## info o sygnale
			sig_info_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('Fs', '>f4'), ('ptspmV', '>f4'),
							('chnl_cnt', '>u2')])
							
		elif ident == 6: ## info o dekompozycji
			dec_info_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('percent', '>f4'), ('maxiterations', '>u4'),
							('dict_size', '>u4'), ('dict_type', '>S1')])
							
		elif ident == 10: #dirac
			return
			atom_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('modulus', '>f4'), ('amplitude', '>f4'),
							('t', '>f4')])
							
		elif ident == 11: #gauss
			atom_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('modulus', '>f4'), ('amplitude', '>f4'),
							('t', '>f4'), ('scale', '>f4')])
							
		elif ident == 12: #sinus
			return
			atom_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('modulus', '>f4'), ('amplitude', '>f4'),
							('f', '>f4'), ('phase', '>f4')])
							
		elif ident == 13: #gabor
			atom_s = np.fromfile(f, '>u1', count=1)[0]
			return np.dtype([('modulus', '>f4'), ('amplitude', '>f4'),
							('t', '>f4'), ('scale', '>f4'),
							('f', '>f4'), ('phase', '>f4')])
							
		else:
			return None
			
	def _get_signal(self, f, epoch_nr, epoch_s):
		'''uruchamiana przy odnalezieniu bloku sygnału'''
		sig_s = np.fromfile(f, '>u4', count=1)[0]
		chnl_nr = np.fromfile(f, '>u2', count=1)[0]
		signal = np.fromfile(f, '>f4', count= epoch_s)
		return chnl_nr, signal
		
	def _get_atoms(self, f):
		'''uruchamiana przy odnalezieniu bloku atomów'''
		atoms = collections.defaultdict(list)
		atoms_s = np.fromfile(f, '>u4', count=1)[0]
		a_chnl_nr = np.fromfile(f, '>u2', count=1)[0]
		ident = np.fromfile(f, '>u1', count=1)
		while ident in [10, 11, 12, 13]:
			atom = np.fromfile(f, self._get_type(ident[0], f), count=1)[0]
			atoms[ident[0]].append(atom)
			ident = np.fromfile(f, '>u1', count=1)
		f.seek(f.tell()-1)
		return atoms, a_chnl_nr
		
	def _gabor(self, amplitude, position, scale, afrequency, phase):
		time = np.linspace(0, self.epoch_s/self.fs, self.epoch_s)
		width = scale
		frequency = afrequency*2*np.pi
		signal = amplitude*np.exp(-np.pi*((time-position)/width)**2)*np.cos(frequency*(time-position) + phase)
		return signal

	def _read_book(self, f):
		try:
			f = open(f, 'rb')
		except Exception:
			f = f
		version = np.fromfile(f, 'S6', count=1)
		print('version: ',version[0])
		data = {} # zachowanie danych z nagłówku
		ident = np.fromfile(f, 'u1', count=1)[0]
		ct = self._get_type(ident, f)
		signals = collections.defaultdict(list)
		atoms = collections.defaultdict(list)
		while ident:
			if ct:
				point = np.fromfile(f,ct, count=1)[0]
				data[ident] = point
			elif ident ==7:
				data_s = np.fromfile(f, '>u4', count=1)[0]
				epoch_nr = np.fromfile(f, '>u2', count=1)[0]
				epoch_s = np.fromfile(f, '>u4', count=1)[0]
				#print('epoch_s', epoch_s)
			elif ident == 8:
				chnl_nr, signal = self._get_signal(f, epoch_nr, epoch_s)
				signals[epoch_nr].append([chnl_nr, signal])
		
			elif ident == 9:
				pl = f.tell()
				atom, a_chnl_nr = self._get_atoms(f)
				atoms[a_chnl_nr] = atom
			
			ident = np.fromfile(f, '>u1', count=1)
			if ident:
				ident = ident[0]
			ct = self._get_type(ident, f)
		return data, signals, atoms, epoch_s

	def _get_atoms_amplitudes(self):
		amps = []
		for trial in self.atoms.keys():
			amps_temp = []
			for i,atom in enumerate(self.atoms[trial][13]):
				position  = atom['t']/self.fs
				width     = atom['scale']/self.fs/2
				frequency = atom['f']*self.fs/2
				amplitude = atom['amplitude']
				phase 	  = atom['phase']
				if self.t_min<position<self.t_max:
					if self.atom_min <= i <= self.atom_max:
						amps_temp.append(amplitude)
			# print(amps_temp)
			amps.append(max(amps_temp))
		return amps	

	def _find_extrema(self, signal):
		deriv = np.gradient(signal)
		deriv[np.where(deriv < 0)] = -1
		deriv[np.where(deriv > 0)] = 1
		d = list(deriv)
		amps = []
		temp = 0
		data = []
		for j,i in enumerate(d):
			try:
				if i == d[j-1] and i == 1:
					if temp == 0:
						start = j-1
					temp += 1
				else:
					if temp != 0:
						data.append((start,temp))
					temp = 0
			except:
				pass
		for frag in data:
			slope_start = frag[0]
			slope_end = frag[0] + frag[1]
			amps.append(abs(signal[slope_end] - signal[slope_start]))
		return max(amps)

	def _get_atoms_peak2peak_amplitudes(self):
		amps = []
		for trial in self.atoms.keys():
			signal_reconstruction = self.reconstruct_signal_freq(self.atoms[trial])
			amps.append(self._find_extrema(signal_reconstruction))
			# py.figure(str(trial))
			# py.plot(signal_reconstruction)
			# py.show()
		return amps	

	def draw_mean_reconstructed_signal(self, atoms, signals):
		N = len(atoms.keys())
		tpr = len(signals[1][0][1])
		t = np.arange(0, tpr/self.fs, 1/self.fs)
		for i,trial in enumerate(atoms.keys()):
			signal_reconstruction = self.reconstruct_signal_freq(atoms[trial])
			signal = signals[1][i][1]
			try:
				signal_a += signal
				signal_reconstruction_a += signal_reconstruction
			except Exception as a:
				print(a, 'objection')
				signal_a = signal
				signal_reconstruction_a = signal_reconstruction
		signal_a /= N
		signal_reconstruction_a /= N
		py.figure('Mean')
		py.plot(t,signal_reconstruction_a,color='g',label='rekonstrukcja')
		py.plot(t,signal_a,color='m',label='sygnal')
		# py.ylim(-30,40)
		py.ylabel('Amplituda [uV]')
		py.xlabel('Czas [s]')
		py.axvline(x=0.8,color='r')
		# py.axvline(x=0.3,color='r')
		py.legend()

	def draw_reconstructed_signal(self):
		tpr = len(self.signals[1][0][1])
		# t = np.arange(-0.5, tpr/self.fs, 1/self.fs)
		t = np.linspace(-0.5, 1, tpr)
		for i,trial in enumerate(self.atoms.keys()):
			if i in xrange(self.id_min, self.id_max+1):
				signal_reconstruction = self.reconstruct_signal_freq(self.atoms[trial])
				signal = self.signals[1][i][1] #- np.mean(self.signals[1][i][1][:0.5*self.fs])
				# amp = np.mean(self.signals[1][i][1][0.25*self.fs:0.35*self.fs])
				# amp_sig.append(amp)
				py.figure('Trial '+str(i))
				py.plot(t,signal_reconstruction,color='g',label='rekonstrukcja')
				py.plot(t,signal,color='m',label='sygnal')
				# py.axvline(x=0.3,color='r')
				py.axvline(x=0.3,color='r')
				py.ylabel('Amplituda [uV]')
				py.xlabel('Czas [s]')
				# py.ylim(-30,40)
				py.xlim(-0.5,1)
				py.legend()


	def reconstruct_signal_freq(self, atoms):
		reconstruction = np.zeros(self.epoch_s)
		for i,atom in enumerate(atoms[13]):
			position  = atom['t']/self.fs
			width     = atom['scale']/self.fs
			frequency = atom['f']*self.fs/2
			amplitude = atom['amplitude']
			phase 	  = atom['phase']
			# if self.t_min < position < self.t_max:
			if self.atom_min <= i < self.atom_max+1:
				reconstruction = reconstruction + self._gabor(amplitude, position, width, frequency, phase)
		return reconstruction

	def perform_linear_regression(self, amps):
		y = np.linspace(1,len(amps),len(amps))
		A = np.array([y, np.ones(len(amps))])
		w = np.linalg.lstsq(A.T,amps)[0]
		line = w[0]*y+w[1]
		py.figure()
		py.plot(y,line,'r-',y,amps,'o')
		py.ylim(-5,22)
		py.ylabel('Amplituda [uV]')
		py.xlabel('Realizacje')

################ rysowanie mapy
# 8 Hz < f < 12 Hz
# a > 5 µV
# s > 1.5 s

	def calculate_mean_map(self,df = 0.05, dt = 1/256., f_a = [0, 64.]):
		N = len(self.atoms.keys())
		tpr = len(self.signals[1][0][1])
		# sigt = np.arange(0, tpr/self.fs, 1/self.fs)
		sigt = np.linspace(-0.5, 1, tpr)
		alpha_mean = []
		for nr, chnl in enumerate(self.atoms.keys()):
			t, f, E, sigt_s,s_s,sr_s = self.calculate_map(self.atoms[chnl],self.signals[1][nr][1], df, dt,  f_a = f_a)
			alpha_recon = self._calculate_alpha_power(self.atoms[chnl],self.signals[1][nr][1], df, dt,  f_a = f_a)
			alpha_mean.append(sum(alpha_recon))
			signal_reconstruction = self._reconstruct_signal(self.atoms[chnl])
			signal = self.signals[1][nr][1]
			try:
				signal_a += signal
				E_a += E
				signal_recontruction_a += signal_reconstruction
			except Exception as a:
				# print(a)
				signal_a = signal
				E_a = E
				signal_recontruction_a = signal_reconstruction
		signal_a /= N
		signal_recontruction_a /= N
		E_a /= N
		return t, f, E_a, sigt, signal_a, signal_recontruction_a,alpha_mean

	def calculate_map(self,atoms,signal, df, dt, contour=0, f_a = [0, 64.]):
		''' atoms - dictionary with trials/channels; for each dictionary with 
		atoms 4 different atom's types (10,11,12,13-Gabors)'''
		tpr = len(signal)
		f = np.arange(f_a[0], f_a[1], df)
		t = np.arange(0, tpr/self.fs, dt)
		lent = len(t)
		lenf = len(f)
		E = np.zeros((lent, lenf)).T
		t, f = np.meshgrid(t, f)
		sigt = np.arange(0, tpr/self.fs, 1/self.fs)
		for atom in atoms[13]:
			exp1 = np.exp(-2*(atom['scale']/self.fs)**(-2)*(t-atom['t']/self.fs)**2)
			exp2 = np.exp(-2*np.pi**2 *(atom['scale']/self.fs)**2*(atom['f']*self.fs/2-f)**2)
			wigners = ((atom['amplitude']/self.ptspmV)**2 * (2*np.pi)**0.5 * 
						atom['scale']/self.fs*exp1*exp2)
			E += atom['modulus']**2 * wigners
		for atom in atoms[12]:
			amp =  atom['modulus']**2*(atom['amplitude']/self.ptspmV)**2
			E[:,int(len(f)*atom['f']/(2*np.pi))] +=  amp
		for atom in atoms[11]:
			exp1 = np.exp(-2*(atom['scale']/self.fs)**(-2)*(t-atom['t']/self.fs)**2)
			exp2 = np.exp(-2*np.pi**2 *(atom['scale']/self.fs)**2*(-f)**2)
			wigners = ((atom['amplitude']/self.ptspmV)**2 * (2*np.pi)**0.5 * 
						atom['scale']/self.fs*exp1*exp2)
			E += atom['modulus']**2 * wigners
		for atom in atoms[10]:
			amp =  atom['modulus']**2*(atom['amplitude']/self.ptspmV)**2
			E[int(lent*atom['t']/tpr)] +=  amp	
		signal_reconstruction = self._reconstruct_signal(atoms)
		return t, f, np.log(E), sigt, signal, signal_reconstruction

	# def _calculate_alpha_power(self,atoms,signal, df, dt, contour=0, f_a = [0, 64.]):
	# 	tpr = len(signal)
	# 	f = np.arange(f_a[0], f_a[1], df)
	# 	t = np.arange(0, tpr/self.fs, dt)
	# 	lent = len(t)
	# 	lenf = len(f)
	# 	E = np.zeros((lent, lenf)).T
	# 	t, f = np.meshgrid(t, f)
	# 	sigt = np.arange(0, tpr/self.fs, 1/self.fs)
	# 	for atom in atoms[13]:
	# 		freq = atom['f']*self.fs/2
	# 		amp = 2*atom['amplitude']/self.ptspmV
	# 		scale = atom['scale']/self.fs
	# 		position = atom['t']/self.fs
	# 		# print(freq,position,scale)
	# 		if 8 < freq < 12: # and position < 0.5:
	# 			print(freq,position,scale)
	# 			exp1 = np.exp(-2*(atom['scale']/self.fs)**(-2)*(t-atom['t']/self.fs)**2)
	# 			exp2 = np.exp(-2*np.pi**2 *(atom['scale']/self.fs)**2*(atom['f']*self.fs/2-f)**2)
	# 			wigners = ((atom['amplitude']/self.ptspmV)**2 * (2*np.pi)**0.5 * 
	# 						atom['scale']/self.fs*exp1*exp2)
	# 			E += atom['modulus']**2 * wigners
	# 	return np.mean(np.mean(E))

	def _calculate_alpha_power(self,atoms,signal, df, dt, contour=0, f_a = [0, 64.]):
		reconstruction = np.zeros(self.epoch_s)
		for atom in atoms[13]:
			amplitude = atom['amplitude']
			position  = atom['t']/self.fs
			width     = atom['scale']/self.fs
			frequency = atom['f']*self.fs/2
			phase 	  = atom['phase']
			if 8 < frequency < 12 and position < 0.5:
			# if 13 < frequency < 30 and position < 0.5: #and amplitude > 5:
				# print(frequency,position,width)
				reconstruction = reconstruction + self._gabor(amplitude,position,width,frequency,phase)
		return np.abs(reconstruction)**2

	def only_draw_map(self, t, f, E_a, sigt, signal_a, signal_recontruction_a, contour=False):
		fig = py.figure()
		gs = gridspec.GridSpec(2,1, height_ratios=[3,1])
		ax1 = fig.add_subplot(gs[0])
		ax1.set_ylabel(u'Czestosc [Hz]')
		ax1.set_title(sys.argv[-1])
		if contour:
			ax1.contour(t, f, E_a)
		else:
			ax1.pcolor(t, f, E_a)
		ax2 = fig.add_subplot(gs[1])
		ax2.plot(sigt, signal_a, 'red')
		ax2.plot(sigt, signal_recontruction_a, 'blue')
		ax2.axvline(x=0,color='r')
		ax2.set_ylabel(u'Amplituda, [$\\mu$V]')
		ax2.set_xlabel(u'Czas [s]')
		ax2.set_xlim(-0.5,1)

	def _reconstruct_signal(self, atoms):
		reconstruction = np.zeros(self.epoch_s)
		for atom in atoms[13]:
			position  = atom['t']/self.fs
			width     = atom['scale']/self.fs
			frequency = atom['f']*self.fs/2
			amplitude = atom['amplitude']
			phase = atom['phase']
			# print(amplitude)
			reconstruction = reconstruction + self._gabor(amplitude,position,width,frequency,phase)
		return reconstruction


if __name__ == '__main__':

	# fstr = './wybrane/newMP/Cz/MMP1_compatible_group_nogo_longer_mmp.b'

	# fstr = './wybrane/MMP1_compatible_ind_go_mmp.b'
	
	# fstr = './wybrane/MMP2_compatible_ind_go_mmp.b'
	# fstr = './wybrane/MMP3_compatible_ind_go_mmp.b'

	# fstr = './wybrane/MMP1_compatible_ind_go_longer_mmp.b'

	fstr = './wybrane/newMP/Cz/MEANS/MMP3_MEANS_GROUP_mmp.b'
	b = BookImporter(fstr, id_min=0, id_max=20, t_min=0, t_max=1.5, atom_min=0, atom_max=0)
	ampsg = b.amps

	fstr = './wybrane/newMP/Cz/MEANS/MMP3_MEANS_INDIVIDUAL_mmp.b'
	b = BookImporter(fstr, id_min=0, id_max=10, t_min=0, t_max=1.5, atom_min=0, atom_max=1)
	ampsi = b.amps

	# fstr = './wybrane/newMP/Cz/MMP1_compatible_ind_nogo_longer_mmp.b'
	# b = BookImporter(fstr, id_min=0, id_max=1, t_min=0.5, t_max=1, atom_min=0, atom_max=2)
	# ind_comp = b.amps
	# ind_c_alpha = b.alpha_mean
	# fstr = './wybrane/newMP/Cz/MMP1_incompatible_ind_nogo_longer_mmp.b'
	# b = BookImporter(fstr, id_min=10, id_max=20, t_min=0.5, t_max=1, atom_min=0, atom_max=2)
	# ind_incomp = b.amps
	# # ind_inc_alpha = b.alpha_mean

	# ind = ind_comp + ind_incomp
	# # ind_alp = ind_c_alpha + ind_inc_alpha

	# fstr = './wybrane/newMP/Cz/MMP1_incompatible_group_nogo_longer_mmp.b'
	# b = BookImporter(fstr, id_min=0, id_max=1, t_min=0.5, t_max=1, atom_min=0, atom_max=2)
	# gr_comp = b.amps
	# # gr_c_alpha = b.alpha_mean
	# fstr = './wybrane/newMP/Cz/MMP1_compatible_group_nogo_longer_mmp.b'
	# b = BookImporter(fstr, id_min=10, id_max=20, t_min=0.5, t_max=1, atom_min=0, atom_max=2)
	# gr_incomp = b.amps
	# # gr_inc_alpha = b.alpha_mean

	# gr = gr_incomp + gr_comp
	# # gr_alp = gr_inc_alpha + gr_c_alpha

	# # print(pearsonr(gr_comp,gr_c_alpha))

	# # print(stats.shapiro(gr)[1],stats.shapiro(ind)[1])
	# print(stats.kruskal(gr,ind)[1])
