#!/usr/bin/env python

import pdb
import sys
import time

import astropy.units as u
import h5py
import matplotlib.pyplot as plt 
import numpy as np

from multiprocessing import Pool

from astropy.units import Quantity
from matplotlib import dates
from sunpy.time import TimeRange

sys.path.insert(0,"/home/obs/Pearse")
from lofar_bst import LOFAR_Raw, LOFAR_Raw_357

def bg(data, amount=0.05):
	#tmp = data - np.mean(data, axis=0) 
	sdevs = np.mean(data, 1) 
	cand = sorted(range(data.shape[0]), key=lambda y: sdevs[y]) 
	realcand = cand[:max(1, int(amount*len(cand)))] 
	bg = np.mean(data[realcand, :], 0) 
	return bg 
start = time.time()
d_file = sys.argv[1]
trange =  TimeRange("2020/06/03 06:35:00", 60*u.s) 
frange = [20,45]*u.MHz
raw_357 = LOFAR_Raw_357(d_file)#, trange, frange)
#raw_357 = LOFAR_Raw(d_file,np.arange(56, 542), 3,  trange)#, frange)
#bst_357 = Lofar_bst_357("20200602_071428_bst_00X.dat")
read_time = time.time() - start
print("Time to read 357: {} seconds".format(read_time))

new_dt = 1e-3
resample_start = time.time()
sum_int = int(np.floor(new_dt/raw_357.dt.sec))
sum_shape = int(np.round(raw_357.data.shape[1]/sum_int))
def resample(i):
	return np.sum(raw_357.data[:,i*sum_int:(i+1)*sum_int], axis=1)
with Pool() as pool:
	#resample_part = functools.partial(resample, data = raw_357.data)
	resample_data = pool.map(resample,range(sum_shape))
resample_data = np.array(resample_data)
resample_time = time.time() - resample_start

bg_data = resample_data/ bg(resample_data,0.05)
print("Time to resample multi: {} seconds".format(resample_time))

fig, ax = plt.subplots(figsize=(10,8))
if raw_357.trange is None:
	ax.imshow(bg_data.T, aspect="auto", 
		extent=[raw_357.obs_start.plot_date, raw_357.obs_end.plot_date,
			raw_357.freqs[-1][-1].value, raw_357.freqs[0][0].value],
		vmin=np.percentile(bg_data, 1), vmax=np.percentile(bg_data,99))
else:
	ax.imshow(bg_data.T, aspect="auto", 
		extent=[raw_357.trange.start.plot_date, raw_357.trange.end.plot_date,
			raw_357.freqs[-1][-1].value, raw_357.freqs[0][0].value],
		vmin=np.percentile(bg_data, 1), vmax=np.percentile(bg_data,99))
date_format = dates.DateFormatter("%H:%M:%S")
ax.xaxis_date()
ax.xaxis.set_major_formatter(date_format)
#plt.savefig("raw_plot_sbursts.png")
end = time.time() - start
print("Time to run: {} seconds".format(end))
plt.show()
