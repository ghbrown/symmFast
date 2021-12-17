#!/usr/bin/env python3
"""
# Created: Thu Dec 16 15:12:54 2021 (-0500)
# @author: jacobfaibussowitsch
"""
import collections
import matplotlib.pyplot as plt
import numpy as np

def plot(x,data,sub,title,xlabel,ylabel):
  for key,val in data.items():
    if key == 'dense':
      label = 'canonical 1-D column'
    else:
      label = 'fast 2-D upper triangular'
    plt.loglog(x,val[sub],label=label,marker='x')
  plt.ylabel(ylabel)
  plt.xlabel(xlabel)
  plt.title(title)
  plt.grid()
  plt.legend()
  ax = plt.gca()
  ax.yaxis.set_tick_params(which='major',size=8,width=2,direction='in',right='on')
  ax.yaxis.set_tick_params(which='minor',size=4,width=1,direction='in',right='on')
  plt.show()
  return

def main(infile):
  data = collections.defaultdict(list)
  mat_size  = []
  num_procs = []
  dim_r     = []
  for line in infile.read_text().splitlines():
    entries = line.split()
    if len(entries) == 2:
      # timings
      data[entries[0]].append(float(entries[1]))
    else:
      mat_size.append(int(entries[0]))
      num_procs.append(int(entries[1]))
      dim_r.append(int(entries[2]))
  mat_size  = np.array(mat_size)
  num_procs = np.array(num_procs)
  dim_r     = np.array(dim_r)
  for k,v in data.items():
    data[k] = {'time' : np.array(v)}
    dk = data[k]
    data[k]['speedup'] = dk['time'][0]/dk['time']
    data[k]['efficiency'] = np.divide(dk['speedup'],num_procs)
  suffix    = infile.stem.split('_')[1].title()
  plot(num_procs,data,'speedup',f'Strong Scaling Speedup - {suffix}','ranks','speedup')
  plot(num_procs,data,'efficiency',f'Strong Scaling Efficiency - {suffix}','ranks','efficiency')
  plot(num_procs,data,'time',f'Execution Time - {suffix}','ranks','time [s]')
  return

if __name__ == '__main__':
  import argparse
  import pathlib

  parser = argparse.ArgumentParser()
  parser.add_argument('data_file',type=pathlib.Path,help='data file containing timings')
  args = parser.parse_args()

  main(args.data_file.resolve())
