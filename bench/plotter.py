#!/usr/bin/env python3
"""
# Created: Thu Dec 16 15:12:54 2021 (-0500)
# @author: jacobfaibussowitsch
"""
import collections
import matplotlib.pyplot as plt

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
  import ipdb; ipdb.set_trace()
  return

if __name__ == '__main__':
  import argparse
  import pathlib

  parser = argparse.ArgumentParser()
  parser.add_argument('data_file',type=pathlib.Path,help='data file containing timings')
  args = parser.parse_args()

  main(args.data_file.resolve())
