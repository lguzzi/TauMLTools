from __future__ import print_function
import sys; PYTHON_MAJOR = int(sys.version_info.major)
dict_iterator = 'items' if PYTHON_MAJOR == 3 else 'iteritems'

import ROOT
import glob
import json
from collections import OrderedDict

import argparse
parser = argparse.ArgumentParser('''
The script runs a binned KS test between different chunks of the same RDataFrame. For simplicity, all chunks are compared to the first one.
If a KS test is below the threshold, a warning message is printed on screen.
NOTE: the binning of each variable must be hard coded in the script (using the BINS dictionary)
NOTE: pvalue = 99 means that one of the two histograms is empty, -99 means never filled.
''')


parser.add_argument('--input'       , required = True, type = str, help = 'input file. Accepts glob patterns (use quotes)')
parser.add_argument('--output'      , required = True, type = str, help = 'output directory name')
parser.add_argument('--nsplit'      , default  = 100 , type = int, help = 'number of chunks per file')
parser.add_argument('--pvthreshold' , default  = .05 , type = str, help = 'threshold of KS test (above = ok)')
parser.add_argument('--n_threads'   , default  = 1   , type = int, help = 'enable ROOT implicit multithreading')

parser.add_argument('--visual', action = 'store_true', help = 'Won\'t run the script in batch mode')
parser.add_argument('--legend', action = 'store_true', help = 'Draw a TLegent on canvases')
args = parser.parse_args()

ROOT.gROOT.SetBatch(not args.visual)
ROOT.gStyle.SetOptStat(0)

import os
pdf_dir = '/'.join([args.output, 'pdf'])
if not os.path.exists(pdf_dir):
  os.makedirs(pdf_dir)

JSON_DICT      = OrderedDict()
OUTPUT_ROOT    = ROOT.TFile.Open('{}/histograms.root'.format(args.output), 'RECREATE')
OUTPUT_JSON    = open('{}/pvalues.json'.format(args.output), 'w')
N_SPLITS       = args.nsplit
PVAL_THRESHOLD = args.pvthreshold

## binning of tested variables (cannot use unbinned distributions with python before root 6.18)
BINS = {
  'tau_pt'    : (50, 0, 5000),
  'tau_eta'   : (5, -2.5, 2.5),
  'lepton_gen_match' : (20, -1, 19),
  'sampleType': (20, -1, 19),      
  'dataset_id': (20, -1, 19),
  'dataset_group_id': (20, -1, 19),
}

class Entry:
  def __init__(self, var, dataframe, bin_dir = None):
    self.var      = var
    self.dframe   = dataframe
    self.tdir     = '/'.join([bin_dir, self.var]) if not bin_dir is None else self.var
    self.hptr     = None
    self.histos   = []
    self.pvalues  = []
    self.nev      = self.dframe.Count()
  
  def load_hpointers(self):
    self.dframe = self.dframe.Define('chunk_id', 'rdfentry_ % {}'.format(N_SPLITS))
    model = (self.var, '') + (N_SPLITS, 0, N_SPLITS) + BINS[self.var]
    self.hptr = self.dframe.Histo2D(model, 'chunk_id', self.var)
  
  def load_histograms(self, norm = True):
    if not self.nev.GetValue():
      return False
    
    histo2D = self.hptr.GetValue()
    self.histos = [histo2D.ProjectionY('chunk_{}'.format(jj), jj+1, jj+1) for jj in range(N_SPLITS)]
    self.histos = [hh for hh in self.histos if hh.Integral()]
    
    if not len(self.histos):
      return False
    
    self.histos[0].SetMarkerStyle(20)
    for jj, hh in enumerate(self.histos):
      hh.SetTitle(self.tdir)
      hh.Sumw2()
      hh.SetLineColor(jj+1)
      if norm:
        hh.Scale(1. / hh.Integral())
    
    return True
      
  def run_KS_test(self):
    self.pvalues = [self.histos[0].KolmogorovTest(hh) if self.histos[0].Integral()*hh.Integral() else 99 for hh in self.histos]
    if not self.histos[0].Integral():
      print ('[WARNING] control histogram is empty inside {}'.format(self.tdir))
    
    if not all([pv >= PVAL_THRESHOLD for pv in self.pvalues]):
      print ('[WARNING] KS test failed for step {}. p-values are:'.format(self.tdir))
      print ('\t', self.pvalues)
    
  def save_data(self):
    OUTPUT_ROOT.cd()
    
    if not OUTPUT_ROOT.GetDirectory(self.tdir):
      OUTPUT_ROOT.mkdir(self.tdir)
    
    OUTPUT_ROOT.cd(self.tdir)
    
    can = ROOT.TCanvas()
    leg = ROOT.TLegend(0.9, 0.1, 1., 0.9, "p-values (KS with the first chunk)")
    for ii, hh in enumerate(self.histos):
      hh.Write()
      hh.Draw('PE'+' SAME'*(ii != 0))
      leg.AddEntry(hh, 'chunk %d - pval = %.3f' %(ii, self.pvalues[ii]), 'lep')
    if args.legend:
      leg.Draw("SAME")
    
    can.SaveAs('{}/pdf/{}.pdf'.format(args.output, self.tdir.replace('/', '_')), 'pdf')
    can.Write()

    OUTPUT_ROOT.cd()

    json_here = JSON_DICT
    for here in self.tdir.split('/'):
      if not here in json_here.keys():
        json_here[here] = OrderedDict()
      json_here = json_here[here]
    json_here['pvalues'] = self.pvalues


def groupby(dataframe, by):
  return {tt: dataframe.Filter('{} == {}'.format(by, tt)) for tt in range(*BINS[by][1:])}

if __name__ == '__main__':
  print ('[INFO] reading files', args.input)
  
  if args.n_threads > 1:
    ROOT.ROOT.EnableImplicitMT(args.n_threads)

  input_files = ROOT.std.vector('std::string')()
  
  for file in glob.glob(args.input):
    input_files.push_back(str(file))

  main_dir = 'KS_test'
  dataframe = ROOT.RDataFrame('taus', input_files)
  dataframe_lgm = groupby(dataframe, 'lepton_gen_match')
  dataframe_st  = groupby(dataframe, 'sampleType')
  dataframe_dgi = groupby(dataframe, 'dataset_group_id')
  dataframe_di  = groupby(dataframe, 'dataset_id')

  main_entries = [
    Entry(var = vv, dataframe = dataframe) for vv in ['lepton_gen_match', 'sampleType', 'dataset_group_id', 'dataset_id']
  ]
  sub_entries_lgm = {
    'tau_pt'  : [Entry(var = 'tau_pt' , dataframe = df, bin_dir = 'lepton_gen_match/{}'.format(bb)) for bb, df in getattr(dataframe_lgm, dict_iterator)()],
    'tau_eta' : [Entry(var = 'tau_eta', dataframe = df, bin_dir = 'lepton_gen_match/{}'.format(bb)) for bb, df in getattr(dataframe_lgm, dict_iterator)()],
  }
  sub_entries_st = {
    'tau_pt'  : [Entry(var = 'tau_pt' , dataframe = df, bin_dir = 'sampleType/{}'.format(bb)) for bb, df in getattr(dataframe_st, dict_iterator)()],
    'tau_eta' : [Entry(var = 'tau_eta', dataframe = df, bin_dir = 'sampleType/{}'.format(bb)) for bb, df in getattr(dataframe_st, dict_iterator)()],
  }
  sub_entries_dgi = {
    'tau_pt'  : [Entry(var = 'tau_pt' , dataframe = df, bin_dir = 'dataset_group_id/{}'.format(bb)) for bb, df in getattr(dataframe_dgi, dict_iterator)()],
    'tau_eta' : [Entry(var = 'tau_eta', dataframe = df, bin_dir = 'dataset_group_id/{}'.format(bb)) for bb, df in getattr(dataframe_dgi, dict_iterator)()],
  }
  sub_entries_di = {
    'tau_pt'  : [Entry(var = 'tau_pt' , dataframe = df, bin_dir = 'dataset_id/{}'.format(bb)) for bb, df in getattr(dataframe_di, dict_iterator)()],
    'tau_eta' : [Entry(var = 'tau_eta', dataframe = df, bin_dir = 'dataset_id/{}'.format(bb)) for bb, df in getattr(dataframe_di, dict_iterator)()],
  }
  entries = main_entries +\
    sub_entries_lgm['tau_pt'] + sub_entries_lgm['tau_eta'] +\
    sub_entries_st ['tau_pt'] + sub_entries_st ['tau_eta'] +\
    sub_entries_di ['tau_pt'] + sub_entries_di ['tau_eta'] +\
    sub_entries_dgi['tau_pt'] + sub_entries_dgi['tau_eta'] 

  for ee in entries:
    ee.load_hpointers()

  for ee in entries:
    if ee.load_histograms():
      ee.run_KS_test()
      ee.save_data()

  OUTPUT_ROOT.Close()
  json.dump(JSON_DICT, OUTPUT_JSON, indent = 4)
  print ('[INFO] all done. Files saved in', args.output)