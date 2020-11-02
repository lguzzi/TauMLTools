from __future__ import print_function
import sys; PYTHON_MAJOR = int(sys.version_info.major)
dict_iterator = 'items' if PYTHON_MAJOR == 3 else 'iteritems'
from itertools import product

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

parser.add_argument('--visual', action = 'store_true', help = 'Won\'t run the script in batch mode')
parser.add_argument('--legend', action = 'store_true', help = 'Draw a TLegent on canvases')
args = parser.parse_args()

ROOT.gROOT.SetBatch(not args.visual)
ROOT.gStyle.SetOptStat(0)
ROOT.ROOT.EnableImplicitMT(4)
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

class Entry_3D:
  def __init__(self, var, pointer, binned_var, tdir = None):
    self.var  = var
    self.ptr  = pointer
    self.tdir = '/'.join([var, binned_var]) if tdir is None else tdir
    self.binned_var = binned_var
    self.entries_2D = []
  
  def get_histograms(self):
    self.histo = self.ptr.GetValue()
    self.entries_2D = [
      Entry_2D(var = self.var, histo2D = self.get_histo_2D(vbin = bb+1))
        for bb in range(*BINS[self.binned_var][1:])
    ]
  
  def get_histo_2D(self, vbin):
    histo2D = ROOT.TH2F(self.var, '', *((N_SPLITS, 0, N_SPLITS)+BINS[self.var]))
    for x, y in product(range(N_SPLITS), range(*BINS[self.var][1:])):
      histo2D.Fill(x+1,y+1,self.histo.GetBinContent(x+1,y+1,vbin))
    
    return histo2D
  
  def run_KS_test(self):
    for ee in self.entries_2D:
      ee.run_KS_test()
  
  def save_data(self):
    for ee in self.entries_2D:
      ee.save_data()


class Entry_2D:
  def __init__(self, var, histo2D = None, pointer = None, tdir = None):
    self.var  = var
    self.ptr  = pointer
    self.tdir = self.var if tdir is None else tdir

    self.histo2D = histo2D if pointer is None else None
    self.histos  = []
  
  def get_histograms(self):
    self.histo2D = self.ptr.GetValue()
  
  def get_chunks(self, norm = True):
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

if __name__ == '__main__':
  print ('[INFO] reading files', args.input)
  
  input_files = ROOT.std.vector('std::string')()
  
  for file in glob.glob(args.input):
    input_files.push_back(str(file))

  dataframe = ROOT.RDataFrame('taus', input_files)
  dataframe = dataframe.Define('chunk_id', 'rdfentry_ % {}'.format(N_SPLITS))
  
  entries = [
    ## unbinned distributions
    Entry_2D(var = 'lepton_gen_match', pointer = dataframe.Histo2D(('lepton_gen_match', '', N_SPLITS, 0, N_SPLITS)+BINS['lepton_gen_match'], 'chunk_id', 'lepton_gen_match')),
    Entry_2D(var = 'sampleType'      , pointer = dataframe.Histo2D(('sampleType'      , '', N_SPLITS, 0, N_SPLITS)+BINS['sampleType']      , 'chunk_id', 'sampleType'))      ,
    Entry_2D(var = 'dataset_group_id', pointer = dataframe.Histo2D(('dataset_group_id', '', N_SPLITS, 0, N_SPLITS)+BINS['dataset_group_id'], 'chunk_id', 'dataset_group_id')),
    Entry_2D(var = 'dataset_id'      , pointer = dataframe.Histo2D(('dataset_id'      , '', N_SPLITS, 0, N_SPLITS)+BINS['dataset_id']      , 'chunk_id', 'dataset_id'))      ,
    ## binned distributions
    Entry_3D(var = 'tau_pt', pointer = dataframe.Histo3D(('tau_pt', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_pt']+BINS['lepton_gen_match'], 'chunk_id', 'tau_pt', 'lepton_gen_match'), binned_var =  'lepton_gen_match'),
    Entry_3D(var = 'tau_pt', pointer = dataframe.Histo3D(('tau_pt', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_pt']+BINS['sampleType']      , 'chunk_id', 'tau_pt', 'sampleType'      ), binned_var =  'sampleType')      ,
    Entry_3D(var = 'tau_pt', pointer = dataframe.Histo3D(('tau_pt', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_pt']+BINS['dataset_group_id'], 'chunk_id', 'tau_pt', 'dataset_group_id'), binned_var =  'dataset_group_id'),
    Entry_3D(var = 'tau_pt', pointer = dataframe.Histo3D(('tau_pt', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_pt']+BINS['dataset_id']      , 'chunk_id', 'tau_pt', 'dataset_id'      ), binned_var =  'dataset_id')      ,
    
    Entry_3D(var = 'tau_eta', pointer =  dataframe.Histo3D(('tau_eta', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_eta']+BINS['lepton_gen_match'], 'chunk_id', 'tau_eta', 'lepton_gen_match'), binned_var =  'lepton_gen_match'),
    Entry_3D(var = 'tau_eta', pointer =  dataframe.Histo3D(('tau_eta', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_eta']+BINS['sampleType']      , 'chunk_id', 'tau_eta', 'sampleType'      ), binned_var =  'sampleType')      ,
    Entry_3D(var = 'tau_eta', pointer =  dataframe.Histo3D(('tau_eta', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_eta']+BINS['dataset_group_id'], 'chunk_id', 'tau_eta', 'dataset_group_id'), binned_var =  'dataset_group_id'),
    Entry_3D(var = 'tau_eta', pointer =  dataframe.Histo3D(('tau_eta', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_eta']+BINS['dataset_id']      , 'chunk_id', 'tau_eta', 'dataset_id'      ), binned_var =  'dataset_id')      ,

    Entry_3D(var = 'dataset_id', pointer = dataframe.Histo3D(('dataset_id', '', N_SPLITS, 0, N_SPLITS)+BINS['tau_eta']+BINS['dataset_group_id'], 'chunk_id', 'dataset_id', 'dataset_group_id'), binned_var =  'dataset_group_id'),
  ]
  for ee in entries:
    if ee.get_histograms():
      ee.run_KS_test()
      ee.save_data()

  '''
  tau_pt_lgm_histo3D = dataframe.Histo3D()
  tau_pt_st_histo3D  = dataframe.Histo3D()
  tau_pt_dgi_histo3D = dataframe.Histo3D()
  tau_pt_di_histo3D  = dataframe.Histo3D()

  tau_eta_lgm_histo3D = dataframe.Histo3D()
  tau_eta_st_histo3D  = dataframe.Histo3D()
  tau_eta_dgi_histo3D = dataframe.Histo3D()
  tau_eta_di_histo3D  = dataframe.Histo3D()
  '''


  OUTPUT_ROOT.Close()
  json.dump(JSON_DICT, OUTPUT_JSON, indent = 4)
  print ('[INFO] all done. Files saved in', args.output)
  #print ('RDataFrame was run', dataframe.GetNRuns(), 'times')