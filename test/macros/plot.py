#!/usr/bin/env python
import os, sys
from ROOT import gStyle,gROOT,std,ROOT,TFile,TTree,TH1D,TH2D,TCanvas,TLegend,kRed,kBlue,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)

#from tex import latex2pdf

def setCosmetics(h,c):
    h.SetLineColor(c)
    h.SetLineWidth(2)
    
Path = '/Users/skhalil/Desktop/Analysis/KUSUSY/Analysis/'

# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()                                                                                                                                         
parser.add_option('--isb', metavar='T', type='string', action='store',
                  default='b',
                  dest='isb',
                  help='options are: b, notb')

parser.add_option('--hist', metavar='T', type='string', action='store',
                  default='hCutflow_',
                  dest='hist',
                  help='name of the histogram from input files')

parser.add_option('--rebin', metavar='T', type='string', action='store',
                  default='1',
                  dest='rebin',
                  help='rebin the histograms')

parser.add_option('--cutflow', action='store_true',
                  default=False,
                  dest='cutflow',
                  help='draw a cutflow plot')

(options,args) = parser.parse_args()
# ==========end: options =============

f_signal    = TFile(Path+'SMS-T2bb_mSbot-1650to2600_out.root')
f_top       = TFile(Path+'TTJets_out.root')


h_cutflow_sig  = f_signal.Get(options.hist + options.isb)
h_cutflow_top  = f_top.Get(options.hist + options.isb)
hname = h_cutflow_top.GetName()
setCosmetics(h_cutflow_sig, kBlue)
setCosmetics(h_cutflow_top, kRed)
    
if options.cutflow:

    bin1_top_value = 0
    bin1_sig_value = 0
    nBins = h_cutflow_sig.GetNbinsX()

    f = open('Plots/'+hname+'_table.tex', 'wb')
    
    print '\n'

    #print '\\documentclass{article}'
    #f.write('\\documentclass{article} \n')
    print '\\documentclass{standalone}'
    f.write('\\documentclass{standalone} \n')
    print '\\begin{document}'
    f.write('\\begin{document} \n')
    #print '\\begin{table}'
    #f.write('\\begin{table} \n')
    #print '\\begin{center}'
    #f.write('\\begin{center} \n')
    print '\\begin{tabular}{l|c|r}'
    f.write('\\begin{tabular}{l|c|r} \n')
    print '\hline'
    f.write('\hline \n')
    print ' & top & signal\\\\ '
    f.write(' & top & signal\\\\ \n')
    print '\hline'
    f.write('\hline \n')
    for ibin in range(1, nBins+1):
        isig = h_cutflow_sig.GetBinContent(ibin)
        itop = h_cutflow_top.GetBinContent(ibin)
        isig_err = h_cutflow_sig.GetBinError(ibin)
        itop_err = h_cutflow_top.GetBinError(ibin)
        if ibin==1:
            bin1_top_value = itop
            bin1_sig_value = isig
        
        print '{0:<5} & {1:<5.0f} $\pm$ {2:<5.1f}  & {3:<5.0f} $\pm$ {4:<5.1f} \\\\ '. format('bin'+str(ibin), itop, itop_err, isig, isig_err)
        f.write('{0:<5} & {1:<5.0f} $\pm$ {2:<5.1f}  & {3:<5.0f} $\pm$ {4:<5.1f} \\\\ \n'. format('bin'+str(ibin), itop, itop_err, isig, isig_err))
    print '\hline'
    f.write('\hline \n')            
    print '\end{tabular} \n'
    f.write('\end{tabular} \n')
    #print '\end{center}'
    #f.write('\end{center} \n')
    #print '\end{table}'
    #f.write('\end{table} \n')
    print '\end{document}'
    f.write('\end{document} \n')
    f.close()
         
    h_cutflow_top.Scale(1/bin1_top_value)
    h_cutflow_sig.Scale(1/bin1_sig_value)
    c1 = TCanvas('c1', 'c1', 800, 600)
    h_cutflow_sig.Draw("Hist")
    h_cutflow_top.Draw("Hist, same")
    #h_cutflow_sig.Draw("Hist, same")


leg = TLegend(0.50,0.92,0.95,0.80)
    #leg.SetY1(1-0.05*Ysize)
leg.SetBorderSize(1)
leg.SetFillColor(10)
leg.AddEntry(h_cutflow_top, "top", 'l')
leg.AddEntry(h_cutflow_sig, "signal", 'l')
leg.Draw()
    
hname = h_cutflow_top.GetName()
if options.cutflow: c1.SaveAs("Plots/"+hname+".pdf")
else:
    c2 = TCanvas('c2', 'c2', 800, 600)
    
    h_cutflow_top.DrawNormalized("Hist")
    h_cutflow_sig.DrawNormalized("Hist, same")
    leg.Draw()
    c2.SaveAs("Plots/"+hname+".pdf")
    
        
#raw_input('...')
