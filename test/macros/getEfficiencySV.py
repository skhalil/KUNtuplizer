#!/usr/bin/env python
import os, sys
from ROOT import gStyle,gROOT,std,ROOT,TFile,TTree,TH1D,TH2D,TCanvas,TLegend,kRed,kBlue,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)

#from tex import latex2pdf

def pullEff(f, dbin, nbin):
    # input: root file and bin with final cuts(numerator), bin with selection (denominator)
    # output: list of eff and eff_err (floats) and list of hist
    h_list = []; eff_list = []; eff_err_list = [];
    
    for c in cat:
        h  = f.Get('hCutflow_' + c)
        nBins = h.GetNbinsX()
        h_list.append(h)
        h_eff     = h.Clone()
        h_eff.Scale(1/h.GetBinContent(dbin))
        eff     = h_eff.GetBinContent(nbin)
        eff_err = h_eff.GetBinError(nbin)
        eff_list.append(eff)
        eff_err_list.append(eff_err)

        if 'l' in c:            
            h_lightOther = h.Clone()
            h_lightOther.SetName('hCutflow_lOther')
        if 'other' in c:
            h_lightOther.Add(h)
            h_list.append(h_lightOther)                        
            h_eff_lo     = h_lightOther.Clone()
            h_eff_lo.Scale(1/h_lightOther.GetBinContent(dbin))
            eff_lo     = h_eff_lo.GetBinContent(nbin)
            eff_err_lo = h_eff_lo.GetBinError(nbin)
            eff_list.append(eff_lo)
            eff_err_list.append(eff_err_lo)
        
    return h_list, eff_list, eff_err_list 
    
def setCosmetics(h,color):
    h.SetLineColor(color)
    h.SetLineWidth(2)
    
Path = '/Users/skhalil/Desktop/Analysis/KUSUSY/Analysis/input/'


f_T2bb    = TFile(Path+'SMS-T2bb_mSbot-1650to2600_out.root')
f_top       = TFile(Path+'TTJets_out.root')

cat = ['unmatched', 'bnotg', 'cnotg', 'l', 'other']

#h_signal = []
h_top_list, eff_top_list, eff_err_top_list = pullEff(f_top, 3, 7)
h_T2bb_list, eff_T2bb_list, eff_err_T2bb_list = pullEff(f_T2bb, 3, 7)

print h_top_list



print 'eff_top: ', eff_top_list, ',\n  eff_err_top: ', eff_err_top_list

print '    '
print 'eff_t2bb: ', eff_T2bb_list, ',\n  eff_err_t2bb: ', eff_err_T2bb_list

    
       
