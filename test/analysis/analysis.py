#!/usr/bin/env python
import os, sys
#from itertools import permutations
from itertools import combinations 
from ROOT import gROOT,std,ROOT,TFile,TTree,TH1D,TH2D,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")

def setAxisBins(h, name):
   ibin=0
   for n in name:
      ibin = ibin+1
      h.GetXaxis().SetBinLabel(ibin, n)
   return h

def overUnderFlow(hist):
    xbins = hist.GetNbinsX()
    hist.SetBinContent(xbins, hist.GetBinContent(xbins)+hist.GetBinContent(xbins+1))
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinError(xbins, TMath.Sqrt(TMath.Power(hist.GetBinError(xbins),2)+TMath.Power(hist.GetBinError(xbins+1),2)))
    hist.SetBinError(1, TMath.Sqrt(TMath.Power(hist.GetBinError(0),2)+TMath.Power(hist.GetBinError(1),2)))
    hist.SetBinContent(xbins+1, 0.)
    hist.SetBinContent(0, 0.)
    hist.SetBinError(xbins+1, 0.)
    hist.SetBinError(0, 0.)
    
def fillSVProperties(drMinSVGen_i, drSVGen_i, dxy_i, dxyerr_i, dxysig_i, d3d_i, d3derr_i, d3dsig_i, cosTheta_i, ndof_i, ndau_i, hbin):
   hdrMinSVGen[hbin].Fill(drMinSVGen_i, evtwt) ;overUnderFlow(hdrMinSVGen[hbin])
   hdrSVGen[hbin].Fill(drSVGen_i, evtwt)       ;overUnderFlow(hdrSVGen[hbin])
   hdxy[hbin].Fill(dxy_i ,evtwt)               ;overUnderFlow(hdxy[hbin])
   hdxyerr[hbin].Fill(dxyerr_i ,evtwt)         ;overUnderFlow(hdxyerr[hbin])
   hdxysig[hbin].Fill(dxysig_i ,evtwt)         ;overUnderFlow(hdxysig[hbin])              
   hd3d[hbin].Fill(d3d_i ,evtwt)               ;overUnderFlow(hd3d[hbin])
   hd3derr[hbin].Fill(d3derr_i ,evtwt)         ;overUnderFlow(hd3derr[hbin])
   hd3dsig[hbin].Fill(d3dsig_i ,evtwt)         ;overUnderFlow(hd3dsig[hbin])             
   hcosTheta[hbin].Fill(cosTheta_i, evtwt)     ;overUnderFlow(hcosTheta[hbin])
   hndof[hbin].Fill(ndof_i, evtwt)             ;overUnderFlow(hndof[hbin])
   hndau[hbin].Fill(ndau_i, evtwt)             ;overUnderFlow(hndau[hbin])


def fillSVPropertiesNm1(drMinSVGen_i, drSVGen_i, dxy_i, dxyerr_i, dxysig_i, d3d_i, d3derr_i, d3dsig_i, cosTheta_i, ndof_i, ndau_i, hbin):
   
   if dxy[i]<3 and cosTheta[i]>0.98 and d3dsig[i]>4:
      hndau[hbin].Fill(ndau_i, evtwt)             ;overUnderFlow(hndau[hbin])
      
   if cosTheta[i]>0.98 and ndau[i] > 2 and d3dsig[i]>4:
      hdxy[hbin].Fill(dxy_i ,evtwt)               ;overUnderFlow(hdxy[hbin])
      
   if dxy[i]<3 and ndau[i] > 2 and d3dsig[i]>4:
      hcosTheta[hbin].Fill(cosTheta_i, evtwt)     ;overUnderFlow(hcosTheta[hbin])
      
   if dxy[i]<3 and cosTheta[i]>0.98 and ndau[i] > 2:
      hd3dsig[hbin].Fill(d3dsig_i ,evtwt)         ;overUnderFlow(hd3dsig[hbin])
            
   if dxy[i]<3 and cosTheta[i]>0.98 and ndau[i] > 2 and d3dsig[i]>4:
      hdrMinSVGen[hbin].Fill(drMinSVGen_i, evtwt) ;overUnderFlow(hdrMinSVGen[hbin])
      hdrSVGen[hbin].Fill(drSVGen_i, evtwt)       ;overUnderFlow(hdrSVGen[hbin])
      hdxyerr[hbin].Fill(dxyerr_i ,evtwt)         ;overUnderFlow(hdxyerr[hbin])
      hdxysig[hbin].Fill(dxysig_i ,evtwt)         ;overUnderFlow(hdxysig[hbin])              
      hd3d[hbin].Fill(d3d_i ,evtwt)               ;overUnderFlow(hd3d[hbin])
      hd3derr[hbin].Fill(d3derr_i ,evtwt)         ;overUnderFlow(hd3derr[hbin])
      hndof[hbin].Fill(ndof_i, evtwt)             ;overUnderFlow(hndof[hbin])

def fillSVP4Values(svPt_i, svEta_i, svPhi_i, svM_i, genPt_i, genEta_i, genPhi_i, genM_i, hbin):
   hSVPt[hbin].Fill(svPt_i, evtwt)         ;overUnderFlow(hSVPt[hbin])
   hSVEta[hbin].Fill(svEta_i, evtwt)       ;overUnderFlow(hSVEta[hbin])
   hSVPhi[hbin].Fill(svPhi_i, evtwt)       ;overUnderFlow(hSVPhi[hbin])
   hSVMass[hbin].Fill(svM_i, evtwt)        ;overUnderFlow(hSVMass[hbin])
   hSVPtGen[hbin].Fill(genPt_i, evtwt)     ;overUnderFlow(hSVPtGen[hbin])
   hSVEtaGen[hbin].Fill(genEta_i, evtwt)   ;overUnderFlow(hSVEtaGen[hbin])
   hSVPhiGen[hbin].Fill(genPhi_i, evtwt)   ;overUnderFlow(hSVPhiGen[hbin])
   hSVMassGen[hbin].Fill(genM_i, evtwt)    ;overUnderFlow(hSVMassGen[hbin])
   
def fillPtEtaBinnedCutFlow (hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt, svEta, binN, cat):
   bg = False; bnotg = False; other = False;
   if   cat == 'bg': bg = True
   elif cat == 'bnotg': bnotg = True
   elif cat == 'other':  other  = True
   else: raise Exception ('give me a resonable category to plot histogram')
        
   if svPt <= 5.0 and abs(svEta) <= 4.0:
      if bg:         hCutflow_pt[0].Fill(binN, evtwt)
      elif bnotg:    hCutflow_pt[3].Fill(binN, evtwt)
      elif other:     hCutflow_pt[6].Fill(binN, evtwt)      
   if svPt > 5.0  and svPt <= 10.0 and abs(svEta) <= 4.0:
      if bg:         hCutflow_pt[1].Fill(binN, evtwt)
      elif bnotg:    hCutflow_pt[4].Fill(binN, evtwt)
      elif other:     hCutflow_pt[7].Fill(binN, evtwt)      
   if svPt > 10.0 and svPt <= 20.0 and abs(svEta) <= 4.0:
      if bg:         hCutflow_pt[2].Fill(binN, evtwt)
      elif bnotg:    hCutflow_pt[5].Fill(binN, evtwt)
      elif other:     hCutflow_pt[8].Fill(binN, evtwt)

   if abs(svEta) <= 1.5 and svPt <= 20.0:
      if bg:         hCutflow_eta[0].Fill(binN, evtwt)
      elif bnotg:    hCutflow_eta[3].Fill(binN, evtwt)
      elif other:     hCutflow_eta[6].Fill(binN, evtwt)
   if abs(svEta) > 1.5 and abs(svEta) <= 2.0 and svPt <= 20.0:
      if bg:         hCutflow_eta[1].Fill(binN, evtwt)
      elif bnotg:    hCutflow_eta[4].Fill(binN, evtwt)
      elif other:     hCutflow_eta[7].Fill(binN, evtwt)
   if abs(svEta) > 2.0 and abs(svEta) <= 2.4 and svPt <= 20.0:
      if bg:         hCutflow_eta[2].Fill(binN, evtwt)
      elif bnotg:    hCutflow_eta[5].Fill(binN, evtwt)
      elif other:     hCutflow_eta[8].Fill(binN, evtwt)


   if svPt <= 5.0:                        ##--> pt block
      if abs(svEta) <= 1.5:                                   #'_pt_5_eta_1p5'
         if bg:         hCutflow_pt_eta[0].Fill(binN, evtwt)  
         elif bnotg:    hCutflow_pt_eta[9].Fill(binN, evtwt)  
         elif other:    hCutflow_pt_eta[18].Fill(binN, evtwt) 
      elif abs(svEta) > 1.5 and abs(svEta) <= 2.0 :           #'_pt_5_eta_2'
         if bg:         hCutflow_pt_eta[1].Fill(binN, evtwt)  
         elif bnotg:    hCutflow_pt_eta[10].Fill(binN, evtwt) 
         elif other:     hCutflow_pt_eta[19].Fill(binN, evtwt) 
      elif abs(svEta) > 2.0 and abs(svEta) <= 2.4 :           #'_pt_5_eta_2p4'
         if bg:         hCutflow_pt_eta[2].Fill(binN, evtwt)  
         elif bnotg:    hCutflow_pt_eta[11].Fill(binN, evtwt) 
         elif other:     hCutflow_pt_eta[20].Fill(binN, evtwt) 
   elif svPt > 5.0  and svPt <= 10.0:     ##--> pt block
      if abs(svEta) <= 1.5:                                   #'_pt_10_eta_1p5'
         if bg:         hCutflow_pt_eta[3].Fill(binN, evtwt) 
         elif bnotg:    hCutflow_pt_eta[12].Fill(binN, evtwt) 
         elif other:     hCutflow_pt_eta[21].Fill(binN, evtwt) 
      elif abs(svEta) > 1.5 and abs(svEta) <= 2.0 :           #'_pt_10_eta_2'
         if bg:         hCutflow_pt_eta[4].Fill(binN, evtwt) 
         elif bnotg:    hCutflow_pt_eta[13].Fill(binN, evtwt) 
         elif other:     hCutflow_pt_eta[22].Fill(binN, evtwt) 
      elif abs(svEta) > 2.0 and abs(svEta) <= 2.4 :           #'_pt_10_eta_2p4'
         if bg:         hCutflow_pt_eta[5].Fill(binN, evtwt) 
         elif bnotg:    hCutflow_pt_eta[14].Fill(binN, evtwt)
         elif other:     hCutflow_pt_eta[23].Fill(binN, evtwt)
   elif svPt > 10.0 and svPt <= 20.0:     ##-->pt block 
      if abs(svEta) <= 1.5:                                   #'_pt_20_eta_1p5'
         if bg:         hCutflow_pt_eta[6].Fill(binN, evtwt)
         elif bnotg: hCutflow_pt_eta[15].Fill(binN, evtwt)
         elif other:  hCutflow_pt_eta[24].Fill(binN, evtwt)
      elif abs(svEta) > 1.5 and abs(svEta) <= 2.0 :           #'_pt_20_eta_2'
         if bg:         hCutflow_pt_eta[7].Fill(binN, evtwt)
         elif bnotg:    hCutflow_pt_eta[16].Fill(binN, evtwt)
         elif other:     hCutflow_pt_eta[25].Fill(binN, evtwt)
      elif abs(svEta) > 2.0 and abs(svEta) <= 2.4 :           #'_pt_20_eta_2p4'
         if bg:         hCutflow_pt_eta[8].Fill(binN, evtwt)
         elif bnotg:    hCutflow_pt_eta[17].Fill(binN, evtwt)
         elif other:     hCutflow_pt_eta[26].Fill(binN, evtwt)


from optparse import OptionParser

# Create a command line option parser
options  = OptionParser()

options.add_option('--inDir', metavar='T', type='string', action='store',
                  default='input/TTJets', 
                  dest='inDir', 
                  help='input data directory name')
options.add_option('-f', '--files',  
                   dest="files", 
                   default="input/TTJets/TTJets.txt", #SMS-T2bb_mSbot-1650to2600.txt, TTJets.txt
                   type="string")
options.add_option('-n', '--maxEvts',  
                   dest="maxEvts", 
                   default=-1,
                   type="int")
options.add_option('-l', '--lowPtCut',  
                   dest="lowPtCut", 
                   default=5,
                   type="int")
options.add_option('--outFile', metavar='P', type='string', action='store',
                  default='out.root', 
                  dest='outFile',
                  help='output file')

(options,args) = options.parse_args()

print options

maxEvts = options.maxEvts
inDir = options.inDir
ptCut = options.lowPtCut
outFile = options.outFile
# Define the output histogram file
fdir = options.files
fname = fdir.rstrip().split('/')[1]
print fname
#fout = TFile(fname.replace('.txt', '_out.root'), 'RECREATE')
fout  = TFile(outFile, "RECREATE")
fout.cd()


cutsName = ['preselection', 'dR(gen,SV) #leq 0.05', 'IsCat', 'ndau > 2',  'cos#theta_{PV,SV} > 0.98', 'SIP3D > 4', 'IP2D < 3']
cutsNameUnMatch = ['preselection', 'dR(gen,SV) > 0.05', 'IsCat', 'IsCat', 'ndau > 2',  'cos#theta_{PV,SV} > 0.98', 'SIP3D > 4', 'IP2D < 3']


cat = ['_all', '_bg', '_bnotg', '_cg', '_cnotg', '_l', '_other', '_unmatched']

hCutflow_eta    = []; eta_range    = ['_eta_1p5', '_eta_2', '_eta_2p4'] 
hCutflow_pt     = []; pt_range     = ['_pt_5', '_pt_10', '_pt_20'] 
hCutflow_pt_eta = []; pt_eta_range = ['_pt_5_eta_1p5', '_pt_5_eta_2', '_pt_5_eta_2p4', '_pt_10_eta_1p5', '_pt_10_eta_2', '_pt_10_eta_2p4', '_pt_20_eta_1p5', '_pt_20_eta_2', '_pt_20_eta_2p4'] 

hSVPt     = []; pt_type       = ['Pt', 'Pt_f']
hSVEta    = []; eta_type      = ['Eta', 'Eta_f']
hSVPhi    = []; phi_type      = ['Phi', 'Phi_f']
hSVMass   = []; mass_type     = ['Mass', 'Mass_f']
hSVPtGen  = []; gen_pt_type   = ['Pt_gen', 'Pt_gen_f']
hSVEtaGen = []; gen_eta_type  = ['Eta_gen', 'Eta_gen_f']
hSVPhiGen = []; gen_phi_type  = ['Phi_gen', 'Phi_gen_f']
hSVMassGen= []; gen_mass_type = ['Mass_gen', 'Mass_gen_f']

hCutflow = [];

hdrMinSVGen = [];  drMinSVGen_type = ['drMinSVGen',  'drMinSVGen_f']
hdrSVGen    = [];  drSVGen_type    = ['drSVGen',     'drSVGen_f']
hcosTheta   = [];  cosTheta_type   = ['cosTheta',    'cosTheta_n1']
hdxy        = [];  dxy_type        = ['dxy',         'dxy_n1']
hdxyerr     = [];  dxyerr_type     = ['dxyerr',      'dxyerr_f']
hdxysig     = [];  dxysig_type     = ['dxysig',      'dxysig_f'] 
hd3d        = [];  d3d_type        = ['d3d',         'd3d_f']
hd3derr     = [];  d3derr_type     = ['d3derr',      'd3derr_f']
hd3dsig     = [];  d3dsig_type     = ['d3dsig',      'd3dsig_n1']
hndof       = [];  ndof_type       = ['ndof',        'ndof_f']
hndau       = [];  ndau_type       = ['ndau',        'ndau_n1']


var = ['d3dsig', 'cosTheta', 'dxy', 'Mass', 'Eta']
l = list(combinations(var, 2))
#print l      

#set 2D plots for all 10 combinations
#exit()

for c in cat:
   
   h = TH1D("hCutflow"+c ,";;Events;" , 7 , 0.5 ,7.5 )
   if 'unmatched' in c:
      setAxisBins(h, cutsNameUnMatch)      
   else:   
      setAxisBins(h, cutsName)
   hCutflow.append(h)

   for s in drMinSVGen_type:
      h = TH1D("h"+s+c,";;Events;" , 50 , 0 , 0.05 )
      hdrMinSVGen.append(h)

   for s in drSVGen_type:
      h = TH1D("h"+s+c,";;Events;" , 100 , 0 , 0.1)
      hdrSVGen.append(h)
      
   for s in cosTheta_type:
      if 'n1' in s:
         h = TH1D("h"+s+c ,";;Events;" , 100 , 0.5 ,1.05 )
      else:   
         h = TH1D("h"+s+c ,";;Events;" , 200 , -1.05 ,1.05 )
      hcosTheta.append(h)
      
   for s in dxy_type:
      h = TH1D("h"+s+c ,";;Events;" , 50 , 0 , 5 )
      hdxy.append(h)
      
   for s in dxyerr_type:
      h = TH1D("h"+s+c ,";;Events;" , 40 , 0 , 4 )
      hdxyerr.append(h)

   for s in dxysig_type:
      h = TH1D("h"+s+c ,";;Events;" , 100 , 0 , 20 )
      hdxysig.append(h)

   for s in d3d_type:   
      h = TH1D("h"+s+c ,";;Events;" , 50 , 0 , 5 )
      hd3d.append(h)

   for s in d3derr_type:   
      h = TH1D("h"+s+c ,";;Events;" , 40 , 0 , 4 )
      hd3derr.append(h)

   for s in d3dsig_type:   
      h = TH1D("h"+s+c ,";;Events;" , 100 , 0 , 50 )
      hd3dsig.append(h)

   for s in ndof_type:
      h = TH1D("h"+s+c ,";;Events;" , 8 , 0 , 8 )
      hndof.append(h)
      
   for s in ndau_type:
      h = TH1D("h"+s+c ,";;Events;" , 6 , 0 , 6 )
      hndau.append(h)
   
   '''
   for s in eta_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 6 , 0.5 ,6.5 )
      setAxisBins(h, cutsName)      
      hCutflow_eta.append (h)
        
   for s in pt_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 6 , 0.5 ,6.5 )
      setAxisBins(h, cutsName)      
      hCutflow_pt.append (h)

   for s in pt_eta_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 6 , 0.5 ,5.5 )
      setAxisBins(h, cutsName)      
      hCutflow_pt_eta.append (h)
   '''
####
   for s in pt_type:
      h =  TH1D("hSV"+s+c, "p_{T}; p_{T}(GeV); Events/1 GeV;", 20, 0, 20)
      hSVPt.append(h)

   for s in eta_type:
      h =  TH1D("hSV"+s+c, "#eta; #eta; Events/10 bins;", 80, -4.0, 4.0)
      hSVEta.append(h)

   for s in phi_type:
      h =  TH1D("hSV"+s+c, "#phi; #phi; Events/10 bins;", 80, -4.0, 4.0)
      hSVPhi.append(h)   

   for s in mass_type:
      h =  TH1D("hSV"+s+c, "mass; mass(GeV); Events GeV;", 100, 0.0, 4)
      hSVMass.append(h)
#
   for s in gen_pt_type:
      h =  TH1D("hSV"+s+c, "p_{T}; p_{T}(GeV); Events/1 GeV;", 20, 0, 20)
      hSVPtGen.append(h)

   for s in gen_eta_type:
      h =  TH1D("hSV"+s+c, "#eta; #eta; Events/10 bins;", 80, -4.0, 4.0)
      hSVEtaGen.append(h)
      
   for s in gen_phi_type:
      h =  TH1D("hSV"+s+c, "#phi; #phi; Events/10 bins;", 80, -4.0, 4.0)
      hSVPhiGen.append(h)
      
   for s in gen_mass_type:
      h =  TH1D("hSV"+s+c, "mass; mass(GeV); Events GeV;", 100, 0.0, 10)
      hSVMassGen.append(h)       

# 2D plots
   

# Open the input ntuples
fnames = [line.strip() for line in open(fdir, 'r')]
ievt = 0

for fname in fnames:
    if maxEvts > 0 and ievt > maxEvts: break
    #if ievt%1000 == 0: print " Processing evt %i" % ievt

    print 'Opening file %s' % fname
    f = TFile.Open(fname)
    print f.ls()

    tree = f.Get("makeNtuples/svtree")
    entries = tree.GetEntriesFast()
    nmatch = 0; nb = 0; total = 0;
    for t in tree:
       
       if maxEvts > 0 and ievt > maxEvts: break
       if ievt%1000 == 0: print " Processing evt %i" % ievt
       ievt += 1
       
       ncut = 0
       
       evtwt = t.GenEvt_genWt
       lhewt = t.GenEvt_lheWts
       
       nGoodPV = t.Evt_nGoodVtx
       if nGoodPV <= 1: continue
            
       
       nSV = t.Evt_nSV
       if nSV <= 1: continue
       
       svPt      = t.Evt_ptSV
       svEta     = t.Evt_etaSV
       svPhi     = t.Evt_phiSV
       svM       = t.Evt_massSV
       
       genPt     = t.Evt_matchedGenPtSV
       genEta    = t.Evt_matchedGenEtaSV
       genPhi    = t.Evt_matchedGenPhiSV
       genM      = t.Evt_matchedGenMassSV

       isSVgen   = t.Evt_isSVgen
       isSVbgen  = t.Evt_isSVb
       isSVb4rmg = t.Evt_isSVbg
       isSVc4rmg = t.Evt_isSVcg
       isSVcgen  = t.Evt_isSVc
       isSVlgen   = t.Evt_isSVl
       isSVogen   = t.Evt_isSVo

       nb        = t.Evt_nbsSV
       nc        = t.Evt_ncsSV
       nl        = t.Evt_nlsSV
       ng        = t.Evt_ngsSV
       no        = t.Evt_nothersSV
             
       dxy       = t.Evt_dxySV
       dxyerr    = t.Evt_dxyerrSV
       dxysig    = t.Evt_dxysigSV
       
       d3d      = t.Evt_d3dSV
       d3derr    = t.Evt_d3derrSV
       d3dsig    = t.Evt_d3dsigSV
       
       cosTheta  = t.Evt_costhetaSvPvSV
       ndau      = t.Evt_ndauSV
       ndof      = t.Evt_ndofSV
       normChi2  = t.Evt_normChi2SV
       chi2      = t.Evt_chi2SV
       
       drSVFlightJet = t.Evt_dRFlightJetSV
       drSVJet       = t.Evt_dRJetSV
       drMinSVGen    = t.Evt_matchedGenMindRSV
       drSVGen       = t.Evt_matchedGendRSV
       svMom         = t.Evt_matchedGenMomSV
       isLastCopyBeforeFSR  = t.Evt_matchedGenIsLastCopyBeforeFSR
       isHardProcess        = t.Evt_matchedGenIsHardProcess
       momID                = t.Evt_matchedGenMomSV

       
       #for i in range (0, len(isSVgen)):
       #   total = total + len(isSVgen)
       #   if isSVgen[i] == 1 : nmatch+1
       #   if isSVbgen[i] == 1: nb+1

        
       #print ('total: ', total, ', gen match: ', nmatch, ', b match: ', nb)
       
       
       for i in range(0, len(t.Evt_ptSV)):
          
          if svPt[i] > ptCut: continue
          
          hCutflow[0].Fill(1, evtwt)
          hCutflow[1].Fill(1, evtwt)
          hCutflow[2].Fill(1, evtwt)
          hCutflow[3].Fill(1, evtwt)
          hCutflow[4].Fill(1, evtwt)
          hCutflow[5].Fill(1, evtwt)
          hCutflow[6].Fill(1, evtwt)
          hCutflow[7].Fill(1, evtwt)
          
          # if dR(gen, SV) < 0.05
          if isSVgen[i] == 1:
             hCutflow[0].Fill(2, evtwt)
             hCutflow[1].Fill(2, evtwt)
             hCutflow[2].Fill(2, evtwt)
             hCutflow[3].Fill(2, evtwt)
             hCutflow[4].Fill(2, evtwt)
             hCutflow[5].Fill(2, evtwt)
             hCutflow[6].Fill(2, evtwt)
             hCutflow[7].Fill(2, evtwt)

             # category 0: fill histograms for all matched gen particles
             hCutflow[0].Fill(3, evtwt)              
             fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i],  0)
             fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 0)           
             # fill N-1 plots:
             fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 1)
             
             if ndau[i] > 2: 
                hCutflow[0].Fill(4, evtwt)
                # Fill all the 2D plots
                if cosTheta[i]>0.98:
                   hCutflow[0].Fill(5, evtwt) 
                   if d3dsig[i] > 4:
                      hCutflow[0].Fill(6, evtwt)
                      if dxy[i] < 3:
                         hCutflow[0].Fill(7, evtwt)
                         fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 1)
             
             # category 1: fill histograms if b is from glon splitting
             if isSVbgen[i] == 1 and isSVb4rmg[i] == 1:             
               hCutflow[1].Fill(3, evtwt)              
               fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i],  2)
               fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 2)           
               # fill N-1 plots:
               fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 3)
               
               if ndau[i] > 2:
                  hCutflow[1].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[1].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[1].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[1].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 3)
                                               
             # category 2: fill histograms if b is not from glon splitting                     
             elif isSVbgen[i] == 1 and isSVb4rmg[i] == 0:
                hCutflow[2].Fill(3, evtwt)
                fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 4)
                fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 4)   
                # fill N-1 plots here:
                fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 5)
                
                if ndau[i] > 2:
                   hCutflow[2].Fill(4, evtwt)               
                   if cosTheta[i]>0.98:
                     hCutflow[2].Fill(5, evtwt)                   
                     if d3dsig[i] > 4:
                        hCutflow[2].Fill(6, evtwt)                     
                        if dxy[i] < 3:
                           hCutflow[2].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 5)

             # category 3: fill histograms if c is from glon splitting
             elif isSVcgen[i] == 1 and isSVc4rmg[i] == 1:             
               hCutflow[3].Fill(3, evtwt)              
               fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i],  6)
               fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 6)           
               # fill N-1 plots:
               fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 7)
               
               if ndau[i] > 2:
                  hCutflow[3].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[3].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[3].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[3].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 7)
                                               
             # category 4: fill histograms if c is not from glon splitting                     
             elif isSVcgen[i] == 1 and isSVc4rmg[i] == 0:
                hCutflow[4].Fill(3, evtwt)
                fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 8)
                fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 8)   
                # fill N-1 plots here:
                fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 9)

                if ndau[i] > 2:
                  hCutflow[4].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[4].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[4].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[4].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 9)

             # category 5: fill histograms if match to light jet                     
             elif isSVlgen[i]:
                hCutflow[5].Fill(3, evtwt)
                fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 10)
                fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 10)   
                # fill N-1 plots here:
                fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 11)

                if ndau[i] > 2:
                  hCutflow[5].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[5].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[5].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[5].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 11)                
                                                    
             # category 6: fill histograms if SV is not matched to b,c,l                
             elif isSVogen[i] == 1:            
                hCutflow[6].Fill(3, evtwt)                
                fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 12)
                fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 12)
                # fill N-1 plots here:
                fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 13)

                if ndau[i] > 2:
                  hCutflow[6].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[6].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[6].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[6].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 13)
                           
          #category 7: fill histograms if SV is not matched any gen particle     
          else: 
             hCutflow[7].Fill(2, evtwt)
             hCutflow[7].Fill(3, evtwt)
             fillSVProperties(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 14)
             fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 14)             
             # fill N-1 plots here:
             fillSVPropertiesNm1(drMinSVGen[i], drSVGen[i], dxy[i], dxyerr[i], dxysig[i], d3d[i], d3derr[i], d3dsig[i], cosTheta[i], ndof[i], ndau[i], 15)

             if ndau[i] > 2:
                  hCutflow[7].Fill(4, evtwt)                               
                  if cosTheta[i]>0.98:
                     hCutflow[7].Fill(5, evtwt) 
                     if d3dsig[i] > 4:
                        hCutflow[7].Fill(6, evtwt)
                        if dxy[i] < 3:
                           hCutflow[7].Fill(7, evtwt)
                           fillSVP4Values(svPt[i], svEta[i], svPhi[i], svM[i], genPt[i], genEta[i], genPhi[i], genM[i], 15)
                           
             #print ('b: ', isSVbgen[i], ', c: ', isSVcgen[i], ', l: ', isSVlgen[i], ', o: ', isSVogen[i])
          #print '..No category to fill....'
fout.Write()
fout.Close()
