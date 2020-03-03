#!/usr/bin/env python
import os, sys
from ROOT import gROOT,std,ROOT,TFile,TTree,TH1D,TH2D,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")

def setAxisBins(h, name):
   ibin=0
   for n in name:
      ibin = ibin+1
      h.GetXaxis().SetBinLabel(ibin, n)
   return h

def fillPtEtaBinnedCutFlow (hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt, svEta, binN, b):
   
   if svPt <= 5.0 and abs(svEta) <= 4.0:
      if b: hCutflow_pt[0].Fill(binN, evtwt)
      else: hCutflow_pt[3].Fill(binN, evtwt)
   if svPt > 5.0  and svPt <= 10.0 and abs(svEta) <= 4.0:
      if b: hCutflow_pt[1].Fill(binN, evtwt)
      else: hCutflow_pt[4].Fill(binN, evtwt)
   if svPt > 10.0 and svPt <= 20.0 and abs(svEta) <= 4.0:
      if b: hCutflow_pt[2].Fill(binN, evtwt)
      else: hCutflow_pt[5].Fill(binN, evtwt)

   if abs(svEta) <= 1.5 and svPt <= 20.0:
      if b: hCutflow_eta[0].Fill(binN, evtwt)
      else: hCutflow_eta[3].Fill(binN, evtwt)
   if abs(svEta) > 1.5 and abs(svEta) <= 2.0 and svPt <= 20.0:
      if b: hCutflow_eta[1].Fill(binN, evtwt)
      else: hCutflow_eta[4].Fill(binN, evtwt)
   if abs(svEta) > 2.0 and abs(svEta) <= 2.4 and svPt <= 20.0:
      if b: hCutflow_eta[2].Fill(binN, evtwt)
      else: hCutflow_eta[5].Fill(binN, evtwt)

   if svPt <= 5.0 and abs(svEta) <= 1.5:
      if b: hCutflow_pt_eta[0].Fill(binN, evtwt)
      else: hCutflow_pt_eta[9].Fill(binN, evtwt)
   if svPt <= 5.0 and abs(svEta) > 1.5 and abs(svEta) <= 2.0 :
      if b: hCutflow_pt_eta[1].Fill(binN, evtwt)
      else: hCutflow_pt_eta[10].Fill(binN, evtwt)
   if svPt <= 5.0 and abs(svEta) > 2.0 and abs(svEta) <= 2.4 :
      if b: hCutflow_pt_eta[2].Fill(binN, evtwt)
      else: hCutflow_pt_eta[11].Fill(binN, evtwt)
   if svPt > 5.0  and svPt <= 10.0 and abs(svEta) <= 1.5:
      if b: hCutflow_pt_eta[3].Fill(binN, evtwt)
      else: hCutflow_pt_eta[12].Fill(binN, evtwt)
   if svPt > 5.0  and svPt <= 10.0 and abs(svEta) > 1.5 and abs(svEta) <= 2.0:
      if b: hCutflow_pt_eta[4].Fill(binN, evtwt)
      else: hCutflow_pt_eta[13].Fill(binN, evtwt)
   if svPt > 5.0  and svPt <= 10.0 and abs(svEta) > 2.0 and abs(svEta) <= 2.4 :
      if b: hCutflow_pt_eta[5].Fill(binN, evtwt)
      else: hCutflow_pt_eta[13].Fill(binN, evtwt)
   if svPt > 10.0 and svPt <= 20.0 and abs(svEta) <= 1.5:
      if b: hCutflow_pt_eta[6].Fill(binN, evtwt)
      else: hCutflow_pt_eta[14].Fill(binN, evtwt)
   if svPt > 10.0 and svPt <= 20.0 and abs(svEta) > 1.5 and abs(svEta) <= 2.0:
      if b: hCutflow_pt_eta[7].Fill(binN, evtwt)
      else: hCutflow_pt_eta[15].Fill(binN, evtwt)
   if svPt > 10.0 and svPt <= 20.0 and abs(svEta) > 2.0 and abs(svEta) <= 2.4 :  
      if b: hCutflow_pt_eta[8].Fill(binN, evtwt)
      else: hCutflow_pt_eta[16].Fill(binN, evtwt)
      
      
     #_pt_5_eta_1p5', '_pt_5_eta_2', '_pt_5_eta_2p4', '_pt_10_eta_1p5', '_pt_10_eta_2', '_pt_10_eta_2p4', '_pt_20_eta_1p5', '_pt_20_eta_2', '_pt_20_eta_2p4']  

from optparse import OptionParser

# Create a command line option parser
options  = OptionParser()

options.add_option('--inDir', metavar='T', type='string', action='store',
                  default='input', 
                  dest='inDir', 
                  help='input data directory name')
options.add_option('-f', '--files',  
                   dest="files", 
                   default="SMS-T2bb_mSbot-1650to2600.txt", #SMS-T2bb_mSbot-1650to2600.txt, TTJets.txt
                   type="string")
options.add_option('-n', '--maxEvts',  
                   dest="maxEvts", 
                   default=-1,
                   type="int")

(options,args) = options.parse_args()

print options

maxEvts = options.maxEvts
inDir = options.inDir

# Define the output histogram file
fdir = inDir+'/'+options.files
fname = fdir.rstrip().split('/')[1]
fout = TFile(fname.replace('.txt', '_out.root'), 'RECREATE')
fout.cd()


cutsName = ['preselection', 'IP2D < 3', 'cos#theta_{PV,SV} > 0.98', 'ndau > 2', 'SIP3D > 4']
#cutsName = ['#geq 1 good PV', '#geq 1 SV', 'preselection', '#DeltaR(SV,b)<0.4', 'DeltaR(SV,j)>0.4', 'IP2D < 3', 'cos#theta_{PV,SV} > 0.98', 'ndau > 2', 'SIP3D > 4']

bNOTb = ['_b', '_notb']

hCutflow_eta = [];    eta_range = ['_eta_1p5', '_eta_2', '_eta_2p4'] 
hCutflow_pt  = [];    pt_range  = ['_pt_5', '_pt_10', '_pt_20'] 
hCutflow_pt_eta = []; pt_eta_range = ['_pt_5_eta_1p5', '_pt_5_eta_2', '_pt_5_eta_2p4', '_pt_10_eta_1p5', '_pt_10_eta_2', '_pt_10_eta_2p4', '_pt_20_eta_1p5', '_pt_20_eta_2', '_pt_20_eta_2p4'] 

hSVPt  = []; pt_type = ['Pt', 'Pt_gen', 'Pt_f', 'Pt_gen_f']
hSVEta = []; eta_type= ['Eta', 'Eta_gen', 'Eta_f', 'Eta_gen_f']
hSVMass= []; mass_type = ['Mass', 'Mass_gen', 'Mass_f', 'Mass_gen_f']

hCutflow = []

for c in bNOTb:
   
   h = TH1D("hCutflow"+c ,";;Events;" , 5 , 0.5 ,5.5 )
   setAxisBins(h, cutsName)
   hCutflow.append(h)
   
   for s in eta_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 5 , 0.5 ,5.5 )
      setAxisBins(h, cutsName)      
      hCutflow_eta.append (h)
        
   for s in pt_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 5 , 0.5 ,5.5 )
      setAxisBins(h, cutsName)      
      hCutflow_pt.append (h)

   for s in pt_eta_range:
      h = TH1D("hCutflow"+s+c ,";;Events;" , 5 , 0.5 ,5.5 )
      setAxisBins(h, cutsName)      
      hCutflow_pt_eta.append (h)

   for s in pt_type:
      h =  TH1D("hSV"+s+c, "p_{T}; p_{T}(GeV); Events/1 GeV;", 20, 0, 20)
      hSVPt.append(h)

   for s in eta_type:
      h =  TH1D("hSV"+s+c, "#eta; #eta; Events/10 bins;", 80, -4.0, 4.0)
      hSVEta.append(h)

   for s in mass_type:
      h =  TH1D("hSV"+s+c, "mass; mass(GeV); Events/0.1 GeV;", 100, 0.0, 10)
      hSVMass.append(h)    
   

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
        
        isSVInJet = t.Evt_isMatchedToJetSV
        isSVb     = t.Evt_isSVb
        isSVc     = t.Evt_isSVc
        isSVg     = t.Evt_isSVgen
        dxySV     = t.Evt_dxySV
        cosTheta  = t.Evt_costhetaSvPvSV
        ndau      = t.Evt_ndauSV
        d3dsigSV  = t.Evt_d3dsigSV
        
        
        for i in range(0, len(t.Evt_ptSV)):

            # genCand matching to svCand is fixed to dR<=0.4, but it can tighten as per need
            svCand  = TLorentzVector(0.0, 0.0, 0.0, 0.0)
            svCand.SetPtEtaPhiM(svPt[i], svEta[i], svM[i], svPhi[i])

            genCand = TLorentzVector(0.0, 0.0, 0.0, 0.0)
            genCand.SetPtEtaPhiM(genPt[i], genEta[i], genM[i], genPhi[i])
                                                       
            if isSVInJet[i]==1: continue #matching to dR<=0.4
                         
            if isSVb[i] == 1:
               
               hSVPt[0].Fill(svPt[i], evtwt)
               hSVEta[0].Fill(svEta[i], evtwt)
               hSVMass[0].Fill(svM[i], evtwt)
               hSVPt[1].Fill(genPt[i], evtwt)
               hSVEta[1].Fill(genEta[i], evtwt)
               hSVMass[1].Fill(genM[i], evtwt)
               
            
               hCutflow[0].Fill(1, evtwt)
               
               fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 1, isSVb[i])
            
               if dxySV[i]<3:
                  hCutflow[0].Fill(2, evtwt)
                  fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 2, isSVb[i])
               
                  if cosTheta[i]<0.98:
                     hCutflow[0].Fill(3, evtwt)
                     fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 3, isSVb[i])
                  
                     if ndau[i] > 2:
                        hCutflow[0].Fill(4, evtwt)
                        fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 4, isSVb[i])
                     
                        if d3dsigSV[i]>4:
                           hCutflow[0].Fill(5, evtwt)
                           fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 5, isSVb[i])

                           hSVPt[2].Fill(svPt[i], evtwt)
                           hSVEta[2].Fill(svEta[i], evtwt)
                           hSVMass[2].Fill(svM[i], evtwt)
                           hSVPt[3].Fill(genPt[i], evtwt)
                           hSVEta[3].Fill(genEta[i], evtwt)
                           hSVMass[3].Fill(genM[i], evtwt)
                     
            if isSVb[i] == 0:
               
               hSVPt[4].Fill(svPt[i], evtwt) # 4, 5, 6, 7 for not b
               hSVEta[4].Fill(svEta[i], evtwt)
               hSVMass[4].Fill(svM[i], evtwt)
               hSVPt[5].Fill(genPt[i], evtwt)
               hSVEta[5].Fill(genEta[i], evtwt)
               hSVMass[5].Fill(genM[i], evtwt)                
               
               hCutflow[1].Fill(1, evtwt)
               fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 1, isSVb[i])
            
               if dxySV[i]<3:
                  hCutflow[1].Fill(2, evtwt)
                  fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 2, isSVb[i])
               
                  if cosTheta[i]<0.98:
                     hCutflow[1].Fill(3, evtwt)
                     fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 3, isSVb[i])
                  
                     if ndau[i] > 2:
                        hCutflow[1].Fill(4, evtwt)
                        fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 4, isSVb[i])
                     
                        if d3dsigSV[i]>4:
                           hCutflow[1].Fill(5, evtwt)
                           fillPtEtaBinnedCutFlow(hCutflow_pt, hCutflow_eta, hCutflow_pt_eta, svPt[i], svEta[i], 5, isSVb[i])
                           
                           hSVPt[6].Fill(svPt[i], evtwt)
                           hSVEta[6].Fill(svEta[i], evtwt)
                           hSVMass[6].Fill(svM[i], evtwt)
                           hSVPt[7].Fill(genPt[i], evtwt)
                           hSVEta[7].Fill(genEta[i], evtwt)
                           hSVMass[7].Fill(genM[i], evtwt)

        #print '....'    
        
fout.Write()
fout.Close()
