from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os


fnames = ["mu-_mom27_bite5_Ne01atm_gausbeam",
          "mu-_mom26_bite5_Ne01atm_gausbeam",
          "mu-_mom25_bite5_Ne01atm_gausbeam",
          "mu-_mom24_bite5_Ne01atm_gausbeam",
          "mu-_mom23_bite5_Ne01atm_gausbeam",
          "mu-_mom22_bite5_Ne01atm_gausbeam",
          "mu-_mom21_bite5_Ne01atm_gausbeam",
          "mu-_mom20_bite5_Ne01atm_gausbeam",
          "mu-_mom19_bite5_Ne01atm_gausbeam",
          "mu-_mom18_bite5_Ne01atm_gausbeam",
          "mu-_mom17_bite5_Ne01atm_gausbeam",
          "mu-_mom16_bite5_Ne01atm_gausbeam"]

ntgstop = np.zeros_like(fnames,dtype=np.float)

#beamI = np.array([6.0,4.0,2.8,2.3,1.7])# 20, 17, 16, 15, 14
#bscale = beamI/beamI[0]
bscale = np.ones_like(fnames,dtype=np.float)

nbin=100
minx=-10.
maxx=990.

fs = [ROOT.TFile.Open("./geant_"+fnames[i]+".root") for i in xrange(len(fnames))]
hs = [ROOT.TH1F("h%d"%(i),"h%d"%(i),nbin,minx,maxx) for i in xrange(len(fnames))]
htgs = [ROOT.TH1F("htg%d"%(i),"htg%d"%(i),nbin,minx,maxx) for i in xrange(len(fnames))]

ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("c1","c1")
ROOT.gPad.SetGrid(0,0);
ROOT.gPad.SetTicks();
legend = ROOT.TLegend(0.6, 0.40, 0.9, 0.80)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.02)
for i,(f,fname,h) in enumerate(zip(fs,fnames,hs)):
    f.tree.Draw("beamlastpos[0][2]>>h%d"%(i),"beampdg[0]==13 && beamproc[0]==60","GOFF")
    h.SetLineColor(30+i*2)
    h.Scale(bscale[i])
    if i==0:
        bw=h.GetBinWidth(1)
        h.SetTitle("mu- stopped position")
        h.GetXaxis().SetTitle("Z pos [mm]")
        h.GetYaxis().SetTitle("Counts / %.1f mm"%bw)
        h.GetYaxis().SetRangeUser(0,200)
        h.Draw("hist")
    else:
        h.Draw("histsame")
    legend.AddEntry(h.GetName(), fname, "l");
legend.Draw("same")
c1.SaveAs("mu-_stoppos_momdep_Ne01atm_gausbeam.pdf")



c2 = ROOT.TCanvas("c2","c2")
ROOT.gPad.SetGrid(0,0);
ROOT.gPad.SetTicks();
legend = ROOT.TLegend(0.5, 0.40, 0.9, 0.80)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.02)
for i,(f,fname,htg) in enumerate(zip(fs,fnames,htgs)):
    f.tree.Draw("tglastpos[0][2]>>htg%d"%(i),"tgpdg[0]==13 && tgproc[0]==60","GOFF")
    htg.SetLineColor(30+i*2)
    htg.Scale(bscale[i])
    ntgstop[i] = htg.GetEntries()
    if i==0:
        bw=htg.GetBinWidth(1)
        htg.SetTitle("mu- stopped position on NeGas target / 10k")
        htg.GetXaxis().SetTitle("Z pos [mm]")
        htg.GetYaxis().SetTitle("Counts / %.1f mm"%bw)
        htg.GetYaxis().SetRangeUser(0,200)
        htg.Draw("hist")
    else:
        htg.Draw("histsame")
    legend.AddEntry(htg.GetName(), fname+" stpr=%.2f"%(ntgstop[i]/1e4), "l");
legend.Draw("same")
c2.SaveAs("mu-_target_stoppos_momdep_Ne01atm_gausbeam.pdf")




