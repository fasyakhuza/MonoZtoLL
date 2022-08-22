#created by Fasya Khuzaimah on 2020.04.17
#ctau distribution

import ROOT
from ROOT import TFile, TTree, TH1F, TCanvas, TLegend, TAxis, TLatex, gStyle
import array as arr

gStyle.SetOptFit(111111)

openfile = TFile("histo_test.root", "read")
h_ctau_lab = openfile.Get("h_ctau_lab")
h_ctau_proper = openfile.Get("h_ctau_proper")


#draw ctau_lab histogram
c1 = TCanvas("c1","c1",900,700) #width-height
#c1.SetLeftMargin(0.15)
#gStyle.SetOptStat(0)

#leg = TLegend(0.65,0.7,0.85,0.87)
#leg.SetBorderSize(0)
#leg.SetTextSize(0.027)

#h_ctau_lab.SetAxisRange(0,5,"X")
h_ctau_lab.Draw()
h_ctau_lab.Fit("expo")


c1.cd()
c1.Modified()
c1.Update()
c1.SaveAs("test_ctau_lab_histo.pdf")


#draw ctau_proper histogram
c2 = TCanvas("c2","c2",900,700) #width-height
#c2.SetLeftMargin(0.15)
#gStyle.SetOptStat(0)

#leg = TLegend(0.65,0.7,0.85,0.87)
#leg.SetBorderSize(0)
#leg.SetTextSize(0.027)

#h_ctau_proper.SetAxisRange(0,3,"X")
h_ctau_proper.Draw()
h_ctau_proper.Fit("expo")


c2.cd()
c2.Modified()
c2.Update()
c2.SaveAs("test_ctau_proper_histo.pdf")