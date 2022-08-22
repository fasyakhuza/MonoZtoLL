#created by Fasya Khuzaimah on 2022.06.28
#plot 3Dsig, 3Dsiglog, chi3D, alpha

import ROOT
from ROOT import TFile, TTree, TH1F, TH1D, TCanvas, TLegend, TAxis, TLatex, gStyle, TLorentzVector
import array as arr

gStyle.SetOptFit(111111)

openfile1 = TFile("Output_TrackSelection_Mchi2-150_Mchi1-1_ctau1p0_test.root", "read")
openfile2 = TFile("Output_TrackSelection_Mchi2-150_Mchi1-1_ctau1p0_withTrackHighPurity_test.root", "read")
openfile3 = TFile("Mx2_150_old.root", "read")

#-----------#
#for h1 tree and Histograms from file1 and 2#
#-----------#
#File 1#
#h_eventWeight_withoutTHP = TH1D("h_eventWeight_withoutTHP", "", 6, -3, 3)

h_track3Dsig_withoutTHP = openfile1.Get("Histograms/h_track3Dsig")
h_track3Dsiglog_withoutTHP = openfile1.Get("Histograms/h_track3Dsiglog")

h_Chi3D_withoutTHP = openfile1.Get("Histograms/h_trackChi3D")
h_trackChi3Dlog_withoutTHP = openfile1.Get("Histograms/h_trackChi3Dlog")

h_alpha_withoutTHP = openfile1.Get("Histograms/h_passJetAlpha")
h_alpha2_withoutTHP = openfile1.Get("Histograms/h_passJetAlpha2")
h_alpha3_withoutTHP = openfile1.Get("Histograms/h_passJetAlpha3")
h_alpha4_withoutTHP = openfile1.Get("Histograms/h_passJetAlpha4")
h_alphaChi_withoutTHP = openfile1.Get("Histograms/h_passJetAlphaChi")

h_tracksdR_withoutTHP = TH1D("h_tracksdR_withoutTHP", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

#File 2#
#h_eventWeight_withTHP = TH1D("h_eventWeight_withTHP", "", 6, -3, 3)

h_track3Dsig_withTHP = openfile2.Get("Histograms/h_track3Dsig")
h_track3Dsiglog_withTHP = openfile2.Get("Histograms/h_track3Dsiglog")

h_Chi3D_withTHP = openfile2.Get("Histograms/h_trackChi3D")
h_trackChi3Dlog_withTHP = openfile2.Get("Histograms/h_trackChi3Dlog")

h_alpha_withTHP = openfile2.Get("Histograms/h_passJetAlpha")
h_alpha2_withTHP = openfile2.Get("Histograms/h_passJetAlpha2")
h_alpha3_withTHP = openfile2.Get("Histograms/h_passJetAlpha3")
h_alpha4_withTHP = openfile2.Get("Histograms/h_passJetAlphaChi")
h_alphaChi_withTHP = openfile2.Get("Histograms/h_passJetAlphaChi")

h_tracksdR_withTHP = TH1D("h_tracksdR_withTHP", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

#File 3#
h_track3Dsig_2016 = TH1D("h_track3Dsig_2016", "", 100, 0, 1000)
h_Chi3D_2016 = TH1D("h_Chi3D_2016", "", 100, 0, 1000)

h_alpha_2016 = TH1D("h_alpha_2016", "", 50, 0, 1)
h_alpha2_2016 = TH1D("h_alpha2_2016", "", 50, 0, 1)
h_alpha3_2016 = TH1D("h_alpha3_2016", "", 50, 0, 1)
h_alpha4_2016 = TH1D("h_alpha4_2016", "", 50, 0, 1)

h_track3Dsiglog_2016 = TH1D("h_track3Dsiglog_2016", "", 50, -5, 5)
h_trackChi3Dlog_2016 = TH1D("h_trackChi3Dlog_2016", "", 50, -5, 5)


#canvas
c1 = TCanvas("c1","c1",900,700) #width-height

#entries file 1
opentree1 = openfile1.Get("h1")
treeEvents1 = opentree1.GetEntries()

for jEntry1 in range(treeEvents1):
    opentree1.GetEntry(jEntry1)
    D_weight1 = getattr(opentree1, 'D_weight')
    #h_eventWeight.Fill(D_weight1)

    v_Trackdr1 = getattr(opentree1, 'v_Trackdr')

    for i in range(len(v_Trackdr1)):
        h_tracksdR_withoutTHP.Fill(v_Trackdr1[i], D_weight1)

#entries file 2
opentree2 = openfile2.Get("h1")
treeEvents2 = opentree2.GetEntries()

for jEntry2 in range(treeEvents2):
    opentree2.GetEntry(jEntry2)
    D_weight2 = getattr(opentree2, 'D_weight')
    #h_eventWeight.Fill(D_weight2)

    v_Trackdr2 = getattr(opentree2, 'v_Trackdr')

    for i in range(len(v_Trackdr2)):
        h_tracksdR_withTHP.Fill(v_Trackdr2[i], D_weight2)

#entries file 3
opentree3 = openfile3.Get("T_tree")
treeEvents3 = opentree3.GetEntries()

for jEntry3 in range(treeEvents3):
    opentree3.GetEntry(jEntry3)
    I_weight3 = getattr(opentree3, 'I_weight')

    v_track3Dsig_2016 = getattr(opentree3, 'v_Chi3D')
    v_track3Dsiglog_2016 = getattr(opentree3, 'v_Chi3Dlog')
    v_Chi3D_2016 = getattr(opentree3, 'v_Chi3DPaper')
    v_trackChi3Dlog_2016 = getattr(opentree3, 'v_Chi3DlogPaper')
    v_alpha_2016 = getattr(opentree3, 'v_fakealpha')
    v_alpha2_2016 = getattr(opentree3, 'v_fakealpha2')
    v_alpha3_2016 = getattr(opentree3, 'v_fakealpha3')
    v_alpha4_2016 = getattr(opentree3, 'v_fakealpha4')

    #Plot 3Dsig, 3Dsiglog, chi3D, alpha
    for i in range(len(v_track3Dsig_2016)):
        h_track3Dsig_2016.Fill(v_track3Dsig_2016[i], I_weight3)
    
    '''h_Chi3D_2016.Fill(v_Chi3D_2016, I_weight3)

    h_alpha_2016.Fill(v_alpha_2016, I_weight3)
    h_alpha2_2016.Fill(v_alpha2_2016, I_weight3)
    h_alpha3_2016.Fill(v_alpha3_2016, I_weight3)
    h_alpha4_2016.Fill(v_alpha4_2016, I_weight3)

    h_track3Dsiglog_2016.Fill(v_track3Dsiglog_2016, I_weight3)
    h_trackChi3Dlog_2016.Fill(v_trackChi3Dlog_2016, I_weight3)'''


#stack each plot
leg = TLegend(0.65,0.7,0.85,0.87)
leg.SetBorderSize(0)
leg.SetTextSize(0.027)

#3Dsig
h_track3Dsig_2016.SetLineColor(809)
h_track3Dsig_2016.GetXaxis().SetTitle("3D Significance")
h_track3Dsig_2016.GetYaxis().SetTitle("Events/Bin")
h_track3Dsig_2016.SetMaximum(50000)
leg.AddEntry(h_track3Dsig_2016, "2016 MC", "lep")

h_track3Dsig_withoutTHP.SetLineColor(821)#923
leg.AddEntry(h_track3Dsig_withoutTHP, "Without THP", "lep")

h_track3Dsig_withTHP.SetLineColor(813)
leg.AddEntry(h_track3Dsig_withTHP, "With THP", "lep")

h_track3Dsig_withoutTHP.Draw("hist")
h_track3Dsig_withTHP.Draw("histsame")
h_track3Dsig_2016.Draw("histsame")
leg.Draw()
c1.Print("trackVar_plots/3Dsig_plot.pdf")

#3Dsiglog
#alpha
#alpa2
#alpha3
#alpha4
#chi3D
#chi3Dlog
#dR1 and dR2