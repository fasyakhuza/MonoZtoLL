#created by Fasya Khuzaimah on 2022.06.28
#plot 3Dsig, 3Dsiglog, chi3D, alpha

import ROOT
from ROOT import TFile, TTree, TH1F, TH1D, TCanvas, TLegend, TAxis, TLatex, gStyle, TLorentzVector
import array as arr

#gStyle.SetOptFit(111111)
gStyle.SetOptStat(111111)

openfile1 = TFile("Output_TrackSelection_Mchi2-150_Mchi1-1_ctau1p0_test.root", "read")
openfile2 = TFile("Output_TrackSelection_Mchi2-150_Mchi1-1_ctau1p0_withTrackHighPurity_test.root", "read")
openfile3 = TFile("Mx2_150_old.root", "read")

#-----------#
#for h1 tree and Histograms from file1 and 2#
#-----------#
#File 1#
#h_eventWeight_withoutTHP = TH1F("h_eventWeight_withoutTHP", "", 6, -3, 3)

h_track3Dsig_withoutTHP = TH1F("h_track3Dsig_withoutTHP", "", 100, 0, 100)
h_Chi3D_withoutTHP = TH1F("h_Chi3D_withoutTHP", "", 100, 0, 100)

h_alpha_withoutTHP = TH1F("h_alpha_withoutTHP", "", 50, 0, 1.5)#1.5)
h_alpha2_withoutTHP = TH1F("h_alpha2_withoutTHP", "", 50, 0, 1.5)#1.5)
h_alpha3_withoutTHP = TH1F("h_alpha3_withoutTHP", "", 50, 0, 1.5)#1.5)
h_alpha4_withoutTHP = TH1F("h_alpha4_withoutTHP", "", 50, 0, 1.5)#1.5)

h_track3Dsiglog_withoutTHP = openfile1.Get("Histograms/h_track3Dsiglog")
h_trackChi3Dlog_withoutTHP = openfile1.Get("Histograms/h_trackChi3Dlog")

h_tracksdR_withoutTHP = TH1F("h_tracksdR_withoutTHP", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

#File 2#
#h_eventWeight_withTHP = TH1F("h_eventWeight_withTHP", "", 6, -3, 3)

h_track3Dsig_withTHP = TH1F("h_track3Dsig_withTHP", "", 100, 0, 100)
h_Chi3D_withTHP = TH1F("h_Chi3D_withTHP", "", 100, 0, 100)

h_alpha_withTHP = TH1F("h_alpha_withTHP", "", 50, 0, 1.5)#1.5)
h_alpha2_withTHP = TH1F("h_alpha2_withTHP", "", 50, 0, 1.5)#1.5)
h_alpha3_withTHP = TH1F("h_alpha3_withTHP", "", 50, 0, 1.5)#1.5)
h_alpha4_withTHP = TH1F("h_alpha4_withTHP", "", 50, 0, 1.5)#1.5)

h_track3Dsiglog_withTHP = openfile2.Get("Histograms/h_track3Dsiglog")
h_trackChi3Dlog_withTHP = openfile2.Get("Histograms/h_trackChi3Dlog")

h_tracksdR_withTHP = TH1F("h_tracksdR_withTHP", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

#File 3#
h_track3Dsig_2016 = TH1F("h_track3Dsig_2016", "", 100, 0, 100)
h_Chi3D_2016 = TH1F("h_Chi3D_2016", "", 100, 0, 100)

h_alpha_2016 = TH1F("h_alpha_2016", "", 50, 0, 1.5)#1.5)
h_alpha2_2016 = TH1F("h_alpha2_2016", "", 50, 0, 1.5)#1.5)
h_alpha3_2016 = TH1F("h_alpha3_2016", "", 50, 0, 1.5)#1.5)
h_alpha4_2016 = TH1F("h_alpha4_2016", "", 50, 0, 1.5)#1.5)

h_track3Dsiglog_2016 = TH1F("h_track3Dsiglog_2016", "", 50, -5, 5)
h_trackChi3Dlog_2016 = TH1F("h_trackChi3Dlog_2016", "", 50, -5, 5)

h_tracksdR_2016 = TH1F("h_tracksdR_2016P", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

#canvas
c1 = TCanvas("c1","c1",900,700) #width-height

#entries file 1
opentree1 = openfile1.Get("h1")
treeEvents1 = opentree1.GetEntries()

for jEntry1 in range(treeEvents1):
    opentree1.GetEntry(jEntry1)
    D_weight1 = getattr(opentree1, 'D_weight')
    #h_eventWeight.Fill(I_weight)

    v_track3Dsig1 = getattr(opentree1, 'v_track3Dsig')
    v_trackChi3D_1 = getattr(opentree1, 'v_trackChi3D')
    v_passjetalpha_1 = getattr(opentree1, 'v_passjetalpha')
    v_passjetalpha2_1 = getattr(opentree1, 'v_passjetalpha2')
    v_passjetalpha3_1 = getattr(opentree1, 'v_passjetalpha3')
    v_passjetalpha4_1 = getattr(opentree1, 'v_passjetalpha4')

    v_Trackdr1 = getattr(opentree1, 'v_Trackdr')

    #Plot 3Dsig, 3Dsiglog, chi3D, alpha
    for i in range(len(v_track3Dsig1)):
        h_track3Dsig_withoutTHP.Fill(v_track3Dsig1[i], D_weight1)

    for i in range(len(v_trackChi3D_1)):
        h_Chi3D_withoutTHP.Fill(v_trackChi3D_1[i], D_weight1)

    for i in range(len(v_passjetalpha_1)):
        h_alpha_withoutTHP.Fill(v_passjetalpha_1[i], D_weight1)

    for i in range(len(v_passjetalpha2_1)):
        h_alpha2_withoutTHP.Fill(v_passjetalpha2_1[i], D_weight1)
    
    for i in range(len(v_passjetalpha3_1)):
        h_alpha3_withoutTHP.Fill(v_passjetalpha3_1[i], D_weight1)

    for i in range(len(v_passjetalpha4_1)):
        h_alpha4_withoutTHP.Fill(v_passjetalpha4_1[i], D_weight1)

    for i in range(len(v_Trackdr1)):
        h_tracksdR_withoutTHP.Fill(v_Trackdr1[i], D_weight1)



#entries file 2
opentree2 = openfile2.Get("h1")
treeEvents2 = opentree2.GetEntries()

for jEntry2 in range(treeEvents2):
    opentree2.GetEntry(jEntry2)
    D_weight2 = getattr(opentree2, 'D_weight')
    #h_eventWeight.Fill(I_weight)

    v_track3Dsig2 = getattr(opentree2, 'v_track3Dsig')
    v_trackChi3D_2 = getattr(opentree2, 'v_trackChi3D')
    v_passjetalpha_2 = getattr(opentree2, 'v_passjetalpha')
    v_passjetalpha2_2 = getattr(opentree2, 'v_passjetalpha2')
    v_passjetalpha3_2 = getattr(opentree2, 'v_passjetalpha3')
    v_passjetalpha4_2 = getattr(opentree2, 'v_passjetalpha4')

    v_Trackdr2 = getattr(opentree2, 'v_Trackdr')

    #Plot 3Dsig, 3Dsiglog, chi3D, alpha
    for j in range(len(v_track3Dsig2)):
        h_track3Dsig_withTHP.Fill(v_track3Dsig2[j], D_weight2)

    for j in range(len(v_trackChi3D_2)):
        h_Chi3D_withTHP.Fill(v_trackChi3D_2[j], D_weight2)

    for j in range(len(v_passjetalpha_2)):
        h_alpha_withTHP.Fill(v_passjetalpha_2[j], D_weight2)

    for j in range(len(v_passjetalpha2_2)):
        h_alpha2_withTHP.Fill(v_passjetalpha2_2[j], D_weight2)
    
    for j in range(len(v_passjetalpha3_2)):
        h_alpha3_withTHP.Fill(v_passjetalpha3_2[j], D_weight2)

    for j in range(len(v_passjetalpha4_2)):
        h_alpha4_withTHP.Fill(v_passjetalpha4_2[j], D_weight2)

    for j in range(len(v_Trackdr2)):
        h_tracksdR_withTHP.Fill(v_Trackdr2[j], D_weight2)



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

    v_Trackdr_2016 = getattr(opentree3, 'v_Trackdr')

    #Plot 3Dsig, 3Dsiglog, chi3D, alpha
    for k in range(len(v_track3Dsig_2016)):
        h_track3Dsig_2016.Fill(v_track3Dsig_2016[k], I_weight3)

    for k in range(len(v_track3Dsiglog_2016)):
        h_track3Dsiglog_2016.Fill(v_track3Dsiglog_2016[k], I_weight3)

    for k in range(len(v_Chi3D_2016)):
        h_Chi3D_2016.Fill(v_Chi3D_2016[k], I_weight3)

    for k in range(len(v_trackChi3Dlog_2016)):
        h_trackChi3Dlog_2016.Fill(v_trackChi3Dlog_2016[k], I_weight3)

    for k in range(len(v_alpha_2016)):
        h_alpha_2016.Fill(v_alpha_2016[k], I_weight3)

    for k in range(len(v_alpha2_2016)):
        h_alpha2_2016.Fill(v_alpha2_2016[k], I_weight3)

    for k in range(len(v_alpha3_2016)):
        h_alpha3_2016.Fill(v_alpha3_2016[k], I_weight3)

    for k in range(len(v_alpha4_2016)):
        h_alpha4_2016.Fill(v_alpha4_2016[k], I_weight3)

    for k in range(len(v_Trackdr_2016)):
        h_tracksdR_2016.Fill(v_Trackdr_2016[k], I_weight3)


    
    

#stack each plot
#3Dsig
leg3Dsig = TLegend(0.4,0.7,0.65,0.87)
leg3Dsig.SetBorderSize(0)
leg3Dsig.SetTextSize(0.027)

h_track3Dsig_2016.SetLineColor(4)
h_track3Dsig_2016.GetXaxis().SetTitle("3D Significance")
h_track3Dsig_2016.GetYaxis().SetTitle("Events/Bin")
h_track3Dsig_2016.SetMaximum(7500)
leg3Dsig.AddEntry(h_track3Dsig_2016, "2016 MC", "l")

h_track3Dsig_withoutTHP.SetLineColor(809)#923
leg3Dsig.AddEntry(h_track3Dsig_withoutTHP, "Without THP", "l")

h_track3Dsig_withTHP.SetLineColor(813)
leg3Dsig.AddEntry(h_track3Dsig_withTHP, "With THP", "l")

h_track3Dsig_2016.Draw("hist")
h_track3Dsig_withoutTHP.Draw("histsame")
h_track3Dsig_withTHP.Draw("histsame")
leg3Dsig.Draw()
c1.Print("plots/3Dsig_plot.pdf")

#3Dsiglog
leg3Dsiglog = TLegend(0.2,0.7,0.45,0.87)
leg3Dsiglog.SetBorderSize(0)
leg3Dsiglog.SetTextSize(0.027)

h_track3Dsiglog_2016.SetLineColor(4)
h_track3Dsiglog_2016.GetXaxis().SetTitle("3D Significance")
h_track3Dsiglog_2016.GetYaxis().SetTitle("Events/Bin")
#h_track3Dsiglog_2016.SetMaximum(7500)
leg3Dsiglog.AddEntry(h_track3Dsiglog_2016, "2016 MC", "l")

h_track3Dsiglog_withoutTHP.SetLineColor(809)#923
leg3Dsiglog.AddEntry(h_track3Dsiglog_withoutTHP, "Without THP", "l")

h_track3Dsiglog_withTHP.SetLineColor(813)
leg3Dsiglog.AddEntry(h_track3Dsiglog_withTHP, "With THP", "l")

h_track3Dsiglog_2016.Draw("hist")
h_track3Dsiglog_withoutTHP.Draw("histsame")
h_track3Dsiglog_withTHP.Draw("histsame")
leg3Dsiglog.Draw()
c1.Print("plots/3Dsiglog_plot.pdf")

#chi3D
legChi3D = TLegend(0.4,0.7,0.65,0.87)
legChi3D.SetBorderSize(0)
legChi3D.SetTextSize(0.027)

h_Chi3D_2016.SetLineColor(4)
h_Chi3D_2016.GetXaxis().SetTitle("3D Significance")
h_Chi3D_2016.GetYaxis().SetTitle("Events/Bin")
h_Chi3D_2016.SetMaximum(8500)
legChi3D.AddEntry(h_Chi3D_2016, "2016 MC", "l")

h_Chi3D_withoutTHP.SetLineColor(809)#923
legChi3D.AddEntry(h_Chi3D_withoutTHP, "Without THP", "l")

h_Chi3D_withTHP.SetLineColor(813)
legChi3D.AddEntry(h_Chi3D_withTHP, "With THP", "l")

h_Chi3D_2016.Draw("hist")
h_Chi3D_withoutTHP.Draw("histsame")
h_Chi3D_withTHP.Draw("histsame")
legChi3D.Draw()
c1.Print("plots/chi3Dsig_plot.pdf")

#chi3Dlog
legChi3Dlog = TLegend(0.2,0.7,0.45,0.87)
legChi3Dlog.SetBorderSize(0)
legChi3Dlog.SetTextSize(0.027)

h_trackChi3Dlog_2016.SetLineColor(4)
h_trackChi3Dlog_2016.GetXaxis().SetTitle("3D Significance")
h_trackChi3Dlog_2016.GetYaxis().SetTitle("Events/Bin")
#h_trackChi3Dlog_2016.SetMaximum(7500)
legChi3Dlog.AddEntry(h_trackChi3Dlog_2016, "2016 MC", "l")

h_trackChi3Dlog_withoutTHP.SetLineColor(809)#923
legChi3Dlog.AddEntry(h_trackChi3Dlog_withoutTHP, "Without THP", "l")

h_trackChi3Dlog_withTHP.SetLineColor(813)
legChi3Dlog.AddEntry(h_trackChi3Dlog_withTHP, "With THP", "l")

h_trackChi3Dlog_2016.Draw("hist")
h_trackChi3Dlog_withoutTHP.Draw("histsame")
h_trackChi3Dlog_withTHP.Draw("histsame")
legChi3Dlog.Draw()
c1.Print("plots/chi3Dsiglog_plot.pdf")

#alpha
legAlpha = TLegend(0.4,0.7,0.65,0.87)
legAlpha.SetBorderSize(0)
legAlpha.SetTextSize(0.027)

h_alpha_2016.SetLineColor(4)
h_alpha_2016.GetXaxis().SetTitle("3D Significance")
h_alpha_2016.GetYaxis().SetTitle("Events/Bin")
#h_alpha_2016.SetMaximum(7500)
legAlpha.AddEntry(h_alpha_2016, "2016 MC", "l")

h_alpha_withoutTHP.SetLineColor(809)#923
legAlpha.AddEntry(h_alpha_withoutTHP, "Without THP", "l")

h_alpha_withTHP.SetLineColor(813)
legAlpha.AddEntry(h_alpha_withTHP, "With THP", "l")

h_alpha_2016.Draw("hist")
h_alpha_withoutTHP.Draw("histsame")
h_alpha_withTHP.Draw("histsame")
legAlpha.Draw()
c1.Print("plots/alpha_plot.pdf")

#alpa2
legAlpha2 = TLegend(0.2,0.7,0.45,0.87)
legAlpha2.SetBorderSize(0)
legAlpha2.SetTextSize(0.027)

h_alpha2_2016.SetLineColor(4)
h_alpha2_2016.GetXaxis().SetTitle("3D Significance")
h_alpha2_2016.GetYaxis().SetTitle("Events/Bin")
h_alpha2_2016.SetMaximum(5000)
legAlpha2.AddEntry(h_alpha2_2016, "2016 MC", "l")

h_alpha2_withoutTHP.SetLineColor(809)#923
legAlpha2.AddEntry(h_alpha2_withoutTHP, "Without THP", "l")

h_alpha2_withTHP.SetLineColor(813)
legAlpha2.AddEntry(h_alpha2_withTHP, "With THP", "l")

h_alpha2_2016.Draw("hist")
h_alpha2_withoutTHP.Draw("histsame")
h_alpha2_withTHP.Draw("histsame")
legAlpha2.Draw()
c1.Print("plots/alpha2_plot.pdf")

#alpha3
legAlpha3 = TLegend(0.2,0.7,0.45,0.87)
legAlpha3.SetBorderSize(0)
legAlpha3.SetTextSize(0.027)

h_alpha3_2016.SetLineColor(4)
h_alpha3_2016.GetXaxis().SetTitle("3D Significance")
h_alpha3_2016.GetYaxis().SetTitle("Events/Bin")
h_alpha3_2016.SetMaximum(8000)
legAlpha3.AddEntry(h_alpha3_2016, "2016 MC", "l")

h_alpha3_withoutTHP.SetLineColor(809)#923
legAlpha3.AddEntry(h_alpha3_withoutTHP, "Without THP", "l")

h_alpha3_withTHP.SetLineColor(813)
legAlpha3.AddEntry(h_alpha3_withTHP, "With THP", "l")

h_alpha3_2016.Draw("hist")
h_alpha3_withoutTHP.Draw("histsame")
h_alpha3_withTHP.Draw("histsame")
legAlpha3.Draw()
c1.Print("plots/alpha3_plot.pdf")

#alpha4
legAlpha4 = TLegend(0.2,0.7,0.45,0.87)
legAlpha4.SetBorderSize(0)
legAlpha4.SetTextSize(0.027)

h_alpha4_2016.SetLineColor(4)
h_alpha4_2016.GetXaxis().SetTitle("3D Significance")
h_alpha4_2016.GetYaxis().SetTitle("Events/Bin")
h_alpha4_2016.SetMaximum(9000)
legAlpha4.AddEntry(h_alpha4_2016, "2016 MC", "l")

h_alpha4_withoutTHP.SetLineColor(809)#923
legAlpha4.AddEntry(h_alpha4_withoutTHP, "Without THP", "l")

h_alpha4_withTHP.SetLineColor(813)
legAlpha4.AddEntry(h_alpha4_withTHP, "With THP", "l")

h_alpha4_2016.Draw("hist")
h_alpha4_withoutTHP.Draw("histsame")
h_alpha4_withTHP.Draw("histsame")
legAlpha4.Draw()
c1.Print("plots/alpha4_plot.pdf")

#dR1 and dR2
legdR = TLegend(0.2,0.7,0.45,0.87)
legdR.SetBorderSize(0)
legdR.SetTextSize(0.027)

h_tracksdR_2016.SetLineColor(4)
h_tracksdR_2016.GetXaxis().SetTitle("3D Significance")
h_tracksdR_2016.GetYaxis().SetTitle("Events/Bin")
h_tracksdR_2016.SetMaximum(3500)
legdR.AddEntry(h_tracksdR_2016, "2016 MC", "l")

h_tracksdR_withoutTHP.SetLineColor(809)#923
legdR.AddEntry(h_tracksdR_withoutTHP, "Without THP", "l")

h_tracksdR_withTHP.SetLineColor(813)
legdR.AddEntry(h_tracksdR_withTHP, "With THP", "l")

h_tracksdR_2016.Draw("hist")
h_tracksdR_withoutTHP.Draw("histsame")
h_tracksdR_withTHP.Draw("histsame")
legdR.Draw()
c1.Print("plots/deltaR_plot.pdf")