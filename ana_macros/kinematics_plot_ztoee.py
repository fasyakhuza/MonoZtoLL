#created by Fasya Khuzaimah on 2022.06.08
#kinematics distributions and efficiency in Z(ee)

import ROOT
from ROOT import TFile, TTree, TH1F, TH1D, TCanvas, TLegend, TAxis, TLatex, gStyle, TLorentzVector
import array as arr

gStyle.SetOptFit(111111)

openfile = TFile("Output_Preselection_Ztoee_Mchi2-150_Mchi1-1_ctau1p0_test.root", "read")

c1 = TCanvas("c1","c1",900,700) #width-height

#----------#
#for T_tree#
#----------#
h_eventWeight = TH1D("h_eventWeight", "", 6, -3, 3)

h_Zboson_Pt = TH1F("h_Zboson_Pt", "Z (ee) pT", 100, 0, 1000)
h_Zboson_Mass = TH1F("h_Zboson_Mass", "Z (ee) Mass", 100, 0, 200)
h_Zboson_Eta = TH1F("h_Zboson_Eta", "Z (ee) Eta", 50, -5, 5)
h_Zboson_Phi = TH1F("h_Zboson_Phi", "Z (ee) Phi", 50, -5, 5)

h_goodEle_Pt = TH1F("h_goodEle_Pt", "Ele pT from Tree", 100, 0, 1000)
h_goodEle_Mass = TH1F("h_goodEle_Mass", "Ele Mass from Tree", 50, 0, 0.5)
h_goodEle_Eta = TH1F("h_goodEle_Eta ", "Ele Eta from Tree", 50, -5, 5)
h_goodEle_Phi = TH1F("h_goodEle_Phi", "Ele Phi from Tree", 50, -5, 5)

h_goodMu_Pt = TH1F("h_goodMu_Pt", "Mu pT from Tree", 100, 0, 1000)
h_goodMu_Mass = TH1F("h_goodMu_Mass", "Mu Mass from Tree", 50, 0, 0.5)
h_goodMu_Eta = TH1F("h_goodMu_Eta", "Mu Eta from Tree", 50, -5, 5)
h_goodMu_Phi = TH1F("h_goodMu_Phi", "Mu Phi from Tree", 50, -5, 5)

h_goodTau_Pt = TH1F("h_goodTau_Pt", "Tau pT from Tree", 100, 0, 1000)
h_goodTau_Mass = TH1F("h_goodTau_Mass", "Tau Mass from Tree", 50, 0, 0.5)
h_goodTau_Eta = TH1F("h_goodTau_Eta", "Tau Eta from Tree", 50, -5, 5)
h_goodTau_Phi = TH1F("h_goodTau_Phi", "Tau Phi from Tree", 50, -5, 5)

#Tracks
h_tracksdR = TH1D("h_tracksdR", "deltaR between passAK4jet and its tracks", 100, -0.5, 0.5)

opentree = openfile.Get("T_tree")
treeEvents = opentree.GetEntries()

for jEntry in range(treeEvents):
    opentree.GetEntry(jEntry)
    I_weight = getattr(opentree, 'I_weight')
    h_eventWeight.Fill(I_weight)

    f_ZbosonPt = getattr(opentree, 'f_ZbosonPt')
    f_ZbosonMass = getattr(opentree, 'f_ZbosonMass')
    f_ZbosonEta = getattr(opentree, 'f_ZbosonEta')
    f_ZbosonPhi = getattr(opentree, 'f_ZbosonPhi')

    f_goodElePt = getattr(opentree, 'f_goodElePt')
    f_goodEleMass = getattr(opentree, 'f_goodEleMass')
    f_goodEleEta = getattr(opentree, 'f_goodEleEta')
    f_goodElePhi = getattr(opentree, 'f_goodElePhi')

    f_goodMuPt = getattr(opentree, 'f_goodMuPt')
    f_goodMuMass = getattr(opentree, 'f_goodMuMass')
    f_goodMuEta = getattr(opentree, 'f_goodMuEta')
    f_goodMuPhi = getattr(opentree, 'f_goodMuPhi')

    f_goodTauPt = getattr(opentree, 'f_goodTauPt')
    f_goodTauMass = getattr(opentree, 'f_goodTauMass')
    f_goodTauEta = getattr(opentree, 'f_goodTauEta')
    f_goodTauPhi = getattr(opentree, 'f_goodTauPhi')

    #tracks
    v_Trackdr = getattr(opentree, 'v_Trackdr')

    #Plot kinematics variable of leptons, Z boson from T_tree
    h_Zboson_Pt.Fill(f_ZbosonPt, I_weight)
    h_Zboson_Mass.Fill(f_ZbosonMass, I_weight)
    h_Zboson_Eta.Fill(f_ZbosonEta, I_weight)
    h_Zboson_Phi.Fill(f_ZbosonPhi, I_weight)

    h_goodEle_Pt.Fill(f_goodElePt, I_weight)
    h_goodEle_Mass.Fill(f_goodEleMass, I_weight)
    h_goodEle_Eta.Fill(f_goodEleEta, I_weight)
    h_goodEle_Phi.Fill(f_goodElePhi, I_weight)

    h_goodMu_Pt.Fill(f_goodMuPt, I_weight)
    h_goodMu_Mass.Fill(f_goodMuMass, I_weight)
    h_goodMu_Eta.Fill(f_goodMuEta, I_weight)
    h_goodMu_Phi.Fill(f_goodMuPhi, I_weight)

    h_goodTau_Pt.Fill(f_goodTauPt, I_weight)
    h_goodTau_Mass.Fill(f_goodTauMass, I_weight)
    h_goodTau_Eta.Fill(f_goodTauEta, I_weight)
    h_goodTau_Phi.Fill(f_goodTauPhi, I_weight)

    #tracks
    #h_tracksdR.Fill(v_Trackdr, I_weight)

h_Zboson_Pt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_Zboson_Pt.pdf")

h_Zboson_Mass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_Zboson_Mass.pdf")

h_Zboson_Eta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_Zboson_Eta.pdf")

h_Zboson_Phi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_Zboson_Phi.pdf")

h_goodEle_Pt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEle_Pt_fromTree.pdf")

h_goodEle_Mass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEle_Mass_fromTree.pdf")

h_goodEle_Eta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEle_Eta_fromTree.pdf")

h_goodEle_Phi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEle_Phi_fromTree.pdf")

h_goodMu_Pt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMu_Pt_fromTree.pdf")

h_goodMu_Mass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMu_Mass_fromTree.pdf")

h_goodMu_Eta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMu_Eta_fromTree.pdf")

h_goodMu_Phi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMu_Phi_fromTree.pdf")

h_goodTau_Pt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTau_Pt_fromTree.pdf")

h_goodTau_Mass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTau_Mass_fromTree.pdf")

h_goodTau_Eta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTau_Eta_fromTree.pdf")

h_goodTau_Phi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTau_Phi_fromTree.pdf")

#tracks
#h_tracksdR.Draw("hist")
#c1.Print("output_plot/recoeeEvent-and-Zboson/h_tracksdR.pdf")

#---------------#
#Event Variables#
#---------------#
h_totevent = openfile.Get("Event_Variable/h_totevent")
h_recoee_event = openfile.Get("Event_Variable/h_recoee_event")
h_ele_n = openfile.Get("Event_Variable/h_ele_n")
h_mu_n = openfile.Get("Event_Variable/h_mu_n")
h_tau_n = openfile.Get("Event_Variable/h_tau_n")
h_Zboson_n = openfile.Get("Event_Variable/h_Zboson_n")
h_ee_npass = openfile.Get("Event_Variable/h_ee_npass")

#Kinematics variables of leptons after each lepton selection
h_goodElePt = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodElePt")
h_goodEleMass = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodEleMass")
h_goodEleEta = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodEleEta")
h_goodElePhi = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodElePhi")
h_goodMuPt = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodMuPt")
h_goodMuMass = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodMuMass")
h_goodMuEta = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodMuEta")
h_goodMuPhi = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodMuPhi")
h_goodTauPt = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodTauPt")
h_goodTauMass = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodTauMass")
h_goodTauEta = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodTauEta")
h_goodTauPhi = openfile.Get("Kinematics_Variable_afterEachLeptonSelection/h_goodTauPhi")


#Plot Event Variables
h_totevent.Draw()
c1.Print("output_plot/recoeeEvent-and-Zboson/h_total_event.pdf")

h_recoee_event.Draw()
c1.Print("output_plot/recoeeEvent-and-Zboson/h_recoee_event.pdf")

h_ele_n.Draw()
c1.Print("output_plot/recoeeEvent-and-Zboson/number_of_good_ele.pdf")

h_Zboson_n.Draw()
c1.Print("output_plot/recoeeEvent-and-Zboson/number_of_Zboson.pdf")

h_ee_npass.Draw()
c1.Print("output_plot/recoeeEvent-and-Zboson/recoee_cutflow.pdf")


#Plot kinematics variables of leptons after each lepton selection
eventWeight = h_eventWeight.Integral()
print("eventWeight = ", eventWeight)

#h_goodElePt = h_goodElePt * h_eventWeight
h_goodElePt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodElePt_aftrEleSelection.pdf")

#h_goodEleMass = h_goodEleMass * h_eventWeight
h_goodEleMass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEleMass_aftrEleSelection.pdf")

#h_goodEleEta = h_goodEleEta * h_eventWeight
h_goodEleEta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodEleEta_aftrEleSelection.pdf")

#h_goodElePhi = h_goodElePhi * h_eventWeight
h_goodElePhi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodElePhi_aftrEleSelection.pdf")

#h_goodMuPt = h_goodMuPt * h_eventWeight
h_goodMuPt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMuPt_aftrMuSelection.pdf")

#h_goodMuPt = h_goodMuPt * h_eventWeight
h_goodMuMass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMuMass_aftrMuSelection.pdf")

#h_goodMuEta = h_goodMuEta * h_eventWeight
h_goodMuEta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMuEta_aftrMuSelection.pdf")

#h_goodMuPhi = h_goodMuPhi * h_eventWeight
h_goodMuPhi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodMuPhi_aftrMuSelection.pdf")

#h_goodTauPt = h_goodTauPt * h_eventWeight
h_goodTauPt.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTauPt_aftrTauSelection.pdf")

#h_goodTauMass = h_goodTauMass * h_eventWeight
h_goodTauMass.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTauMass_aftrTauSelection.pdf")

#h_goodTauEta = h_goodTauEta * h_eventWeight
h_goodTauEta.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTauEta_aftrTauSelection.pdf")

#h_goodTauPhi = h_goodTauPhi * h_eventWeight
h_goodTauPhi.Draw("hist")
c1.Print("output_plot/recoeeEvent-and-Zboson/h_goodTauPhi_aftrTauSelection.pdf")


