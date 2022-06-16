#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>
#include <TH1D.h>
#include <TH1F.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TAxis.h>
#include <math.h>
#include <string>
using namespace std;

void efferr(float nsig, float ntotal, float factor = 1)
{
    float eff = nsig / ntotal;
    float err = sqrt((1 - eff) * eff / ntotal);
    cout << "efficiency = " << eff * factor << " +- " << err * factor << endl;
}

bool pt_greater(const TLorentzVector a, const TLorentzVector b)
{
    double A = a.Pt();
    double B = b.Pt();
    return (A > B);
}

float median_value(vector<float> tmpvector)
{
    float med_value = 0.0;
    sort(tmpvector.begin(), tmpvector.end());
    
    if (tmpvector.size() % 2 == 0) //even
    {
        med_value = (tmpvector[tmpvector.size() / 2 - 1] + tmpvector[tmpvector.size() / 2]) / 2;
    }

    else //odd
    {
        med_value = tmpvector[tmpvector.size() / 2];
    }

    return (med_value);
}

float mean_value(const vector<float> &ttmpvector)
{
    float sum_value = 0.0;
    for (float x : ttmpvector)
    {
        sum_value += x;
    }
    return (sum_value / ttmpvector.size());
}

float cal_dphi(float phi1, float phi2)
{
    float dphi = phi1 - phi2;
    while (dphi >= TMath::Pi()) dphi -= 2 * TMath::Pi();
    while (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
    return TMath::Abs(dphi);
}

void xAna_ztoee_test(string inputfile = "ExoPieElementTuples_MonoZLL_NLO_Mphi-500_Mchi2-150_Mchi1-1_gSM-0p25_gDM-1p0_ctau1p0.root", string outputfile = "Output_Preselection_Ztoee_Mchi2-150_Mchi1-1_ctau1p0_test.root")
{

    //------------------
    // Create histrogram
    //------------------
    TH1D *h_genee_event = new TH1D("h_genee_event", "gen events", 5, 0, 5);
    h_genee_event->Sumw2();

    TH1D *h_recoee_event = new TH1D("h_recoee_event", "reco events", 5, 0, 5);
    h_recoee_event->Sumw2();

    TH1D *h_totevent = new TH1D("h_totevent", "total events", 5, 0, 5);
    h_totevent->Sumw2();

    TH1F *h_HT_eventCout = new TH1F("h_HT_eventCout", "", 10, 0, 10);
    h_HT_eventCout->SetYTitle("N event");
    h_HT_eventCout->Sumw2();

    TH1F *h_ee_npass = new TH1F("h_ee_npass", "", 10, 0, 10);
    h_ee_npass->SetXTitle("npass");
    h_ee_npass->Sumw2();

    TH1D *gen_chi2numb = new TH1D("gen_chi2numb", "Chi2", 5, 0, 5);
    TH1D *gen_dquarknumb = new TH1D("gen_dquarknumb", "d quark", 10, 0, 10);
    TH1D *gen_eenumber = new TH1D("gen_eenumber", "e", 10, 0, 10);
    TH1D *match_dquarknumb = new TH1D("match_dquarknumb", "d quark", 10, 0, 10);

    TH1D *h_ele_n = new TH1D("h_ele_n", "Nr of electron", 10, 0, 10);
    h_ele_n->GetXaxis()->SetTitle("Number of electron");
    h_ele_n->GetYaxis()->SetTitle("Number of Events");
    h_ele_n->Sumw2();

    TH1D *h_tau_n = new TH1D("h_tau_n", "Nr of tau", 10, 0, 10);
    h_tau_n->GetXaxis()->SetTitle("Number of tau");
    h_tau_n->GetYaxis()->SetTitle("Number of Events");
    h_tau_n->Sumw2();

    TH1D *h_mu_n = new TH1D("h_mu_n", "Nr of Muon", 10, 0, 10);
    h_mu_n->GetXaxis()->SetTitle("Number of muon");
    h_mu_n->GetYaxis()->SetTitle("Number of Events");
    h_mu_n->Sumw2();

    TH1D *h_Zboson_n = new TH1D("h_Zboson_n", "Nr of Z boson", 10, 0, 10);
    h_Zboson_n->GetXaxis()->SetTitle("Number of Z boson");
    h_Zboson_n->GetYaxis()->SetTitle("Number of Events");
    h_Zboson_n->Sumw2();

    TH1D *Z_eemass = new TH1D("Z_eemass", "Z->ee", 150, 0, 150);
    Z_eemass->GetXaxis()->SetTitle("Z mass");
    Z_eemass->GetYaxis()->SetTitle("Number of Events");
    Z_eemass->Sumw2();

    TH1D *h_jet_n = new TH1D("h_jet_n", "Nr of Matched jets", 15, 0, 15);
    h_jet_n->GetXaxis()->SetTitle("Number of Jets");
    h_jet_n->GetYaxis()->SetTitle("Number of Events");
    h_jet_n->Sumw2();

    TH1D *h_jet_rank = new TH1D("h_jet_rank", "Nr of rank Matched jets", 15, 0, 15);
    h_jet_rank->GetXaxis()->SetTitle("Jet_rank");
    h_jet_rank->GetYaxis()->SetTitle("Number of Events");
    h_jet_rank->Sumw2();

    TH1D *ratioTrackInferror = new TH1D("ratioTrackInferror", "", 10, 0, 1);
    ratioTrackInferror->Sumw2();

    //for goodleptons kinematics before reco-level
    TH1D *h_goodElePt = new TH1D("h_goodElePt", "Electrons pT after ele selection", 100, 0, 1000);
    h_goodElePt->GetXaxis()->SetTitle("Ele pT (GeV)");
    h_goodElePt->GetYaxis()->SetTitle("Number of Events");
    h_goodElePt->Sumw2();

    TH1D *h_goodEleMass = new TH1D("h_goodEleMass", "Electrons Mass after ele selection", 100, -0.5, 0.5);
    h_goodEleMass->GetXaxis()->SetTitle("Ele Mass (GeV)");
    h_goodEleMass->GetYaxis()->SetTitle("Number of Events");
    h_goodEleMass->Sumw2();

    TH1D *h_goodEleEta = new TH1D("h_goodEleEta", "Electrons Eta after ele selection", 100, -5.0, 5.0);
    h_goodEleEta->GetXaxis()->SetTitle("Ele Eta");
    h_goodEleEta->GetYaxis()->SetTitle("Number of Events");
    h_goodEleEta->Sumw2();

    TH1D *h_goodElePhi = new TH1D("h_goodElePhi", "Electrons Phi after ele selection", 100, -5.0, 5.0);
    h_goodElePhi->GetXaxis()->SetTitle("Ele Phi");
    h_goodElePhi->GetYaxis()->SetTitle("Number of Events");
    h_goodElePhi->Sumw2();

    TH1D *h_goodMuPt = new TH1D("h_goodMuPt", "Muons pT bfr after Mu selection", 100, 0, 1000);
    h_goodMuPt->GetXaxis()->SetTitle("Mu pT (GeV)");
    h_goodMuPt->GetYaxis()->SetTitle("Number of Events");
    h_goodMuPt->Sumw2();

    TH1D *h_goodMuMass = new TH1D("h_goodMuMass", "Muons Mass after Mu selection", 100, -0.5, 0.5);
    h_goodMuMass->GetXaxis()->SetTitle("Mu Mass (GeV)");
    h_goodMuMass->GetYaxis()->SetTitle("Number of Events");
    h_goodMuMass->Sumw2();

    TH1D *h_goodMuEta = new TH1D("h_goodMuEta", "Muons Eta after Mu selection", 100, -5.0, 5.0);
    h_goodMuEta->GetXaxis()->SetTitle("Mu Eta");
    h_goodMuEta->GetYaxis()->SetTitle("Number of Events");
    h_goodMuEta->Sumw2();

    TH1D *h_goodMuPhi = new TH1D("h_goodMuPhi", "Muons Phi after Mu selection", 100, -5.0, 5.0);
    h_goodMuPhi->GetXaxis()->SetTitle("Mu Phi");
    h_goodMuPhi->GetYaxis()->SetTitle("Number of Events");
    h_goodMuPhi->Sumw2();

    TH1D *h_goodTauPt = new TH1D("h_goodTauPt", "Good Taus pT after Tau selection", 100, 0, 1000);
    h_goodTauPt->GetXaxis()->SetTitle("Tau pT (GeV)");
    h_goodTauPt->GetYaxis()->SetTitle("Number of Events");
    h_goodTauPt->Sumw2();

    TH1D *h_goodTauMass = new TH1D("h_goodTauMass", "Good Taus Mass after Tau selection", 100, -0.5, 0.5);
    h_goodTauMass->GetXaxis()->SetTitle("Tau Mass (GeV)");
    h_goodTauMass->GetYaxis()->SetTitle("Number of Events");
    h_goodTauMass->Sumw2();

    TH1D *h_goodTauEta = new TH1D("h_goodTauEta", "Good Taus Eta after Tau selection", 100, -5.0, 5.0);
    h_goodTauEta->GetXaxis()->SetTitle("Tau Eta");
    h_goodTauEta->GetYaxis()->SetTitle("Number of Events");
    h_goodTauEta->Sumw2();

    TH1D *h_goodTauPhi = new TH1D("h_goodTauPhi", "Good Taus Phi after Tau selection", 100, -5.0, 5.0);
    h_goodTauPhi->GetXaxis()->SetTitle("Tau Phi");
    h_goodTauPhi->GetYaxis()->SetTitle("Number of Events");
    h_goodTauPhi->Sumw2();

    //for calculating efficiency
    TH1D *h_recoee_Vtxpass = new TH1D("h_recoee_Vtxpass", "Nr of Reco ee Events after nVtx selection", 10, 0, 10);
    h_recoee_Vtxpass->GetXaxis()->SetTitle("");
    h_recoee_Vtxpass->GetYaxis()->SetTitle("Number of Events");
    h_recoee_Vtxpass->Sumw2();

    TH1D *h_recoee_vetoTau = new TH1D("h_recoee_vetoTau", "Nr of Reco ee Events after vetoing Tau", 10, 0, 10);
    h_recoee_vetoTau->GetXaxis()->SetTitle("");
    h_recoee_vetoTau->GetYaxis()->SetTitle("Number of Events");
    h_recoee_vetoTau->Sumw2();

    TH1D *h_recoee_eePtpass = new TH1D("h_recoee_eePtpass", "Nr of Reco ee Events after ee pT selection", 10, 0, 10);
    h_recoee_eePtpass->GetXaxis()->SetTitle("");
    h_recoee_eePtpass->GetYaxis()->SetTitle("Number of Events");
    h_recoee_eePtpass->Sumw2();

    TH1D *h_recoee_deltaMasspass = new TH1D("h_recoee_deltaMasspass", "Nr of Reco ee Events after deltaMass selection", 10, 0, 10);
    h_recoee_deltaMasspass->GetXaxis()->SetTitle("");
    h_recoee_deltaMasspass->GetYaxis()->SetTitle("Number of Events");
    h_recoee_deltaMasspass->Sumw2();

    TH1D *h_recoee_nAK4pass = new TH1D("h_recoee_nAK4pass", "Nr of Reco ee Events with nAK4 gte 1", 10, 0, 10);
    h_recoee_nAK4pass->GetXaxis()->SetTitle("");
    h_recoee_nAK4pass->GetYaxis()->SetTitle("Number of Events");
    h_recoee_nAK4pass->Sumw2();

    //----------------------
    // Void Tree variable
    //----------------------
    Int_t I_event;
    Int_t I_weight;
    ULong64_t I_eventID;
    Float_t f_Met;
    //Float_t f_HT;
    //Float_t f_dileptonmass;

    Float_t f_goodElePt;
    Float_t f_goodEleMass;
    Float_t f_goodEleEta;
    Float_t f_goodElePhi;
    Float_t f_goodMuPt;
    Float_t f_goodMuMass;
    Float_t f_goodMuEta;
    Float_t f_goodMuPhi;
    Float_t f_goodTauPt;
    Float_t f_goodTauMass;
    Float_t f_goodTauEta;
    Float_t f_goodTauPhi;

    Float_t f_dileptonPT;
    Float_t f_dileptonMass;
    Float_t f_ZbosonPt;
    Float_t f_ZbosonMass;
    Float_t f_ZbosonEta;
    Float_t f_ZbosonPhi;
    vector<float> v_met_lepdeltaPhi;
    Int_t I_nThinJets;
    //vector<float> v_passJetindex;
    vector<float> v_passJetEta;
    vector<float> v_passJetPt;
    vector<float> v_passJetCSV;
    vector<int> v_passjethadronflavor;
    //vector<int> v_passjetpartonflavor;
    vector<int> v_nTrack;
    vector<float> v_TrackPT;
    vector<float> v_TrackEta;
    vector<float> v_Trackdr;
    //vector<float> v_Trackindex;
    vector<float> v_Trackdz;
    vector<float> v_Trackdzerror;
    vector<float> v_Trackdxy;
    vector<float> v_Trackdxyerror;
    
    vector<float> v_goodElePx;
    vector<float> v_goodElePy;
    vector<float> v_goodElePz;
    vector<float> v_goodEleE;
    vector<float> v_goodMuPx;
    vector<float> v_goodMuPy;
    vector<float> v_goodMuPz;
    vector<float> v_goodMuE;
    vector<float> v_goodTauPx;
    vector<float> v_goodTauPy;
    vector<float> v_goodTauPz;
    vector<float> v_goodTauE;

    TTree *T_tree = new TTree("T_tree", "Tree");
    T_tree->Branch("I_event", &I_event);
    T_tree->Branch("I_weight", &I_weight);
    T_tree->Branch("I_eventID", &I_eventID);
    T_tree->Branch("f_Met", &f_Met);
    //T_tree->Branch("f_HT", &f_HT);
    //T_tree->Branch("f_dileptonmass", &f_dileptonmass);
    T_tree->Branch("f_dileptonPT", &f_dileptonPT);
    T_tree->Branch("f_dileptonMass", &f_dileptonMass);
    T_tree->Branch("f_ZbosonPt", &f_ZbosonPt);
    T_tree->Branch("f_ZbosonMass", &f_ZbosonMass);
    T_tree->Branch("f_ZbosonEta", &f_ZbosonEta);
    T_tree->Branch("f_ZbosonPhi", &f_ZbosonPhi);
    
    T_tree->Branch("f_goodElePt", &f_goodElePt);
    T_tree->Branch("f_goodEleMass", &f_goodEleMass);
    T_tree->Branch("f_goodEleEta", &f_goodEleEta);
    T_tree->Branch("f_goodElePhi", &f_goodElePhi);
    T_tree->Branch("f_goodMuPt", &f_goodMuPt);
    T_tree->Branch("f_goodMuMass", &f_goodMuMass);
    T_tree->Branch("f_goodMuEta", &f_goodMuEta);
    T_tree->Branch("f_goodMuPhi", &f_goodMuPhi);
    T_tree->Branch("f_goodTauPt", &f_goodTauPt);
    T_tree->Branch("f_goodTauMass", &f_goodTauMass);
    T_tree->Branch("f_goodTauEta", &f_goodTauEta);
    T_tree->Branch("f_goodTauPhi", &f_goodTauPhi);
    
    T_tree->Branch("v_met_lepdeltaPhi", &v_met_lepdeltaPhi);
    T_tree->Branch("I_nThinJets", &I_nThinJets);
    //T_tree->Branch("v_passJetindex", &v_passJetindex);
    T_tree->Branch("v_passJetEta", &v_passJetEta);
    T_tree->Branch("v_passJetPt", &v_passJetPt);
    T_tree->Branch("v_passJetCSV", &v_passJetCSV);
    T_tree->Branch("v_passjethadronflavor", &v_passjethadronflavor);
    //T_tree->Branch("v_passjetpartonflavor", &v_passjetpartonflavor);
    T_tree->Branch("v_nTrack", &v_nTrack);
    T_tree->Branch("v_TrackPT", &v_TrackPT);
    T_tree->Branch("v_TrackEta", &v_TrackEta);
    T_tree->Branch("v_Trackdr", &v_Trackdr);
    //T_tree->Branch("v_Trackindex", &v_Trackindex);
    T_tree->Branch("v_Trackdz", &v_Trackdz);
    T_tree->Branch("v_Trackdzerror", &v_Trackdzerror);
    T_tree->Branch("v_Trackdxy", &v_Trackdxy);
    T_tree->Branch("v_Trackdxyerror", &v_Trackdxyerror);
    
    T_tree->Branch("v_goodElePx", &v_goodElePx);
    T_tree->Branch("v_goodElePy", &v_goodElePy);
    T_tree->Branch("v_goodElePz", &v_goodElePz);
    T_tree->Branch("v_goodEleE", &v_goodEleE);
    T_tree->Branch("v_goodMuPx", &v_goodMuPx);
    T_tree->Branch("v_goodMuPy", &v_goodMuPy);
    T_tree->Branch("v_goodMuPz", &v_goodMuPz);
    T_tree->Branch("v_goodMuE", &v_goodMuE);
    T_tree->Branch("v_goodTauPx", &v_goodTauPx);
    T_tree->Branch("v_goodTauPy", &v_goodTauPy);
    T_tree->Branch("v_goodTauPz", &v_goodTauPz);
    T_tree->Branch("v_goodTauE", &v_goodTauE);

    //get TTree from file ...
    TreeReader data(inputfile.data());

    Long64_t nTotal = 0;
    Long64_t nZboson = 0;
    int nEleBefore = 0;
    int nEleAfter = 0;
    int nEventAfterAllPreselection = 0;

    //Event Loop
    for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++)
    {
        if (jEntry % 1000 == 0)
        {
            fprintf(stderr, "Processing event %lli of %lli\n", jEntry+1, data.GetEntriesFast());
        }

        //-------------------
        // void some variable
        // clear Tree vector for each event
        //-------------------
        v_met_lepdeltaPhi.clear();
        //v_passJetindex.clear();
        v_passJetEta.clear();
        v_passJetPt.clear();
        v_passJetCSV.clear();
        v_passjethadronflavor.clear();
        //v_passjetpartonflavor.clear();
        v_nTrack.clear();
        v_TrackPT.clear();
        v_TrackEta.clear();
        v_Trackdr.clear();
        //v_Trackindex.clear();
        v_Trackdz.clear();
        v_Trackdzerror.clear();
        v_Trackdxy.clear();
        v_Trackdxyerror.clear();
        
        v_goodElePx.clear();
        v_goodElePy.clear();
        v_goodElePz.clear();
        v_goodEleE.clear();
        v_goodMuPx.clear();
        v_goodMuPy.clear();
        v_goodMuPz.clear();
        v_goodMuE.clear();
        v_goodTauPx.clear();
        v_goodTauPy.clear();
        v_goodTauPz.clear();
        v_goodTauE.clear();
        
        data.GetEntry(jEntry);
        nTotal ++;

        Float_t mcWeight = data.GetFloat("mcWeight");
        Double_t eventWeight = mcWeight;
        if (eventWeight > 0)
            {
                eventWeight = 1;
            }
        else if (eventWeight < 0)
            {
                eventWeight = -1;
            }
        else
            {
                eventWeight = 1;
            }

        //---------------------------
        // Store Total event number
        //---------------------------
        h_totevent->Fill(1.0, eventWeight);

        // Generator-Level
        bool matchee = false;
        vector<TLorentzVector> myEles;
        vector<TLorentzVector> dquark;
        vector<TLorentzVector> chi2s;
        dquark.clear();
        myEles.clear();
        chi2s.clear();
        // 0. check the generator-level information and make sure there is a Z->e+e-
        Int_t nGenPar = data.GetInt("nGenPar");
        Float_t* genParPx = data.GetPtrFloat("genParPx");
        Float_t* genParPy = data.GetPtrFloat("genParPy");
        Float_t* genParPz = data.GetPtrFloat("genParPz");
        Float_t* genParE = data.GetPtrFloat("genParE");
        Int_t *genParId = data.GetPtrInt("genParId");
        Int_t *genParSt = data.GetPtrInt("genParSt");
        Int_t *genMomParId = data.GetPtrInt("genMomParId");
        for (int ig = 0; ig < nGenPar; ig++)
        {
            TLorentzVector *thisGen = new TLorentzVector(genParPx[ig], genParPy[ig], genParPz[ig], genParE[ig]);
            int pid = genParId[ig];
            int mompid = genMomParId[ig];
            int status = genParSt[ig];
            if (abs(pid) == 11 && mompid == 23) //ele = 11 and Z is 23
            {
                matchee = true;
                myEles.push_back(*thisGen);
            }
            // chi2Id 18
            if (abs(pid) == 1 && abs(mompid) == 18) //dquark is 1 and chi2 is 18
            {
                dquark.push_back(*thisGen);
            }
            if (abs(pid) == 18)
            {
                chi2s.push_back(*thisGen);
            }
        }
        gen_dquarknumb->Fill(dquark.size(), eventWeight);
        gen_chi2numb->Fill(chi2s.size(), eventWeight);
        gen_eenumber->Fill(myEles.size(), eventWeight);

        if (matchee)
        {
            h_genee_event->Fill(1.0, eventWeight);

            // Starting reco-level
            // 1. electron
            Int_t nEle = data.GetInt("nEle");
            Float_t* elePx = data.GetPtrFloat("elePx");
            Float_t* elePy = data.GetPtrFloat("elePy");
            Float_t* elePz = data.GetPtrFloat("elePz");
            Float_t* eleEnergy = data.GetPtrFloat("eleEnergy");
            vector<bool>& eleIsPassLoose = *((vector<bool>*)data.GetPtr("eleIsPassLoose"));
            vector<bool>& eleIsPassMedium = *((vector<bool>*)data.GetPtr("eleIsPassMedium"));
            vector<bool>& eleIsPassVeto = *((vector<bool>*)data.GetPtr("eleIsPassVeto"));
            //Float_t *eleChHadIso = data.GetPtrFloat("eleChHadIso"); //can't find this information in the Tuples
            //Float_t *eleNeHadIso = data.GetPtrFloat("eleNeHadIso");
            //Float_t *eleGamIso = data.GetPtrFloat("eleGamIso");
            //Float_t *elePUPt = data.GetPtrFloat("elePUPt");

            vector<TLorentzVector> goodElectrons;
            goodElectrons.clear();
            
            for (int ie = 0; ie < nEle; ie++)
            {
                nEleBefore ++;
                TLorentzVector* myEle = new TLorentzVector(elePx[ie], elePy[ie], elePz[ie], eleEnergy[ie]);
                double myElePt = myEle->Pt();
                double myEleEta = myEle->Eta();
                if (myElePt < 20)
                {
                    continue;
                }
                if (fabs(myEleEta) > 2.4)
                {
                    continue;
                }
                if (!eleIsPassMedium[ie])
                {
                    continue;
                }
                goodElectrons.push_back(*myEle);
                
                h_goodElePt->Fill(myEle->Pt(), eventWeight);
                h_goodEleMass->Fill(myEle->M(), eventWeight);
                h_goodEleEta->Fill(myEle->Eta(), eventWeight);
                h_goodElePhi->Fill(myEle->Phi(), eventWeight);

                f_goodElePt = myEle->Pt();
                f_goodEleMass = myEle->M();
                f_goodEleEta = myEle->Eta();
                f_goodElePhi = myEle->Phi();
                v_goodElePx.push_back(elePx[ie]);
                v_goodElePy.push_back(elePy[ie]);
                v_goodElePz.push_back(elePz[ie]);
                v_goodEleE.push_back(eleEnergy[ie]);
                
                nEleAfter ++;
            } //End of Ele Loop

            h_ele_n->Fill(goodElectrons.size(), eventWeight);

            sort(goodElectrons.begin(), goodElectrons.end(), pt_greater);
            
            // 2. Muon
            Int_t nMu = data.GetInt("nMu");
            Float_t* muPx = data.GetPtrFloat("muPx");
            Float_t* muPy = data.GetPtrFloat("muPy");
            Float_t* muPz = data.GetPtrFloat("muPz");
            Float_t* muEnergy = data.GetPtrFloat("muEnergy");
            vector<bool>& isTightMuon = *((vector<bool>*)data.GetPtr("isTightMuon"));
            vector<bool>& isSoftMuon = *((vector<bool>*)data.GetPtr("isSoftMuon"));
            Float_t* muChHadIso = data.GetPtrFloat("muChHadIso");
            Float_t* muNeHadIso = data.GetPtrFloat("muNeHadIso");
            Float_t* muGamIso = data.GetPtrFloat("muGamIso");
            Float_t* muPUPt = data.GetPtrFloat("muPUPt");
            //Int_t *muTrkLayers = data.GetPtrInt("muTrkLayers"); //can't find this information in the Tuples

            vector<TLorentzVector> goodmuons;
            goodmuons.clear();

            for (int imu = 0; imu < nMu; imu++)
            {
                TLorentzVector* myMu = new TLorentzVector(muPx[imu], muPy[imu], muPz[imu], muEnergy[imu]);
                double myMuPt = myMu->Pt();
                double myMuEta = myMu->Eta();
                if (myMuPt < 20)
                {
                    continue;
                }
                if (fabs(myMuEta) > 2.4)
                {
                    continue;
                }

                double myMuIso = (muChHadIso[imu] + max(0., muNeHadIso[imu] + muGamIso[imu] - 0.5 * muPUPt[imu])) / myMuPt;
                if (myMuIso > 0.15)
                {
                    continue;
                }

                //if (muTrkLayers[imu] < 5) continue;
                if (!isSoftMuon[imu])
                //if (!isTightMuon[imu])
                {
                    continue;
                }

                goodmuons.push_back(*myMu);
                
                h_goodMuPt->Fill(myMu->Pt(), eventWeight);
                h_goodMuMass->Fill(myMu->M(), eventWeight);
                h_goodMuEta->Fill(myMu->Eta(), eventWeight);
                h_goodMuPhi->Fill(myMu->Phi(), eventWeight);

                f_goodMuPt = myMu->Pt();
                f_goodMuMass = myMu->M();
                f_goodMuEta = myMu->Eta();
                f_goodMuPhi = myMu->Phi();
                v_goodMuPx.push_back(muPx[imu]);
                v_goodMuPy.push_back(muPy[imu]);
                v_goodMuPz.push_back(muPz[imu]);
                v_goodMuE.push_back(muEnergy[imu]);
            } //End of Muon Loop

            h_mu_n->Fill(goodmuons.size(), eventWeight);

            bool recoeeEvent = false;
            bool recouuEvent = false;
            if (goodmuons.size() == goodElectrons.size())
            {
                continue;
            }
            if (goodElectrons.size() >= 2 && goodmuons.size() < 2)
            {
                recoeeEvent = true;
            }
            if (goodmuons.size() >= 2 && goodElectrons.size() < 2)
            {
                recouuEvent = true;
            }

            if (recoeeEvent)
            {
                h_recoee_event->Fill(1, eventWeight);
                h_ee_npass->Fill(1, eventWeight);

                // 3. Good Vertex
                int nVtx = data.GetInt("nVtx");
                if (nVtx < 1)
                    continue;
                h_recoee_Vtxpass->Fill(1, eventWeight);
                h_ee_npass->Fill(2, eventWeight);

                // 4. Tau Veto
                int nTau = data.GetInt("HPSTau_n");
                Float_t* tauPx = data.GetPtrFloat("HPSTau_Px");
                Float_t* tauPy = data.GetPtrFloat("HPSTau_Py");
                Float_t* tauPz = data.GetPtrFloat("HPSTau_Pz");
                Float_t* tauEnergy = data.GetPtrFloat("HPSTau_Energy");

                vector<bool> &disc_decayModeFinding = *((vector<bool>*)data.GetPtr("disc_decayModeFinding")); // DecayModeFinding method?
                vector<bool> &disc_decayModeFindingNewDMs = *((vector<bool>*)data.GetPtr("disc_decayModeFindingNewDMs"));
                //vector<bool> &disc_byVVTightIsolationMVA3newDMwLT = *((vector<bool> *)data.GetPtr("disc_byVVTightIsolationMVA3newDMwLT"));
                //vector<bool> &disc_byVTightIsolationMVA3newDMwLT = *((vector<bool> *)data.GetPtr("disc_byVTightIsolationMVA3newDMwLT"));
                //vector<bool> &disc_byTightIsolationMVA3newDMwLT = *((vector<bool> *)data.GetPtr("disc_byTightIsolationMVA3newDMwLT"));
                vector<bool>& disc_byVVTightDeepTau2017v2p1VSjet = *((vector<bool>*)data.GetPtr("disc_byVVTightDeepTau2017v2p1VSjet"));
                vector<bool>& disc_byVTightDeepTau2017v2p1VSjet = *((vector<bool>*)data.GetPtr("disc_byVTightDeepTau2017v2p1VSjet"));
                vector<bool>& disc_byTightDeepTau2017v2p1VSjet = *((vector<bool>*)data.GetPtr("disc_byTightDeepTau2017v2p1VSjet"));
            
                vector<TLorentzVector> goodtau;
                goodtau.clear();
                for (int it = 0; it < nTau; it++)
                {
                    TLorentzVector* myTau = new TLorentzVector(tauPx[it], tauPy[it], tauPz[it], tauEnergy[it]);
                    if (myTau->Pt() < 18)
                        continue;
                    if (fabs(myTau->Eta()) > 2.5)
                        continue;
                    if (!disc_decayModeFindingNewDMs[it])
                        continue;
                    if (!disc_byVTightDeepTau2017v2p1VSjet[it])
                        continue;
                    goodtau.push_back(*myTau);
                    
                    h_goodTauPt->Fill(myTau->Pt(), eventWeight);
                    h_goodTauMass->Fill(myTau->M(), eventWeight);
                    h_goodTauEta->Fill(myTau->Eta(), eventWeight);
                    h_goodTauPhi->Fill(myTau->Phi(), eventWeight);

                    f_goodTauPt = myTau->Pt();
                    f_goodTauMass = myTau->M();
                    f_goodTauEta = myTau->Eta();
                    f_goodTauPhi = myTau->Phi();
                    v_goodTauPx.push_back(tauPx[it]);
                    v_goodTauPy.push_back(tauPy[it]);
                    v_goodTauPz.push_back(tauPz[it]);
                    v_goodTauE.push_back(tauEnergy[it]);
                } //End loop of Tau
                h_tau_n->Fill(goodtau.size(), eventWeight);

                bool tauee = false;
                if (goodtau.size() > 0)
                {
                    for (int i = 0; i < goodtau.size(); i++)
                    {
                        float thisdr = goodtau[i].DeltaR(goodElectrons[0]);
                        float thatdr = goodtau[i].DeltaR(goodElectrons[1]);

                        if (thisdr > 0.4 && thatdr > 0.4)
                        {
                            tauee = true;
                            break;
                        }
                    } //End of Tau Veto Loop
                }
                if (tauee)
                    continue;
                h_recoee_vetoTau->Fill(1, eventWeight);
                h_ee_npass->Fill(3, eventWeight);

                // 5. MET
                Float_t met = data.GetFloat("pfMetCorrPt");
                Float_t metPhi = data.GetFloat("pfMetCorrPhi");

                for (auto i : goodElectrons)
                {
                    float lepPhi = i.Phi();
                    v_met_lepdeltaPhi.push_back(cal_dphi(metPhi, lepPhi)); 
                }

                // 6. Z boson
                float PDGZmass = 91.1876;
                TLorentzVector dilep = goodElectrons[0] + goodElectrons[1];
                float dilepPt = dilep.Pt();
                float dilepMass = dilep.M();
                if (goodElectrons[0].Pt() < 25 && goodElectrons[1].Pt() < 20)
                {
                    continue;
                }
                h_ee_npass->Fill(4, eventWeight);
                h_recoee_eePtpass->Fill(1, eventWeight);

                //TLorentzVector Z_boson_ee = goodElectrons[0] + goodElectrons[1];
                //float Z_boson_ee_pt = Z_boson_ee.Pt();
                float dilep_withZbosonMass = dilep.M(); //Z_boson_ee.M();
                float deltaMass = abs(PDGZmass - dilep_withZbosonMass);
                if (deltaMass > 15)
                {
                    continue;
                }
                h_ee_npass->Fill(5, eventWeight);
                h_recoee_deltaMasspass->Fill(1, eventWeight);
        
                //----------------------
                // To reduce diboson case (veto extra leptons)
                //----------------------
                if (goodElectrons.size() > 2 || goodmuons.size() > 2)
                {
                    continue;
                }
                TLorentzVector Zboson = goodElectrons[0] + goodElectrons[1];
                float Zboson_pt = Zboson.Pt();
                float Zboson_mass = Zboson.M();
                float Zboson_eta = Zboson.Eta();
                float Zboson_phi = Zboson.Phi();
                h_Zboson_n->Fill(1, eventWeight);
                h_ee_npass->Fill(6, eventWeight);
                nZboson ++;

                // 6. Thin Jet
                int nTHINjets = data.GetInt("THINnJet");
                Float_t* THINjetPx = data.GetPtrFloat("THINjetPx");
                Float_t* THINjetPy = data.GetPtrFloat("THINjetPy");
                Float_t* THINjetPz = data.GetPtrFloat("THINjetPz");
                Float_t* THINjetEnergy = data.GetPtrFloat("THINjetEnergy");
                Int_t* THINjetHadronFlavor = data.GetPtrInt("THINjetHadronFlavor");
                Float_t* THINjetCSV = data.GetPtrFloat("THINjetCISVV2");
                vector<int> indexForPassAK4;
                indexForPassAK4.clear();
                vector<TLorentzVector> goodTHINjets;
                goodTHINjets.clear();

                for (int itj = 0; itj < nTHINjets; itj++)
                {
                    TLorentzVector* thisTHINjet = new TLorentzVector(THINjetPx[itj], THINjetPy[itj], THINjetPz[itj], THINjetEnergy[itj]);
                    if (thisTHINjet->Pt() < 30)
                    {
                        continue;
                    }
                    if (fabs(thisTHINjet->Eta()) > 2.5)
                    {
                        continue;
                    }

                    bool thinJetee = false;
                    for (int j = 0; j < goodElectrons.size(); j++)
                    {
                        if (thisTHINjet->DeltaR(goodElectrons[j]) < 0.4)
                        {
                            thinJetee = true;
                            break;
                        }
                    }
                    if (thinJetee)
                    {
                        continue;
                    }

                    v_passJetPt.push_back(thisTHINjet->Pt());
                    v_passJetEta.push_back(thisTHINjet->Eta());
                    v_passJetCSV.push_back(THINjetCSV[itj]);
                    v_passjethadronflavor.push_back(THINjetHadronFlavor[itj]);
                    indexForPassAK4.push_back(itj);
                    h_jet_rank->Fill(itj, eventWeight);
                    goodTHINjets.push_back(*thisTHINjet);
                } // End of nTHINjets loop

                // Sort goodTHINjets by pT
                sort(goodTHINjets.begin(), goodTHINjets.end(), pt_greater);
                if (indexForPassAK4.size() == 0)
                    continue;
                h_recoee_nAK4pass->Fill(1, eventWeight);
                h_ee_npass->Fill(7, eventWeight);
                nEventAfterAllPreselection ++;

                //---------------------------
                // Study match Jet's tracks of Thin Jets
                //---------------------------
                Int_t* THINjetNTracks = data.GetPtrInt("THINjetNTracks");
                vector<float>* THINjetTrackImpdz = data.GetPtrVectorFloat("THINjetTrackImpdz", nTHINjets);
                vector<float>* THINjetTrackImpdzError = data.GetPtrVectorFloat("THINjetTrackImpdzError", nTHINjets);
                vector<float>* THINjetTrackImpdxy = data.GetPtrVectorFloat("THINjetTrackImpdxy", nTHINjets);
                vector<float>* THINjetTrackImpdxyError = data.GetPtrVectorFloat("THINjetTrackImpdxyError", nTHINjets);
                vector<float>* THINjetTrackPt = data.GetPtrVectorFloat("THINjetTrackPt", nTHINjets);
                vector<float>* THINjetTrackEta = data.GetPtrVectorFloat("THINjetTrackEta", nTHINjets);
                vector<float>* THINjetTrackPhi = data.GetPtrVectorFloat("THINjetTrackPhi", nTHINjets);
                for (int i = 0; i < indexForPassAK4.size(); i++)
                {
                    int correctjet_index = indexForPassAK4[i];
                    int Ntrack = 0;
                    for (int j = 0; j < THINjetNTracks[correctjet_index]; j++)
                    {
                        float dz = THINjetTrackImpdz[correctjet_index][j];
                        float dzError = THINjetTrackImpdzError[correctjet_index][j];
                        float dxy = THINjetTrackImpdxy[correctjet_index][j];
                        float dxyError = THINjetTrackImpdxyError[correctjet_index][j];
                        float trackPt = THINjetTrackPt[correctjet_index][j];
                        //----------------------
                        // Calculate DeltaR
                        //----------------------
                        float trackEta = THINjetTrackEta[correctjet_index][j];
                        float trackPhi = THINjetTrackPhi[correctjet_index][j];
                        TLorentzVector* thisJet = new TLorentzVector(THINjetPx[correctjet_index], THINjetPy[correctjet_index], THINjetPz[correctjet_index], THINjetEnergy[correctjet_index]);
                        float thinJetEta = thisJet->Eta();
                        float thinJetPhi = thisJet->Phi();
                        float deta = thinJetEta - trackEta;
                        float dphi = thinJetPhi - trackPhi;
                        while (dphi > TMath::Pi())
                            dphi -= 2 * TMath::Pi();
                        while (dphi <= -TMath::Pi())
                            dphi += 2 * TMath::Pi(); 
                        float dR = sqrt(deta * deta + dphi * dphi);
                        Ntrack++;
                        v_Trackdr.push_back(dR);
                        v_TrackPT.push_back(trackPt);
                        v_TrackEta.push_back(trackEta);
                        v_Trackdz.push_back(dz);
                        v_Trackdzerror.push_back(dzError);
                        v_Trackdxy.push_back(dxy);
                        v_Trackdxyerror.push_back(dxyError);
                    } // End of THINjetNTracks loop
                    v_nTrack.push_back(Ntrack);
                } // End of index of passed AK4 loop

                //Store the result in outputfile (root)
                //---------------------------
                //  Fill Tree event variable
                //---------------------------
                I_weight = eventWeight;
                ULong64_t eventId = data.GetLong64("eventId");
                I_eventID = eventId;
                f_dileptonPT = dilepPt;
                f_dileptonMass = dilepMass;
                f_ZbosonPt = Zboson_pt;
                f_ZbosonMass = Zboson_mass;
                f_ZbosonEta = Zboson_eta;
                f_ZbosonPhi = Zboson_phi;
                //f_HT = HT;
                f_Met = met;
                I_nThinJets = indexForPassAK4.size();
                T_tree->Fill();
            } // End of recoeeEvent Loop
        }     // End of matchee (gen ee)
    }         // End of Event Entries loop

    cout << "nEleBefore = " << nEleBefore << "\n";
    cout << "nEleAfter = " << nEleAfter << "\n";
    cout << "nZboson = " << nZboson << "\n";
    cout << "nEventAfterAllPreselection = " << nEventAfterAllPreselection << "\n";

    //--------------------//
    //Calculate efficiency//
    //-------------------//
    float ntotevent = h_totevent->Integral();
    cout << "number of weighted totEvents = " << ntotevent << "\n";
    float ngeneeEvent = h_genee_event->Integral();
    cout << "number of weighted gen-level ee Events = " << ngeneeEvent << "\n";
    cout << "efficiency of weighted gen-lev ee compared with totEvents:" << "\n";
    efferr(ngeneeEvent, ntotevent);
    
    float nrecoee_event = h_recoee_event->Integral();
    cout << "number of weighted reco-level ee Events = " << nrecoee_event << "\n";
    //cout << "efficiency of weighted reco-lev ee compared to totEvents:" << "\n";
    //efferr(nrecoee_event, ntotevent);
    cout << "efficiency of weighted reco-lev ee compared to gen-lev ee:" << "\n";
    efferr(nrecoee_event, ngeneeEvent);
    
    float nrecoee_event_Vtxpass = h_recoee_Vtxpass->Integral();
    cout << "number of weighted Events with nVtx > 1 = " << nrecoee_event_Vtxpass << "\n";
    cout << "efficiency of weighted Events passing nVtx compared to gen-lev ee:" << "\n";
    efferr(nrecoee_event_Vtxpass, ngeneeEvent);
    //cout << "efficiency of weighted Events passing nVtx compared to reco-lev ee:" << "\n";
    //efferr(nrecoee_event_Vtxpass, nrecoee_event);

    float nrecoee_event_tauVeto = h_recoee_vetoTau->Integral();
    cout << "number of weighted Events after vetoing tau = " << nrecoee_event_tauVeto << "\n";
    cout << "efficiency of weighted Events after vetoing tau compared to gen-lev ee:" << "\n";
    efferr(nrecoee_event_tauVeto, ngeneeEvent);
    //cout << "efficiency of weighted Events after vetoing tau compared to reco-lev ee:" << "\n";
    //efferr(nrecoee_event_tauVeto, nrecoee_event);

    float nrecoee_event_eePtpass = h_recoee_eePtpass->Integral();
    cout << "number of weighted Events after Pt selection = " << nrecoee_event_eePtpass << "\n";
    cout << "efficiency of weighted Events passing Pt selection compared to gen-lev ee:" << "\n";
    efferr(nrecoee_event_eePtpass, ngeneeEvent);
    //cout << "efficiency of weighted Events passing Pt selection compared to reco-lev ee:" << "\n";
    //efferr(nrecoee_event_eePtpass, nrecoee_event);

    float nrecoee_event_deltaMasspass = h_recoee_deltaMasspass->Integral();
    cout << "number of weighted Events after deltaMassZ selection = " << nrecoee_event_deltaMasspass << "\n";
    cout << "efficiency of weighted Events passing deltaMassZ selection compared to gen-lev ee:" << "\n";
    efferr(nrecoee_event_deltaMasspass, ngeneeEvent);
    //cout << "efficiency of weighted Events passing deltaMassZ selection compared to reco-lev ee:" << "\n";
    //efferr(nrecoee_event_deltaMasspass, nrecoee_event);

    float nZbosonEvent = h_Zboson_n->Integral();
    cout << "number of weighted Events after vetoing extra leptons (nZboson Events) = " << nZbosonEvent << "\n";
    //efficiency of Z boson
    cout << "efficiency of weighted Events after vetoing extra leptons (nZboson Events) compared to gen-lev ee:" << "\n";
    efferr(nZbosonEvent, ngeneeEvent);
    //cout << "efficiency of weighted Events after vetoing extra leptons (nZboson Events) compared to reco-lev ee:" << "\n";
    //efferr(nZbosonEvent, nrecoee_event);

    float nrecoee_events_nAK4pass = h_recoee_nAK4pass->Integral();
    cout << "number of weighted Events with AK4 >= 1 = " << nrecoee_events_nAK4pass << "\n";
    cout << "efficiency of weighted Events with AK4 >= 1 compared to gen-lev ee:" << "\n";
    efferr(nrecoee_events_nAK4pass, ngeneeEvent);
    //cout << "efficiency of weighted Events with AK4 >= 1 compared to reco-lev ee:" << "\n";
    //efferr(nrecoee_events_nAK4pass, nrecoee_event);

    // out Tree branches
    TFile *outFile = new TFile(outputfile.c_str(), "RECREATE");
    outFile->cd();
    T_tree->Write();
    outFile->mkdir("Event_Variable", "Event_Variable")->cd();
    h_totevent->Write();
    h_genee_event->Write();
    h_recoee_event->Write();
    h_ele_n->Write();
    h_mu_n->Write();
    h_tau_n->Write();
    h_Zboson_n->Write();
    //Z_eemass->Write();
    h_ee_npass->Write();
    gen_chi2numb->Write();
    gen_dquarknumb->Write();
    gen_eenumber->Write();
    //match_dquarknumb->Write();
    //h_HT_eventCout->Write();
    outFile->cd("/");
    outFile->mkdir("Jet_Variable", "Jet_Variable")->cd();
    h_jet_rank->Write();
    //h_jet_n->Write();
    outFile->cd("/");
    //outFile->mkdir("Track_Variable", "Track_Variable")->cd();
    //ratioTrackInferror->Write();
    //outFile->cd("/");
    outFile->mkdir("Kinematics_Variable_afterEachLeptonSelection", "Kinematics_Variable_afterEachLeptonSelection")->cd();
    h_goodElePt->Write();
    h_goodEleMass->Write();
    h_goodEleEta->Write();
    h_goodElePhi->Write();
    h_goodMuPt->Write();
    h_goodMuMass->Write();
    h_goodMuEta->Write();
    h_goodMuPhi->Write();
    h_goodTauPt->Write();
    h_goodTauMass->Write();
    h_goodTauEta->Write();
    h_goodTauPhi->Write();
    outFile->cd("/");
    outFile->Close();
} // big end
//Store the result in outputfile (root)