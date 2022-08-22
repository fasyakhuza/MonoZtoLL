#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <algorithm>
#include <TH1D.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
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

float median_value(vector<float> tmpvector)
{
    float med_value = 0.0;
    sort(tmpvector.begin(), tmpvector.end());

    if (tmpvector.size() % 2 == 0) // even
    {
        med_value = (tmpvector[tmpvector.size() / 2 - 1] + tmpvector[tmpvector.size() / 2]) / 2;
    }
    else // odd
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
//void xAna_checkDY_forJetBins(string inputtxtFilename = "DYlist_forCheck.txt")//, TString outputfile)
void xAna_checkDY(string inputtxtFilename = "DYlist_forCheckNew.txt")//, TString outputfile)
{

    Double_t lumi = 40500.0; // integrated luminosity; unit: pb^{-1}
    Double_t XS_DYincl = 6077.22; // unit: pb
    Double_t XS_DY1Jets = 1112.14; // unit: pb
    Double_t XS_DY2Jets = 344.987; // unit: pb
    Double_t XS_DY3Jets = 102.468; // unit: pb
    Double_t XS_DY4Jets = 48.8496; // unit: pb

    // For Initial
    TH1F* h_totevent = new TH1F("h_totevent", "", 5, 0, 5);
    h_totevent->Sumw2();

    TH1F* h_unnormalized_nAK4pass = new TH1F("h_unnormalized_nAK4pass", "", 5, 0, 5);
    h_unnormalized_nAK4pass->Sumw2();   

    TH1F *h_Numerator = new TH1F("h_Numerator", "", 5, 0, 5);
    h_Numerator->Sumw2();

    TH1F *h_ZbosonPt = new TH1F("h_ZbosonPt", "", 100, 0, 400);
    h_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_ZbosonPt->GetYaxis()->SetTitle("");
    h_ZbosonPt->Sumw2();

    TH1F *h_ZbosonEta = new TH1F("h_ZbosonEta", "", 50, -5, 5);
    h_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_ZbosonEta->GetYaxis()->SetTitle("");
    h_ZbosonEta->Sumw2();

    TH1F *h_jetPt = new TH1F("h_jetPt", "", 100, 0, 400);
    h_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_jetPt->GetYaxis()->SetTitle("");
    h_jetPt->Sumw2();

    TH1F *h_jetEta = new TH1F("h_jetEta", "", 50, -5, 5);
    h_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_jetEta->GetYaxis()->SetTitle("");
    h_jetEta->Sumw2();


    //DY Inclusive
    TH1F* h_DYinclDenominator = new TH1F("h_DYinclDenominator", "", 5, 0, 5);
    h_DYinclDenominator->Sumw2();

    TH1F* h_DYinclNumerator = new TH1F("h_DYinclNumerator", "", 5, 0, 5);
    h_DYinclNumerator->Sumw2();

    TH1F* h_DYincl_ZbosonPt = new TH1F("h_DYincl_ZbosonPt", "", 100, 0, 400);
    h_DYincl_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_DYincl_ZbosonPt->GetYaxis()->SetTitle("");
    h_DYincl_ZbosonPt->Sumw2();

    TH1F* h_DYincl_ZbosonEta = new TH1F("h_DYincl_ZbosonEta", "", 50, -5, 5);
    h_DYincl_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_DYincl_ZbosonEta->GetYaxis()->SetTitle("");
    h_DYincl_ZbosonEta->Sumw2();

    TH1F* h_DYincl_jetPt = new TH1F("h_DYincl_jetPt", "", 100, 0, 400);
    h_DYincl_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_DYincl_jetPt->GetYaxis()->SetTitle("");
    h_DYincl_jetPt->Sumw2();

    TH1F* h_DYincl_jetEta = new TH1F("h_DYincl_jetEta", "", 50, -5, 5);
    h_DYincl_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_DYincl_jetEta->GetYaxis()->SetTitle("");
    h_DYincl_jetEta->Sumw2();


    //DY1Jets
    TH1F* h_DY1JetsDenominator = new TH1F("h_DY1JetsDenominator", "", 5, 0, 5);
    h_DY1JetsDenominator->Sumw2();

    TH1F* h_DY1JetsNumerator = new TH1F("h_DY1JetsNumerator", "", 5, 0, 5);
    h_DY1JetsNumerator->Sumw2();

    TH1F* h_DY1Jets_ZbosonPt = new TH1F("h_DY1Jets_ZbosonPt", "", 100, 0, 400);
    h_DY1Jets_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_DY1Jets_ZbosonPt->GetYaxis()->SetTitle("");
    h_DY1Jets_ZbosonPt->Sumw2();

    TH1F* h_DY1Jets_ZbosonEta = new TH1F("h_DY1Jets_ZbosonEta", "", 50, -5, 5);
    h_DY1Jets_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_DY1Jets_ZbosonEta->GetYaxis()->SetTitle("");
    h_DY1Jets_ZbosonEta->Sumw2();

    TH1F* h_DY1Jets_jetPt = new TH1F("h_DY1Jets_jetPt", "", 100, 0, 400);
    h_DY1Jets_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_DY1Jets_jetPt->GetYaxis()->SetTitle("");
    h_DY1Jets_jetPt->Sumw2();

    TH1F* h_DY1Jets_jetEta = new TH1F("h_DY1Jets_jetEta", "", 50, -5, 5);
    h_DY1Jets_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_DY1Jets_jetEta->GetYaxis()->SetTitle("");
    h_DY1Jets_jetEta->Sumw2();


    //DY2Jets
    TH1F* h_DY2JetsDenominator = new TH1F("h_DY2JetsDenominator", "", 5, 0, 5);
    h_DY2JetsDenominator->Sumw2();

    TH1F* h_DY2JetsNumerator = new TH1F("h_DY2JetsNumerator", "", 5, 0, 5);
    h_DY2JetsNumerator->Sumw2();

    TH1F* h_DY2Jets_ZbosonPt = new TH1F("h_DY2Jets_ZbosonPt", "", 100, 0, 400);
    h_DY2Jets_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_DY2Jets_ZbosonPt->GetYaxis()->SetTitle("");
    h_DY2Jets_ZbosonPt->Sumw2();

    TH1F* h_DY2Jets_ZbosonEta = new TH1F("h_DY2Jets_ZbosonEta", "", 50, -5, 5);
    h_DY2Jets_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_DY2Jets_ZbosonEta->GetYaxis()->SetTitle("");
    h_DY2Jets_ZbosonEta->Sumw2();

    TH1F* h_DY2Jets_jetPt = new TH1F("h_DY2Jets_jetPt", "", 100, 0, 400);
    h_DY2Jets_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_DY2Jets_jetPt->GetYaxis()->SetTitle("");
    h_DY2Jets_jetPt->Sumw2();

    TH1F* h_DY2Jets_jetEta = new TH1F("h_DY2Jets_jetEta", "", 50, -5, 5);
    h_DY2Jets_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_DY2Jets_jetEta->GetYaxis()->SetTitle("");
    h_DY2Jets_jetEta->Sumw2();


    //DY3Jets
    TH1F* h_DY3JetsDenominator = new TH1F("h_DY3JetsDenominator", "", 5, 0, 5);
    h_DY3JetsDenominator->Sumw2();

    TH1F* h_DY3JetsNumerator = new TH1F("h_DY3JetsNumerator", "", 5, 0, 5);
    h_DY3JetsNumerator->Sumw2();

    TH1F* h_DY3Jets_ZbosonPt = new TH1F("h_DY3Jets_ZbosonPt", "", 100, 0, 400);
    h_DY3Jets_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_DY3Jets_ZbosonPt->GetYaxis()->SetTitle("");
    h_DY3Jets_ZbosonPt->Sumw2();

    TH1F* h_DY3Jets_ZbosonEta = new TH1F("h_DY3Jets_ZbosonEta", "", 50, -5, 5);
    h_DY3Jets_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_DY3Jets_ZbosonEta->GetYaxis()->SetTitle("");
    h_DY3Jets_ZbosonEta->Sumw2();

    TH1F* h_DY3Jets_jetPt = new TH1F("h_DY3Jets_jetPt", "", 100, 0, 400);
    h_DY3Jets_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_DY3Jets_jetPt->GetYaxis()->SetTitle("");
    h_DY3Jets_jetPt->Sumw2();

    TH1F* h_DY3Jets_jetEta = new TH1F("h_DY3Jets_jetEta", "", 50, -5, 5);
    h_DY3Jets_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_DY3Jets_jetEta->GetYaxis()->SetTitle("");
    h_DY3Jets_jetEta->Sumw2();


    //DY4Jets
    TH1F* h_DY4JetsDenominator = new TH1F("h_DY4JetsDenominator", "", 5, 0, 5);
    h_DY4JetsDenominator->Sumw2();

    TH1F* h_DY4JetsNumerator = new TH1F("h_DY4JetsNumerator", "", 5, 0, 5);
    h_DY4JetsNumerator->Sumw2();

    TH1F* h_DY4Jets_ZbosonPt = new TH1F("h_DY4Jets_ZbosonPt", "", 100, 0, 400);
    h_DY4Jets_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_DY4Jets_ZbosonPt->GetYaxis()->SetTitle("");
    h_DY4Jets_ZbosonPt->Sumw2();

    TH1F* h_DY4Jets_ZbosonEta = new TH1F("h_DY4Jets_ZbosonEta", "", 50, -5, 5);
    h_DY4Jets_ZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_DY4Jets_ZbosonEta->GetYaxis()->SetTitle("");
    h_DY4Jets_ZbosonEta->Sumw2();

    TH1F* h_DY4Jets_jetPt = new TH1F("h_DY4Jets_jetPt", "", 100, 0, 400);
    h_DY4Jets_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_DY4Jets_jetPt->GetYaxis()->SetTitle("");
    h_DY4Jets_jetPt->Sumw2();

    TH1F* h_DY4Jets_jetEta = new TH1F("h_DY4Jets_jetEta", "", 50, -5, 5);
    h_DY4Jets_jetEta->GetXaxis()->SetTitle("Jet Eta");
    h_DY4Jets_jetEta->GetYaxis()->SetTitle("");
    h_DY4Jets_jetEta->Sumw2();


    //summation
    TH1F* h_sumZbosonPt = new TH1F("h_sumZbosonPt", "", 100, 0, 400);
    h_sumZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    h_sumZbosonPt->GetYaxis()->SetTitle("");
    h_sumZbosonPt->Sumw2();

    TH1F* h_sumZbosonEta = new TH1F("h_sumZbosonEta", "", 50, -5, 5);
    h_sumZbosonEta->GetXaxis()->SetTitle("Z Eta");
    h_sumZbosonEta->GetYaxis()->SetTitle("");
    h_sumZbosonEta->Sumw2();

    TH1F* h_sumJetPt = new TH1F("h_sumJetPt", "", 100, 0, 400);
    h_sumJetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    h_sumJetPt->GetYaxis()->SetTitle("");
    h_sumJetPt->Sumw2();

    TH1F* h_sumJetEta = new TH1F("h_sumJetEta", "", 50, -5, 5);
    h_sumJetEta->GetXaxis()->SetTitle("Jet Eta");
    h_sumJetEta->GetYaxis()->SetTitle("");
    h_sumJetEta->Sumw2();


    cout << "inputtxtFilename = " << inputtxtFilename << endl;
    ifstream flist(inputtxtFilename.data());
    string inputFile;
    int filenumber = 0;
    while (getline(flist, inputFile))
    {
        h_totevent->Reset();
        h_unnormalized_nAK4pass->Reset();
        h_Numerator->Reset();
        h_ZbosonEta->Reset();
        h_ZbosonPt->Reset();
        h_jetPt->Reset();
        h_jetEta->Reset();

        filenumber++;
        //cout << inputFile << endl;
        
        TString myFile = inputFile;
        cout << "\n" << myFile << endl;

        TFile* file = TFile::Open(myFile);
        //cout << "successfully opened" << endl;

        // For get branches in T_tree
        Int_t I_weight;
        Int_t I_nThinJets;
        Float_t f_ZbosonPt;
        Float_t f_ZbosonMass;
        Float_t f_ZbosonEta;
        vector<int> *v_passJetIndex = new vector<int>();
        vector<float> *v_met = new vector<float>();
        vector<float> *v_passJetPt = new vector<float>();
        vector<float> *v_passJetEta = new vector<float>();

        TTree *Ttree;
        file->GetObject("T_tree", Ttree);
        Ttree->SetBranchAddress("I_weight", &I_weight);
        Ttree->SetBranchAddress("I_nThinJets", &I_nThinJets);
        Ttree->SetBranchAddress("f_ZbosonPt", &f_ZbosonPt);
        Ttree->SetBranchAddress("f_ZbosonMass", &f_ZbosonMass);
        Ttree->SetBranchAddress("f_ZbosonEta", &f_ZbosonEta);
        Ttree->SetBranchAddress("v_passJetIndex", &v_passJetIndex);
        Ttree->SetBranchAddress("v_met", &v_met);
        Ttree->SetBranchAddress("v_passJetPt", &v_passJetPt);
        Ttree->SetBranchAddress("v_passJetEta", &v_passJetEta);


        //TH1F* h_totevent = static_cast<TH1F*>(file->Get("Event_Variable/h_totevent"));
        //TH1F* h_unnormalized_nAK4pass = static_cast<TH1F*>(file->Get("Event_Variable/h_recoee_nAK4pass"));
        h_totevent = static_cast<TH1F*>(file->Get("Event_Variable/h_totevent"));
        h_unnormalized_nAK4pass = static_cast<TH1F*>(file->Get("Event_Variable/h_recoee_nAK4pass"));

        Double_t nProd = h_totevent->Integral();
        cout << "nProduction: " << nProd << endl;
        cout << "nAK4jets before jet cut and normalization: " << h_unnormalized_nAK4pass->Integral() << endl;


        //get the numerator and kinematics
        for (int evt = 0; evt < Ttree->GetEntries(); evt++)
        {
            Ttree->GetEntry(evt);

            Int_t mcWeight = I_weight;
            Double_t eventWeight = mcWeight;
            float* MET = v_met->data();

            if (I_nThinJets < 2)
                continue;
            //if (*MET < 140.0)
            //    continue;
            h_Numerator->Fill(1, eventWeight);
            for (int ij = 0; ij < I_nThinJets; ij++)
            {
                h_ZbosonPt->Fill(f_ZbosonPt, eventWeight);
                h_ZbosonEta->Fill(f_ZbosonEta, eventWeight);
                h_jetPt->Fill((*v_passJetPt)[ij], eventWeight);
                h_jetEta->Fill((*v_passJetEta)[ij], eventWeight);
            }
        } //end of getEntries

        if (inputFile.find("DYJetsToLL") != string::npos)
        {
            cout << "Processing " << inputFile << endl;

            //calculate the denominator
            Double_t norm = lumi * XS_DYincl / nProd;
            cout << "DYincl normalization: " << norm << endl;
            h_DYinclDenominator = (TH1F*)h_unnormalized_nAK4pass->Clone("h_DYinclDenominator");
            h_DYinclDenominator->Scale(norm);

            cout << "DYincl numerator before normalization: " << h_Numerator->Integral() << endl;
            h_DYinclNumerator = (TH1F*)h_Numerator->Clone("h_DYinclNumerator");
            h_DYinclNumerator->Scale(norm);

            h_DYincl_ZbosonPt = (TH1F*)h_ZbosonPt->Clone("h_DYincl_ZbosonPt");
            h_DYincl_ZbosonPt->Scale(norm);

            h_DYincl_ZbosonEta = (TH1F*)h_ZbosonEta->Clone("h_DYincl_ZbosonEta");
            h_DYincl_ZbosonEta->Scale(norm);

            h_DYincl_jetPt = (TH1F*)h_jetPt->Clone("h_DYincl_jetPt");
            h_DYincl_jetPt->Scale(norm);

            h_DYincl_jetEta = (TH1F*)h_jetEta->Clone("h_DYincl_jetEta");
            h_DYincl_jetEta->Scale(norm);

        }
        if (inputFile.find("DY1JetsToLL") != string::npos)
        {
            continue;
            /*cout << "Processing " << inputFile << endl;

            //calculate the denominator
            Double_t norm = lumi * XS_DY1Jets / nProd;
            cout << "DY1Jets normalization: " << norm << endl;
            h_DY1JetsDenominator = (TH1F*)h_unnormalized_nAK4pass->Clone("h_DY1JetsDenominator");
            h_DY1JetsDenominator->Scale(norm);

            h_DY1JetsNumerator = (TH1F*)h_Numerator->Clone("h_DY1JetsNumerator");
            h_DY1JetsNumerator->Scale(norm);

            h_DY1Jets_ZbosonPt = (TH1F*)h_ZbosonPt->Clone("h_DY1Jets_ZbosonPt");
            h_DY1Jets_ZbosonPt->Scale(norm);

            h_DY1Jets_ZbosonEta = (TH1F*)h_ZbosonEta->Clone("h_DY1Jets_ZbosonEta");
            h_DY1Jets_ZbosonEta->Scale(norm);

            h_DY1Jets_jetPt = (TH1F*)h_jetPt->Clone("h_DY1Jets_jetPt");
            h_DY1Jets_jetPt->Scale(norm);

            h_DY1Jets_jetEta = (TH1F*)h_jetEta->Clone("h_DY1Jets_jetEta");
            h_DY1Jets_jetEta->Scale(norm);*/

        } // close DY1JetsToLL

        if (inputFile.find("DY2JetsToLL") != string::npos)
        {
            cout << "Processing " << inputFile << endl;

            //calculate the denominator
            Double_t norm = lumi * XS_DY2Jets / nProd;
            cout << "DY2Jets normalization: " << norm << endl;
            //TH1F h_DY2JetsDenominator = norm * (*h_unnormalized_nAK4pass);
            h_DY2JetsDenominator = (TH1F*)h_unnormalized_nAK4pass->Clone("h_DY2JetsDenominator");
            h_DY2JetsDenominator->Scale(norm);

            cout << "DY2Jets numerator before normalization: " << h_Numerator->Integral() << endl;
            h_DY2JetsNumerator = (TH1F*)h_Numerator->Clone("h_DY2JetsNumerator");
            h_DY2JetsNumerator->Scale(norm);

            h_DY2Jets_ZbosonPt = (TH1F*)h_ZbosonPt->Clone("h_DY2Jets_ZbosonPt");
            h_DY2Jets_ZbosonPt->Scale(norm);

            h_DY2Jets_ZbosonEta = (TH1F*)h_ZbosonEta->Clone("h_DY2Jets_ZbosonEta");
            h_DY2Jets_ZbosonEta->Scale(norm);

            h_DY2Jets_jetPt = (TH1F*)h_jetPt->Clone("h_DY2Jets_jetPt");
            h_DY2Jets_jetPt->Scale(norm);

            h_DY2Jets_jetEta = (TH1F*)h_jetEta->Clone("h_DY2Jets_jetEta");
            h_DY2Jets_jetEta->Scale(norm);

        } // close DY2JetsToLL 

        if (inputFile.find("DY3JetsToLL") != string::npos)
        {
            cout << "Processing " << inputFile << endl;

            //calculate the denominator
            Double_t norm = lumi * XS_DY3Jets / nProd;
            cout << "DY3Jets normalization: " << norm << endl;
            //TH1F h_DY3JetsDenominator = norm * (*h_unnormalized_nAK4pass);
            h_DY3JetsDenominator = (TH1F*)h_unnormalized_nAK4pass->Clone("h_DY3JetsDenominator");
            h_DY3JetsDenominator->Scale(norm);

            cout << "DY3Jets numerator before normalization: " << h_Numerator->Integral() << endl;
            h_DY3JetsNumerator = (TH1F*)h_Numerator->Clone("h_DY3JetsNumerator");
            h_DY3JetsNumerator->Scale(norm);

            h_DY3Jets_ZbosonPt = (TH1F*)h_ZbosonPt->Clone("h_DY3Jets_ZbosonPt");
            h_DY3Jets_ZbosonPt->Scale(norm);

            h_DY3Jets_ZbosonEta = (TH1F*)h_ZbosonEta->Clone("h_DY3Jets_ZbosonEta");
            h_DY3Jets_ZbosonEta->Scale(norm);

            h_DY3Jets_jetPt = (TH1F*)h_jetPt->Clone("h_DY3Jets_jetPt");
            h_DY3Jets_jetPt->Scale(norm);

            h_DY3Jets_jetEta = (TH1F*)h_jetEta->Clone("h_DY3Jets_jetEta");
            h_DY3Jets_jetEta->Scale(norm);

        } // close DY3JetsToLL

        if (inputFile.find("DY4JetsToLL") != string::npos)
        {
            cout << "Processing " << inputFile << endl;

            //calculate the denominator
            Double_t norm = lumi * XS_DY4Jets / nProd;
            cout << "DY4Jets normalization: " << norm << endl;
            //TH1F h_DY4JetsDenominator = norm * (*h_unnormalized_nAK4pass);
            h_DY4JetsDenominator = (TH1F*)h_unnormalized_nAK4pass->Clone("h_DY4JetsDenominator");
            h_DY4JetsDenominator->Scale(norm);

            cout << "DY4Jets numerator before normalization: " << h_Numerator->Integral() << endl;
            h_DY4JetsNumerator = (TH1F*)h_Numerator->Clone("h_DY4JetsNumerator");
            h_DY4JetsNumerator->Scale(norm);

            h_DY4Jets_ZbosonPt = (TH1F*)h_ZbosonPt->Clone("h_DY4Jets_ZbosonPt");
            h_DY4Jets_ZbosonPt->Scale(norm);

            h_DY4Jets_ZbosonEta = (TH1F*)h_ZbosonEta->Clone("h_DY4Jets_ZbosonEta");
            h_DY4Jets_ZbosonEta->Scale(norm);

            h_DY4Jets_jetPt = (TH1F*)h_jetPt->Clone("h_DY4Jets_jetPt");
            h_DY4Jets_jetPt->Scale(norm);

            h_DY4Jets_jetEta = (TH1F*)h_jetEta->Clone("h_DY4Jets_jetEta");
            h_DY4Jets_jetEta->Scale(norm);

        } // close DY4JetsToLL
        
    }// end of flist
    
    //Inclusive DY efficiency
    float inclNumerator = h_DYinclNumerator->Integral();
    cout << "\ninclusive DY numerator: " << inclNumerator << endl;
    float inclDenominator = h_DYinclDenominator->Integral();
    cout << "inclusive DY denominator: " << inclDenominator << endl;

    cout << "efficiency (inclusive DY): " << endl;
    efferr(inclNumerator, inclDenominator);


    //DY+2-4jets efficiency
    //float sumNumerator = (h_DY1JetsNumerator->Integral()) + (h_DY2JetsNumerator->Integral()) + (h_DY3JetsNumerator->Integral()) + (h_DY4JetsNumerator->Integral());
    float sumNumerator = (h_DY2JetsNumerator->Integral()) + (h_DY3JetsNumerator->Integral()) + (h_DY4JetsNumerator->Integral());
    cout << "\n\nsum of all numerator: " << sumNumerator << endl;
    //float sumDenominator = (h_DY1JetsDenominator->Integral()) + (h_DY2JetsDenominator->Integral()) + (h_DY3JetsDenominator->Integral()) + (h_DY4JetsDenominator->Integral());
    float sumDenominator = (h_DY2JetsDenominator->Integral()) + (h_DY3JetsDenominator->Integral()) + (h_DY4JetsDenominator->Integral());
    cout << "sum of all denominator: " << sumDenominator << endl;

    cout << "efficiency (DY+2-4jets): " << endl;
    efferr(sumNumerator, sumDenominator);


    //canvas
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1","c1",900,700); //width-height

    //legend
    auto leg = new TLegend(0.75,0.6,0.85,0.7); //x1,y1,x2,y2
    leg->SetBorderSize(0);
    leg->SetTextSize(0.027);

    // sum the kinematics of all DYXJets samples
    h_sumZbosonPt = (TH1F*)h_DY2Jets_ZbosonPt->Clone("h_sumZbosonPt");
    h_sumZbosonPt->Add(h_DY3Jets_ZbosonPt);
    h_sumZbosonPt->Add(h_DY4Jets_ZbosonPt);
    
    h_sumZbosonEta = (TH1F*)h_DY2Jets_ZbosonEta->Clone("h_sumZbosonEta");
    h_sumZbosonEta->Add(h_DY3Jets_ZbosonEta);
    h_sumZbosonEta->Add(h_DY4Jets_ZbosonEta);

    h_sumJetPt = (TH1F*)h_DY2Jets_jetPt->Clone("h_sumJetPt");
    h_sumJetPt->Add(h_DY3Jets_jetPt);
    h_sumJetPt->Add(h_DY4Jets_jetPt);

    h_sumJetEta = (TH1F*)h_DY2Jets_jetEta->Clone("h_sumJetEta");
    h_sumJetEta->Add(h_DY3Jets_jetEta);
    h_sumJetEta->Add(h_DY4Jets_jetEta);


    //stack and draw histograms
    //Z boson pT
    leg->Clear();
    h_DYincl_ZbosonPt->SetLineColor(4); //blue
    h_DYincl_ZbosonPt->GetYaxis()->SetTitle("Events/Bin");
    h_DYincl_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    leg->AddEntry(h_DYincl_ZbosonPt, "Inclusive DY", "l");
    
    h_sumZbosonPt->SetLineColor(809); //red
    leg->AddEntry(h_sumZbosonPt, "DY+2-4jets", "l");
    
    h_DYincl_ZbosonPt->Draw("hist");
    h_sumZbosonPt->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Stack_ZbosonPt.pdf");

    //Z boson Eta
    leg->Clear();
    h_DYincl_ZbosonEta->SetLineColor(4); //blue
    h_DYincl_ZbosonEta->SetMaximum(70000);
    h_DYincl_ZbosonEta->GetYaxis()->SetTitle("Events/Bin");
    h_DYincl_ZbosonEta->GetXaxis()->SetTitle("Z #eta");
    leg->AddEntry(h_DYincl_ZbosonEta, "Inclusive DY", "l");

    h_sumZbosonEta->SetLineColor(809); //red
    leg->AddEntry(h_sumZbosonEta, "DY+2-4jets", "l");

    h_DYincl_ZbosonEta->Draw("hist");
    h_sumZbosonEta->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Stack_ZbosonEta.pdf");

    //Jet pT
    leg->Clear();
    h_DYincl_jetPt->SetLineColor(4); //blue
    h_DYincl_jetPt->GetYaxis()->SetTitle("Events/Bin");
    h_DYincl_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    leg->AddEntry(h_DYincl_jetPt, "Inclusive DY", "l");

    h_sumJetPt->SetLineColor(809); //red
    leg->AddEntry(h_sumJetPt, "DY+2-4jets", "l");

    h_DYincl_jetPt->Draw("hist");
    h_sumJetPt->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Stack_JetPt.pdf");    

    //Jet Eta
    leg->Clear();
    h_DYincl_jetEta->SetLineColor(4); //blue
    h_DYincl_jetEta->GetYaxis()->SetTitle("Events/Bin");
    h_DYincl_jetEta->GetXaxis()->SetTitle("Jet #eta");
    leg->AddEntry(h_DYincl_jetEta, "Inclusive DY", "l");

    h_sumJetEta->SetLineColor(809); //red
    leg->AddEntry(h_sumJetEta, "DY+2-4jets", "l");

    h_DYincl_jetEta->Draw("hist");
    h_sumJetEta->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Stack_JetEta.pdf");


    //scale each histo to 1 and stack the histo
    //Scaled Z boson pT
    leg->Clear();
    //cout << "integral of h_DYincl_ZbosonPt: " << h_DYincl_ZbosonPt->Integral() << endl;
    h_DYincl_ZbosonPt->Scale(1/(h_DYincl_ZbosonPt->Integral()));
    h_DYincl_ZbosonPt->SetLineColor(4); //blue
    h_DYincl_ZbosonPt->GetYaxis()->SetTitle("Scaling");
    h_DYincl_ZbosonPt->GetXaxis()->SetTitle("Z pT (GeV)");
    leg->AddEntry(h_DYincl_ZbosonPt, "Inclusive DY", "l");
    
    h_sumZbosonPt->Scale(1/(h_sumZbosonPt->Integral()));
    h_sumZbosonPt->SetLineColor(809); //red
    leg->AddEntry(h_sumZbosonPt, "DY+2-4jets", "l");
    
    h_DYincl_ZbosonPt->Draw("hist");
    h_sumZbosonPt->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Scaled_ZbosonPt.pdf");

    //Scaled Z boson Eta
    leg->Clear();
    h_DYincl_ZbosonEta->Scale(1/(h_DYincl_ZbosonEta->Integral()));
    h_DYincl_ZbosonEta->SetLineColor(4); //blue
    h_DYincl_ZbosonEta->SetMaximum(0.055);
    h_DYincl_ZbosonEta->GetYaxis()->SetTitle("Scaling");
    h_DYincl_ZbosonEta->GetXaxis()->SetTitle("Z #eta");
    leg->AddEntry(h_DYincl_ZbosonEta, "Inclusive DY", "l");

    h_sumZbosonEta->Scale(1/(h_sumZbosonEta->Integral()));
    h_sumZbosonEta->SetLineColor(809); //red
    leg->AddEntry(h_sumZbosonEta, "DY+2-4jets", "l");

    h_DYincl_ZbosonEta->Draw("hist");
    h_sumZbosonEta->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Scaled_ZbosonEta.pdf");

    //Scaled Jet pT
    leg->Clear();
    h_DYincl_jetPt->Scale(1/(h_DYincl_jetPt->Integral()));
    h_DYincl_jetPt->SetLineColor(4); //blue
    h_DYincl_jetPt->GetYaxis()->SetTitle("Scaling");
    h_DYincl_jetPt->GetXaxis()->SetTitle("Jet pT (GeV)");
    leg->AddEntry(h_DYincl_jetPt, "Inclusive DY", "l");

    h_sumJetPt->Scale(1/(h_sumJetPt->Integral()));
    h_sumJetPt->SetLineColor(809); //red
    leg->AddEntry(h_sumJetPt, "DY+2-4jets", "l");

    h_DYincl_jetPt->Draw("hist");
    h_sumJetPt->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Scaled_JetPt.pdf"); 

    //Scaled Jet Eta
    leg->Clear();
    h_DYincl_jetEta->Scale(1/(h_DYincl_jetEta->Integral()));
    h_DYincl_jetEta->SetLineColor(4); //blue
    h_DYincl_jetEta->GetYaxis()->SetTitle("Scaling");
    h_DYincl_jetEta->GetXaxis()->SetTitle("Jet #eta");
    leg->AddEntry(h_DYincl_jetEta, "Inclusive DY", "l");

    h_sumJetEta->Scale(1/(h_sumJetEta->Integral()));
    h_sumJetEta->SetLineColor(809); //red
    leg->AddEntry(h_sumJetEta, "DY+2-4jets", "l");

    h_DYincl_jetEta->Draw("hist");
    h_sumJetEta->Draw("histsame");
    leg->Draw();
    c1->Print("plot_JetsAndZkinematicsGE2jets/Scaled_JetEta.pdf");

} // end of main loop
