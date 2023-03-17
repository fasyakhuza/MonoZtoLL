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
#include <TGraph.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <math.h>
#include "Cross_section.h"
#include <string>
using namespace std;

void efferr(float nsig, float ntotal, float factor = 1)
{
    float eff = nsig / ntotal;
    float err = sqrt((1 - eff) * eff / ntotal);
    cout << "efficiency = " << eff * factor << " +- " << err * factor << endl;
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


void xPlot_background_eachType(string inputtxtFilename = "inputListBkg.txt", string outputtxtFilename = "allBkg_wDYincl_test.root")
{

    Double_t lumi = 41500.0; // integrated luminosity; unit: pb^{-1}
    Double_t XS_DYincl = 6077.22; // unit: pb


    // For Initial
    TH1D* h_total_mcweight = new TH1D("h_total_mcweight", "", 5, 0, 5);
    h_total_mcweight->Sumw2();

    TH1D* h_total_mcweight_new = new TH1D("h_total_mcweight_new", "", 5, 0, 5);
    h_total_mcweight_new->Sumw2();

    TH1D *h_HT_eventCount = new TH1D("h_HT_eventCount", "", 10, 0, 10);
    h_HT_eventCount->SetYTitle("N event");
    h_HT_eventCount->Sumw2();

    TH1D *h_HT_eventCount_new = new TH1D("h_HT_eventCount_new", "", 10, 0, 10);
    h_HT_eventCount_new->SetYTitle("N event");
    h_HT_eventCount_new->Sumw2();

    TH1D* h_DilepPt_afterRecoee = new TH1D("h_DilepPt_afterRecoee", "Dilepton pT after recoee", 1000, 0, 1000);
    h_DilepPt_afterRecoee->Sumw2();

    TH1D* h_DilepPt_afterNvtx = new TH1D("h_DilepPt_afterNvtx", "Dilepton pT after nVtx", 1000, 0, 1000);
    h_DilepPt_afterNvtx->Sumw2();

    TH1D* h_DilepPt_afterTauVeto = new TH1D("h_DilepPt_afterTauVeto", "Dilepton pT after Tau Veto", 1000, 0, 1000);
    h_DilepPt_afterTauVeto->Sumw2();

    TH1D* h_DilepPt_afterElePairPt = new TH1D("h_DilepPt_afterElePairPt", "Dilepton pT after Ele Pair Pt", 1000, 0, 1000);
    h_DilepPt_afterElePairPt->Sumw2();

    TH1D* h_DilepPt_afterZmass = new TH1D("h_DilepPt_afterZmass", "Dilepton pT after Z and Dilepton mass delta", 1000, 0, 1000);
    h_DilepPt_afterZmass->Sumw2();

    TH1D* h_DilepPt_afterExtraLepVeto = new TH1D("h_DilepPt_afterExtraLepVeto", "Dilepton pT after Extra Lepton Veto", 1000, 0, 1000);
    h_DilepPt_afterExtraLepVeto->Sumw2();

    TH1D* h_DilepPt_afterNthinjets = new TH1D("h_DilepPt_afterNthinjets", "Dilepton pT after nThinJets", 1000, 0, 1000);
    h_DilepPt_afterNthinjets->Sumw2();

    TH1D* h_DilepPt_afterDilepPt = new TH1D("h_DilepPt_afterDilepPt", "Dilepton pT after Dilepton Pt", 1000, 0, 1000);
    h_DilepPt_afterDilepPt->Sumw2();




    //histo arrays for DY incl for intial

    string suffix[9] = {"HT0to70", "HT70to100", "HT100to200", "HT200to400", "HT400to600", "HT600to800", "HT800to1200", "HT1200to2500", "HT2500toInf"};
    
    TH1D *h_DilepPt_afterRecoee_DYincl[9];
    char name_h_DilepPt_afterRecoee_DYincl[50];
    char title_h_DilepPt_afterRecoee_DYincl[50];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_DYincl, "h_DilepPt_afterRecoee_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_DYincl, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_DYincl[i] = new TH1D(name_h_DilepPt_afterRecoee_DYincl, title_h_DilepPt_afterRecoee_DYincl, 1000, 0, 1000);
        h_DilepPt_afterRecoee_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_DYincl[9];
    char name_h_DilepPt_afterNvtx_DYincl[50];
    char title_h_DilepPt_afterNvtx_DYincl[50];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_DYincl, "h_DilepPt_afterNvtx_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_DYincl, "Dilepton pT after nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_DYincl[i] = new TH1D(name_h_DilepPt_afterNvtx_DYincl, title_h_DilepPt_afterNvtx_DYincl, 1000, 0, 1000);
        h_DilepPt_afterNvtx_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_DYincl[9];
    char name_h_DilepPt_afterTauVeto_DYincl[50];
    char title_h_DilepPt_afterTauVeto_DYincl[50];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_DYincl, "h_DilepPt_afterTauVeto_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_DYincl, "Dilepton pT after Tau Veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_DYincl[i] = new TH1D(name_h_DilepPt_afterTauVeto_DYincl, title_h_DilepPt_afterTauVeto_DYincl, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_DYincl[9];
    char name_h_DilepPt_afterElePairPt_DYincl[50];
    char title_h_DilepPt_afterElePairPt_DYincl[50];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_DYincl, "h_DilepPt_afterElePairPt_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_DYincl, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_DYincl[i] = new TH1D(name_h_DilepPt_afterElePairPt_DYincl, title_h_DilepPt_afterElePairPt_DYincl, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_DYincl[i]->Sumw2();
    }




    //-------------------------------------//
    // NORMALIZED TO INTEGRATED LUMINOSITY //
    //-------------------------------------//
    //-------------------------------//
    //Histo array DY-incl XS weighted//
    //-------------------------------//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DYincl[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DYincl[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DYincl[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DYincl, "h_DilepPt_afterRecoee_xsWeighted_DYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DYincl, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DYincl[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DYincl, title_h_DilepPt_afterRecoee_xsWeighted_DYincl, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DYincl[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DYincl[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DYincl[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DYincl, "h_DilepPt_afterNvtx_xsWeighted_DYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DYincl, "Dilepton pT after nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DYincl[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DYincl, title_h_DilepPt_afterNvtx_xsWeighted_DYincl, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DYincl[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DYincl[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DYincl[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl, "h_DilepPt_afterTauVeto_xsWeighted_DYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DYincl, "Dilepton pT after Tau veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DYincl[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl, title_h_DilepPt_afterTauVeto_xsWeighted_DYincl, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DYincl[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DYincl[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DYincl[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DYincl[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl, "h_DilepPt_afterElePairPt_xsWeighted_DYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DYincl, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DYincl[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl, title_h_DilepPt_afterElePairPt_xsWeighted_DYincl, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DYincl[i]->Sumw2();
    }


    //For Each DY HT//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach, "h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach, title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach, "h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach, "Dilepton pT after Nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach, title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach, "h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach, "Dilepton pT after Tau Veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach, title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach, "h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach, title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[i]->Sumw2();
    }


    //For All DY incl and each HT plots//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT, "h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT, "Dilepton pT after recoee for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT, title_h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT, "h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT, "Dilepton pT after Nvtx for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT, title_h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT, "h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT, "Dilepton pT after Tau Veto for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT, title_h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT, "h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT, "Dilepton pT after Ele Pair Pt for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT, title_h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[i]->Sumw2();
    }


    //-----------------------------//
    //Histo array DY-HT XS weighted//
    //-----------------------------//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DY[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DY[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DY, "h_DilepPt_afterRecoee_xsWeighted_DY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DY, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DY[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DY, title_h_DilepPt_afterRecoee_xsWeighted_DY, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DY[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DY[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DY, "h_DilepPt_afterNvtx_xsWeighted_DY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DY, "Dilepton pT after Nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DY[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DY, title_h_DilepPt_afterNvtx_xsWeighted_DY, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DY[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DY[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DY, "h_DilepPt_afterTauVeto_xsWeighted_DY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DY, "Dilepton pT after Tau Veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DY[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DY, title_h_DilepPt_afterTauVeto_xsWeighted_DY, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DY[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DY[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DY, "h_DilepPt_afterElePairPt_xsWeighted_DY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DY, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DY[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DY, title_h_DilepPt_afterElePairPt_xsWeighted_DY, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DY[i]->Sumw2();
    }


    //For Each//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DY_forEach[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DY_forEach[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DY_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DY_forEach, "h_DilepPt_afterRecoee_xsWeighted_DY_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DY_forEach, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DY_forEach[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DY_forEach, title_h_DilepPt_afterRecoee_xsWeighted_DY_forEach, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DY_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DY_forEach[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DY_forEach[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DY_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DY_forEach, "h_DilepPt_afterNvtx_xsWeighted_DY_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DY_forEach, "Dilepton pT after Nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DY_forEach[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DY_forEach, title_h_DilepPt_afterNvtx_xsWeighted_DY_forEach, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DY_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach, "h_DilepPt_afterTauVeto_xsWeighted_DY_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach, "Dilepton pT after Tau Veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach, title_h_DilepPt_afterTauVeto_xsWeighted_DY_forEach, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach, "h_DilepPt_afterElePairPt_xsWeighted_DY_forEach_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach, title_h_DilepPt_afterElePairPt_xsWeighted_DY_forEach, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[i]->Sumw2();
    }


    //For All DY incl and each HT plots//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT, "h_DilepPt_afterRecoee_xsWeighted_DY_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT, "Dilepton pT after recoee for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT, title_h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT, "h_DilepPt_afterNvtx_xsWeighted_DY_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT, "Dilepton pT after Nvtx for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT, title_h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT, "h_DilepPt_afterTauVeto_xsWeighted_DY_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT, "Dilepton pT after Tau Veto for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT, title_h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT, "h_DilepPt_afterElePairPt_xsWeighted_DY_forDYincland%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT, "Dilepton pT after Ele Pair Pt for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT, title_h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[i]->Sumw2();
    }
    


    //--------------------------------------------------//
    // Sum After Normalization to Integrated Luminosity //
    //--------------------------------------------------//
    //for saving to output root file
    TH1D* h_sumDY_DilepPt_afterRecoee_xsWeighted = new TH1D("h_sumDY_DilepPt_afterRecoee_xsWeighted", "", 1000, 0, 1000);
    //TH1D* h_sum_DilepPt_afterRecoee_xsWeighted = new TH1D("h_sum_DilepPt_afterRecoee_xsWeighted", "Sum of Zpt distribution of all bkg before Zpt cut and after XS weight", 500, 0, 500);
    h_sumDY_DilepPt_afterRecoee_xsWeighted->Sumw2();

    TH1D* h_sumDY_DilepPt_afterNvtx_xsWeighted = new TH1D("h_sumDY_DilepPt_afterNvtx_xsWeighted", "", 1000, 0, 1000);
    h_sumDY_DilepPt_afterNvtx_xsWeighted->Sumw2();

    TH1D* h_sumDY_DilepPt_afterTauVeto_xsWeighted = new TH1D("h_sumDY_DilepPt_afterTauVeto_xsWeighted", "", 1000, 0, 1000);
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->Sumw2();

    TH1D* h_sumDY_DilepPt_afterElePairPt_xsWeighted = new TH1D("h_sumDY_DilepPt_afterElePairPt_xsWeighted", "", 1000, 0, 1000);
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->Sumw2();


    //for saving DY categorized by HT
    TH1D *h_DilepPt_afterRecoee_xsWeighted_sumDY[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_sumDY[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_sumDY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_sumDY, "h_DilepPt_afterRecoee_xsWeighted_sumDY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_sumDY, "Dilepton pT after recoee for %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_sumDY[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_sumDY, title_h_DilepPt_afterRecoee_xsWeighted_sumDY, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_sumDY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_sumDY[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_sumDY[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_sumDY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_sumDY, "h_DilepPt_afterNvtx_xsWeighted_sumDY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_sumDY, "Dilepton pT after Nvtx for %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_sumDY[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_sumDY, title_h_DilepPt_afterNvtx_xsWeighted_sumDY, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_sumDY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_sumDY[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_sumDY[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_sumDY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_sumDY, "h_DilepPt_afterTauVeto_xsWeighted_sumDY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_sumDY, "Dilepton pT after Tau Veto for %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_sumDY, title_h_DilepPt_afterTauVeto_xsWeighted_sumDY, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_sumDY[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_sumDY[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_sumDY[70];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_sumDY, "h_DilepPt_afterElePairPt_xsWeighted_sumDY_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_sumDY, "Dilepton pT after Ele Pair Pt for %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_sumDY, title_h_DilepPt_afterElePairPt_xsWeighted_sumDY, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[i]->Sumw2();
    }


    //For saving all DY and add each HT to it//
    TH1D *h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[9];
    char name_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[70];
    char title_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT, "h_DilepPt_afterRecoee_xsWeighted_sumDYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT, "Dilepton pT after recoee for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[i] = new TH1D(name_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT, title_h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT, 1000, 0, 1000);
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[9];
    char name_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[70];
    char title_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT, "h_DilepPt_afterNvtx_xsWeighted_sumDYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT, "Dilepton pT after Nvtx for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[i] = new TH1D(name_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT, title_h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT, 1000, 0, 1000);
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[9];
    char name_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[70];
    char title_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT, "h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT, "Dilepton pT after Tau Veto for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[i] = new TH1D(name_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT, title_h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT, 1000, 0, 1000);
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[i]->Sumw2();
    }

    TH1D *h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[9];
    char name_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[70];
    char title_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[100];
    for (int i = 0; i < 9; i++)
    {
        sprintf(name_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT, "h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_%s", suffix[i].c_str());
        sprintf(title_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT, "Dilepton pT after Ele Pair Pt for DY incl and %s", suffix[i].c_str());
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[i] = new TH1D(name_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT, title_h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT, 1000, 0, 1000);
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[i]->Sumw2();
    }



    
    
    cout << "inputtxtFilename = " << inputtxtFilename << endl;
    ifstream flist(inputtxtFilename.data());
    string inputFile;
    while (getline(flist, inputFile))
    {
        //h_total_mcweight->Reset();
        h_total_mcweight_new->Reset();
        h_HT_eventCount_new->Reset();
        if (inputFile.find("DYJetsToLL_M-50_TuneCP5") == string::npos)
        {
            h_DilepPt_afterRecoee->Reset();
            h_DilepPt_afterNvtx->Reset();
            h_DilepPt_afterTauVeto->Reset();
            h_DilepPt_afterElePairPt->Reset();
        }

        TString myFile = inputFile;
        cout << "\n" << myFile << endl;

        TFile* file = TFile::Open(myFile);
        //cout << "successfully opened" << endl;

        //TH1F* h_totevent = static_cast<TH1F*>(file->Get("Event_Variable/h_totevent"));
        h_total_mcweight = static_cast<TH1D*>(file->Get("h_total_mcweight"));
        h_total_mcweight_new = static_cast<TH1D*>(file->Get("h_total_mcweight_new"));
        if (inputFile.find("DYJetsToLL_M-50_TuneCP5") == string::npos)
        {
            h_DilepPt_afterRecoee = static_cast<TH1D*>(file->Get("h_DilepPt_afterRecoee"));
            h_DilepPt_afterNvtx = static_cast<TH1D*>(file->Get("h_DilepPt_afterNvtx"));
            h_DilepPt_afterTauVeto = static_cast<TH1D*>(file->Get("h_DilepPt_afterTauVeto"));
            h_DilepPt_afterElePairPt = static_cast<TH1D*>(file->Get("h_DilepPt_afterElePairPt"));
        }
            

        Double_t nTotalBeforePreselection = h_total_mcweight_new->Integral();

        Double_t nTotalBfrPreselection_HT0to70;
        Double_t nTotalBfrPreselection_HT70to100;
        Double_t nTotalBfrPreselection_HT100to200;
        Double_t nTotalBfrPreselection_HT200to400;
        Double_t nTotalBfrPreselection_HT400to600;
        Double_t nTotalBfrPreselection_HT600to800;
        Double_t nTotalBfrPreselection_HT800to1200;
        Double_t nTotalBfrPreselection_HT1200to2500;
        Double_t nTotalBfrPreselection_HT2500toInf;
        
        if (inputFile.find("DYJetsToLL_M-50_TuneCP5") != string::npos)
        {
            h_HT_eventCount = static_cast<TH1D*>(file->Get("h_HT_eventCount"));
            nTotalBfrPreselection_HT0to70 = h_HT_eventCount->GetBinContent(2);
            //cout << "nTotalBfrPreselection_HT0to70: " << nTotalBfrPreselection_HT0to70 << endl; 
            nTotalBfrPreselection_HT70to100 = h_HT_eventCount->GetBinContent(3);
            //cout << "nTotalBfrPreselection_HT70to100: " << nTotalBfrPreselection_HT70to100 << endl;
            nTotalBfrPreselection_HT100to200 = h_HT_eventCount->GetBinContent(4);
            nTotalBfrPreselection_HT200to400 = h_HT_eventCount->GetBinContent(5);
            nTotalBfrPreselection_HT400to600 = h_HT_eventCount->GetBinContent(6);
            nTotalBfrPreselection_HT600to800 = h_HT_eventCount->GetBinContent(7);
            nTotalBfrPreselection_HT800to1200 = h_HT_eventCount->GetBinContent(8);
            nTotalBfrPreselection_HT1200to2500 = h_HT_eventCount->GetBinContent(9);
            nTotalBfrPreselection_HT2500toInf = h_HT_eventCount->GetBinContent(10);


            char inputname_h_DilepPt_afterRecoee_DYincl[70];
            for (int i = 0; i < 9; i++)
            {
                sprintf(inputname_h_DilepPt_afterRecoee_DYincl, "h_DilepPt_afterRecoee_DYincl_%s", suffix[i].c_str());
                //cout << inputname_h_DilepPt_afterRecoee_DYincl << endl;
                h_DilepPt_afterRecoee_DYincl[i] = static_cast<TH1D*>(file->Get(inputname_h_DilepPt_afterRecoee_DYincl));

                cout << h_DilepPt_afterRecoee_DYincl[i]->Integral() << endl;
            }

            char inputname_h_DilepPt_afterNvtx_DYincl[70];
            for (int i = 0; i < 9; i++)
            {
                sprintf(inputname_h_DilepPt_afterNvtx_DYincl, "h_DilepPt_afterNvtx_DYincl_%s", suffix[i].c_str());
                h_DilepPt_afterNvtx_DYincl[i] = static_cast<TH1D*>(file->Get(inputname_h_DilepPt_afterNvtx_DYincl));

                cout << h_DilepPt_afterNvtx_DYincl[i]->Integral() << endl;
            }

            char inputname_h_DilepPt_afterTauVeto_DYincl[70];
            for (int i = 0; i < 9; i++)
            {
                sprintf(inputname_h_DilepPt_afterTauVeto_DYincl, "h_DilepPt_afterTauVeto_DYincl_%s", suffix[i].c_str());
                h_DilepPt_afterTauVeto_DYincl[i] = static_cast<TH1D*>(file->Get(inputname_h_DilepPt_afterTauVeto_DYincl));

                cout << h_DilepPt_afterTauVeto_DYincl[i]->Integral() << endl;
            }

            char inputname_h_DilepPt_afterElePairPt_DYincl[70];
            for (int i = 0; i < 9; i++)
            {
                sprintf(inputname_h_DilepPt_afterElePairPt_DYincl, "h_DilepPt_afterElePairPt_DYincl_%s", suffix[i].c_str());
                h_DilepPt_afterElePairPt_DYincl[i] = static_cast<TH1D*>(file->Get(inputname_h_DilepPt_afterElePairPt_DYincl));

                cout << h_DilepPt_afterElePairPt_DYincl[i]->Integral() << endl;
            }


            //For DY Inclusive HT0to70//
            Double_t HT0to70weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT0to70CS) / nTotalBfrPreselection_HT0to70;
            
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[0] = (TH1D*)h_DilepPt_afterRecoee_DYincl[0]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[0]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[0]->Scale(HT0to70weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0] = (TH1D*)h_DilepPt_afterRecoee_DYincl[0]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0]->Scale(HT0to70weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[0] = (TH1D*)h_DilepPt_afterRecoee_DYincl[0]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[0]->Scale(HT0to70weight);
            cout << "DY HT0to70: " << h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0]->Integral() << endl;

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[0] = (TH1D*)h_DilepPt_afterNvtx_DYincl[0]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[0]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[0]->Scale(HT0to70weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0] = (TH1D*)h_DilepPt_afterNvtx_DYincl[0]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0]->Scale(HT0to70weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[0] = (TH1D*)h_DilepPt_afterNvtx_DYincl[0]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[0]->Scale(HT0to70weight);
            cout << "DY HT0to70: " << h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0]->Integral() << endl;

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[0] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[0]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[0]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[0]->Scale(HT0to70weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[0]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0]->Scale(HT0to70weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[0] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[0]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[0]->Scale(HT0to70weight);
            cout << "DY HT0to70: " << h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0]->Integral() << endl;

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[0] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[0]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[0]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[0]->Scale(HT0to70weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[0]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0]->Scale(HT0to70weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[0] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[0]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[0]->Scale(HT0to70weight);
            cout << "DY HT0to70: " << h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0]->Integral() << endl;
        }

        if (inputFile.find("HT-70to100") != string::npos)
        {
            Double_t HT70to100weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT70to100CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT70to100);

            //For DY HT-70to100//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[1] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[1]);
            h_DilepPt_afterRecoee_xsWeighted_DY[1]->Scale(HT70to100weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[1] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[1]->Scale(HT70to100weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[1] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[1]);
            h_DilepPt_afterNvtx_xsWeighted_DY[1]->Scale(HT70to100weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[1] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[1]->Scale(HT70to100weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[1] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[1]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[1]->Scale(HT70to100weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[1] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[1]->Scale(HT70to100weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[1] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[1]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[1]->Scale(HT70to100weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[1] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[1]->Scale(HT70to100weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[1]->Scale(HT70to100weight);


            //For DY Inclusive HT70to100
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[1] = (TH1D*)h_DilepPt_afterRecoee_DYincl[1]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[1]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[1]->Scale(HT70to100weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[1] = (TH1D*)h_DilepPt_afterRecoee_DYincl[1]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[1]->Scale(HT70to100weight);
            cout << "DY HT70to100 from DY Inclusive: " << h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[1]->Integral() << endl;

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterRecoee_DYincl[1]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[1] = (TH1D*)h_DilepPt_afterNvtx_DYincl[1]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[1]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[1]->Scale(HT70to100weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[1] = (TH1D*)h_DilepPt_afterNvtx_DYincl[1]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[1]->Scale(HT70to100weight);
            cout << "DY HT70to100 from DY Inclusive: " << h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[1]->Integral() << endl;

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterNvtx_DYincl[1]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[1] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[1]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[1]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[1]->Scale(HT70to100weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[1] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[1]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[1]->Scale(HT70to100weight);
            cout << "DY HT70to100 from DY Inclusive: " << h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[1]->Integral() << endl;

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[1]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[1]->Scale(HT70to100weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[1] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[1]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[1]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[1]->Scale(HT70to100weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[1] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[1]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[1]->Scale(HT70to100weight);
            cout << "DY HT70to100 from DY Inclusive: " << h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[1]->Integral() << endl;

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[1] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[1]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[1]->Scale(HT70to100weight);
        }

        if (inputFile.find("HT-100to200") != string::npos)
        {
            Double_t HT100to200weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT100to200CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT100to200);
            
            //For DY HT-100to200//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[2] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[2]);
            h_DilepPt_afterRecoee_xsWeighted_DY[2]->Scale(HT100to200weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[2] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[2] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[2]);
            h_DilepPt_afterNvtx_xsWeighted_DY[2]->Scale(HT100to200weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[2] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[2] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[2]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[2]->Scale(HT100to200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[2] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[2] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[2]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[2]->Scale(HT100to200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[2] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[2]->Scale(HT100to200weight);


            //For DY Inclusive HT100to200
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[2] = (TH1D*)h_DilepPt_afterRecoee_DYincl[2]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[2]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[2]->Scale(HT100to200weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[2] = (TH1D*)h_DilepPt_afterRecoee_DYincl[2]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterRecoee_DYincl[2]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[2] = (TH1D*)h_DilepPt_afterNvtx_DYincl[2]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[2]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[2]->Scale(HT100to200weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[2] = (TH1D*)h_DilepPt_afterNvtx_DYincl[2]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterNvtx_DYincl[2]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[2] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[2]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[2]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[2]->Scale(HT100to200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[2] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[2]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[2]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[2]->Scale(HT100to200weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[2] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[2]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[2]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[2]->Scale(HT100to200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[2] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[2]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[2]->Scale(HT100to200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[2] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[2]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[2]->Scale(HT100to200weight);
        }

        if (inputFile.find("HT-200to400") != string::npos)
        {
            Double_t HT200to400weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT200to400CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT200to400);
            
            //For DY HT-200to400//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[3] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[3]);
            h_DilepPt_afterRecoee_xsWeighted_DY[3]->Scale(HT200to400weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[3] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[3] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[3]);
            h_DilepPt_afterNvtx_xsWeighted_DY[3]->Scale(HT200to400weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[3] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[3] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[3]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[3]->Scale(HT200to400weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[3] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[3] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[3]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[3]->Scale(HT200to400weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[3] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[3]->Scale(HT200to400weight);


            //For DY Inclusive HT200to400//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[3] = (TH1D*)h_DilepPt_afterRecoee_DYincl[3]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[3]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[3]->Scale(HT200to400weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[3] = (TH1D*)h_DilepPt_afterRecoee_DYincl[3]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterRecoee_DYincl[3]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[3] = (TH1D*)h_DilepPt_afterNvtx_DYincl[3]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[3]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[3]->Scale(HT200to400weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[3] = (TH1D*)h_DilepPt_afterNvtx_DYincl[3]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterNvtx_DYincl[3]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[3] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[3]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[3]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[3]->Scale(HT200to400weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[3] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[3]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[3]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[3]->Scale(HT200to400weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[3] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[3]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[3]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[3]->Scale(HT200to400weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[3] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[3]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[3]->Scale(HT200to400weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[3] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[3]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[3]->Scale(HT200to400weight);
        }

        if (inputFile.find("HT-400to600") != string::npos)
        {
            Double_t HT400to600weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT400to600CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT400to600);
            
            //For DY HT-400to600//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[4] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[4]);
            h_DilepPt_afterRecoee_xsWeighted_DY[4]->Scale(HT400to600weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[4] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[4] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[4]);
            h_DilepPt_afterNvtx_xsWeighted_DY[4]->Scale(HT400to600weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[4] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[4] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[4]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[4]->Scale(HT400to600weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[4] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[4] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[4]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[4]->Scale(HT400to600weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[4] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[4]->Scale(HT400to600weight);


            //For DY Inclusive HT400to600//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[4] = (TH1D*)h_DilepPt_afterRecoee_DYincl[4]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[4]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[4]->Scale(HT400to600weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[4] = (TH1D*)h_DilepPt_afterRecoee_DYincl[4]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterRecoee_DYincl[4]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[4] = (TH1D*)h_DilepPt_afterNvtx_DYincl[4]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[4]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[4]->Scale(HT400to600weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[4] = (TH1D*)h_DilepPt_afterNvtx_DYincl[4]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterNvtx_DYincl[4]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[4] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[4]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[4]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[4]->Scale(HT400to600weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[4] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[4]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[4]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[4]->Scale(HT400to600weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[4] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[4]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[4]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[4]->Scale(HT400to600weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[4] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[4]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[4]->Scale(HT400to600weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[4] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[4]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[4]->Scale(HT400to600weight);
        }

        if (inputFile.find("HT-600to800") != string::npos)
        {
            Double_t HT600to800weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT600to800CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT600to800);
            
            //For DY HT-600to800//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[5] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[5]);
            h_DilepPt_afterRecoee_xsWeighted_DY[5]->Scale(HT600to800weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[5] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[5] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[5]);
            h_DilepPt_afterNvtx_xsWeighted_DY[5]->Scale(HT600to800weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[5] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[5] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[5]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[5]->Scale(HT600to800weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[5] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[5] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[5]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[5]->Scale(HT600to800weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[5] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[5]->Scale(HT600to800weight);


            //For DY Inclusive HT600to800//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[5] = (TH1D*)h_DilepPt_afterRecoee_DYincl[5]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[5]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[5]->Scale(HT600to800weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[5] = (TH1D*)h_DilepPt_afterRecoee_DYincl[5]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterRecoee_DYincl[5]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[5] = (TH1D*)h_DilepPt_afterNvtx_DYincl[5]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[5]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[5]->Scale(HT600to800weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[5] = (TH1D*)h_DilepPt_afterNvtx_DYincl[5]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterNvtx_DYincl[5]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[5] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[5]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[5]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[5]->Scale(HT600to800weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[5] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[5]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[5]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[5]->Scale(HT600to800weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[5] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[5]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[5]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[5]->Scale(HT600to800weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[5] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[5]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[5]->Scale(HT600to800weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[5] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[5]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[5]->Scale(HT600to800weight);
        }
        
        if (inputFile.find("HT-800to1200") != string::npos)
        {
            Double_t HT800to1200weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT800to1200CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT800to1200);
            
            //For DY HT-800to1200//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[6] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[6]);
            h_DilepPt_afterRecoee_xsWeighted_DY[6]->Scale(HT800to1200weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[6] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[6] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[6]);
            h_DilepPt_afterNvtx_xsWeighted_DY[6]->Scale(HT800to1200weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[6] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[6] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[6]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[6]->Scale(HT800to1200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[6] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[6] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[6]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[6]->Scale(HT800to1200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[6] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[6]->Scale(HT800to1200weight);


            //For DY Inclusive HT800to1200//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[6] = (TH1D*)h_DilepPt_afterRecoee_DYincl[6]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[6]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[6]->Scale(HT800to1200weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[6] = (TH1D*)h_DilepPt_afterRecoee_DYincl[6]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterRecoee_DYincl[6]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[6] = (TH1D*)h_DilepPt_afterNvtx_DYincl[6]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[6]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[6]->Scale(HT800to1200weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[6] = (TH1D*)h_DilepPt_afterNvtx_DYincl[6]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterNvtx_DYincl[6]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[6] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[6]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[6]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[6]->Scale(HT800to1200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[6] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[6]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[6]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[6]->Scale(HT800to1200weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[6] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[6]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[6]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[6]->Scale(HT800to1200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[6] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[6]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[6]->Scale(HT800to1200weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[6] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[6]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[6]->Scale(HT800to1200weight);
        }

        if (inputFile.find("HT-1200to2500") != string::npos)
        {
            Double_t HT1200to2500weight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT1200to2500CS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT1200to2500);
            
            //For DY HT-1200to2500//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[7] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[7]);
            h_DilepPt_afterRecoee_xsWeighted_DY[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[7] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[7] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[7]);
            h_DilepPt_afterNvtx_xsWeighted_DY[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[7] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[7] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[7]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[7] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[7] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[7]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[7] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[7]->Scale(HT1200to2500weight);


            //For DY Inclusive HT1200to2500//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[7] = (TH1D*)h_DilepPt_afterRecoee_DYincl[7]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[7]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[7] = (TH1D*)h_DilepPt_afterRecoee_DYincl[7]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterRecoee_DYincl[7]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[7] = (TH1D*)h_DilepPt_afterNvtx_DYincl[7]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[7]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[7] = (TH1D*)h_DilepPt_afterNvtx_DYincl[7]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterNvtx_DYincl[7]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[7] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[7]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[7]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[7] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[7]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[7]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[7]->Scale(HT1200to2500weight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[7] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[7]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[7]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[7] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[7]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[7]->Scale(HT1200to2500weight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[7] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[7]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[7]->Scale(HT1200to2500weight);
        }

        if (inputFile.find("HT-2500toInf") != string::npos)
        {
            Double_t HT2500toInfweight = (GlobalConstants::Lumi2017) * 1000 * (GlobalConstants::HT2500toInfCS) / (nTotalBeforePreselection + nTotalBfrPreselection_HT2500toInf);
            
            //For DY HT-2500toInf//
            ////after Recoee
            h_DilepPt_afterRecoee_xsWeighted_DY[8] = (TH1D*)h_DilepPt_afterRecoee->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DY[8]);
            h_DilepPt_afterRecoee_xsWeighted_DY[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[8] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterRecoee->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            ////after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DY[8] = (TH1D*)h_DilepPt_afterNvtx->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DY[8]);
            h_DilepPt_afterNvtx_xsWeighted_DY[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[8] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterNvtx->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            ////after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DY[8] = (TH1D*)h_DilepPt_afterTauVeto->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DY[8]);
            h_DilepPt_afterTauVeto_xsWeighted_DY[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[8] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterTauVeto->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            ////after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DY[8] = (TH1D*)h_DilepPt_afterElePairPt->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DY[8]);
            h_DilepPt_afterElePairPt_xsWeighted_DY[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[8] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterElePairPt->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[8]->Scale(HT2500toInfweight);


            //For DY Inclusive HT2500toInf//
            //after Recoeee
            h_DilepPt_afterRecoee_xsWeighted_DYincl[8] = (TH1D*)h_DilepPt_afterRecoee_DYincl[8]->Clone();//(h_DilepPt_afterRecoee_xsWeighted_DYincl[8]);
            h_DilepPt_afterRecoee_xsWeighted_DYincl[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[8] = (TH1D*)h_DilepPt_afterRecoee_DYincl[8]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterRecoee_DYincl[8]->Clone();
            h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            //after Nvtx
            h_DilepPt_afterNvtx_xsWeighted_DYincl[8] = (TH1D*)h_DilepPt_afterNvtx_DYincl[8]->Clone();//(h_DilepPt_afterNvtx_xsWeighted_DYincl[8]);
            h_DilepPt_afterNvtx_xsWeighted_DYincl[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[8] = (TH1D*)h_DilepPt_afterNvtx_DYincl[8]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterNvtx_DYincl[8]->Clone();
            h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            //after Tau Veto
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[8] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[8]->Clone();//(h_DilepPt_afterTauVeto_xsWeighted_DYincl[8]);
            h_DilepPt_afterTauVeto_xsWeighted_DYincl[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[8] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[8]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterTauVeto_DYincl[8]->Clone();
            h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[8]->Scale(HT2500toInfweight);

            //after Ele Pair Pt
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[8] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[8]->Clone();//(h_DilepPt_afterElePairPt_xsWeighted_DYincl[8]);
            h_DilepPt_afterElePairPt_xsWeighted_DYincl[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[8] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[8]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[8]->Scale(HT2500toInfweight);

            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[8] = (TH1D*)h_DilepPt_afterElePairPt_DYincl[8]->Clone();
            h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[8]->Scale(HT2500toInfweight);
        }
    }// end of flist
    

    //Sum each type of background samples//
    for (int idyin = 0; idyin < 9; idyin++)
    {
        //after Recoee
        h_sumDY_DilepPt_afterRecoee_xsWeighted->Add(h_DilepPt_afterRecoee_xsWeighted_DYincl[idyin]);

        for (int idyex = 0; idyex < 9; idyex++)
        {
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->Add(h_DilepPt_afterRecoee_xsWeighted_DYincl_forDYinclandEachHT[idyex]);
        }
        cout << "\n(Recoee) DY incl only: " << h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->Integral() << endl;

        //after Nvtx
        h_sumDY_DilepPt_afterNvtx_xsWeighted->Add(h_DilepPt_afterNvtx_xsWeighted_DYincl[idyin]);

        for (int idyex = 0; idyex < 9; idyex++)
        {
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->Add(h_DilepPt_afterNvtx_xsWeighted_DYincl_forDYinclandEachHT[idyex]);
        }
        cout << "(Nvtx) DY incl only: " << h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->Integral() << endl;

        //after Tau Veto
        h_sumDY_DilepPt_afterTauVeto_xsWeighted->Add(h_DilepPt_afterTauVeto_xsWeighted_DYincl[idyin]);

        for (int idyex = 0; idyex < 9; idyex++)
        {
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->Add(h_DilepPt_afterTauVeto_xsWeighted_DYincl_forDYinclandEachHT[idyex]);
        }
        cout << "(Tau Veto) DY incl only: " << h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->Integral() << endl;

        //after Ele Pair Pt
        h_sumDY_DilepPt_afterElePairPt_xsWeighted->Add(h_DilepPt_afterElePairPt_xsWeighted_DYincl[idyin]);

        for (int idyex = 0; idyex < 9; idyex++)
        {
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->Add(h_DilepPt_afterElePairPt_xsWeighted_DYincl_forDYinclandEachHT[idyex]);
        }
        cout << "(Ele Pair Pt) DY incl only: " << h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->Integral() << endl;
    }

    for (int idyht = 1; idyht < 9; idyht++)
    {
        //after Recoee
        h_sumDY_DilepPt_afterRecoee_xsWeighted->Add(h_DilepPt_afterRecoee_xsWeighted_DY[idyht]);

        ////For each DY categorized by HT////
        h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[idyht]);
        h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterRecoee_xsWeighted_DY_forEach[idyht]);
        cout << "\n(Recoee) DY events from all samples with " << suffix[idyht] << ": " << h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->Integral() << endl;


        ////for all DY incl and each DY HT////
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyht]->Add(h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[idyht]);
        cout << suffix[idyht] << " sample only: " << h_DilepPt_afterRecoee_xsWeighted_DY_forDYinclandEachHT[idyht]->Integral() << endl;
        cout << "(Recoee) DY incl and DY " << suffix[idyht] << " samples: " << h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyht]->Integral() << endl;


        //after Nvtx
        h_sumDY_DilepPt_afterNvtx_xsWeighted->Add(h_DilepPt_afterNvtx_xsWeighted_DY[idyht]);

        ////For each DY categorized by HT////
        h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[idyht]);
        h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterNvtx_xsWeighted_DY_forEach[idyht]);
        cout << "\n(Nvtx) DY events from all samples with " << suffix[idyht] << ": " << h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->Integral() << endl;


        ////for all DY incl and each DY HT////
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyht]->Add(h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[idyht]);
        cout << suffix[idyht] << " sample only: " << h_DilepPt_afterNvtx_xsWeighted_DY_forDYinclandEachHT[idyht]->Integral() << endl;
        cout << "(Nvtx) DY incl and DY " << suffix[idyht] << " samples: " << h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyht]->Integral() << endl;


        //after Tau Veto
        h_sumDY_DilepPt_afterTauVeto_xsWeighted->Add(h_DilepPt_afterTauVeto_xsWeighted_DY[idyht]);

        ////For each DY categorized by HT////
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[idyht]);
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterTauVeto_xsWeighted_DY_forEach[idyht]);
        cout << "\n(Tau Veto) DY events from all samples with " << suffix[idyht] << ": " << h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->Integral() << endl;


        ////for all DY incl and each DY HT////
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyht]->Add(h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[idyht]);
        cout << suffix[idyht] << " sample only: " << h_DilepPt_afterTauVeto_xsWeighted_DY_forDYinclandEachHT[idyht]->Integral() << endl;
        cout << "(Tau Veto) DY incl and DY " << suffix[idyht] << " samples: " << h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyht]->Integral() << endl;


        //after Ele Pair Pt
        h_sumDY_DilepPt_afterElePairPt_xsWeighted->Add(h_DilepPt_afterElePairPt_xsWeighted_DY[idyht]);

        ////For each DY categorized by HT////
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[idyht]);
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->Add(h_DilepPt_afterElePairPt_xsWeighted_DY_forEach[idyht]);
        cout << "\n(Ele Pair Pt) DY events from all samples with " << suffix[idyht] << ": " << h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->Integral() << endl;


        ////for all DY incl and each DY HT////
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyht]->Add(h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[idyht]);
        cout << suffix[idyht] << " sample only: " << h_DilepPt_afterElePairPt_xsWeighted_DY_forDYinclandEachHT[idyht]->Integral() << endl;
        cout << "(Ele Pair Pt) DY incl and DY " << suffix[idyht] << " samples: " << h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyht]->Integral() << endl;
    }

    


    //Try to calculate number of events
    /*cout << "\nNumber of Events of ZpT distribution before ZpT cut: " << h_sum_DilepPt_afterRecoee_xsWeighted->Integral() << endl;
    cout << "Number of Events of ZpT distribution after ZpT cut > 110 GeV: " << h_sum_Zpt_afterZpTcut_xsWeighted->Integral() << endl;*/



    TString outputfile(outputtxtFilename);

    TFile *outFile = TFile::Open(outputfile, "RECREATE");
    outFile->cd();
    h_sumDY_DilepPt_afterRecoee_xsWeighted->Write();
    h_sumDY_DilepPt_afterNvtx_xsWeighted->Write();
    /*h_HT_eventCount->Write();*/
    outFile->Close();

    //check the integral of each background
    cout << "\n" << endl;
    cout << "\n(Recoee) Drell-Yan : " << h_sumDY_DilepPt_afterRecoee_xsWeighted->Integral() << endl;
    cout << "\n(Nvtx) Drell-Yan : " << h_sumDY_DilepPt_afterNvtx_xsWeighted->Integral() << endl;
    cout << "\n(Tau Veto) Drell-Yan : " << h_sumDY_DilepPt_afterTauVeto_xsWeighted->Integral() << endl;
    cout << "\n(Ele Pair Pt) Drell-Yan : " << h_sumDY_DilepPt_afterElePairPt_xsWeighted->Integral() << endl;

    //canvas
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1","c1",1200,700); //width-height

    //legend
    auto leg = new TLegend(0.44,0.65,0.75,0.85); //x1,y1,x2,y2
    leg->SetBorderSize(0);
    leg->SetTextSize(0.027);

    auto leg2 = new TLegend(0.44,0.45,0.75,0.85); //x1,y1,x2,y2
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.027);

    char legtitle[70];
    int dyht_color[8] = {880, 632, 800, 416, 922, 623, 829, 432}; //kviolet, kred, korange, kgreen, gray, pink, spring, cyan

    char legtitle2[70];
    int dyht_color2[9] = {860, 880, 632, 800, 416, 922, 623, 829, 432}; //kazure, kviolet, kred, korange, kgreen, gray, pink, spring, cyan

    //---------//
    //Drell-Yan//
    //---------//
    //------//
    //Recoee//
    //------//
    //h_sum_DilepPt_afterRecoee_xsWeighted->SetTitle("Sum of weighted events of Zpt distributions before ZpT cut");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg->AddEntry(h_sumDY_DilepPt_afterRecoee_xsWeighted, "All Drell Yan", "l");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->Draw("hist");
    leg->Draw();
    c1->Print("outputplot/afterRecoee_Histo_DY_DileptonPt.pdf");
    c1->Print("outputplot/afterRecoee_Histo_DY_DileptonPt.png");


    ///Overlay total Drell-Yan and each Drell-Yan category///
    leg2->Clear();
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterRecoee_xsWeighted, "All Drell Yan", "l");

    h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0]->SetLineColor(860); //azure
    h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0], "Drell-Yan HT0to70", "l");

    for (int idyht = 1; idyht < 9; idyht++)
    {
        sprintf(legtitle, "Drell-Yan %s", suffix[idyht].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->SetLineColor(dyht_color[idyht-1]);
        h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->SetLineWidth(2);
        leg2->AddEntry(h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht], legtitle, "l");
    }

    h_sumDY_DilepPt_afterRecoee_xsWeighted->Draw("hist");
    h_DilepPt_afterRecoee_xsWeighted_DYincl_forEach[0]->Draw("ehistsame");
    for (int idyht = 1; idyht < 9; idyht++)
    {
        h_DilepPt_afterRecoee_xsWeighted_sumDY[idyht]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterRecoee_Histo_eachDY_DileptonPt.pdf");
    c1->Print("outputplot/afterRecoee_Histo_eachDY_DileptonPt.png");


    ///Overlay all DY and all DY incl and each HT plots///
    leg2->Clear();
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterRecoee_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    //h_sumDY_DilepPt_afterRecoee_xsWeighted->SetMaximum(11000);
    h_sumDY_DilepPt_afterRecoee_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterRecoee_xsWeighted, "All Drell Yan", "l");
    
    h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[0]->SetLineColor(860);
    h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[0], "DY inclusive only", "l");

    for (int idyin = 1; idyin < 9; idyin++)
    {
        sprintf(legtitle2, "DY Inclusive and DY %s", suffix[idyin].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineColor(dyht_color2[idyin]);
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineWidth(2);
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->GetXaxis()->SetRangeUser(0,400);
        leg2->AddEntry(h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin], legtitle2, "l");
    }

    h_sumDY_DilepPt_afterRecoee_xsWeighted->Draw("hist");
    for (int idyin = 0; idyin < 9; idyin++)
    {
        h_DilepPt_afterRecoee_xsWeighted_sumDYincl_and_eachHT[idyin]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterRecoee_Histo_DYinclAndEachHT_DileptonPt.pdf");
    c1->Print("outputplot/afterRecoee_Histo_DYinclAndEachHT_DileptonPt.png");



    //-----//
    //Nvtx//
    //-----//
    //h_sum_DilepPt_afterNvtx_xsWeighted->SetTitle("Sum of weighted events of Zpt distributions before ZpT cut");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg->AddEntry(h_sumDY_DilepPt_afterNvtx_xsWeighted, "All Drell Yan", "l");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->Draw("hist");
    leg->Draw();
    c1->Print("outputplot/afterNvtx_Histo_DY_DileptonPt.pdf");
    c1->Print("outputplot/afterNvtx_Histo_DY_DileptonPt.png");

    ///Overlay total Drell-Yan and each Drell-Yan category///
    leg2->Clear();
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterNvtx_xsWeighted, "All Drell Yan", "l");

    h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0]->SetLineColor(860); //azure
    h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0], "Drell-Yan HT0to70", "l");

    for (int idyht = 1; idyht < 9; idyht++)
    {
        sprintf(legtitle, "Drell-Yan %s", suffix[idyht].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->SetLineColor(dyht_color[idyht-1]);
        h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->SetLineWidth(2);
        leg2->AddEntry(h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht], legtitle, "l");
    }

    h_sumDY_DilepPt_afterNvtx_xsWeighted->Draw("hist");
    h_DilepPt_afterNvtx_xsWeighted_DYincl_forEach[0]->Draw("ehistsame");
    for (int idyht = 1; idyht < 9; idyht++)
    {
        h_DilepPt_afterNvtx_xsWeighted_sumDY[idyht]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterNvtx_Histo_eachDY_DileptonPt.pdf");
    c1->Print("outputplot/afterNvtx_Histo_eachDY_DileptonPt.png");

    ///Overlay all DY and all DY incl and each HT plots///
    leg2->Clear();
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterNvtx_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    //h_sumDY_DilepPt_afterNvtx_xsWeighted->SetMaximum(11000);
    h_sumDY_DilepPt_afterNvtx_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterNvtx_xsWeighted, "All Drell Yan", "l");
    
    h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[0]->SetLineColor(860);
    h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[0], "DY inclusive only", "l");

    for (int idyin = 1; idyin < 9; idyin++)
    {
        sprintf(legtitle2, "DY Inclusive and DY %s", suffix[idyin].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineColor(dyht_color2[idyin]);
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineWidth(2);
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->GetXaxis()->SetRangeUser(0,400);
        leg2->AddEntry(h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin], legtitle2, "l");
    }

    h_sumDY_DilepPt_afterNvtx_xsWeighted->Draw("hist");
    for (int idyin = 0; idyin < 9; idyin++)
    {
        h_DilepPt_afterNvtx_xsWeighted_sumDYincl_and_eachHT[idyin]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterNvtx_Histo_DYinclAndEachHT_DileptonPt.pdf");
    c1->Print("outputplot/afterNvtx_Histo_DYinclAndEachHT_DileptonPt.png");



    //--------//
    //Tau Veto//
    //--------//
    //h_sum_DilepPt_afterTauVeto_xsWeighted->SetTitle("Sum of weighted events of Zpt distributions before ZpT cut");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg->AddEntry(h_sumDY_DilepPt_afterTauVeto_xsWeighted, "All Drell Yan", "l");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->Draw("hist");
    leg->Draw();
    c1->Print("outputplot/afterTauVeto_Histo_DY_DileptonPt.pdf");
    c1->Print("outputplot/afterTauVeto_Histo_DY_DileptonPt.png");

    ///Overlay total Drell-Yan and each Drell-Yan category///
    leg2->Clear();
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterTauVeto_xsWeighted, "All Drell Yan", "l");

    h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0]->SetLineColor(860); //azure
    h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0], "Drell-Yan HT0to70", "l");

    for (int idyht = 1; idyht < 9; idyht++)
    {
        sprintf(legtitle, "Drell-Yan %s", suffix[idyht].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->SetLineColor(dyht_color[idyht-1]);
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->SetLineWidth(2);
        leg2->AddEntry(h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht], legtitle, "l");
    }

    h_sumDY_DilepPt_afterTauVeto_xsWeighted->Draw("hist");
    h_DilepPt_afterTauVeto_xsWeighted_DYincl_forEach[0]->Draw("ehistsame");
    for (int idyht = 1; idyht < 9; idyht++)
    {
        h_DilepPt_afterTauVeto_xsWeighted_sumDY[idyht]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterTauVeto_Histo_eachDY_DileptonPt.pdf");
    c1->Print("outputplot/afterTauVeto_Histo_eachDY_DileptonPt.png");

    ///Overlay all DY and all DY incl and each HT plots///
    leg2->Clear();
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    //h_sumDY_DilepPt_afterTauVeto_xsWeighted->SetMaximum(11000);
    h_sumDY_DilepPt_afterTauVeto_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterTauVeto_xsWeighted, "All Drell Yan", "l");
    
    h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[0]->SetLineColor(860);
    h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[0], "DY inclusive only", "l");

    for (int idyin = 1; idyin < 9; idyin++)
    {
        sprintf(legtitle2, "DY Inclusive and DY %s", suffix[idyin].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineColor(dyht_color2[idyin]);
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineWidth(2);
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->GetXaxis()->SetRangeUser(0,400);
        leg2->AddEntry(h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin], legtitle2, "l");
    }

    h_sumDY_DilepPt_afterTauVeto_xsWeighted->Draw("hist");
    for (int idyin = 0; idyin < 9; idyin++)
    {
        h_DilepPt_afterTauVeto_xsWeighted_sumDYincl_and_eachHT[idyin]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterTauVeto_Histo_DYinclAndEachHT_DileptonPt.pdf");
    c1->Print("outputplot/afterTauVeto_Histo_DYinclAndEachHT_DileptonPt.png");



    //-----------//
    //Ele Pair Pt//
    //-----------//
    //h_sum_DilepPt_afterElePairPt_xsWeighted->SetTitle("Sum of weighted events of Zpt distributions before ZpT cut");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg->AddEntry(h_sumDY_DilepPt_afterElePairPt_xsWeighted, "All Drell Yan", "l");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->Draw("hist");
    leg->Draw();
    c1->Print("outputplot/afterElePairPt_Histo_DY_DileptonPt.pdf");
    c1->Print("outputplot/afterElePairPt_Histo_DY_DileptonPt.png");

    ///Overlay total Drell-Yan and each Drell-Yan category///
    leg2->Clear();
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterElePairPt_xsWeighted, "All Drell Yan", "l");

    h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0]->SetLineColor(860); //azure
    h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0], "Drell-Yan HT0to70", "l");

    for (int idyht = 1; idyht < 9; idyht++)
    {
        sprintf(legtitle, "Drell-Yan %s", suffix[idyht].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->SetLineColor(dyht_color[idyht-1]);
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->SetLineWidth(2);
        leg2->AddEntry(h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht], legtitle, "l");
    }

    h_sumDY_DilepPt_afterElePairPt_xsWeighted->Draw("hist");
    h_DilepPt_afterElePairPt_xsWeighted_DYincl_forEach[0]->Draw("ehistsame");
    for (int idyht = 1; idyht < 9; idyht++)
    {
        h_DilepPt_afterElePairPt_xsWeighted_sumDY[idyht]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterElePairPt_Histo_eachDY_DileptonPt.pdf");
    c1->Print("outputplot/afterElePairPt_Histo_eachDY_DileptonPt.png");

    ///Overlay all DY and all DY incl and each HT plots///
    leg2->Clear();
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetTitle("");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineColor(804); //brown
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetLineWidth(2);
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetYaxis()->SetTitle("Number of Events");
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetTitle("p_{T,ll} (GeV)");
    //h_sumDY_DilepPt_afterElePairPt_xsWeighted->SetMaximum(11000);
    h_sumDY_DilepPt_afterElePairPt_xsWeighted->GetXaxis()->SetRangeUser(0,400);
    leg2->AddEntry(h_sumDY_DilepPt_afterElePairPt_xsWeighted, "All Drell Yan", "l");
    
    h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[0]->SetLineColor(860);
    h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[0]->SetLineWidth(2);
    leg2->AddEntry(h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[0], "DY inclusive only", "l");

    for (int idyin = 1; idyin < 9; idyin++)
    {
        sprintf(legtitle2, "DY Inclusive and DY %s", suffix[idyin].c_str());

        //For each DY categorized by HT//
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineColor(dyht_color2[idyin]);
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->SetLineWidth(2);
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->GetXaxis()->SetRangeUser(0,400);
        leg2->AddEntry(h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin], legtitle2, "l");
    }

    h_sumDY_DilepPt_afterElePairPt_xsWeighted->Draw("hist");
    for (int idyin = 0; idyin < 9; idyin++)
    {
        h_DilepPt_afterElePairPt_xsWeighted_sumDYincl_and_eachHT[idyin]->Draw("ehistsame");
    }
    leg2->Draw();
    c1->Print("outputplot/afterElePairPt_Histo_DYinclAndEachHT_DileptonPt.pdf");
    c1->Print("outputplot/afterElePairPt_Histo_DYinclAndEachHT_DileptonPt.png");
} // end of main loop
