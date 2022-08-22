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

void xAna_ztoee_genlevelonly_test_2(string inputfile = "ExoPieElementTuples_MonoZLL_NLO_Mphi-500_Mchi2-150_Mchi1-1_gSM-0p25_gDM-1p0_ctau1p0.root", string outputfile = "Kinematics_Mchi2-150_Mchi1-1_ctau1p0_test.root")
{
    //get TTree from file ...
    TreeReader data(inputfile.data());

    Long64_t nTotal = 0;

    TH1F* h_genParPx_e = new TH1F("h_genParPx_e", "Px of Eles", 100, 0, 1000);
    TH1F* h_genParPx_mu = new TH1F("h_genParPx_mu", "Px of Muons", 100, 0, 1000);
    TH1F* h_genParPx_Z = new TH1F("h_genParPx_Z", "Px of Z bosons", 100, 0, 1000);

    TH1F* h_genParPy_e = new TH1F("h_genParPy_e", "Py of Eles", 100, 0, 1000);
    TH1F* h_genParPy_mu = new TH1F("h_genParPy_mu", "Py of Muons", 100, 0, 1000);
    TH1F* h_genParPy_Z = new TH1F("h_genParPy_Z", "Py of Z bosons", 100, 0, 1000);

    TH1F* h_genParPz_e = new TH1F("h_genParPz_e", "Pz of Eles", 100, 0, 1000);
    TH1F* h_genParPz_mu = new TH1F("h_genParPz_mu", "Pz of Muons", 100, 0, 1000);
    TH1F* h_genParPz_Z = new TH1F("h_genParPz_Z", "Pz of Z bosons", 100, 0, 1000);

    TH1F* h_genParE_e = new TH1F("h_genParE_e", "Energy of Eles", 100, 0, 1000);
    TH1F* h_genParE_mu = new TH1F("h_genParE_mu", "Energy of Muons", 100, 0, 1000);
    TH1F* h_genParE_Z = new TH1F("h_genParE_Z", "Energy of Z bosons", 100, 0, 1000);

    TH1F* h_pt_e = new TH1F("h_pt_e", "Pt of Eles", 100, 0, 1000);
    TH1F* h_pt_mu = new TH1F("h_pt_mu", "Pt of Muons", 100, 0, 1000);
    TH1F* h_pt_Z = new TH1F("h_pt_Z", "Pt of Z bosons", 100, 0, 1000);

    TH1F* h_phi_e = new TH1F("h_phi_e", "Phi of Eles", 100, -5, 5);
    TH1F* h_phi_mu = new TH1F("h_phi_mu", "Phi of Muons", 100, -5, 5);
    TH1F* h_phi_Z = new TH1F("h_phi_Z", "Phi of Z bosons", 100, -5, 5);

    TH1F* h_eta_e = new TH1F("h_eta_e", "Eta of Eles", 100, -10, 10);
    TH1F* h_eta_mu = new TH1F("h_eta_mu", "Eta of Muons", 100, -10, 10);
    TH1F* h_eta_Z = new TH1F("h_eta_Z", "Eta of Z bosons", 100, -10, 10);

    TH1F* h_mass_e = new TH1F("h_mass_e", "Mass of Eles", 100, 0, 1000);
    TH1F* h_mass_mu = new TH1F("h_mass_mu", "Mass of Muons", 100, 0, 1000);
    TH1F* h_mass_Z = new TH1F("h_mass_Z", "Mass of Z bosons", 100, 0, 200);

    //Event Loop
    for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++)
    {
        if (jEntry % 1000 == 0)
        {
            fprintf(stderr, "Processing event %lli of %lli\n", jEntry+1, data.GetEntriesFast());
        }

        data.GetEntry(jEntry);
        nTotal ++;

        // 0. check the generator-level information
        Int_t nGenPar = data.GetInt("nGenPar");
        Int_t* genParId = data.GetPtrInt("genParId");
        Int_t* genMomParId = data.GetPtrInt("genMomParId");

        Float_t* genParPx = data.GetPtrFloat("genParPx");
        Float_t* genParPy = data.GetPtrFloat("genParPy");
        Float_t* genParPz = data.GetPtrFloat("genParPz");
        Float_t* genParE = data.GetPtrFloat("genParE");
        
        // store the kinematics of generator particles
        for (int ig=0; ig < nGenPar; ig++)
        {
            int pid = genParId[ig];
            int mompid = genMomParId[ig];
            //h_genParPx->Fill(genParPx[ig]);
            //h_genParPy->Fill(genParPy[ig]);
            //h_genParPz->Fill(genParPz[ig]);
            //h_genParE->Fill(genParE[ig]);

            //kinematics that independent to z axis
            TLorentzVector* genPar_p4 = new TLorentzVector(genParPx[ig], genParPy[ig], genParPz[ig], genParE[ig]);
            //float pt = sqrt(genParPx[ig] * genParPx[ig] + genParPy[ig] * genParPy[ig]);
            //cout << "pt_ig = " << pt << "\n";
            //h_pt->Fill(pt);

            //float phi = atan(genParPy[ig] / genParPx[ig]);
            //h_phi->Fill(phi);

            //float p = sqrt(genParPx[ig] * genParPx[ig] + genParPy[ig] * genParPy[ig] + genParPz[ig] * genParPz[ig]);
            //float theta = asin(pt / p);

            //float eta = -log(tan(theta /2));
            //float eta = 0.5 * log((p + genParPz[ig]) / (p - genParPz[ig]));
            //h_eta->Fill(eta);

            //float mass = sqrt(genParE[ig] * genParE[ig] - p * p);
            //h_mass->Fill(mass);


            // get eles coming from Z
            if (abs(pid) == 11 && mompid == 23)
                {
                    h_genParPx_e->Fill(genParPx[ig]);
                    h_genParPy_e->Fill(genParPy[ig]);
                    h_genParPz_e->Fill(genParPz[ig]);
                    h_genParE_e->Fill(genParE[ig]);

                    h_pt_e->Fill(genPar_p4->Pt());
                    h_phi_e->Fill(genPar_p4->Phi());
                    h_eta_e->Fill(genPar_p4->Eta());
                    h_mass_e->Fill(genPar_p4->M());
                }

            //get muons from Z
            if (abs(pid) == 13 && mompid == 23)
                {
                    h_genParPx_mu->Fill(genParPx[ig]);
                    h_genParPy_mu->Fill(genParPy[ig]);
                    h_genParPz_mu->Fill(genParPz[ig]);
                    h_genParE_mu->Fill(genParE[ig]);

                    h_pt_mu->Fill(genPar_p4->Pt());
                    h_phi_mu->Fill(genPar_p4->Phi());
                    h_eta_mu->Fill(genPar_p4->Eta());
                    h_mass_mu->Fill(genPar_p4->M());
                }
            
            //get Z boson
            if (abs(pid) == 23)
                {
                    h_genParPx_Z->Fill(genParPx[ig]);
                    h_genParPy_Z->Fill(genParPy[ig]);
                    h_genParPz_Z->Fill(genParPz[ig]);
                    h_genParE_Z->Fill(genParE[ig]);

                    h_pt_Z->Fill(genPar_p4->Pt());
                    h_phi_Z->Fill(genPar_p4->Phi());
                    h_eta_Z->Fill(genPar_p4->Eta());
                    h_mass_Z->Fill(genPar_p4->M());
                }

        } // end Generatot Particle loop
    } // end Event Loop

    TFile* outFile = new TFile(outputfile.data(),"recreate");
    h_genParPx_e->Write();
    h_genParPy_e->Write();
    h_genParPz_e->Write();
    h_genParE_e->Write();

    h_pt_e->Write();
    h_phi_e->Write();
    h_eta_e->Write();
    h_mass_e->Write();

    h_genParPx_mu->Write();
    h_genParPy_mu->Write();
    h_genParPz_mu->Write();
    h_genParE_mu->Write();

    h_pt_mu->Write();
    h_phi_mu->Write();
    h_eta_mu->Write();
    h_mass_mu->Write();

    h_genParPx_Z->Write();
    h_genParPy_Z->Write();
    h_genParPz_Z->Write();
    h_genParE_Z->Write();

    h_pt_Z->Write();
    h_phi_Z->Write();
    h_eta_Z->Write();
    h_mass_Z->Write();

    outFile->Close();
    cout << "nTotal = " << nTotal << endl;

} // big end