#include <iostream>
#include <vector>
#include <TTree.h>
#include <TFile.h>
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
void xAna_trk_withTrackHighPurity_test(TString inputfile = "Output_Preselection_Ztoee_Mchi2-150_Mchi1-1_ctau1p0_withTrackHighPurity_test.root", TString outputfile = "Output_TrackSelection_Mchi2-150_Mchi1-1_ctau1p0_withTrackHighPurity_test.root")
{

    TFile *file = TFile::Open(inputfile);

    TH1F *h_trackIP2D = new TH1F("h_trackIP2D", "", 100, -5, 5);
    h_trackIP2D->GetXaxis()->SetTitle("");
    h_trackIP2D->GetYaxis()->SetTitle("");
    h_trackIP2D->Sumw2();
    
    TH1F *h_trackIP2Dlog = new TH1F("h_trackIP2Dlog", "", 100, -5, 5);
    h_trackIP2Dlog->GetXaxis()->SetTitle("");
    h_trackIP2Dlog->GetYaxis()->SetTitle("");
    h_trackIP2Dlog->Sumw2();

    TH1F *h_mean_IP2D = new TH1F("h_mean_IP2D", "", 100, -5, 5);
    h_mean_IP2D->GetXaxis()->SetTitle("");
    h_mean_IP2D->GetYaxis()->SetTitle("");
    h_mean_IP2D->Sumw2();

    TH1F *h_absmedian_IP2D = new TH1F("h_absmedian_IP2D", "", 100, -5, 5);
    h_absmedian_IP2D->GetXaxis()->SetTitle("");
    h_absmedian_IP2D->GetYaxis()->SetTitle("");
    h_absmedian_IP2D->Sumw2();

    TH1F *h_track3Dsig = new TH1F("h_track3Dsig", "", 100, 0, 100);
    h_track3Dsig->GetXaxis()->SetTitle("");
    h_track3Dsig->GetYaxis()->SetTitle("");
    h_track3Dsig->Sumw2();

    TH1F *h_track3Dsiglog = new TH1F("h_track3Dsiglog", "", 100, -5, 5);
    h_track3Dsiglog->GetXaxis()->SetTitle("");
    h_track3Dsiglog->GetYaxis()->SetTitle("");
    h_track3Dsiglog->Sumw2();

    TH1F *h_trackChi3D = new TH1F("h_trackChi3D", "", 100, 0, 100);
    h_trackChi3D->GetXaxis()->SetTitle("");
    h_trackChi3D->GetYaxis()->SetTitle("");
    h_trackChi3D->Sumw2();

    TH1F *h_trackChi3Dlog = new TH1F("h_trackChi3Dlog", "", 100, -5, 5);
    h_trackChi3Dlog->GetXaxis()->SetTitle("");
    h_trackChi3Dlog->GetYaxis()->SetTitle("");
    h_trackChi3Dlog->Sumw2();

    //plot alpha value
    TH1F *h_passJetAlpha = new TH1F("h_passJetAlpha", "", 100, 0, 1);
    h_passJetAlpha->GetXaxis()->SetTitle("");
    h_passJetAlpha->GetYaxis()->SetTitle("");
    h_passJetAlpha->Sumw2();

    TH1F *h_passJetAlpha2 = new TH1F("h_passJetAlpha2", "", 100, 0, 1);
    h_passJetAlpha2->GetXaxis()->SetTitle("");
    h_passJetAlpha2->GetYaxis()->SetTitle("");
    h_passJetAlpha2->Sumw2();

    TH1F *h_passJetAlpha3 = new TH1F("h_passJetAlpha3", "", 100, 0, 1);
    h_passJetAlpha3->GetXaxis()->SetTitle("");
    h_passJetAlpha3->GetYaxis()->SetTitle("");
    h_passJetAlpha3->Sumw2();

    TH1F *h_passJetAlpha4 = new TH1F("h_passJetAlpha4", "", 100, 0, 1);
    h_passJetAlpha4->GetXaxis()->SetTitle("");
    h_passJetAlpha4->GetYaxis()->SetTitle("");
    h_passJetAlpha4->Sumw2();

    TH1F *h_passJetAlphaChi = new TH1F("h_passJetAlphaChi", "", 100, 0, 1);
    h_passJetAlphaChi->GetXaxis()->SetTitle("");
    h_passJetAlphaChi->GetYaxis()->SetTitle("");
    h_passJetAlphaChi->Sumw2();





    //-----------------------------------
    // Original Tree Variable
    //-----------------------------------
    Int_t I_event;
    Int_t I_weight;
    ULong64_t I_eventID;
    Float_t f_Met;
    //Float_t f_HT;
    //Float_t f_dileptonmass;
    Float_t f_dileptonPT;
    Float_t f_dileptonMass;
    Float_t f_ZbosonPt;
    Float_t f_ZbosonMass;
    Float_t f_ZbosonEta;
    Float_t f_ZbosonPhi;
    vector<float> *v_met_lepdeltaPhi = new vector<float>();
    Int_t I_nThinJets;
    //vector<float> *v_passJetindex = new vector<float>();
    vector<float> *v_passJetEta = new vector<float>();
    vector<float> *v_passJetPt = new vector<float>();
    vector<float> *v_passJetCSV = new vector<float>();
    vector<int> *v_passJetHadronFlavor = new vector<int>();
    //vector<int> *v_passjetPartonFlavor = new vector<int>();
    vector<int> *v_nTrack = new vector<int>();
    vector<float> *v_TrackPT = new vector<float>();
    vector<float> *v_TrackEta = new vector<float>();
    vector<float> *v_Trackdr = new vector<float>();
    //vector<float> *v_Trackindex = new vector<float>();
    vector<float> *v_Trackdz = new vector<float>();
    vector<float> *v_Trackdzerror = new vector<float>();
    vector<float> *v_Trackdxy = new vector<float>();
    vector<float> *v_Trackdxyerror = new vector<float>();
    vector<int> *v_TrackStatus = new vector<int>();
    vector<int> *v_TrackHighPurity = new vector<int>();

    v_met_lepdeltaPhi->clear();
    //v_passJetindex->clear();
    v_passJetEta->clear();
    v_passJetPt->clear();
    v_passJetCSV->clear();
    v_passJetHadronFlavor->clear();
    //v_passjetPartonFlavor->clear();
    v_nTrack->clear();
    v_TrackPT->clear();
    v_TrackEta->clear();
    v_Trackdr->clear();
    //v_Trackindex->clear();
    v_Trackdz->clear();
    v_Trackdzerror->clear();
    v_Trackdxy->clear();
    v_Trackdxyerror->clear();
    v_TrackStatus->clear();
    v_TrackHighPurity->clear();

    TTree *T_tree;
    file->GetObject("T_tree", T_tree);
    T_tree->SetBranchAddress("I_event", &I_event);
    T_tree->SetBranchAddress("I_weight", &I_weight);
    T_tree->SetBranchAddress("I_eventID", &I_eventID);
    T_tree->SetBranchAddress("f_Met", &f_Met);
    //T_tree->SetBranchAddress("f_HT", &f_HT);
    //T_tree->SetBranchAddress("f_dileptonmass", &f_dileptonmass);
    T_tree->SetBranchAddress("f_dileptonPT", &f_dileptonPT);
    T_tree->SetBranchAddress("f_dileptonMass", &f_dileptonMass);
    T_tree->SetBranchAddress("f_ZbosonPt", &f_ZbosonPt);
    T_tree->SetBranchAddress("f_ZbosonMass", &f_ZbosonMass);
    T_tree->SetBranchAddress("f_ZbosonEta", &f_ZbosonEta);
    T_tree->SetBranchAddress("f_ZbosonPhi", &f_ZbosonPhi);
    T_tree->SetBranchAddress("v_met_lepdeltaPhi", &v_met_lepdeltaPhi);
    T_tree->SetBranchAddress("I_nThinJets", &I_nThinJets);
    //T_tree->SetBranchAddress("v_passJetindex", &v_passJetindex);
    T_tree->SetBranchAddress("v_passJetEta", &v_passJetEta);
    T_tree->SetBranchAddress("v_passJetPt", &v_passJetPt);
    T_tree->SetBranchAddress("v_passJetCSV", &v_passJetCSV);
    T_tree->SetBranchAddress("v_passjethadronflavor", &v_passJetHadronFlavor);
    //T_tree->SetBranchAddress("v_passjetpartonflavor", &v_passjetPartonFlavor);
    T_tree->SetBranchAddress("v_nTrack", &v_nTrack);
    T_tree->SetBranchAddress("v_TrackPT", &v_TrackPT);
    T_tree->SetBranchAddress("v_TrackEta", &v_TrackEta);
    T_tree->SetBranchAddress("v_Trackdr", &v_Trackdr);
    //T_tree->SetBranchAddress("v_Trackindex", &v_Trackindex);
    T_tree->SetBranchAddress("v_Trackdz", &v_Trackdz);
    T_tree->SetBranchAddress("v_Trackdzerror", &v_Trackdzerror);
    T_tree->SetBranchAddress("v_Trackdxy", &v_Trackdxy);
    T_tree->SetBranchAddress("v_Trackdxyerror", &v_Trackdxyerror);
    T_tree->SetBranchAddress("v_TrackStatus", &v_TrackStatus);
    T_tree->SetBranchAddress("v_TrackHighPurity", &v_TrackHighPurity);

    //----------------------------
    // New Tree Variables
    //----------------------------
    Int_t I_event_new;
    Double_t D_weight_new; //Int_t I_weight_new;
    ULong64_t I_eventID_new;
    Float_t f_Met_new;
    //Float_t f_HT_new;
    //Float_t f_dileptonmass_new;
    Float_t f_dileptonPT_new;
    Float_t f_dileptonMass_new;
    Float_t f_ZbosonPt_new;
    Float_t f_ZbosonMass_new;
    Float_t f_ZbosonEta_new;
    Float_t f_ZbosonPhi_new;
    Int_t I_nThinJets_new;
    vector<double> v_met_lepdeltaPhi_new;
    //vector<float> v_passJetindex_new;
    vector<double> v_passJetEta_new;
    vector<double> v_passJetPt_new;
    vector<double> v_passJetCSV_new;
    vector<int> v_passJetHadronFlavor_new;
    //vector<int> v_passjetPartonFlavor_new;
    
    vector<int> v_nTrack_new;
    vector<double> v_TrackPT_new;
    vector<double> v_TrackEta_new;
    vector<double> v_Trackdr_new;
    //vector<float> v_Trackindex_new;
    //vector<float> v_Trackdz_new;
    //vector<float> v_Trackdzerror_new;
    //vector<float> v_Trackdxy_new;
    //vector<float> v_Trackdxyerror_new;
    vector<int> v_TrackStatus_new;
    vector<int> v_TrackHighPurity_new;

    vector<double> v_passjetalpha;
    vector<double> v_passjetalpha2;
    vector<double> v_passjetalpha3;
    vector<double> v_passjetalpha4;
    vector<double> v_passjetalphaChiTrack;
    vector<double> v_passjetmedian_IP2D;
    vector<double> v_passjetmean_IP2D;
    vector<double> v_track3Dsig;
    vector<double> v_trackChi3D;

    

    //---------------
    // Create New Tree
    //----------------
    TFile *output = TFile::Open(outputfile, "RECREATE");
    TTree *h1 = new TTree("h1", "NewTree");
    h1->Branch("I_event", &I_event_new);
    h1->Branch("D_weight", &D_weight_new); //h1->Branch("I_weight", &I_weight_new);
    h1->Branch("I_eventID", &I_eventID_new);
    h1->Branch("f_Met", &f_Met_new);
    //h1->Branch("f_HT", &f_HT_new);
    //h1->Branch("f_dileptonmass", &f_dileptonmass_new);
    h1->Branch("f_dileptonPT", &f_dileptonPT_new);
    h1->Branch("f_dileptonMass", &f_dileptonMass_new);
    h1->Branch("f_ZbosonPt", &f_ZbosonPt_new);
    h1->Branch("f_ZbosonMass", &f_ZbosonMass_new);
    h1->Branch("f_ZbosonEta", &f_ZbosonEta_new);
    h1->Branch("f_ZbosonPhi", &f_ZbosonPhi_new);
    h1->Branch("v_met_lepdeltaPhi", &v_met_lepdeltaPhi);
    h1->Branch("I_nThinJets", &I_nThinJets_new);
    //h1->Branch("v_passJetindex", &v_passJetindex_new);
    h1->Branch("v_passJetEta", &v_passJetEta_new);
    h1->Branch("v_passJetPt", &v_passJetPt_new);
    h1->Branch("v_passJetCSV", &v_passJetCSV_new);
    h1->Branch("v_passJetHadronFlavor", &v_passJetHadronFlavor_new);
    //h1->Branch("v_passjetpartonflavor", &v_passjetPartonFlavor_new);
    
    h1->Branch("v_nTrack", &v_nTrack_new);
    h1->Branch("v_TrackPT", &v_TrackPT_new);
    h1->Branch("v_TrackEta", &v_TrackEta_new);
    h1->Branch("v_Trackdr", &v_Trackdr_new);
    //h1->Branch("v_Trackindex", &v_Trackindex_new);
    //h1->Branch("v_Trackdz", &v_Trackdz_new);
    //h1->Branch("v_Trackdzerror", &v_Trackdzerror_new);
    //h1->Branch("v_Trackdxy", &v_Trackdxy_new);
    //h1->Branch("v_Trackdxyerror", &v_Trackdxyerror_new);
    h1->Branch("v_TrackStatus", &v_TrackStatus_new);
    h1->Branch("v_TrackHighPurity", &v_TrackHighPurity_new);

    h1->Branch("v_passjetalpha", &v_passjetalpha);
    h1->Branch("v_passjetalpha2", &v_passjetalpha2);
    h1->Branch("v_passjetalpha3", &v_passjetalpha3);
    h1->Branch("v_passjetalpha4", &v_passjetalpha4);
    h1->Branch("v_passjetalphaChiTrack", &v_passjetalphaChiTrack);
    h1->Branch("v_passjetmedian_IP2D", &v_passjetmedian_IP2D);
    h1->Branch("v_passjetmean_IP2D", &v_passjetmean_IP2D);
    h1->Branch("v_track3Dsig", &v_track3Dsig);
    h1->Branch("v_trackChi3D", &v_trackChi3D);

    int nTrackBefore = 0;
    int nTrackAfter = 0;
    int nJetBefore = 0;
    int nJetAfter= 0;
    int nTrkdzErrInf = 0;
    int nTrkdxyErrInf = 0;
    int nTrkdzErrZero = 0;
    int nTrkdxyErrZero = 0;
    int n3DsigInf = 0;
    int nNoSumPt = 0;

    for (int evt = 0; evt < T_tree->GetEntries(); evt++)
    {
        //v_met_lepdeltaPhi_new.clear();
        //v_passJetindex_new.clear();
        v_passJetEta_new.clear();
        v_passJetPt_new.clear();
        v_passJetCSV_new.clear();
        v_passJetHadronFlavor_new.clear();
        //v_passjetpartonflavor_new.clear();
        
        v_nTrack_new.clear();
        v_TrackPT_new.clear();
        v_TrackEta_new.clear();
        v_Trackdr_new.clear();
        //v_Trackindex_new.clear();
        //v_Trackdz_new.clear();
        //v_Trackdzerror_new.clear();
        //v_Trackdxy_new.clear();
        //v_Trackdxyerror_new.clear();
        v_TrackStatus_new.clear();
        v_TrackHighPurity_new.clear();

        v_passjetalpha.clear();
        v_passjetalpha2.clear();
        v_passjetalpha3.clear();
        v_passjetalpha4.clear();
        v_passjetalphaChiTrack.clear();
        v_passjetmedian_IP2D.clear();
        v_passjetmean_IP2D.clear();
        v_track3Dsig.clear();
        v_trackChi3D.clear();
        

        T_tree->GetEntry(evt);

        Int_t mcWeight = I_weight;
        Double_t eventWeight = mcWeight;

        // Track Selection
        int iTrack = 0;
        int nJet = 0;
        for (size_t iJet = 0; iJet < v_nTrack->size(); iJet++)
        {
            nJetBefore++;
            int nnTrack = 0;
            vector<float> v_trackIP2D;
            v_trackIP2D.clear();
            float SumTrackPT = 0;
            float SumTrackPTcut = 0;

            float SumTrackPTcut2 = 0;
            float SumTrackPTcut3 = 0;
            float SumTrackPTcut4 = 0;
            float SumTrackPTcutChiTrack = 0;

            for (int j = iTrack; j < iTrack + (*v_nTrack)[iJet]; j++)
            {
                nTrackBefore++;
                nnTrack++;
                float dz = (*v_Trackdz)[j];
                float dzerror = (*v_Trackdzerror)[j];
                float dxy = (*v_Trackdxy)[j];
                float dxyerror = (*v_Trackdxyerror)[j];
                if (isinf(dzerror) == 1)
                {
                    nTrkdzErrInf++;
                    continue;
                }
                if (isinf(dxyerror) == 1)
                {
                    nTrkdxyErrInf++;
                    continue;
                }
                if (dzerror == 0.0)
                    nTrkdzErrZero++;
                if (dxyerror == 0.0)
                    nTrkdxyErrZero++;
                
                if ((*v_TrackPT)[j] < 1)
                    continue;
                if ((*v_TrackHighPurity)[j] != 1)
                    continue;
                
                nTrackAfter++;
                
                double Sig3D = sqrt(pow(dz / dzerror, 2) + pow(dxy / dxyerror, 2));
                if (isinf(Sig3D))
                {
                    //cout << "bug" << "\n";
                    //cout << "3Dsig value is infinite" << "\n";
                    n3DsigInf++;
                }
                v_track3Dsig.push_back(Sig3D);
                h_track3Dsig->Fill(Sig3D, I_weight);
                h_track3Dsiglog->Fill(log10(Sig3D), I_weight);

                float IP2D = dxy; //float IP2D = (*v_Trackdxy)[j];
                v_trackIP2D.push_back(IP2D);
                h_trackIP2D->Fill(IP2D, I_weight);
                h_trackIP2Dlog->Fill(log10(abs(IP2D)), I_weight);
                
                double ChiTrack = sqrt(pow(dz / 0.01, 2) + pow(dxy / dxyerror, 2));
                v_trackChi3D.push_back(ChiTrack);
                h_trackChi3D->Fill(ChiTrack, I_weight);
                h_trackChi3Dlog->Fill(log10(ChiTrack), I_weight);
    
                //save track pt, eta, dr in track_new
                //v_nTrack_new.push_back();
                v_TrackPT_new.push_back((*v_TrackPT)[j]);
                v_TrackEta_new.push_back((*v_TrackEta)[j]);
                v_Trackdr_new.push_back((*v_Trackdr)[j]);

                //----------------------------------------------------
                // Choose the log10(Sig3D) value for alpha calculation
                //----------------------------------------------------
                SumTrackPT += (*v_TrackPT)[j];
                if (log10(Sig3D) < 1.0)
                {
                    SumTrackPTcut += (*v_TrackPT)[j];
                }

                //==========================
                //Studying the Track's alpha
                //==========================
                if (log10(Sig3D) < 2.0)
                {
                    SumTrackPTcut2 += (*v_TrackPT)[j];
                }
                if (log10(Sig3D) < 3.0)
                {
                    SumTrackPTcut3 += (*v_TrackPT)[j];
                }
                if (log10(Sig3D) < 4.0)
                {
                    SumTrackPTcut4 += (*v_TrackPT)[j];
                }
                if (log10(ChiTrack) < 1.0)
                {
                    SumTrackPTcutChiTrack += (*v_TrackPT)[j];
                }

            } // End of Track loop

            //save passjet eta, csv, pt [ijet] in passjet_new

            v_passJetEta_new.push_back((*v_passJetEta)[iJet]);
            v_passJetPt_new.push_back((*v_passJetPt)[iJet]);
            v_passJetCSV_new.push_back((*v_passJetCSV)[iJet]);
            v_passJetHadronFlavor_new.push_back((*v_passJetHadronFlavor)[iJet]);
            
            if (v_trackIP2D.size() == 0 && SumTrackPT != 0)
            {
                for (int j = iTrack; j < iTrack + (*v_nTrack)[iJet]; j++)
                {
                    cout << "Track PT = " << (*v_TrackPT)[j] << endl;
                    cout << "dzerror = " << (*v_Trackdzerror)[j] << endl;
                    cout << "dxyerror = " << (*v_Trackdxyerror)[j] << endl;
                    cout << "SumTrackPT = " << SumTrackPT << endl;
                }
            }
            iTrack += (*v_nTrack)[iJet];
            if (nnTrack != (*v_nTrack)[iJet])
            {
                cout << "bug" << endl;
            }
            //------------------------------------
            // Cut off jet if Jet no track pass cut
            //-------------------------------------
            if (SumTrackPT == 0)
            {
                continue;
                nNoSumPt++;
            }
            //----------------
            // Calculate alpha
            //----------------
            nJet++;
            nJetAfter++;
            double alpha = (SumTrackPTcut / SumTrackPT);

            //==========================
            //Studying alpha calculation
            //==========================
            double alpha2 = (SumTrackPTcut2 / SumTrackPT);
            double alpha3 = (SumTrackPTcut3 / SumTrackPT);
            double alpha4 = (SumTrackPTcut4 / SumTrackPT);
            double alphaChiTrack = (SumTrackPTcutChiTrack / SumTrackPT);

            v_passjetmean_IP2D.push_back(mean_value(v_trackIP2D));
            v_passjetmedian_IP2D.push_back(median_value(v_trackIP2D));
            h_mean_IP2D->Fill(log10(abs(mean_value(v_trackIP2D))), I_weight);
            h_absmedian_IP2D->Fill(log10(abs(median_value(v_trackIP2D))), I_weight);
            
            v_passjetalpha.push_back(alpha);
            h_passJetAlpha->Fill(alpha, I_weight);
            
            v_passjetalpha2.push_back(alpha2);
            h_passJetAlpha2->Fill(alpha2, I_weight);

            v_passjetalpha3.push_back(alpha3);
            h_passJetAlpha3->Fill(alpha3, I_weight);

            v_passjetalpha4.push_back(alpha4);
            h_passJetAlpha4->Fill(alpha4, I_weight);

            v_passjetalphaChiTrack.push_back(alphaChiTrack);
            h_passJetAlphaChi->Fill(alphaChiTrack, I_weight);

        } // End of Jet loop

        //v_met_lepdeltaPhi_new.push_back((*v_met_lepdeltaPhi)[evt]);

        // Fill all variable
        I_event_new = I_event;
        D_weight_new = eventWeight; //I_weight_new = I_weight;
        I_eventID_new = I_eventID;
        f_Met_new = f_Met;
        //f_HT_new = f_HT;
        //f_dileptonmass_new = f_dileptonmass;
        f_dileptonPT_new = f_dileptonPT;
        f_dileptonMass_new = f_dileptonMass;
        f_ZbosonPt_new = f_ZbosonPt;
        f_ZbosonMass_new = f_ZbosonMass;
        f_ZbosonEta_new = f_ZbosonEta;
        f_ZbosonPhi_new = f_ZbosonPhi;
        I_nThinJets_new = nJet;
        h1->Fill();
    }

    //Check Track Efficiency
    cout << "number of tracks before selection (with JetTrackHighPurity) = " << nTrackBefore << "\n";
    cout << "number of tracks after selection (with JetTrackHighPurity) = " <<  nTrackAfter << "\n";
    efferr(nTrackAfter, nTrackBefore);

    cout << "nTrkdzErrInf = " << nTrkdzErrInf << "\n";
    cout << "nTrkdxyErrInf = " << nTrkdxyErrInf << "\n";
    cout << "nTrkdzErrZero = " << nTrkdzErrZero << "\n";
    cout << "nTrkdxyErrZero = " << nTrkdxyErrZero << "\n";
    cout << "n3DsigInf = " << n3DsigInf << "\n";
    cout << "nNoSumPt = " << nNoSumPt << "\n";

    cout << "number of jets before track selection (with JetTrackHighPurity) = " << nJetBefore << "\n";
    cout << "number of jets after track selection (with JetTrackHighPurity) = " <<  nJetAfter << "\n";
    efferr(nJetAfter, nJetBefore);

    // out Tree branches
    //TFile *outFile = new TFile(outputfile.c_str(), "RECREATE");
    output->cd();
    h1->Write();
    output->mkdir("Histograms", "Histograms")->cd();
    h_trackIP2D->Write();
    h_trackIP2Dlog->Write();
    h_mean_IP2D->Write();
    h_absmedian_IP2D->Write();
    h_track3Dsig->Write();
    h_track3Dsiglog->Write();
    h_trackChi3D->Write();
    h_trackChi3Dlog->Write();
    h_passJetAlpha->Write();
    h_passJetAlpha2->Write();
    h_passJetAlpha3->Write();
    h_passJetAlpha4->Write();
    h_passJetAlphaChi->Write();
    output->cd("/");
    //output->cd();
    //h1->Write();
    output->Close();
}
