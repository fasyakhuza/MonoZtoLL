// example code to check the proper decay length of the chi2 in the long-lived
// extension of monoZ

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;
void xAna_monoZ_ctau(std::string inputFile="ExoPieElementTuples.root",std::string outputFile="histo_test.root",bool debug=false){

  //get TTree from file ...
  TreeReader data(inputFile.data());
  
  Long64_t nTotal=0;
  //Long64_t nPass[20]={0};

  //  TH1F* h_ctau=new TH1F("h_ctau","",10000,0,100);
  //  TH1F* h_ctau=new TH1F("h_ctau","",10000,0,10);
  TH1F* h_ctau1=new TH1F("h_ctau1","",100,0,5);
  TH1F* h_ctau2=new TH1F("h_ctau2","",100,0,2);
  TH1F* h_ctau_lab  = (TH1F*)h_ctau1->Clone("h_ctau_lab");
  TH1F* h_ctau_proper  = (TH1F*)h_ctau2->Clone("h_ctau_proper");

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    std::cout << "\nEvent-" << jEntry+1 << "\n";

    data.GetEntry(jEntry);
    nTotal ++;

    // 0. check the generator-level information and make sure there is a Z->e+e-
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genMomParId      = data.GetPtrInt("genMomParId");
    //Int_t* genMo1      = data.GetPtrInt("genMo1");
    
    Float_t* genParPx  = data.GetPtrFloat("genParPx");
    Float_t* genParPy  = data.GetPtrFloat("genParPy");
    Float_t* genParPz  = data.GetPtrFloat("genParPz");
    Float_t* genParE  = data.GetPtrFloat("genParE");

    //    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
    //TClonesArray* genParVtx = (TClonesArray*) data.GetPtrTObject("genParVtx");
    Float_t* genParVtxX = data.GetPtrFloat("genParVtxX");
    Float_t* genParVtxY = data.GetPtrFloat("genParVtxY");
    Float_t* genParVtxZ = data.GetPtrFloat("genParVtxZ");


    // 1. find chi2 and chi1

    for(int ig=0; ig < nGenPar; ig++){

      int momIndex;
      if(abs(genParId[ig])==18) //chi2 = 18
        {
          std::cout << "ig : " << ig << "\n";
          momIndex = ig;
          std::cout << "momIndex : " << momIndex << "\n";
        }

      TLorentzVector* chi2_p4 = new TLorentzVector(
	  					   genParPx[momIndex], 
	 					   genParPy[momIndex],
	 					   genParPz[momIndex],
	 					   genParE[momIndex]
	 		   );

      if(genParId[ig]!=5000522)continue; // chi1 added by Pythia8
      //int mo1 = genMo1[ig];
      //if(mo1<0)continue; //what is mother index in genMo1? why if genMo1 is negative, we ignore it?

      int momId=genMomParId[ig];
      if(abs(momId)!=18)continue; //just for double check

      float chi1_VtxX = genParVtxX[ig];
      float chi1_VtxY = genParVtxY[ig];
      float chi1_VtxZ = genParVtxZ[ig];

      float chi2_VtxX = genParVtxX[momIndex];
      float chi2_VtxY = genParVtxY[momIndex];
      float chi2_VtxZ = genParVtxZ[momIndex];

      float distX = chi1_VtxX - chi2_VtxX;
      float distY = chi1_VtxY - chi2_VtxY;
      float distZ = chi1_VtxZ - chi2_VtxZ;


      TVector3 v_chi1_vtx(genParVtxX[ig], genParVtxY[ig], genParVtxZ[ig]);

      TVector3 v_chi2_vtx(genParVtxX[momIndex], genParVtxY[momIndex], genParVtxZ[momIndex]);

      TVector3 v_dist =  v_chi1_vtx - v_chi2_vtx;


      double ctau_lab = sqrt(distX * distX + distY * distY + distZ * distZ);
      std::cout << "ctau_lab = " << ctau_lab << "\n";
      
      double v_ctau_lab = v_dist.Mag();
      std::cout << "v_ctau_lab = " << v_ctau_lab << "\n";

      //      TLorentzVector* chi1_p4 = (TLorentzVector*)genParP4->At(ig);      
      
      
      // note, the position saved in CMSSW has a unit of cm
      // while in Pythia8 the setting has a unit of mm
      
      //TVector3* chi1_vtx = (TVector3*)genParVtx->At(ig);
      //TVector3* chi2_vtx = (TVector3*)genParVtx->At(mo1);


      double ctau_proper = ctau_lab/chi2_p4->Beta()/chi2_p4->Gamma();
      std::cout << "ctau_proper = " << ctau_proper << "\n";

      double v_ctau_proper = v_ctau_lab/chi2_p4->Beta()/chi2_p4->Gamma();
      std::cout << "v_ctau_proper = " << v_ctau_proper << "\n";

      if(debug){
	cout << "ctau_lab = " << ctau_lab;
	cout << "\t ctau_proper = " << ctau_proper << endl;
      }
      
      h_ctau_lab -> Fill(ctau_lab);
      h_ctau_proper->Fill(ctau_proper);
    }
  }
    


  TFile* outFile = new TFile(outputFile.data(),"recreate");
  h_ctau_lab->Write();
  h_ctau_proper->Write();

  outFile->Close();
  std::cout << "nTotal    = " << nTotal << std::endl;
  //for(int i=0;i<20;i++)
  //  if(nPass[i]>0)
  //    std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    
  
  

}
