#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMatrixD.h"
#include "TMath.h"
#include <iostream>

#ifndef FUNCTIONS
#include "functions.C"
#endif

using namespace std;

void vtxSeparation() {

	TCanvas *canvas = new TCanvas("JFEffCanvas","Graph2D example",0,0,800,800);
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
  
  int totalEvs = 0;
  int usedEvs = 0;
  float jet_pt_jet_i, jet_eta_jet_i, jet_JVT_jet_i;
  
  /////////////////////////////
  // HISTOs

  TH1F * ratioDist = new TH1F("", "", 100, -1.0, 2.0);

  ////////////////////////////
  // TREE READER
  
  TFile *myStdFile = TFile::Open("../ntuples/flav_Akt4EMTo_newVars.root");
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  TTreeReaderValue<vector<float> > jetJVT(myStdReader, "jet_JVT");

  TTreeReaderValue<vector<float> > bHadX(myStdReader, "bH_x");
	TTreeReaderValue<vector<float> > bHadY(myStdReader, "bH_y");
	TTreeReaderValue<vector<float> > bHadZ(myStdReader, "bH_z");
	TTreeReaderValue<vector<float> > cHadX(myStdReader, "cH_x");
	TTreeReaderValue<vector<float> > cHadY(myStdReader, "cH_y");
	TTreeReaderValue<vector<float> > cHadZ(myStdReader, "cH_z");

	TTreeReaderValue<vector<vector<float> > > recoVtxX(myStdReader, "jet_jf_vtx_x");
	TTreeReaderValue<vector<vector<float> > > recoVtxY(myStdReader, "jet_jf_vtx_y");
	TTreeReaderValue<vector<vector<float> > > recoVtxZ(myStdReader, "jet_jf_vtx_z");
	TTreeReaderValue<vector<vector<float> > > recoVtxL3D(myStdReader, "jet_jf_vtx_L3D");
	
	TTreeReaderValue<double> tPVX(myStdReader, "truth_PVx");
	TTreeReaderValue<double> tPVY(myStdReader, "truth_PVx");
	TTreeReaderValue<double> tPVZ(myStdReader, "truth_PVx");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {

    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
    jet_pt_jet_i = jetPt->at(ij);
    jet_eta_jet_i = fabs(jetEta->at(ij));
		jet_JVT_jet_i = jetJVT->at(ij);
    
			if(jetFlav->at(ij) == 5 && 
					(jet_JVT_jet_i > 0.59 || jet_pt_jet_i >60000.0 || fabs(jet_eta_jet_i)>2.4)) {
				
				bool vertexCheck = false;
				int numOfGoodVtx = 0;
				int goodVtx = 0;
			
				if(recoVtxL3D->at(ij).size() == 1 && recoVtxL3D->at(ij).at(0) >= 0) {
				
					vertexCheck = true;
					goodVtx = 0;
					
				} else {
				
					for(int iv = 0; iv < recoVtxL3D->at(ij).size(); iv++) {
						
						if(recoVtxL3D->at(ij).at(iv) >= 0) {
							numOfGoodVtx ++;
							goodVtx = iv;
						}
					}
					if(numOfGoodVtx == 1) {
						vertexCheck = true;
					}
				}
			
				if (vertexCheck) {

					float dReco = std::sqrt(recoVtxX->at(ij).at(goodVtx)*recoVtxX->at(ij).at(goodVtx) + recoVtxY->at(ij).at(goodVtx)*recoVtxY->at(ij).at(goodVtx) + recoVtxZ->at(ij).at(goodVtx)*recoVtxZ->at(ij).at(goodVtx));
					float dBt = std::sqrt(bHadX->at(ij)*bHadX->at(ij) + bHadY->at(ij)*bHadY->at(ij) + bHadZ->at(ij)*bHadZ->at(ij)); 
					float dCt = std::sqrt(cHadX->at(ij)*cHadX->at(ij) + cHadY->at(ij)*cHadY->at(ij) + cHadZ->at(ij)*cHadZ->at(ij)); 

					float ratio = (dReco - dBt) / (dCt - dBt) ;
					ratioDist -> Fill(ratio);
					usedEvs++;
				}
			}
			totalEvs++; 
		}
  }
  
float P = (float)usedEvs / (float)totalEvs;
cout << usedEvs << endl << totalEvs << endl << P << endl;
  
ratioDist -> Draw();

}

