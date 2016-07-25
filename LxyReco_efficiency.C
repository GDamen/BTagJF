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


void LxyReco_efficiency() {

	TCanvas *canvas = new TCanvas("Reconstruction of Lxy","Graph2D example",0,0,800,800);
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;

	TH2F * MultRecoEffSV = new TH2F("", "", 200, 0., 50, 32, -16, 16);
  	MultRecoEffSV->SetTitle("Track Reconstruction per Vertex by SV1 efficiency");
  	MultRecoEffSV->GetXaxis()->SetTitle("Transverse distance PV - SV [mm]");
  	MultRecoEffSV->GetYaxis()->SetTitle("Reconstructed - Generated vertex multiplicity");
  	
  TH2F * MultRecoEffJF = new TH2F("", "", 200, 0., 50, 32, -16, 16);
  	MultRecoEffJF->SetTitle("Track Reconstruction per Vertex by JF efficiency");
  	MultRecoEffJF->GetXaxis()->SetTitle("Transverse distance PV - SV [mm]");
  	MultRecoEffJF->GetYaxis()->SetTitle("Reconstructed - Generated vertex multiplicity");
  	
  TGraph * MultRecoEff_value = new TGraph(7);


	
	////////////////////////////
  // TREE READER
  
  TFile *myStdFile = TFile::Open("ntuples/flav_Akt4EMTo_tt_Std.root");
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  TTreeReaderValue<vector<float> > secVertB(myStdReader, "bH_Lxy");
  TTreeReaderValue<vector<float> > secVertC(myStdReader, "cH_Lxy");
	
  TTreeReaderValue<vector<float> > bHX(myStdReader, "bH_x");
  TTreeReaderValue<vector<float> > bHY(myStdReader, "bH_y");
  TTreeReaderValue<vector<float> > bHZ(myStdReader, "bH_z");
  
  TTreeReaderValue<vector<int> > bHnumBtracks(myStdReader, "bH_nBtracks");
  TTreeReaderValue<vector<int> > bHnumCtracks(myStdReader, "bH_nCtracks");
  TTreeReaderValue<vector<int> > cHnumCtracks(myStdReader, "cH_nCtracks");
  
  TTreeReaderValue<vector<int> > sv1NumTracks(myStdReader, "jet_sv1_ntrk");
  TTreeReaderValue<vector<int> > jfNumTracks(myStdReader, "jet_jf_ntrk");
  
  
  		float diff3 = 0;
  		float diff6 = 0;
  		float diff10 = 0;
  		float diff18 = 0;
  		float diff25 = 0;
  		float diff35 = 0;
  		float diff50 = 0;
  		
  		int iter3 = 0;
  		int iter6 = 0;
  		int iter10 = 0;
  		int iter18 = 0;
  		int iter25 = 0;
  		int iter35 = 0;
  		int iter50 = 0;
  
  while (myStdReader.Next()) {
  
  		
  
  	for(int ij = 0; ij < jetPt->size(); ij++){
  	
  		
  		if(secVertB->at(ij) >= 0 && jetFlav->at(ij) == 5 ) {
  		
  			Float_t diff = ((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij));
  			MultRecoEffSV->Fill(secVertB->at(ij), diff);
  			
  		}
  		
  		if(secVertB->at(ij) >= 0 && jetFlav->at(ij) == 5 && bHnumCtracks->at(ij) == 0) {
  		 	Float_t diff = ((float)jfNumTracks->at(ij) - (float)bHnumBtracks->at(ij));
  			MultRecoEffJF->Fill(secVertB->at(ij), diff);

  		}
  		
  	}
  }

  MultRecoEffSV->Draw("colz");
  //MultRecoEffJF->Draw("colz");

}
