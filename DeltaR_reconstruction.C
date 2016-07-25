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

void DeltaR_reconstruction() {

	TCanvas *canvas = new TCanvas("Reconstruction of Lxy","Graph2D example",0,0,800,800);
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;

	TH1F * positionDiff = new TH1F("", "", 1000, -20, 20);
  	positionDiff->SetTitle("Diff");
  	positionDiff->GetXaxis()->SetTitle("distance");
  	positionDiff->GetYaxis()->SetTitle("lol");

	
	////////////////////////////
  // TREE READER
  
  TFile *myStdFile = TFile::Open("flav_Akt4EMTo_tt_Std.root");
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  TTreeReaderValue<vector<float> > secVertB(myStdReader, "bH_Lxy");
  TTreeReaderValue<vector<float> > secVertC(myStdReader, "cH_Lxy");
	TTreeReaderValue<vector<vector<float> > >  jetSV0X(myStdReader, "jet_sv0_vtx_x");
  TTreeReaderValue<vector<vector<float> > >  jetSV0Y(myStdReader, "jet_sv0_vtx_y");
  TTreeReaderValue<vector<vector<float> > >  jetSV0Z(myStdReader, "jet_sv0_vtx_z");
  TTreeReaderValue<vector<float> > bHX(myStdReader, "bH_x");
  TTreeReaderValue<vector<float> > bHY(myStdReader, "bH_y");
  TTreeReaderValue<vector<float> > bHZ(myStdReader, "bH_z");
  
  while (myStdReader.Next()) {
  
  	for(int ij = 0; ij < jetPt->size(); ij++){
  		
  		if(secVertB->at(ij) >= 0 && secVertB->at(ij) < 30) {
  			for(int iv = 0; iv < jetSV0X->at(ij).size(); iv++) {
  				Float_t distance = sqrt(jetSV0X->at(ij).at(iv)*jetSV0X->at(ij).at(iv) + jetSV0Y->at(ij).at(iv)*jetSV0Y->at(ij).at(iv)) - secVertB->at(ij);
  				positionDiff->Fill(distance);
  			}
  		}
  	}
  }
	positionDiff->Fit("gaus");
  positionDiff->Draw("colz");
}
