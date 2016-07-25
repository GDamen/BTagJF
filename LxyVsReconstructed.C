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

void LxyVsReconstructed() {

	TCanvas *canvas = new TCanvas("Reconstruction of Lxy","Graph2D example",0,0,800,800);
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;

	TH2F * LxVsReco = new TH2F("", "", 1000, 0., 100, 1000, 0., 100);
  	LxVsReco->SetTitle("Vertex reconstruction on X axis");
  	LxVsReco->GetXaxis()->SetTitle("Truth level distance from IP on x [mm]");
  	LxVsReco->GetYaxis()->SetTitle("Reconstructed/Truth level distance from the IP");
  	
  TH2F * LyVsReco = new TH2F("", "", 1000, 0., 100, 1000, 0., 100);
  	LyVsReco->SetTitle("Vertex reconstruction on Y axis");
  	LyVsReco->GetXaxis()->SetTitle("Truth level distance from IP on y [mm]");
  	LyVsReco->GetYaxis()->SetTitle("Reconstructed/Truth level distance from the IP");
  	
  TH2F * LxyVsReco = new TH2F("", "", 1000, 0., 50, 1000, 0., 50);
  	LxyVsReco->SetTitle("Vertex reconstruction on XY plane");
  	LxyVsReco->GetXaxis()->SetTitle("Truth level distance from IP on transverse plane [mm]");
  	LxyVsReco->GetYaxis()->SetTitle("Reconstructed level distance from IP on transverse plane [mm]");
	
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
  		
  		if(bHX->at(ij) > -600) {
  			for(int iv = 0; iv < jetSV0X->at(ij).size(); iv++) {
  				LxVsReco->Fill(bHX->at(ij), jetSV0X->at(ij).at(iv));
  			}
  		}
  			
  		if(bHY->at(ij) > -600) {
  			for(int iv = 0; iv < jetSV0Y->at(ij).size(); iv++) {
  				LyVsReco->Fill(bHY->at(ij), jetSV0Y->at(ij).at(iv));
  			}
  		}
  			
  		if(secVertB->at(ij) >= 0) {
  			for(int iv = 0; iv < jetSV0X->at(ij).size(); iv++) {
  				Float_t distance = sqrt(jetSV0X->at(ij).at(iv)*jetSV0X->at(ij).at(iv) + jetSV0Y->at(ij).at(iv)*jetSV0Y->at(ij).at(iv));
  				LxyVsReco->Fill(bHX->at(ij), distance);
  			}
  		}
  	}
  }

  //LxVsReco->Draw("colz");
  //LyVsReco->Draw("colz");
  LxyVsReco->Draw("colz");
}
