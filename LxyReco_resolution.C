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

void LxyReco_resolution() {

	TCanvas *canvas = new TCanvas("Reconstruction of Lxy","Graph2D example",0,0,800,800);
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.5;
  Float_t bEff = 0.70;
	
	// B decay SV
	TH1F * dLxySV1 = new TH1F("", "", 200, -4, 4);
  	dLxySV1->SetTitle("SV reconstruction resolution");
  	dLxySV1->GetXaxis()->SetTitle("Transverse distance (reconstructed - generated) from the IP [mm]");
  	dLxySV1->GetYaxis()->SetTitle("events");
  	
	TH1F * dLxySV1_bHCuts = new TH1F("", "", 100, -4, 4);
	
	TH1F * dLxyJF = new TH1F("", "", 200, -10, 10);
  	dLxyJF->SetTitle("JF reconstruction resolution");
  	dLxyJF->GetXaxis()->SetTitle("Transverse distance (reconstructed - generated) from the IP [mm]");
  	dLxyJF->GetYaxis()->SetTitle("events");
	
	TH1F * dLxyJF_bHCuts = new TH1F("", "", 100, -10, 10);
	
	////////////////////////////
  // TREE READER
  
  //TFile *myStdFile = TFile::Open("ntuples/flav_Akt4EMTo_MC15ZPrime1000_tt_Std_01.root");
  
  TFile *myStdFile = TFile::Open("ntuples/flav_Akt4EMTo_tt_Std.root");
  
  //TFile *myStdFile = TFile::Open("ntuples/group.det-ibl.6666608.BTAGSTREAM._000001.root");
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
  
  TTreeReaderValue<vector<vector<float> > >  jetSV1X(myStdReader, "jet_sv1_vtx_x");
  TTreeReaderValue<vector<vector<float> > >  jetSV1Y(myStdReader, "jet_sv1_vtx_y");
  TTreeReaderValue<vector<vector<float> > >  jetSV1Z(myStdReader, "jet_sv1_vtx_z");
  
  TTreeReaderValue<vector<float> > bHX(myStdReader, "bH_x");
  TTreeReaderValue<vector<float> > bHY(myStdReader, "bH_y");
  TTreeReaderValue<vector<float> > bHZ(myStdReader, "bH_z");
  
  TTreeReaderValue<vector<int> > bHnumBtracks(myStdReader, "bH_nBtracks");
  TTreeReaderValue<vector<int> > bHnumCtracks(myStdReader, "bH_nCtracks");
  //TTreeReaderValue<vector<int> > cHnumCtracks(myStdReader, "cH_nCtracks");
  
  TTreeReaderValue<vector<float> > JFtheta(myStdReader, "jet_jf_theta");
  TTreeReaderValue<vector<vector<float> > >  jetJFL3D(myStdReader, "jet_jf_vtx_L3D");
  
  
  while (myStdReader.Next()) {
  	
  	// SV1
  	
  	for(int ij = 0; ij < jetPt->size(); ij++){
  	
  	if(secVertB->at(ij) >= 0 && secVertB->at(ij) < 30 && jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta) {
  			for(int iv = 0; iv < jetJFL3D->at(ij).size(); iv++) {
  			if(jetJFL3D->at(ij).at(iv) > 0) {
  				Float_t distance = jetJFL3D->at(ij).at(iv)*TMath::Sin(JFtheta->at(ij)) - secVertB->at(ij);
  				dLxyJF->Fill(distance);
  				}
  			}
  		}
  		
  		if(secVertB->at(ij) >= 0 && secVertB->at(ij) < 30 && jetFlav->at(ij) == 5 && bHnumCtracks->at(ij) == 0 && jetEta->at(ij) < maxEta) {
  			for(int iv = 0; iv < jetJFL3D->at(ij).size(); iv++) {
  			if(jetJFL3D->at(ij).at(iv) > 0) {
  					Float_t distance = jetJFL3D->at(ij).at(iv)*TMath::Sin(JFtheta->at(ij)) - secVertB->at(ij);
  				 	dLxyJF_bHCuts->Fill(distance);
  				 }
  			}
  		}
  		
  		if(secVertB->at(ij) >= 0 && secVertB->at(ij) < 30 && jetFlav->at(ij) == 5 && jetEta->at(ij) < maxEta) {
  			for(int iv = 0; iv < jetSV1X->at(ij).size(); iv++) {
  				Float_t distance = sqrt(jetSV1X->at(ij).at(iv)*jetSV1X->at(ij).at(iv) + jetSV1Y->at(ij).at(iv)*jetSV1Y->at(ij).at(iv)) - secVertB->at(ij);
  				dLxySV1->Fill(distance);
  			}
  		}
  		
  		if(secVertB->at(ij) >= 0 && secVertB->at(ij) < 30 && jetFlav->at(ij) == 5 && bHnumCtracks->at(ij) == 0 && jetEta->at(ij) < maxEta) {
  			for(int iv = 0; iv < jetSV1X->at(ij).size(); iv++) {
  				Float_t distance = sqrt(jetSV1X->at(ij).at(iv)*jetSV1X->at(ij).at(iv) + jetSV1Y->at(ij).at(iv)*jetSV1Y->at(ij).at(iv)) - secVertB->at(ij);
  				dLxySV1_bHCuts->Fill(distance);
  			}
  		}
  	}
  }
  // Gaussian Fits

	dLxySV1->Fit("gaus");
	dLxySV1_bHCuts->Fit("gaus");
	dLxySV1_bHCuts->GetFunction("gaus")->SetLineColor(kBlue);
	//dLxyJF->Fit("gaus");
	//dLxyJF_bHCuts->Fit("gaus");

  TLine *line = new TLine(0,0,0,500);

	// Draw Calls
  //dLxySV1->Draw("colz");
  // dLxyJF->Draw("colz");
  // dLxyJF_bHCuts->Draw("SAME");
  //dLxySV1_bHCuts->Draw("SAME");

  
  line -> Draw("SAME");

}
