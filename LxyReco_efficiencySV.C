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


void LxyReco_efficiencySV() {
  
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
  	
  TGraph * MultRecoEff_valueSV = new TGraph(7);


	
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
  	
  		if(secVertB->at(ij) >= 0  && secVertB->at(ij) < 3 && jetFlav->at(ij) == 5) {
  			diff3 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter3++;
  		}
  		if(secVertB->at(ij) >= 3  && secVertB->at(ij) < 6 && jetFlav->at(ij) == 5) {
  			diff6 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter6++;
  		}
  		if(secVertB->at(ij) >= 6  && secVertB->at(ij) < 10 && jetFlav->at(ij) == 5) {
  			diff10 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter10++;
  		}
  		if(secVertB->at(ij) >= 10  && secVertB->at(ij) < 18 && jetFlav->at(ij) == 5) {
  			diff18 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter18++;
  		}
  		if(secVertB->at(ij) >= 18  && secVertB->at(ij) < 25 && jetFlav->at(ij) == 5) {
  			diff25 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter25++;
  		}
  		if(secVertB->at(ij) >= 25  && secVertB->at(ij) < 35 && jetFlav->at(ij) == 5) {
  			diff35 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter35++;
  		}
  		if(secVertB->at(ij) >= 35  && secVertB->at(ij) < 50 && jetFlav->at(ij) == 5) {
  			diff50 += sqrt((((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))) * (((float)sv1NumTracks->at(ij)) - ((float)bHnumBtracks->at(ij))));
				iter50++;
  		}
  	}
  }
  
  diff3 /= iter3;
  diff6 /= iter6;
  diff10 /= iter10;
  diff18 /= iter18;
  diff25 /= iter25;
  diff35 /= iter35;
  diff50 /= iter50;
  MultRecoEff_valueSV -> SetPoint(0, 3, diff3);
  MultRecoEff_valueSV -> SetPoint(1, 6, diff6);
  MultRecoEff_valueSV -> SetPoint(2, 10, diff10);
  MultRecoEff_valueSV -> SetPoint(3, 18, diff18);
  MultRecoEff_valueSV -> SetPoint(4, 25, diff25);
  MultRecoEff_valueSV -> SetPoint(5, 35, diff35);
  MultRecoEff_valueSV -> SetPoint(6, 50, diff50);
  //MultRecoEffSV->Draw("colz");
  //MultRecoEffJF->Draw("colz");
  MultRecoEff_valueSV->SetTitle("SV multiplicity reconstruction");
  MultRecoEff_valueSV->SetLineColor(1);
  MultRecoEff_valueSV->SetLineWidth(2);

  MultRecoEff_valueSV->GetXaxis()->SetTitle("Distance from IP [mm]");
  MultRecoEff_valueSV->GetYaxis()->SetTitle("Average |MultReco - MultGen|");
  MultRecoEff_valueSV->Draw("alp");
  TText * lineName = new TText(50, diff50, "SV");
     lineName->SetTextColor(1);
  lineName -> Draw("SAME");
}
