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

#include "functions.C"
using namespace std;

void EffVsLxyJF_c() {
  
  TCanvas *canvas = new TCanvas("Histocanvas","Graph2D example",0,0,800,800);
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 2.1;
  Float_t bEff = 0.70;
  
  Float_t llrJFC;
  
  /////////////////////////////
  // HISTOs
  TH1F * JF_llr = new TH1F("llr", "Log Likelihood Ratio for c with JF", 1000, -20, 50);
  TH1F * JF_llr_beampipe_03 = new TH1F("", "llr for c with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_06 = new TH1F("", "llr for c with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_10 = new TH1F("", "llr for c with JF in [3,10[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_18 = new TH1F("", "llr for c with JF in [0,3[ mm range", 1000, -20, 50);
  TH1F * JF_llr_beampipe_25 = new TH1F("", "llr for c with JF in [10,25[ cm range", 1000, -20, 50);
  TH1F * JF_llr_IBL = new TH1F("", "llr for c with JF in [25,35[ cm range (around IBL)", 1000, -20, 50);
  TH1F * JF_llr_B = new TH1F("", "llr for c with JF after 35 mm (IBL to B layer)", 1000, -20, 50);
  TH1F * JF_llr_afterB = new TH1F("", "llr for c with JF after 35 mm (IBL to B layer)", 1000, -20, 50);

  TGraph * JF_efficiency_Lxy = new TGraph(8);
  	JF_efficiency_Lxy->SetLineColor(2);
  	JF_efficiency_Lxy->SetLineWidth(4);
  	JF_efficiency_Lxy->SetMarkerColor(4);
  	JF_efficiency_Lxy->SetMarkerSize(1.5);
  	JF_efficiency_Lxy->SetMarkerStyle(21);
  	JF_efficiency_Lxy->SetTitle("Efficiency degradation for JF reconstruction");
  	JF_efficiency_Lxy->GetXaxis()->SetTitle("Distance from IP [mm]");
  	JF_efficiency_Lxy->GetYaxis()->SetTitle("eff(Lxy)/eff(TOT)");

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

  //JF values
  TTreeReaderValue<vector<float> > jetJFProbQ(myStdReader, "jet_jf_pu");
  TTreeReaderValue<vector<float> > jetJFProbB(myStdReader, "jet_jf_pb");
  TTreeReaderValue<vector<float> > jetJFProbC(myStdReader, "jet_jf_pc");
  TTreeReaderValue<vector<float> > llrJFB(myStdReader, "jet_jf_llr");

  ///////////////////////////////
  // FILL HISTOs
  // Loop over TTree entries
  while (myStdReader.Next()) {

    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
    if(jetFlav->at(ij) == 4 && jetJFProbQ->at(ij) > -1 && jetJFProbC->at(ij) > -1 && jetJFProbB->at(ij) > -1){	//charm flavour jet
    
    	llrJFC = jetJFProbC->at(ij) / (jetJFProbB->at(ij) + jetJFProbQ->at(ij));
    
    	JF_llr->Fill(llrJFC);
    	
    	if(secVertC->at(ij) < 3) {
    	  JF_llr_beampipe_03->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 6) {
    	  JF_llr_beampipe_06->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 10) {
    	  JF_llr_beampipe_10->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 18) {
    	  JF_llr_beampipe_18->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 25) {
    	  JF_llr_beampipe_25->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 35) {
    	  JF_llr_IBL->Fill(llrJFC);
    	}
    	else if(secVertC->at(ij) < 50){
    	  JF_llr_B->Fill(llrJFC);
    	}
    	else {
    	  JF_llr_afterB->Fill(llrJFC);
    	}
      }
    }
  }
  //////////////////////////////////////
  // Working Point evaluation
  
  float working_point_JF = IntegrateEff(JF_llr, bEff);
  TLine * line = new TLine(working_point_JF, 0, working_point_JF, 250);
  //JF_llr->Draw();
  float eff_03 = EvaluateEff(JF_llr_beampipe_03, bEff, working_point_JF);
  float eff_06 = EvaluateEff(JF_llr_beampipe_06, bEff, working_point_JF);
  float eff_10 = EvaluateEff(JF_llr_beampipe_10, bEff, working_point_JF);
  float eff_18 = EvaluateEff(JF_llr_beampipe_18, bEff, working_point_JF);
  float eff_25 = EvaluateEff(JF_llr_beampipe_25, bEff, working_point_JF);
  float eff_IBL = EvaluateEff(JF_llr_IBL, bEff, working_point_JF);
  float eff_B = EvaluateEff(JF_llr_B, bEff, working_point_JF);
  float eff_afterB = EvaluateEff(JF_llr_afterB, bEff, working_point_JF);
  
  JF_efficiency_Lxy -> SetPoint(0, 3, eff_03);
  JF_efficiency_Lxy -> SetPoint(1, 6, eff_06);
  JF_efficiency_Lxy -> SetPoint(2, 10, eff_10);
  JF_efficiency_Lxy -> SetPoint(3, 18, eff_18);
  JF_efficiency_Lxy -> SetPoint(4, 25, eff_25);
  JF_efficiency_Lxy -> SetPoint(5, 35, eff_IBL);
  JF_efficiency_Lxy -> SetPoint(6, 50, eff_B);
  JF_efficiency_Lxy -> SetPoint(7, 60, eff_afterB);
  JF_efficiency_Lxy -> Draw("ALP");
  //JF_llr -> Draw();
}


