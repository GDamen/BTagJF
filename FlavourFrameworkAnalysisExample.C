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
using namespace std;

float IntegrateEff(TH1F * histo, float efficiency) {

  Double_t tot_integral = 0;
  Double_t range_integral = 0;
  Double_t integral_ratio = 0;
  float bin_position = 0;
  Int_t last_bin = 0;
  Int_t current_bin = 0;
  tot_integral = histo->Integral();
  
  for(int i = 0; integral_ratio < efficiency; i++) {
    last_bin = histo->GetXaxis()->GetLast();
    current_bin = last_bin - i;
    range_integral = histo->Integral(current_bin, last_bin);
    integral_ratio = range_integral / tot_integral;
    std::cout << integral_ratio << std::endl;
  }
  std::cout << "Current Bin = " << current_bin-1 << std::endl;
  bin_position = histo->GetXaxis()->GetBinCenter(current_bin-1);
  std::cout << "Current Bin Position = " << bin_position << std::endl;
  return bin_position;
}

float EvaluateEff(TH1F * histo, float efficiency, float working_point) {
  Double_t bin_width = histo->GetXaxis()->GetBinWidth(0);
  Double_t bin_center = histo->GetXaxis()->GetBinCenter(0);
  Double_t working_bin = 0;
  
  working_point -= bin_center;
  working_bin = working_point/bin_width;
  
  Double_t tot_integral = histo->Integral();
  Double_t range_integral = histo->Integral(working_bin, histo->GetXaxis()->GetLast());
  
  float efficiencyRange = range_integral/tot_integral;
  
  std::cout << "eff = " << efficiencyRange << " ratio = " << efficiencyRange/efficiency << std::endl;
  return true;
}


void FlavourFrameworkAnalysisExample() {

  TCanvas *canvas = new TCanvas("Histocanvas","Graph2D example",0,0,800,800);
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);


  Float_t minPt = 350.;	//why these min and max pts?
  Float_t maxPt = 1000.;
  Float_t minTrkPt = 0.9;
  Float_t maxEta = 1.95;	//what ends in this eta?
  Float_t bEff = 0.70;		//70% expected
  
  ////////////////////////////////////////////////////////////////
  // TEST HISTOS:creating 1-dimensional histos with float variables
  TH1F * hPtB = new TH1F("hPtB","",100,0.,1000.);		//plot of High-Pt B
  TH2F * BPTvsJPT = new TH2F("momenta", "B Pt versus JetPt with flavour = 5", 200, 0., 1000., 200, 0., 1000.);
  TH2F * bHVTX = new TH2F("bH", "XY coordinates of b secondary vertex", 400, -200, 200, 200, -200, 200);
  ////////////////////////////////////////////////////////////////
  // Working Point finder Histo
  TH1F * IP3D_llr = new TH1F("llr", "Log Likelihood Ratio for b with IP3D", 1000, -20, 50);
  IP3D_llr->SetLineColor(kRed);
  TH1F * IP3D_llr_beampipe_01 = new TH1F("", "llr for b with IP3D in [0,1[ cm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_12 = new TH1F("", "llr for b with IP3D in [1,2[ cm range", 1000, -20, 50);
  TH1F * IP3D_llr_beampipe_23 = new TH1F("", "llr for b with IP3D in [2,3[ cm range", 1000, -20, 50);
  TH1F * IP3D_llr_IBL = new TH1F("", "llr for b with IP3D in [3,50[ cm range (IBL to B layer)", 1000, -20, 50);
  TH1F * IP3D_llr_B = new TH1F("", "llr for b with IP3D after 50 mm (after B layer)", 1000, -20, 50);
  
  TH1F * SV_llr = new TH1F("llr", "Log Likelihood Ratio for b with SV", 1000, -20, 50);
  SV_llr->SetLineColor(kBlue);
  TH1F * JF_llr = new TH1F("llr", "Log Likelihood Ratio for b with JF", 1000, -20, 50);
  JF_llr->SetLineColor(kGreen);
  ////////////////////////////////////////////////////////////////
  // Reconstruction efficiency histos
  TH2F * hBSV1Eff_pb = new TH2F("", "SV1 Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 0.6);
  hBSV1Eff_pb->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBSV1Eff_pb->GetYaxis()->SetTitle("reconstructed b probability");
  hBSV1Eff_pb->GetXaxis()->CenterTitle();
  hBSV1Eff_pb->GetXaxis()->CenterTitle();
  TH2F * hBJFEff_pb = new TH2F("", "JF Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 1.);
  hBJFEff_pb->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBJFEff_pb->GetYaxis()->SetTitle("reconstructed b probability");
  hBJFEff_pb->GetXaxis()->CenterTitle();
  hBJFEff_pb->GetXaxis()->CenterTitle();
  TH2F * hBSV1Eff_pc = new TH2F("", "SV1 Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 0.6);
  hBSV1Eff_pc->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBSV1Eff_pc->GetYaxis()->SetTitle("reconstructed c probability");
  hBSV1Eff_pc->GetXaxis()->CenterTitle();
  hBSV1Eff_pc->GetXaxis()->CenterTitle();
  TH2F * hBJFEff_pc = new TH2F("", "JF Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 1.);
  hBJFEff_pc->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBJFEff_pc->GetYaxis()->SetTitle("reconstructed c probability");
  hBJFEff_pc->GetXaxis()->CenterTitle();
  hBJFEff_pc->GetXaxis()->CenterTitle();
  TH2F * hBSV1Eff_pu = new TH2F("", "SV1 Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 0.6);
  hBSV1Eff_pu->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBSV1Eff_pu->GetYaxis()->SetTitle("reconstructed u probability");
  hBSV1Eff_pu->GetXaxis()->CenterTitle();
  hBSV1Eff_pu->GetXaxis()->CenterTitle();
  TH2F * hBJFEff_pu = new TH2F("", "JF Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 1.);
  hBJFEff_pu->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBJFEff_pu->GetYaxis()->SetTitle("reconstructed u probability");
  hBJFEff_pu->GetXaxis()->CenterTitle();
  hBJFEff_pu->GetXaxis()->CenterTitle();
  TH2F * hBSV1Eff_ptot = new TH2F("", "SV1 Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0., 1.);
  hBSV1Eff_ptot->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBSV1Eff_ptot->GetYaxis()->SetTitle("reconstructed total (u + b + c) probability");
  hBSV1Eff_ptot->GetXaxis()->CenterTitle();
  hBSV1Eff_ptot->GetXaxis()->CenterTitle();
  TH2F * hBJFEff_ptot = new TH2F("", "JF Reconstruction efficiency for hB versus distance from primary vertex", 100, 0., 150., 100, 0.4, 1.4);
  hBJFEff_ptot->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBJFEff_ptot->GetYaxis()->SetTitle("reconstructed total (u + b + c) probability");
  hBJFEff_ptot->GetXaxis()->CenterTitle();
  hBJFEff_ptot->GetXaxis()->CenterTitle();
  //Transverse distance
  TH1F * hBtranverseDistance = new TH1F("bH_Lxy","B transverse distance",100,0.,100.);
  hBtranverseDistance->GetXaxis()->SetTitle("transverse distance traveled by a B+- [GeV]");
  hBtranverseDistance->GetYaxis()->SetTitle("");
  hBtranverseDistance->GetXaxis()->CenterTitle();
  hBtranverseDistance->GetXaxis()->CenterTitle();
  
  TH1F * hPtJetQ = new TH1F("hPtJetQ","",20,0.,2000.);
  TH1F * hPtJetB = new TH1F("hPtJetB","",200,0.,2000.);
  TH1F * hEtaJetQ = new TH1F("hEtaJetQ","",10,0.,2.5);
  TH1F * hEtaJetB = new TH1F("hEtaJetB","",10,0.,2.5);
  
  TH1F * hNVtxB = new TH1F("hNVtxB","",4,-0.5,3.5);
  TH1F * hNVtxQ = new TH1F("hNVtxQ","",4,-0.5,3.5);
  TH1F * hNTkVtxQ = new TH1F("hNTkVtxQ","",9,1.5,10.5);
  TH1F * hNTkVtxB = new TH1F("hNTkVtxB","",9,1.5,10.5);
  TH1F * hMVtxQ = new TH1F("hMVtxQ","",14,0.01,7.01);
  TH1F * hMVtxB = new TH1F("hMVtxB","",14,0.01,7.01);
  TH1F * hProbIP3DQ = new TH1F("hProbIP3DQ","",100,-10,15);
  TH1F * hProbIP3DB = new TH1F("hProbIP3DB","",100,-10,15);
  TH1F * hQProbJFQ = new TH1F("hQProbJFQ","",15,0,1.01);
  TH1F * hBProbJFQ = new TH1F("hBProbJFQ","",15,0,1.01);
  TH1F * hQProbJFB = new TH1F("hQProbJFB","",15,0,1.01);
  TH1F * hBProbJFB = new TH1F("hBProbJFB","",15,0,1.01);
  TH1F * hProbJFQ = new TH1F("hProbJFQ","",100,-4,4);
  TH1F * hProbJFB = new TH1F("hProbJFB","",100,-4,4);
  
  TH1F * hQProbMV2Q = new TH1F("hQProbMV2Q","",15,0,1.01);
  TH1F * hBProbMV2Q = new TH1F("hBProbMV2Q","",15,0,1.01);
  TH1F * hQProbMV2B = new TH1F("hQProbMV2B","",15,0,1.01);
  TH1F * hBProbMV2B = new TH1F("hBProbMV2B","",15,0,1.01);
  TH1F * hProbMV2Q = new TH1F("hProbMV2Q","",100,-1.01,1.01);
  TH1F * hProbMV2B = new TH1F("hProbMV2B","",100,-1.01,1.01);
  TH1F * hPulld0Q = new TH1F("hPulld0Q","",50,-10,10);
  TH1F * hPullz0Q = new TH1F("hPullz0Q","",50,-10,10);
  TH1F * hPulld0B = new TH1F("hPulld0B","",50,-50,50);
  TH1F * hPullz0B = new TH1F("hPullz0B","",50,-50,50);
  TH1F * hProbSV1Q = new TH1F("hProbSV1Q","",35,-1,12);
  TH1F * hProbSV1B = new TH1F("hProbSV1B","",35,-1,12);
  
  TFile *myStdFile = TFile::Open("flav_Akt4EMTo_tt_Std.root");
  
  TTreeReader myStdReader("bTag_AntiKt4EMTopoJets", myStdFile);

  TTreeReaderValue<double> PVx(myStdReader, "PVx");		//DAMEN to plot non vectorized objects, pass the pointer e.g.*PVx
  TTreeReaderValue<double> PVy(myStdReader, "PVy");
  TTreeReaderValue<double> PVz(myStdReader, "PVz");
  TTreeReaderValue<double> truth_PVx(myStdReader, "truth_PVx");
  TTreeReaderValue<double> truth_PVy(myStdReader, "truth_PVy");
  TTreeReaderValue<double> truth_PVz(myStdReader, "truth_PVz");
  
  TTreeReaderValue<vector<float> > jetPt(myStdReader, "jet_pt");	//DAMEN to plot vectorized objects, pass the array element e.g. jetPt->at(i)
  TTreeReaderValue<vector<float> > jetEta(myStdReader, "jet_eta");
  TTreeReaderValue<vector<float> > jetPhi(myStdReader, "jet_phi");
  TTreeReaderValue<vector<int> >   jetFlav(myStdReader, "jet_truthflav");
  //IP3D values
  TTreeReaderValue<vector<float> > jetIP3DProbQ(myStdReader, "jet_ip3d_pu");
  TTreeReaderValue<vector<float> > jetIP3DProbB(myStdReader, "jet_ip3d_pb");
  TTreeReaderValue<vector<float> > jetIP3DProbC(myStdReader, "jet_ip3d_pc");
  TTreeReaderValue<vector<float> > jetIP3DProb(myStdReader, "jet_ip3d_llr");
  //JetFitter Values
  TTreeReaderValue<vector<float> > jetFitterProbQ(myStdReader, "jet_jf_pu");
  TTreeReaderValue<vector<float> > jetFitterProbB(myStdReader, "jet_jf_pb");
  TTreeReaderValue<vector<float> > jetFitterProbC(myStdReader, "jet_jf_pc");
  TTreeReaderValue<vector<float> > jetFitterProb(myStdReader, "jet_jf_llr");
  TTreeReaderValue<vector<float> > jetFitterCombProbQ(myStdReader, "jet_jfcombnn_pu");
  TTreeReaderValue<vector<float> > jetFitterCombProbB(myStdReader, "jet_jfcombnn_pb");
  TTreeReaderValue<vector<int> >   jetFitterNVtx(myStdReader, "jet_jf_nvtx");
  TTreeReaderValue<vector<float> > jetFitterNTrkAtVtx(myStdReader, "jet_jf_ntrkAtVx");
  TTreeReaderValue<vector<float> > jetFitterMass(myStdReader, "jet_jf_m");
  //SV1 values
  TTreeReaderValue<vector<float> > SV1ProbQ(myStdReader, "jet_sv1_pu");
  TTreeReaderValue<vector<float> > SV1ProbB(myStdReader, "jet_sv1_pb");
  TTreeReaderValue<vector<float> > SV1ProbC(myStdReader, "jet_sv1_pc");
  TTreeReaderValue<vector<float> > SV1Prob(myStdReader, "jet_sv1_llr");
  //MV2M Values
  TTreeReaderValue<vector<double> > mv2ProbQ(myStdReader, "jet_mv2m_pu");
  TTreeReaderValue<vector<double> > mv2ProbB(myStdReader, "jet_mv2m_pb");
  TTreeReaderValue<vector<double> > mv2ProbC00(myStdReader, "jet_mv2c00");
  //Secondary Vertex B values
  TTreeReaderValue<vector<float> > bHPt(myStdReader, "bH_pt");		//DAMEN added to check how we read variables
  TTreeReaderValue<vector<float> > bHX(myStdReader, "bH_x");		//DAMEN the 3 coordinates of the decay point
  TTreeReaderValue<vector<float> > bHY(myStdReader, "bH_y");
  TTreeReaderValue<vector<float> > bHZ(myStdReader, "bH_z");
  TTreeReaderValue<vector<float> > secVertB(myStdReader, "bH_Lxy");
  TTreeReaderValue<vector<float> > secVertC(myStdReader, "cH_Lxy");

  TTreeReaderValue<vector<vector<float> > > trkPt(myStdReader, "jet_trk_pt");
  TTreeReaderValue<vector<vector<float> > > trkEta(myStdReader, "jet_trk_eta");
  TTreeReaderValue<vector<vector<float> > > trkPhi(myStdReader, "jet_trk_phi");
  TTreeReaderValue<vector<vector<int> > >  trkNIBLHits(myStdReader, "jet_trk_nBLHits");
  TTreeReaderValue<vector<vector<int> > >  trkNSharedIBLHits(myStdReader, "jet_trk_nsharedBLHits");
  TTreeReaderValue<vector<vector<int> > >  trkNSplitIBLHits(myStdReader, "jet_trk_nsplitBLHits");
  TTreeReaderValue<vector<vector<int> > >  trkNPixelHits(myStdReader, "jet_trk_nPixHits");
  TTreeReaderValue<vector<vector<float> > >  trkd0(myStdReader, "jet_trk_d0");
  TTreeReaderValue<vector<vector<float> > >  trkz0(myStdReader, "jet_trk_z0");
  TTreeReaderValue<vector<vector<float> > >  trkMCd0(myStdReader, "jet_trk_d0_truth");
  TTreeReaderValue<vector<vector<float> > >  trkMCz0(myStdReader, "jet_trk_z0_truth");
  TTreeReaderValue<vector<vector<float> > >  trkIP3Dd0(myStdReader, "jet_trk_ip3d_d0");
  TTreeReaderValue<vector<vector<float> > >  trkIP3DPulld0(myStdReader, "jet_trk_ip3d_d0sig");
  TTreeReaderValue<vector<vector<float> > >  trkIP3Dz0(myStdReader, "jet_trk_ip3d_z0");
  TTreeReaderValue<vector<vector<float> > >  trkIP3DPullz0(myStdReader, "jet_trk_ip3d_z0sig");
  TTreeReaderValue<vector<vector<int> > >  trkOrigin(myStdReader, "jet_trk_orig");
  TTreeReaderValue<vector<vector<float> > >  jetSV0X(myStdReader, "jet_sv0_vtx_x");
  TTreeReaderValue<vector<vector<float> > >  jetSV0Y(myStdReader, "jet_sv0_vtx_y");
  TTreeReaderValue<vector<vector<float> > >  jetSV0Z(myStdReader, "jet_sv0_vtx_z");
  
  
  Int_t nJQ     = 0;	//number of light quark (jFlav = 0) jets
  Int_t nJB     = 0;	//number of bottom quark (jFlav = 5) jets
  Int_t nJIBLSQ = 0;
  Int_t nJIBLSB = 0;


    // Loop over TTree entries
  while (myStdReader.Next()) {

    //loop over jets
    for(int ij = 0; ij < jetPt->size(); ij++){
    
    if(jetFlav->at(ij) == 5){	//bottom flavour jet
    
    	IP3D_llr->Fill(jetIP3DProb->at(ij));
    	SV_llr->Fill(SV1Prob->at(ij));
    	JF_llr->Fill(jetFitterProb->at(ij));
    	
    	if(secVertB->at(ij) < 35) {
    	  IP3D_llr_beampipe_01->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 50) {
    	  IP3D_llr_beampipe_12->Fill(jetIP3DProb->at(ij));
    	}/*
    	else if(secVertB->at(ij) < 20) {
    	  IP3D_llr_beampipe_23->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 50) {
    	  IP3D_llr_IBL->Fill(jetIP3DProb->at(ij));
    	}
    	else if(secVertB->at(ij) < 20) {
    	  IP3D_llr_B->Fill(jetIP3DProb->at(ij));
    	}*/
    	
    	hBSV1Eff_pb->Fill(secVertB->at(ij), SV1ProbB->at(ij));
	hBJFEff_pb->Fill(secVertB->at(ij), jetFitterProbB->at(ij));
	
	hBSV1Eff_pc->Fill(secVertB->at(ij), SV1ProbC->at(ij));
	hBJFEff_pc->Fill(secVertB->at(ij), jetFitterProbC->at(ij));
	
	hBSV1Eff_pu->Fill(secVertB->at(ij), SV1ProbQ->at(ij));
	hBJFEff_pu->Fill(secVertB->at(ij), jetFitterProbQ->at(ij));
	
	hBSV1Eff_ptot->Fill(secVertB->at(ij), SV1ProbQ->at(ij) + SV1ProbC->at(ij) + SV1ProbB->at(ij));
	hBJFEff_ptot->Fill(secVertB->at(ij), jetFitterProbQ->at(ij) + jetFitterProbC->at(ij) + jetFitterProbB->at(ij));
	
	
	hBtranverseDistance->Fill(secVertB->at(ij));
    }
    
    
      if(jetPt->at(ij)/1000. > minPt && jetPt->at(ij)/1000. < maxPt && fabs(jetEta->at(ij)) < maxEta){

	int nIBLS=0;		//what's that?


	// jet flavour
	if(jetFlav->at(ij) == 0){	//light flavour jet
	  nJQ++;
	}else if(jetFlav->at(ij) == 5){	//bottom flavour jet
	  nJB++;
	}

	TVector3 pVtx(*PVx,*PVy,*PVz);
	
	TVector3 pTruthVtx(*truth_PVx,*truth_PVy,*truth_PVz);

	TVector3 thisJet;
	thisJet.SetPtEtaPhi(jetPt->at(ij)/1000.,jetEta->at(ij),jetPhi->at(ij));


	//DAMEN: B Pt versus jet Pt with flavour=5
	if(jetFlav->at(ij) == 5) {
	  //BPTvsJPT->Fill(bHPt->at(ij)/1000, jetPt->at(ij)/1000);
	  //bHVTX->Fill(bHX->at(ij), bHY->at(ij));
	}
	//#

	//loop over tracks in jet	
	for(int it = 0; it < trkPt->at(ij).size(); it++){
	  if(trkPt->at(ij).at(it)/1000. > minTrkPt && trkPt->at(ij).at(it) < 100000. && trkNPixelHits->at(ij).at(it)>1){
	    nIBLS++;
	    if(trkOrigin->at(ij).at(it) == 2){
	      hPulld0Q->Fill(trkIP3DPulld0->at(ij).at(it));
	      hPullz0Q->Fill(trkIP3DPullz0->at(ij).at(it));
	      
	    }else if(trkOrigin->at(ij).at(it) == 0 || trkOrigin->at(ij).at(it) == 1){
	      hPulld0B->Fill(trkIP3DPulld0->at(ij).at(it));
	      hPullz0B->Fill(trkIP3DPullz0->at(ij).at(it));
	    }
	  }
	}

	//JetFitter output for this jet
	if(jetFitterNVtx->size()){
	  if(jetFlav->at(ij) == 0 ){
	    nJIBLSQ++;
	    hPtJetQ->Fill(jetPt->at(ij)/1000.);
	    hEtaJetQ->Fill(jetEta->at(ij));
	    if(jetFitterNVtx->at(ij)<0){
	      hNVtxQ->Fill(jetFitterNVtx->at(ij)+1);
	    }else{
	      hNVtxQ->Fill(jetFitterNVtx->at(ij));
	    }
	    hMVtxQ->Fill(jetFitterMass->at(ij)/1000.);
	    hNTkVtxQ->Fill(jetFitterNTrkAtVtx->at(ij));
	    hProbIP3DQ->Fill(TMath::Log10(jetIP3DProbB->at(ij)/jetIP3DProbQ->at(ij)));
	    if(jetFitterProbQ->at(ij)>-1){
	      hQProbJFQ->Fill(jetFitterProbQ->at(ij));
	      hProbJFQ->Fill(TMath::Log10(jetFitterProbB->at(ij)/jetFitterProbQ->at(ij)));
	    }
	    if(jetFitterProbB->at(ij)>-1)hBProbJFQ->Fill(jetFitterProbB->at(ij));
	    hQProbMV2Q->Fill(mv2ProbQ->at(ij));
	    hBProbMV2Q->Fill(mv2ProbB->at(ij));
	    if(mv2ProbQ->at(ij)>-1)hProbMV2Q->Fill(mv2ProbC00->at(ij));
	    if(ij<SV1Prob->size())hProbSV1Q->Fill(SV1Prob->at(ij));
	  }else if(jetFlav->at(ij) ==5 ){
	    nJIBLSB++;
	    hPtJetB->Fill(jetPt->at(ij)/1000.);
	    hEtaJetB->Fill(jetEta->at(ij));
	    if(jetFitterNVtx->at(ij)<0){
	      hNVtxB->Fill(jetFitterNVtx->at(ij)+1);
	    }else{
	      hNVtxB->Fill(jetFitterNVtx->at(ij));
	    }
	    hNTkVtxB->Fill(jetFitterNTrkAtVtx->at(ij));
	    hMVtxB->Fill(jetFitterMass->at(ij)/1000.);
	    // Jet Fitter probabilities for this jet
	    if(jetFitterProbQ->at(ij)>-1){
	      hQProbJFB->Fill(jetFitterProbQ->at(ij));
	      //b prob is ratio of b/q 
	      hProbJFB->Fill(TMath::Log10(jetFitterProbB->at(ij)/jetFitterProbQ->at(ij)));
	    }
	    if(jetFitterProbB->at(ij)>-1)hBProbJFB->Fill(jetFitterProbB->at(ij));

	    //MV2, SV1 and IP3D probabilites for this jet
	    hQProbMV2B->Fill(mv2ProbQ->at(ij));
	    hBProbMV2B->Fill(mv2ProbB->at(ij));
	    if(mv2ProbB->at(ij)>-1)hProbMV2B->Fill(mv2ProbC00->at(ij));
	    if(ij<SV1Prob->size())hProbSV1B->Fill(SV1Prob->at(ij));
	    //b prob is ratio of b/q 
	    hProbIP3DB->Fill(TMath::Log10(jetIP3DProbB->at(ij)/jetIP3DProbQ->at(ij)));
	  }
	}
	
      }
    }
  }	//added a Draw function, otherwise don't work
  //hPtJetB->Draw("SAME");
  //hBSV1Eff_pb -> Draw("colz");
  //hBJFEff_pb -> Draw("colz");
  //hBSV1Eff_pc -> Draw("colz");
  //hBJFEff_pc -> Draw("colz");
  //hBSV1Eff_pu -> Draw("colz");
  //hBJFEff_pu -> Draw("colz");
  //hBSV1Eff_ptot -> Draw("colz");
  //hBJFEff_ptot -> Draw("colz");
  //hBtranverseDistance->Draw();
  //JF_llr->Draw();
  //SV_llr->Draw("SAME");
  float working_point_IP3D = IntegrateEff(IP3D_llr, 0.7);
  TLine * line = new TLine(working_point_IP3D, 0, working_point_IP3D, 250);
  //IP3D_llr->Draw();
  EvaluateEff(IP3D_llr_beampipe_12, 0.7, working_point_IP3D);
  //line->Draw("SAME");
  //canvas->SetLogy();	//log scale on y
}
