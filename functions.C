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

#define FUNCTIONS

// Given an Histogram and a value of efficiency, integrates from right to left until it finds the last bin in which the efficiency is less than the required one. Returns a float with the position of this bin.
float IntegrateEff(TH1F * histo, Float_t efficiency) {
  
  std::cout << "FUNCTIONS: Integrate Efficiency (IntegrateEff)" << std::endl;

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
  }
  std::cout << "Bin at " << efficiency << " efficiency = " << current_bin-1 << std::endl;
  bin_position = histo->GetXaxis()->GetBinCenter(current_bin-1);
  std::cout << "In position = " << bin_position << std::endl;
  return bin_position;
}


// Given an Histogram, a working point and a value of efficiency, returns the ratio between the efficiency given an the one in the histo in the selected working point
float EvaluateEff(TH1F * histo, Float_t efficiency, float working_point) {

  std::cout << std::endl;
  std::cout << "FUNCTIONS: Evaluate Efficiency (EvaluateEff)" << std::endl;

  Double_t bin_width = histo->GetXaxis()->GetBinWidth(0);
  Double_t bin_center = histo->GetXaxis()->GetBinCenter(0);
  Double_t working_bin = 0;
  
  working_point -= bin_center;
  working_bin = working_point/bin_width;
  
  Double_t tot_integral = histo->Integral();
  Double_t range_integral = histo->Integral(working_bin, histo->GetXaxis()->GetLast());
  
  float efficiency_range = range_integral/tot_integral;
  
  std::cout << "Efficiency in the working point " << working_point << " = " << efficiency_range << std::endl << "Ratio with total efficiency = " << efficiency_range/efficiency << std::endl;
  return efficiency_range/efficiency;
}

double FindEff(TH1F * histoA, TH1F * histoB) {

  std::cout << std::endl;
  std::cout << "FUNCTIONS: Find Efficiency (FindEff)" << std::endl;


  Double_t bin_width = histoA->GetXaxis()->GetBinWidth(0);
  Int_t last_binA = histoA->GetXaxis()->GetLast();
  Int_t current_bin = 0;
  Double_t histoA_content = 0;
  Double_t histoB_content = 0;
  Double_t range_integral = 0;
  Double_t tot_integral = histoA->Integral();
  cout << tot_integral << std::endl;
  /*
  for(int i = 0; done == true; i++) {
  	current_bin = last_binA - i;
  
  	histoA_content = histoA->GetBinContent(current_bin-500);
  	histoB_content = histoB->GetBinContent(current_bin-500);
  	
  	if(histoA_content < histoB_content) done = true;
  }
  */
  range_integral = histoA->Integral(344, last_binA);

return range_integral/tot_integral;
}





