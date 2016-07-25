#include "LxyReco_efficiencyJF.C"
#include "LxyReco_efficiencySV.C"

void AverageMultiplicityComparison() {

#ifndef CANVAS
TCanvas *canvas = new TCanvas("AllEffCanvas","Graph2D example",0,0,800,800);
#define CANVAS
#endif

leg = new TLegend(0.1,0.7,0.48,0.9);
leg->SetHeader("Reconstruction Algorithms");

LxyReco_efficiencySV();
LxyReco_efficiencyJF();

}
