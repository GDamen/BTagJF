#include "EffVsLxyJF_b.C"
#include "EffVsLxyIP3D_b.C"
#include "EffVsLxySV_b.C"
#include "EffVsLxyMV2_b.C"

void EffVsLxy_comparison() {

#ifndef CANVAS
TCanvas *canvas = new TCanvas("AllEffCanvas","Graph2D example",0,0,800,800);
#define CANVAS
#endif

leg = new TLegend(0.1,0.7,0.48,0.9);
leg->SetHeader("Reconstruction Algorithms");


EffVsLxyIP3D_b();
EffVsLxyJF_b();
EffVsLxySV_b();
EffVsLxyMV2_b();
}
