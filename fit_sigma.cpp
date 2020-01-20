#include "sim_double_gaussian.cpp"

void fit_sigma()
{
    TFile * f = new TFile("exp65_4S_merged.root");
    TTree * tree = (TTree*) f->Get("Lambda");
    TH1F * hist1 = new TH1F("hist1", "hist1", 200, 1.105, 1.125);
    TH1F * hist2 = new TH1F("hist2", "hist2", 200, 1.105, 1.125);
    tree->Draw("M >> hist1", "mva >= .1");
    tree->Draw("M >> hist2", "mva < .1");

    sim_double_gaussian(hist1, hist2);
    
}