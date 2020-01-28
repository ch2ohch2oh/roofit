#include "sim_double_gaussian.cpp"

void fit_lambda()
{
    TFile * f = new TFile("exp55_run30_mva_v1.root");
    TTree * tree = (TTree*) f->Get("lambda");
    TH1F * hist1 = new TH1F("hist1", "hist1", 200, 1.105, 1.125);
    TH1F * hist2 = new TH1F("hist2", "hist2", 200, 1.105, 1.125);
    double mva = 0.999;

    tree->Draw("M >> hist1", TString::Format("mva > %f", mva), "goff");
    tree->Draw("M >> hist2", TString::Format("mva <= %f", mva), "goff");

    sim_double_gaussian(hist1, hist2, mva);
}