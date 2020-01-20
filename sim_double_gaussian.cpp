using namespace RooFit;

void sim_double_gaussian(TH1F * hist1, TH1F * hist2)
{
    int bins = 200;

    RooRealVar x("x", "x", 1.105, 1.125);

    RooRealVar mean("mean", "mean", 1.115, 1.113, 1.117);
    RooRealVar sigma_narrow("sigma_narrow", "narrow sigma", 2e-4, 0, 0.1);
    RooRealVar sigma_ratio("sigma_ratio", "ratio of sigmas", 2., 1, 5);
    RooFormulaVar sigma_wide("sigma_wide", "wide sigma", "@0 * @1", 
        RooArgList(sigma_narrow, sigma_ratio));
    RooGaussian narrow_gauss("narrow_gauss", "narrow gaussian", x, mean, sigma_narrow);
    RooGaussian wide_gauss("wide_gauss", "wide gaussian", x, mean, sigma_wide);
    RooRealVar narrow_frac("narrow_frac", "fraction of narrow gaussian", 0.8, 0.1, 1);

    // sigma_ratio.setConstant(kTRUE);
    // narrow_frac.setConstant(kTRUE);

    RooAddPdf double_gauss("double_gauss", "double gaussian", 
        RooArgList(narrow_gauss, wide_gauss), RooArgList(narrow_frac));
    auto sig_pdf = double_gauss;

    // Background one
    RooRealVar a0("a0", "a0", -0.1, -5, 5);
    RooRealVar a1("a1", "a1", -0.01, -5, 5);
    RooRealVar a2("a2", "a2", -0.001, -5, 0.);
    RooChebychev bkg1_pdf("bkg1", "bkg1", x, RooArgList(a0, a1, a2));
    
    // Background two
    RooRealVar b0("b0", "b0", -0.1, -5, 5);
    RooRealVar b1("b1", "b1", -0.01, -5, 5);
    RooRealVar b2("b2", "b2", -0.001, -5, 5);
    RooChebychev bkg2_pdf("bkg2", "bkg2", x, RooArgList(b0, b1, b2));
    // auto bkg2_pdf = bkg1_pdf;

    // combine signal and background pdf
    RooRealVar nsig1("nsig1", "nsig1", 2e5, 0, 1e9);
    RooRealVar nbkg1("nbkg1", "nbkg1", 1e5, 0, 1e9);
    RooRealVar nsig2("nsig2", "nsig2", 1e2, 0, 1e9);
    RooRealVar nbkg2("nbkg2", "nbkg2", 2e6, 0, 1e9);
    
    // signal efficiency and background rejection
    RooFormulaVar sig_eff("sig_eff", "signal efficiency", "@0 / (@0 + @1)", 
        RooArgList(nsig1, nsig2));
    RooFormulaVar bkg_rej("bkg_rej", "background rejection", "@1 / (@0 + @1)", 
        RooArgList(nbkg1, nbkg2));
    
    RooAddPdf model1("model1", "model1", RooArgList(sig_pdf, bkg1_pdf),
        RooArgList(nsig1, nbkg1));
    RooAddPdf model2("model2", "model2", RooArgList(sig_pdf, bkg2_pdf),
        RooArgList(nsig2, nbkg2));
    
    RooCategory sample("sample", "sample");
    sample.defineType("one");
    sample.defineType("two");

    // load dataset
    RooDataHist data1("data1", "data1", x, Import(*hist1));
    RooDataHist data2("data2", "data2", x, Import(*hist2));
    
    RooDataHist comb_data("comb_data", "comb_data", x, 
        Index(sample),
        Import("one", data1), Import("two", data2));
    
    RooSimultaneous sim_pdf("sim_pdf", "simultaneous pdf", sample);
    sim_pdf.addPdf(model1, "one");
    sim_pdf.addPdf(model2, "two");

    // Need to save otherwise 0 is returned
    RooFitResult * fit_result = sim_pdf.fitTo(comb_data, Save());
    cout << "===== fit_result = " << fit_result << endl;

    // plotting
    RooPlot* frame1 = x.frame(Bins(bins), Title("data one"));

    comb_data.plotOn(frame1, Cut("sample == sample::one"));

    
    sim_pdf.plotOn(frame1, Slice(sample, "one"), 
        Components("double_gauss"),
        ProjWData(sample, comb_data),
        LineColor(kRed));
    sim_pdf.plotOn(frame1, Slice(sample, "one"), 
        Components("narrow_gauss"),
        ProjWData(sample, comb_data),
        LineColor(kGreen));
    sim_pdf.plotOn(frame1, Slice(sample, "one"), 
        Components("wide_gauss"),
        ProjWData(sample, comb_data),
        LineColor(kCyan));
    sim_pdf.plotOn(frame1, Slice(sample, "one"), 
        Components("bkg1"),
        ProjWData(sample, comb_data));
    // IMPORTANT: the pull distribution in relative the last plot so 
    // we have to plot the full distribution at last.
    sim_pdf.plotOn(frame1, Slice(sample, "one"), ProjWData(sample, comb_data));

    RooPlot* frame2 = x.frame(Bins(bins), Title("data two"));

    comb_data.plotOn(frame2, Cut("sample == sample::two"));
    
    
    sim_pdf.plotOn(frame2, Slice(sample, "two"), 
        Components("double_gauss"),
        ProjWData(sample, comb_data),
        LineColor(kRed));
    sim_pdf.plotOn(frame2, Slice(sample, "two"), 
        Components("bkg2"),
        ProjWData(sample, comb_data));
    sim_pdf.plotOn(frame2, Slice(sample, "two"), ProjWData(sample, comb_data));

    RooPlot* frame3 = x.frame(Bins(bins), Title("pull one"));
    RooHist * pull1 = frame1->pullHist();
    frame3->addPlotable(pull1, "P");

    RooPlot* frame4 = x.frame(Bins(bins), Title("pull two"));
    RooHist * pull2 = frame2->pullHist();
    frame4->addPlotable(pull2, "P");

    // create the canvas
    TCanvas * can = new TCanvas("sim_fit", "sim_fit", 800, 400);
    can->Divide(2, 2);
    can->cd(1);
    frame1->Draw();
    can->cd(2);
    frame2->Draw();
    can->cd(3);
    frame3->Draw();
    can->cd(4);
    frame4->Draw();

    // print fit result
    cout << "signal efficiency: " << sig_eff.getVal() << " +- " << 
        sig_eff.getPropagatedError(*fit_result) << endl;
    cout << "background rejection: " << bkg_rej.getVal() << " +- " << 
        bkg_rej.getPropagatedError(*fit_result) << endl;
}