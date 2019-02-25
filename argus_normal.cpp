#include "RooRealVar.h"

using namespace RooFit;

void argus_normal()
{
    RooRealVar mes("mes", "m_{ES} (GeV)", 5.2, 5.3);

    RooRealVar sigmean("sigmean", "B^{#pm} mass", 5.28, 5.20, 5.30);
    RooRealVar sigwidth("sigwidth", "B^{#pm} width", 0.0027, 0.001, 1.);

    RooGaussian signal("signal", "signal PDF", mes, sigmean, sigwidth);

    // Background PDF
    RooRealVar argpar("argpar", "Argus shape parameter", -20., -100., -1.);
    RooArgusBG background("background", "Argus PDF", mes, RooConst(5.291), argpar);

    // Signal + Background
    RooRealVar nsig("nsig", "#signal events", 200, 0, 10000);
    RooRealVar nbkg("nbkg", "#background events", 800, 0, 10000);
    RooAddPdf model("model", "signal + background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    RooDataSet * data = model.generate(mes, 2000);

    model.fitTo(*data);

    RooPlot * mesframe = mes.frame();
    data -> plotOn(mesframe);
    model.plotOn(mesframe, LineColor(kRed));
    model.plotOn(mesframe, Components(background), LineStyle(kDashed));

    mesframe -> Draw();
}