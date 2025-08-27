using namespace RooFit;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPad.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TColor.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include "ACCSEL.h"

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h> 
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>  
#include <RooChebychev.h>
#include <RooProdPdf.h>
#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooCurve.h>
#include <RooExtendPdf.h>

#include <iomanip>
#include <RooStats/SPlot.h>
#include <TPaveText.h>

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

static std::pair<double,double> autoRange(TTree* t, const char* br, double padFrac=0.05) {
  // robust-ish range from tree min/max with a small padding
  double xmin = t->GetMinimum(br);
  double xmax = t->GetMaximum(br);
  if (!(xmax > xmin)) { xmin = 0; xmax = 1; }
  double span = xmax - xmin;
  if (span <= 0) span = std::max(std::abs(xmin), 1.0);
  xmin -= padFrac * span;
  xmax += padFrac * span;
  return {xmin, xmax};
}


// Aux: read y-value from the drawn curve at a given mass (for vertical line heights)
double getYatMass(RooPlot* frame, double mass) {
    for (int i = 0; i < frame->numItems(); ++i) {
        RooCurve* curve = dynamic_cast<RooCurve*>(frame->getObject(i));
        if (!curve) continue;
        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();
        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    return 0.0;
}



// Bs Particle
void total_data_fit_Bs() {
    gStyle->SetOptStat(0);
    using std::string;

    const int nbins_plot = 150; // Number of bins for the plotted data

    double mc_sigma1 = 0.02539;
    double mc_sigma2 = 0.01039;
    double mc_c1 = 0.3723;

    double min_signal = 5.314800259;
    double max_signal = 5.420059741;

    double xlow = 5.2;
    double xhigh = 5.6;
    const double bin_width_plot = (xhigh - xlow) / nbins_plot;  // Used in y-axis label

    const char* mcTreePath = "Bfinder/ntphi";
    const char* bLabel       = "B^{s}";
    const char* mcFile       = "/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bs_phat5_Bfinder.root";
    int nbins                = 150; // Number of bins for the sPlot histogram

    // Load real data tree
    TFile* file = TFile::Open("data_unbinned_Phi_FirstCut.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree not found in file." << std::endl;
        return;
    }

    TFile fmc(mcFile);
    if (fmc.IsZombie()) { Error("splot_multi_withratio","cannot open MC"); return; }
    TTree* tmc = (TTree*)fmc.Get(mcTreePath);
    if (!tmc) { Error("splot_multi_withratio","cannot get MC tree"); return; }

    // ------------------ OBSERVABLE MAP --------------------
    // display label (for axis/title), DATA branch (underscored), MC branch (original)
    struct VarMap { string label, dataBr, mcBr; };
    std::vector<VarMap> vars = {
      {"Bpt",                    "B_pt",                   "Bpt"},
      {"By",                     "B_y",                    "By"},
      {"Balpha",                 "B_alpha",                "Balpha"},
      {"Bnorm_svpvDistance",     "B_norm_svpvDistance",    "Bnorm_svpvDistance"},
      {"Bnorm_svpvDistance_2D",  "B_norm_svpvDistance_2D", "Bnorm_svpvDistance_2D"},
      {"Btrk1dR",                "B_trk1dR",               "Btrk1dR"},
      {"Btrk2dR",                "B_trk2dR",               "Btrk2dR"},
      {"Bcos_dtheta",            "B_cos_dtheta",           "Bcos_dtheta"},
      {"Bchi2cl",                "B_chi2cl",               "Bchi2cl"},
      {"Btrk1Pt",                "B_trk1Pt",               "Btrk1Pt"},
      {"Btrk2Pt",                "B_trk2Pt",               "Btrk2Pt"},
    };

   // User-defined plotting ranges (only used when drawing histograms)
   std::map<std::string, std::pair<double,double>> plotRanges = {
     {"Bpt", {0.0, 50}},
     {"By", {-2.4, 2.4}},
     {"Balpha", {0.0, 3.14}},
     {"Bnorm_svpvDistance", {0.0, 85}},
     {"Bnorm_svpvDistance_2D", {0.0, 85}},
     {"Btrk1dR", {0.0, 3.5}},
     {"Btrk2dR", {0.0, 3.5}},
     {"Bcos_dtheta", {0.995, 1.0}},
     {"Bchi2cl", {0.0, 1.0}},
     {"Btrk1Pt", {0, 10.0}},
     {"Btrk2Pt", {0, 10.0}},
     // ... add others as you like
   };


    // RooFit: variable and data
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);
    B_mass.setRange("gaussRange", min_signal, max_signal);

    RooBinning mainBins(nbins_plot, xlow, xhigh);
    B_mass.setBinning(mainBins, "mainBins");

    // ------------------ build RooDataSet (DATA) -----------
    // DATA uses underscored branch names + B_mass
    RooArgSet dataVars(B_mass);

    // store the final [xmin,xmax] chosen per LABEL so MC uses the same range
    std::map<std::string, std::pair<double,double>> ranges;

    //auto-range from DATA tree using DATA branch name
    for(const auto& v : vars){
         double xmin = 0.0, xmax = 0.0;
         auto r = autoRange(tree, v.dataBr.c_str());
         xmin = r.first;
         xmax = r.second;
    
        ranges[v.label] = {xmin, xmax};

        // RooVar name must match the DATA branch name for RooDataSet attachment
        auto* rv = new RooRealVar(v.dataBr.c_str(), v.dataBr.c_str(), xmin, xmax);
        dataVars.add(*rv);
  }


   RooDataSet dataset("dataset" , "Unbinned dataset from TTree", tree, dataVars);

    // Signal: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.36743, 5.3, 5.4);

    // MC-derived widths (FIXED constants — put your values here)
    // Common positive scale (fit parameter) and effective widths
    RooRealVar Cs("Cs", "Resolution scale", 1.0, 0.2, 3.0);

    RooRealVar sigma1("sigma1", "MC Sigma1", mc_sigma1); 
    sigma1.setConstant(kTRUE);

    RooRealVar sigma2("sigma2", "MC Sigma2", mc_sigma2); 
    sigma2.setConstant(kTRUE);

    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(mc_sigma1, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(mc_sigma2, Cs));

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);

    RooGaussian gaus1("gaus1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gaus2("gaus2", "Gaussian 2", B_mass, mean, sigma2_eff);

    // Signal shape becomes a sum of the two Gaussians
    RooAddPdf signal("signal", "Double Gaussian", RooArgList(gaus1, gaus2), RooArgList(c1));


    RooRealVar Nsig("Nsig", "Signal Yield", 3000, 0, 300000);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    // Background: Exponential model
    RooRealVar lambda("lambda", "Lambda", -0.6, -6.0, -0.01);
    RooExponential expo("expo", "Background", B_mass, lambda);

    // Background: Polynomial model
    //RooRealVar b0("b0", "b0", 1.0, -1e5, 1e5);
    //RooRealVar b1("b1", "b1", 0.0, -1e5, 1e5);
    //RooRealVar b2("b2", "b2", 0.0, -1e5, 1e5);
    //RooRealVar b3("b3", "b3", 0.0, -1e5, 1e5);      
    //RooRealVar b4("b4", "b4", 0.0, -1e5, 1e5);
    //RooPolynomial poly("poly", "4th-degree polynomial background", B_mass, RooArgList(b0, b1, b2, b3, b4));

    // Background: Exponential + Poly model
    //RooRealVar m0("m0","m0",5.5); m0.setConstant(kTRUE);
    //RooRealVar b1("b1", "b1", -2, -50, 50);
    //RooRealVar b2("b2", "b2", 0, -50, 50);
    //RooRealVar b3("b3", "b3", 0, -50, 50);
    //RooRealVar b4("b4", "b4", 0, -50, 50);
    //RooGenericPdf bkg("bkg","exp(b1*(B_mass-m0) + b2*pow(B_mass-m0,2)+ b3*pow(B_mass-m0,3)+b4*pow(B_mass-m0,4))", RooArgList(B_mass,b1,b2,b3,b4,m0));

    RooRealVar Nbkg("Nbkg", "Background Yield", 50000, 0, 999999999);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential Background", expo, Nbkg);
    //RooExtendPdf poly_ext("poly_ext", "Extended Polynomial Background", poly, Nbkg);
    //RooExtendPdf bkg_ext("bkg_ext", "Extended Background", bkg, Nbkg);

    // Combined model
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, expo_ext));
    //RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, poly_ext));
    //RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, bkg_ext));

    // Fit (Extended Maximum Likelihood Method)
    RooFitResult* result = model.fitTo(dataset, Save());

    
    // Compute background-only integrals in signal and sideband regions
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    double frac_bkg_signal = expo.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Background in Signal Region
    double frac_bkg_low    = expo.createIntegral(B_mass, NormSet(B_mass), Range("lowSideband"))->getVal(); // Background in Left Noise Region
    double frac_bkg_high   = expo.createIntegral(B_mass, NormSet(B_mass), Range("highSideband"))->getVal(); // Background in Right Noise Region

    double total_bkg_yield = Nbkg.getVal(); // Total background

    double bkg_in_signal = total_bkg_yield * frac_bkg_signal; // Ammount of Noise in Signal Region 
    double bkg_out_signal = total_bkg_yield * (frac_bkg_low + frac_bkg_high); // Ammount of Noise in Sidebands
    double f_b = bkg_in_signal / bkg_out_signal; // Calculating F_b

    double frac_sig_in_signal = signal.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Signal in Signal Region
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal; // Ammount of Signal in Signal Region



    // Opening and Checking MC File
    TFile *file_mc = TFile::Open("/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bs_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }

    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntphi", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        return;
    }

    // Apply same cuts (Change the cuts accordingly)
    TString cut_mc = Form("Bnorm_svpvDistance>2 && abs(Btktkmass - 1.019455) < 0.015 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
                        ACCcuts_ppRef.Data(),
                        SELcuts_ppRef.Data(),
                        TRGmatching.Data());

    int nbins_mc = 150;

    TH1F *hist_mc = new TH1F("hist_mc", "MC Bmass in Signal Region;Bmass [GeV/c^{2}];Entries", nbins_mc, min_signal, max_signal);

    treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");

    double mc_yield_in_signal = hist_mc->Integral();

    delete hist_mc;
    file_mc->Close();
    // --------------------------------------------------------------------------------------------------
    
    double f_s = sig_yield_in_region / mc_yield_in_signal; // Calculating F_s



    // ---------- Canvas with two pads ----------
    TCanvas* cfit = new TCanvas("cfit", "Bs Fit with Pulls", 800, 800);
    cfit->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)cfit->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(B_mass.getBinning("mainBins")), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));  // Total model
    model.plotOn(frame, Components(expo_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    //model.plotOn(frame, Components(poly_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    //model.plotOn(frame, Components(bkg_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));
    model.plotOn(frame, Components(signal), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("signal"));

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();

    // Vertical dashed lines at signal-region edges (heights taken from drawn curve)
    double y_low  = getYatMass(frame, min_signal);
    double y_high = getYatMass(frame, max_signal);

    TLine* line_low  = new TLine(min_signal, 0, min_signal, y_low);
    TLine* line_high = new TLine(max_signal, 0, max_signal, y_high);
    for (TLine* l : {line_low, line_high}) {
        l->SetLineColor(kBlack);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw("same");
    }

    // Chi2 after plotting on frame
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);


    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.58, 0.66, 0.88, 0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(0);
    legend->AddEntry(frame->findObject("data"), "Data (B_{s} )", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    //legend->AddEntry(frame->findObject("background"), "Background Fit (4degreePoly)", "l");
    //legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential + Poly)", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Signal + Background Fit", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.63, 0.38, 0.88, 0.66, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));

    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2.getVal()));
    pave->AddText(Form("c1 (fixed) = %.4f", c1.getVal()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));

    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    pave->AddText(Form("#lambda = %.5f #pm %.5f", lambda.getVal(), lambda.getError()));
    
    // Background
    /*pave->AddText(Form("b_{0} = %.4f #pm %.4f", b0.getVal(), b0.getError()));
    pave->AddText(Form("b_{1} = %.4f #pm %.4f", b1.getVal(), b1.getError()));
    pave->AddText(Form("b_{2} = %.4f #pm %.4f", b2.getVal(), b2.getError()));
    pave->AddText(Form("b_{3} = %.4f #pm %.4f", b3.getVal(), b3.getError()));
    pave->AddText(Form("b_{4} = %.4f #pm %.4f", b4.getVal(), b4.getError()));*/

    // Background: Exponential + Poly
    //pave->AddText(Form("b_{0} = %.4f #pm %.4f", b0.getVal(), b0.getError()));
    //pave->AddText(Form("b_{1} = %.4f #pm %.4f", b1.getVal(), b1.getError()));
    //pave->AddText(Form("b_{2} = %.4f #pm %.4f", b2.getVal(), b2.getError()));
    //pave->AddText(Form("b_{3} = %.4f #pm %.4f", b3.getVal(), b3.getError()));
    //pave->AddText(Form("b_{4} = %.4f #pm %.4f", b4.getVal(), b4.getError()));

    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    pave->AddText(Form("#chi^{2}/ndf = %.5f", chi2));
    pave->Draw();

    // ---------- f_b / f_s box (same place), on TOP pad ----------
    /*p1->cd();
    TPaveText* pave_fb_fs = new TPaveText(0.46, 0.77, 0.58, 0.88, "NDC");
    pave_fb_fs->SetTextAlign(12);
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);
    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s} = %.3f", f_s));
    pave_fb_fs->Draw();*/

    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)cfit->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame();
    RooHist* pullHist = frame->pullHist("data", "global");   // names must match
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.40);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#Psi #Phi} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");


    TString name = "Bs_First_Fit_splot.pdf";
    cfit->SaveAs(name);

    // ------------------- AFTER THE FIT, BEFORE sPlot -------------------
/* We want sPlot to use a *fixed* (constrained) model:
   - Freeze ALL float (non-constant) parameters from the fit,
     EXCEPT the yield parameters we will pass to SPlot (Nsig, Nbkg).
   - This leaves Nsig/Nbkg free conceptually for sWeight construction,
     but we do not refit; we just freeze the shapes.
*/

RooArgSet* allPars = model.getParameters(dataset);  // params in this pdf given 'data'
RooArgList floatPars = result->floatParsFinal(); // floating params returned by the fit

for (int i = 0; i < floatPars.getSize(); ++i) {
  auto* pFit = dynamic_cast<RooRealVar*>(floatPars.at(i));
  if (!pFit) continue;

  // map back to the *live* parameter inside 'model'
  RooAbsArg* pLiveAbs = allPars->find(pFit->GetName());
  auto* pLive = dynamic_cast<RooRealVar*>(pLiveAbs);
  if (!pLive) continue;

  // Don't freeze the yields we are going to pass to SPlot
  if (std::string(pLive->GetName()) == "Nsig") continue;
  if (std::string(pLive->GetName()) == "Nbkg") continue;
  
  // Freeze everything else (mean, lambda, resolution scale Cs, etc.)
  pLive->setConstant(kTRUE);
}

// Extra safety: these were already const, but fine to keep:
sigma1.setConstant(kTRUE);
sigma2.setConstant(kTRUE);
c1       .setConstant(kTRUE);
// mean / lambda / Cs just got frozen by the loop above
// -------------------------------------------------------

    // ------------------ sWeights --------------------------
    RooStats::SPlot sData("sData","sData", dataset, &model, RooArgList(Nsig, Nbkg));
    std::cout << "sWeights built. Fitted Nsig=" << Nsig.getVal()
              << " Nbkg=" << Nbkg.getVal() << "\n";
    double sum_wsig = 0, sum_wbkg = 0;
    int nEntries = dataset.numEntries();

    for (int i=0;i<nEntries;++i) {
      sum_wsig += sData.GetSWeight(i, "Nsig");
      sum_wbkg += sData.GetSWeight(i, "Nbkg");
    }
    std::cout << "Sum sW(sig) = " << sum_wsig << "  (Nsig)\n";
    std::cout << "Sum sW(bkg) = " << sum_wbkg << "  (Nbkg)\n";

    // ------------------ sPlot: Data vs MC w/ Ratio --------
    TString pdfOut = "Bs_splot_comparison_Signal.pdf";
    TCanvas c_splot("c_splot","sPlot comparisons",900,700);
    c_splot.Print(pdfOut+"(");

for (auto& v : vars) {
  double xmin, xmax;
  auto it = plotRanges.find(v.label);
  if (it != plotRanges.end()) {
    xmin = it->second.first;
    xmax = it->second.second;
  } else {
    // fallback: auto-range from data
    auto r = autoRange(tree, v.dataBr.c_str());
    xmin = r.first; xmax = r.second;
  }



      TH1F hDataS(("hDataS_"+v.label).c_str(),
                  (v.label+";"+v.label+";Normalized entries").c_str(),
                  nbins, xmin, xmax);
      hDataS.Sumw2();

      // fill DATA (sWeighted, using DATA branch name)
      for (int i=0;i<nEntries;++i) {
        const RooArgSet* row = dataset.get(i);
        const double xv = ((RooRealVar*)row->find(v.dataBr.c_str()))->getVal();
        const double wS = sData.GetSWeight(i, "Nsig"); // or "Nsig" for signal sWeights
        hDataS.Fill(xv, wS);
     
      }



      // fill MC (unit weight; use MC branch name)
      // Create the MC hist with the right binning
      TH1F hMCS(("hMCS_"+v.label).c_str(),
          (v.label+";"+v.label+";Normalized entries").c_str(),
          nbins, xmin, xmax);
        hMCS.Sumw2();
        // Attach it to the current directory so Project can find it by name
        hMCS.SetDirectory(gDirectory);

      // Build MC weight/cut string
      //TString mcW = "1";  // or use your MC event weight if you have one
      //TString drawCut = mcW + "*(" + cut_mc + ")";
      //treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");

      // Reset histogram before filling
      hMCS.Reset();

      // Fill the histogram using its actual name
      
     // tmc->Project(hMCS.GetName(), v.mcBr.c_str(), cut_mc);// ORIGINAL CODE
     TString varSpecificCut = cut_mc;
    // Remove the Bmass restriction for individual variable plots
    varSpecificCut.ReplaceAll(Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "");
    tmc->Project(hMCS.GetName(), v.mcBr.c_str(), varSpecificCut);




      // Debug print
      std::cout << "MC entries passing cuts for " << v.label << ": ";
      std::cout << tmc->GetEntries(cut_mc) << std::endl;



      // normalize shapes (no width-aware)
      auto safeScale = [](TH1& h){ double I=h.Integral(); if (I>0) h.Scale(1.0/I); };
      safeScale(hDataS); safeScale(hMCS);

      // style
      const Color_t kDataCol = kRed+1, kMCCol = kGreen+2;
      hDataS.SetLineColor(kDataCol); hDataS.SetMarkerColor(kDataCol);
      hDataS.SetMarkerStyle(20); hDataS.SetMarkerSize(0.9);
      hMCS.SetLineColor(kMCCol); hMCS.SetMarkerColor(kMCCol);
      hMCS.SetMarkerStyle(20); hMCS.SetMarkerSize(0.9);

// ---------------- Pull: (Data - MC) / sqrt(sigma_D^2 + sigma_M^2) ----------------
// (after you've already filled and normalized hDataS and hMCS)
TH1F hPull(("hPull_"+v.label).c_str(),
           (v.label+";"+v.label+";(Data - MC)/#sigma").c_str(),
           nbins, xmin, xmax);

// compute pulls bin-by-bin
for (int b = 1; b <= nbins; ++b) {
    const double d  = hDataS.GetBinContent(b);
    const double m  = hMCS  .GetBinContent(b);
    const double ed = hDataS.GetBinError(b);   // sPlot: sqrt(sum w^2), scaled by your safeScale
    const double em = hMCS  .GetBinError(b);   // MC statistical error, also scaled
    const double den = std::sqrt(ed*ed + em*em);

    hPull.SetBinContent(b, (den > 0.0) ? (d - m)/den : 0.0);
    //hPull.SetBinError(b, 0.0); // no error on pull itself
}

 
      // two pads
      c_splot.Clear();
      TPad *pTop = new TPad(("pTop_"+v.label).c_str(),"",0,0.30,1,1);
      TPad *pBot = new TPad(("pBot_"+v.label).c_str(),"",0,0, 1,0.30);
      pTop->SetBottomMargin(0.03);
      pBot->SetTopMargin(0.06);
      pBot->SetBottomMargin(0.32);
      pBot->SetGridy();
      pTop->Draw(); pBot->Draw();

      // top overlay
      pTop->cd();
      double ymax = std::max(hDataS.GetMaximum(), hMCS.GetMaximum());
      if(v.label == "By" || v.label == "Btrk1Dxy") {
        if(ymax*1.1>1.0){
          hMCS.SetMaximum(1.0);
          hMCS.SetMinimum(0);
       }else{
          hMCS.SetMaximum(ymax*1.1); 
          hMCS.SetMinimum(01);
    }}else{ 
        if(ymax*1.1>1.0){
          hMCS.SetMaximum(1.0);
          hMCS.SetMinimum(0.0);
       }else{
          hMCS.SetMaximum(ymax*1.1); 
          hMCS.SetMinimum(0.0);
      }}
      hMCS.GetYaxis()->SetTitle("Normalized entries");
      hMCS.GetXaxis()->SetLabelSize(0);
      hMCS.Draw("E1"); hDataS.Draw("E1 SAME");
      TLegend leg(0.70,0.70,0.90,0.90);
      leg.SetBorderSize(0); leg.SetFillStyle(0);
      leg.AddEntry(&hMCS,   "Monte Carlo", "lep");
      leg.AddEntry(&hDataS, "sPlot",       "lep");
      leg.Draw();
      TLatex lab; lab.SetNDC(true); lab.SetTextFont(62); lab.SetTextSize(0.06);lab.SetTextColor(kBlack);
      lab.DrawLatex(0.16, 0.8, bLabel);
/*
       TH1F* hRatio = (TH1F*)hMCS->Clone("hRatio");
      hRatio->Add(hDataS, hMCS,1.0, 1.0);
*/


pBot->cd();
hPull.SetTitle("");
hPull.SetMarkerStyle(20);
hPull.SetMarkerSize(0.8);
hPull.SetLineColor(kRed+1);
hPull.SetMarkerColor(kRed+1);
hPull.GetYaxis()->SetTitle("(Data-MC)/#sigma");
hPull.GetYaxis()->CenterTitle(true);
hPull.GetYaxis()->SetNdivisions(505);
hPull.GetYaxis()->SetTitleSize(0.12);
hPull.GetYaxis()->SetTitleOffset(0.40);
hPull.GetYaxis()->SetLabelSize(0.10);
hPull.GetXaxis()->SetTitle(v.label.c_str());
hPull.GetXaxis()->SetTitleSize(0.14);
hPull.GetXaxis()->SetLabelSize(0.12);

// choose a sensible range
hPull.SetMinimum(-3.5);
hPull.SetMaximum( 3.5);

hPull.Draw("E1");

// guide lines at 0, ±1, ±2
const double gxmin = hPull.GetXaxis()->GetXmin();
const double gxmax = hPull.GetXaxis()->GetXmax();
TLine l0(gxmin, 0.0, gxmax, 0.0);
TLine l1p(gxmin, 1.0, gxmax, 1.0), l1m(gxmin,-1.0, gxmax,-1.0);
TLine l2p(gxmin, 2.0, gxmax, 2.0), l2m(gxmin,-2.0, gxmax,-2.0);
l0.SetLineStyle(1);
for (TLine* L : std::vector<TLine*>{&l1p,&l1m,&l2p,&l2m}) { L->SetLineStyle(2); }
l0.SetLineColor(kGray+2);
l1p.SetLineColor(kGray+1); l1m.SetLineColor(kGray+1);
l2p.SetLineColor(kGray+1); l2m.SetLineColor(kGray+1);
l0.Draw("SAME"); l1p.Draw("SAME"); l1m.Draw("SAME"); l2p.Draw("SAME"); l2m.Draw("SAME");

      // save
      c_splot.cd(); c_splot.Print(pdfOut);
      c_splot.SaveAs((v.label+"_splot_vs_mc_ratio.png").c_str());
    }
    c_splot.Print(pdfOut+")");

    // f_b, f_s 

    std::cout << "sPlot comparison saved to " << pdfOut << std::endl;


    std::cout << std::fixed << std::setprecision(2);
    std::cout << " " << std::endl;
    std::cout << "R3 = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1 + R2 = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_{b} = " << f_b << std::endl;

    std::cout << " " << std::endl;

    std::cout << "S_{data} " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_{MC} = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_{s} = " << f_s << std::endl;
    std::cout << " " << std::endl; 

    std::cout << "Output saved to " << name << std::endl;

    delete cfit;
    delete tree;
    delete line_low;
    delete line_high;
}



void data_fit_Bs_splot() {
    total_data_fit_Bs();
}