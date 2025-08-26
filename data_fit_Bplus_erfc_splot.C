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


// Aux: take the y-value from the drawn curve at a given mass (for vertical line heights)
double getYatMass(RooPlot* frame, double mass) {
    if (auto* curve = dynamic_cast<RooCurve*>(frame->findObject("global"))) {
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
    // fallback: scan any curve (your original)
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



// B+ Particle with Error Function
void total_data_fit_erfc_Bu() {
    gStyle->SetOptStat(0);
    using std::string;

    double min_signal = 5.178948768;
    double max_signal = 5.380091232;

    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;

    double xlow = 5.0;
    double xhigh = 5.8;

    const int nbins_plot = 80; // Number of bins for the plot
    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    const char* mcTreePath = "Bfinder/ntKp";
    const char* bLabel       = "B^{+}";
    const char* mcFile       = "/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root";
    int nbins                = 150; // Number of bins for the sPlot histogram

    // Load ROOT file and TTree
    TFile *file = TFile::Open("data_unbinned_Bu_LekaiCut.root");
    // data_unbinned_Bu_first.root
    // data_Rsideband_Bu_afterChi_FinalCutted.root
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Tree 'tree' not found in file." << std::endl;
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
      {"Bcos_dtheta",            "B_cos_dtheta",           "Bcos_dtheta"},
      {"Bchi2cl",                "B_chi2cl",               "Bchi2cl"},
      {"Btrk1Pt",                "B_trk1Pt",               "Btrk1Pt"},
    };

   // User-defined plotting ranges (only used when drawing histograms)
   std::map<std::string, std::pair<double,double>> plotRanges = {
     {"Bpt", {0.0, 50}},
     {"By", {-2.4, 2.4}},
     {"Balpha", {0.0, 3.14}},
     {"Bnorm_svpvDistance", {0.0, 85}},
     {"Bnorm_svpvDistance_2D", {0.0, 85}},
     {"Btrk1dR", {0.0, 1.75}},
     {"Bcos_dtheta", {0.995, 1.0}},
     {"Bchi2cl", {0.0, 1.0}},
     {"Btrk1Pt", {0, 10.0}},
     // ... add others as you like
   };



    // Define the mass variable and dataset
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);

    // Create unbinned RooDataSet from TTree
    B_mass.setRange("gaussRange", min_signal, max_signal);

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

    // Define regions for integrals
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);


    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", 5.278, 5.276, 5.28);


    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);


    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", 1, 0.5, 1.5);

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    RooRealVar Nsig("Nsig", "Signal Yield", 40000, 0, 96300000);
    

    // Background model
    RooRealVar lambda("lambda", "Lambda", -0.2, -2, -0.01);
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", 210000, 0, 968000000);


    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", 5.134, 5.08, 5.16);
    RooRealVar csc("csc", "Scaling Constant", 0.04, 0.001, 1);

    // integral form implemented via erf: 1 - erf(x)
    RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));


    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", 8000, 0, 50000000);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));


    // Fit the model to data (Extended Maximum Likelihood)
    RooFitResult* result = model.fitTo(dataset, Save(), Range(xlow, xhigh), Extended(kTRUE));


    // Compute background-only integrals in signal and BOTH sidebands (still exponential-only, as in your original)
    double frac_bkg_signal = expo.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion")) ->getVal();
    double frac_bkg_low    = expo.createIntegral(B_mass, NormSet(B_mass), Range("lowSideband"))  ->getVal();
    double frac_bkg_high   = expo.createIntegral(B_mass, NormSet(B_mass), Range("highSideband")) ->getVal();
    double frac_bkg_low_erfc = erfc_bkg.createIntegral(B_mass, NormSet(B_mass), Range("lowSideband")) -> getVal();

    double total_bkg_yield = Nbkg.getVal();
    double total_erfc_yield = Nerfc.getVal();
    double bkg_in_signal   = total_bkg_yield * frac_bkg_signal; // R3
    double bkg_out_signal  = total_bkg_yield * (frac_bkg_low + frac_bkg_high) + total_erfc_yield * frac_bkg_low_erfc; // R1 + R2
    double f_b             = bkg_in_signal / bkg_out_signal;


    // Compute signal yield in signal region
    double frac_sig_in_signal = signal.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal();
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal;  // Signal in signal region (S_data)

    // Open and process the MC file for signal region yield
    TFile *file_mc = TFile::Open("/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }
    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntKp", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        file_mc->Close();
        return;
    }

    // Apply the same cuts as data
    TString cut_mc = Form("Bnorm_svpvDistance>2 && Bnorm_svpvDistance_2D>4 && Bchi2cl>0.003 && (%s) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
                        ACCcuts_ppRef_Bu.Data(),
                        SELcuts_ppRef_Bu.Data(),
                        TRGmatching.Data());
    int nbins_mc = int((max_signal - min_signal) / 0.01);
    TH1F *hist_mc = new TH1F("hist_mc", "MC Bmass in Signal Region; Bmass [GeV/c^{2}]; Entries", nbins_mc, min_signal, max_signal);
    treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");
    double mc_yield_in_signal = hist_mc->Integral();  // S_MC
    delete hist_mc;
    file_mc->Close();

    double f_s = sig_yield_in_region / mc_yield_in_signal;  // f_s calculation

    // ---------- Canvas with two pads ----------
    TCanvas* cfit = new TCanvas("cfit", "Bmass Fit with Pulls (ERFC model)", 800, 800);
    cfit->Divide(1, 2);

    // ---------- Top pad (fit) ----------
    TPad* p1 = (TPad*)cfit->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    // Plot data + model with the same naming/styles used for pulls
    dataset.plotOn(frame, Binning(nbins_plot), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global")); // total
    model.plotOn(frame, Components(expo), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background")); // exponential part
    model.plotOn(frame, Components(erfc_bkg), LineColor(kMagenta), LineStyle(kDotted), LineWidth(2), Name("erfc_bkg")); // erfc part
    model.plotOn(frame, Components(signal),   LineColor(kGreen+2), LineStyle(kDashed), LineWidth(2), Name("signal")); // signal

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
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

    // Calculate chi2/ndf for the fit
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);

    // ---------- Legend (same place), on TOP pad ----------
    p1->cd();
    TLegend* legend = new TLegend(0.56, 0.66, 0.88, 0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(0);
    legend->AddEntry(frame->findObject("data"), "Data (B^{+}) Unbinned", "lep");
    legend->AddEntry(frame->findObject("background"), "Background Fit (Exponential)", "l");
    legend->AddEntry(frame->findObject("erfc_bkg"), "Background Fit (ERFC, Left Sideband)", "l");
    legend->AddEntry(frame->findObject("signal"), "Signal Fit (Double Gaussian)", "l");
    legend->AddEntry(frame->findObject("global"), "Total Fit (Signal + Background)", "l");
    legend->Draw();

    // ---------- TPaveText (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave = new TPaveText(0.64, 0.30, 0.88, 0.66, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    // Signal: Double Gaussian
    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} (fixed) = %.5f", sigma1_mc.getVal()));
    pave->AddText(Form("#sigma_{2} (fixed) = %.5f", sigma2_mc.getVal()));
    pave->AddText(Form("c_{1} (fixed) = %.4f", c1.getVal()));
    pave->AddText(Form("C_{s} = %.5f #pm %.5f", Cs.getVal(), Cs.getError()));
    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    // Exponential background
    pave->AddText(Form("#lambda = %.4f #pm %.4f", lambda.getVal(), lambda.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    // ERFC background
    pave->AddText(Form("c_{sf} = %.5f #pm %.5f", csf.getVal(), csf.getError()));
    pave->AddText(Form("c_{sc} = %.5f #pm %.5f", csc.getVal(), csc.getError()));
    pave->AddText(Form("N_{erfc} = %.1f #pm %.1f", Nerfc.getVal(), Nerfc.getError()));
    // Chi2
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();

    
    // ---------- f_b / f_s box (same place), on TOP pad ----------
    p1->cd();
    TPaveText* pave_fb_fs = new TPaveText(0.44, 0.77, 0.56, 0.88, "NDC");
    pave_fb_fs->SetTextAlign(12);
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);
    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s} = %.3f", f_s));
    pave_fb_fs->Draw();
    

    // ---------- Bottom pad (pulls) ----------
    TPad* p2 = (TPad*)cfit->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));
    RooHist* pullHist = frame->pullHist("data", "global");  // names must match Name("data") and Name("global")
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.40);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#Psi K^{+}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    // Zero line
    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");

    // Save the canvas to a file
    TString name_file = "Bu_Total_Fit_Erfc_with_Pulls.pdf";
    cfit->SaveAs(name_file);

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
sigma1_mc.setConstant(kTRUE);
sigma2_mc.setConstant(kTRUE);
c1       .setConstant(kTRUE);
// mean / lambda / Cs just got frozen by the loop above
// -------------------------------------------------------

    // ------------------ sWeights --------------------------
    RooStats::SPlot sData("sData","sData", dataset, &model, RooArgList(Nsig, Nbkg, Nerfc));
    std::cout << "sWeights built. Fitted Nsig=" << Nsig.getVal()
              << " Nbkg=" << Nbkg.getVal() << "\n"
              << " Nerfc=" << Nerfc.getVal() << "\n"; // Added Nerfc output
    double sum_wsig = 0, sum_wbkg = 0, sum_werfc = 0;
    int nEntries = dataset.numEntries();

    for (int i=0;i<nEntries;++i) {
      sum_wsig += sData.GetSWeight(i, "Nsig");
      sum_wbkg += sData.GetSWeight(i, "Nbkg");
      sum_werfc += sData.GetSWeight(i, "Nerfc");
    }
    std::cout << "Sum sW(sig) = " << sum_wsig << "  (Nsig)\n";
    std::cout << "Sum sW(bkg) = " << sum_wbkg << "  (Nbkg)\n";
    std::cout << "Sum sW(erfc) = " << sum_werfc << "  (Nerfc)\n";

    // ------------------ sPlot: Data vs MC w/ Ratio --------
    TString pdfOut = "Bu_splot_comparison_Noise_variant.pdf";
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
        const double wS = sData.GetSWeight(i, "Nbkg"); // or "Nsig" for signal sWeights
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


    // Console output summary
    std::cout << "Double Gaussian + Exponential fit complete. Output saved to " << name_file << std::endl;
    std::cout << std::fixed << std::setprecision(2) << std::endl;
    std::cout << "R3 (bkg in signal region) = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1+R2 (bkg in sidebands) = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_b = " << f_b << std::endl << std::endl;
    std::cout << "S_data (signal in region) = " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_MC = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_s = " << f_s << std::endl << std::endl;

    // Clean up
    delete line_low;
    delete line_high;
    delete zeroLine;
    delete cfit;
    //delete c_splot;
}



void data_fit_Bplus_erfc_splot() {
    total_data_fit_erfc_Bu();
}
