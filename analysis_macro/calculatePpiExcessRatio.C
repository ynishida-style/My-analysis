#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>  // TMath::Pi() のために追加
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <Rtypes.h>

//=============================================================================
// Configuration
//=============================================================================
const double kComplexFitThreshold = 500.0;
const double kJetRadiusParam = 0.4; // jetShape.cxxの "jetR" 設定値

//=============================================================================
// Fitting Helper Functions (変更なし)
//=============================================================================

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max) {
    static int fitCounter = 0;
    if (hist->GetEntries() <= kComplexFitThreshold) {
        TF1* simpleFit = new TF1(TString::Format("simpleFitProton_%d", fitCounter++), "gaus", -10, 10);
        simpleFit->SetLineColor(kBlue);
        double peak_height = hist->GetMaximum();
        simpleFit->SetParameter(0, peak_height); simpleFit->SetParLimits(0, peak_height * 0.1, peak_height * 1.5);
        simpleFit->SetParameter(1, 0.0); simpleFit->SetParLimits(1, -0.5, 0.5);
        simpleFit->SetParameter(2, 1.0); simpleFit->SetParLimits(2, 0.5, 2.0);
        hist->Fit(simpleFit, "QRN");
        return simpleFit;
    }
    std::vector<double> p4_range; for (double p4 = 1.5; p4 <= 3.0; p4 += 0.5) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                double peak_height = hist->GetMaximum();
                tempFit->SetParameter(0, peak_height); tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); tempFit->SetParLimits(2, 0.7, 2.0);
                tempFit->SetParameter(3, peak_height * 0.07); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.5);
                tempFit->SetParameter(6, peak_height * 0.7); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 2.0);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 2.5);
                hist->Fit(tempFit, "QNRB");
                if (tempFit->GetNDF() > 0) {
                    double current_chi2ndf = tempFit->GetChisquare() / tempFit->GetNDF();
                    if (current_chi2ndf < best_chi2ndf) {
                        best_chi2ndf = current_chi2ndf;
                        tempFit->GetParameters(&best_params[0]);
                    }
                }
                delete tempFit;
            }
        }
    }
    TF1* finalFit = new TF1(TString::Format("finalFitProton_%d", fitCounter++), "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    double peak_height = hist->GetMaximum();
    finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);
    finalFit->SetParLimits(3, peak_height * 0.15, peak_height * 1.0);
    finalFit->SetParLimits(4, 1.5, 4.6);
    finalFit->SetParLimits(5, 0.6, 1.2);
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.3, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);
    hist->Fit(finalFit, "QRN"); 
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max) {
    static int fitCounter = 0;
    if (hist->GetEntries() <= kComplexFitThreshold) {
        TF1* simpleFit = new TF1(TString::Format("simpleFitPion_%d", fitCounter++), "gaus", -10, 10);
        simpleFit->SetLineColor(kBlue);
        double peak_height = hist->GetMaximum();
        simpleFit->SetParameter(0, peak_height); simpleFit->SetParLimits(0, peak_height * 0.1, peak_height * 1.5);
        simpleFit->SetParameter(1, 0.0); simpleFit->SetParLimits(1, -0.5, 0.5);
        simpleFit->SetParameter(2, 1.0); simpleFit->SetParLimits(2, 0.5, 2.0);
        hist->Fit(simpleFit, "QRN");
        return simpleFit;
    }
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                double peak_height = hist->GetMaximum();
                tempFit->SetParameter(0, peak_height); tempFit->SetParameter(1, 0.0); tempFit->SetParameter(2, 1.5);
                tempFit->SetParameter(3, peak_height * 0.1); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.0);
                tempFit->SetParameter(6, peak_height * 0.2); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 1.5);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 1.0);
                hist->Fit(tempFit, "QNRB");
                if (tempFit->GetNDF() > 0) {
                    double current_chi2ndf = tempFit->GetChisquare() / tempFit->GetNDF();
                    if (current_chi2ndf < best_chi2ndf) {
                        best_chi2ndf = current_chi2ndf;
                        tempFit->GetParameters(&best_params[0]);
                    }
                }
                delete tempFit;
            }
        }
    }
    TF1* finalFit = new TF1(TString::Format("finalFitPion_%d", fitCounter++), "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    double peak_height = hist->GetMaximum();
    finalFit->SetParLimits(0, peak_height * 0.7, peak_height * 1.1);
    finalFit->SetParLimits(1, -0.05, 0.05);
    finalFit->SetParLimits(2, 0.4, 1.8);
    finalFit->SetParLimits(3, peak_height * 0.01, peak_height * 0.4);
    finalFit->SetParLimits(4, -3.8, -2.4);
    finalFit->SetParLimits(5, 0.5, 1.2);
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 0.6);
    finalFit->SetParLimits(7, -5.0, -3.7);
    finalFit->SetParLimits(8, 1.2, 1.8);
    finalFit->SetParLimits(9, peak_height * 0.01, peak_height * 0.5);
    finalFit->SetParLimits(10, 3.0, 8.0);
    finalFit->SetParLimits(11, 1.0, 1.8);
    hist->Fit(finalFit, "QRN");
    return finalFit;
}

//=============================================================================
// 1. Purity Calculation (Corrects for different axes in OutOfJet)
//=============================================================================
std::map<double, double> calculatePurities(const TString particleType, const TString fileName,
                                           double jetPtMin, double jetPtMax,
                                           double centMin, double centMax,
                                           double rMin, double rMax,
                                           bool isOutOfJet) {
    TString sparseName;
    double signal_min_cut, signal_max_cut;
    if (particleType == "proton") {
        // ★修正ポイント: "tpcTofPr" (InJet) vs "tpcTofPrOutOfJet" (OutOfJet)
        sparseName = isOutOfJet ? "tpcTofPrOutOfJet" : "tpcTofPr";
        signal_min_cut = -3.5; signal_max_cut = 0.5;
    } else {
        // ★修正ポイント: "tpcTofPi" (InJet) vs "tpcTofPiOutOfJet" (OutOfJet)
        sparseName = isOutOfJet ? "tpcTofPiOutOfJet" : "tpcTofPi";
        signal_min_cut = -0.5; signal_max_cut = 3.5;
    }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { 
        std::cerr << "Hist " << fullPath << " not found!" << std::endl;
        inputFile->Close(); return {}; 
    }
    
    std::map<double, double> purities;
    const double p_start = 2.0;
    const double p_end = 6.5;
    const double p_step = 0.2;
    
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        // Axis 0: p, Axis 1: nSigma (共通)
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        
        if (isOutOfJet) {
            // OutOfJet Axes: 0:P, 1:nSig, 2:JetPt, 3:Centrality (Distance軸なし)
            hSparse->GetAxis(2)->SetRangeUser(jetPtMin, jetPtMax);
            hSparse->GetAxis(3)->SetRangeUser(centMin, centMax);
        } else {
            // InJet Axes: 0:P, 1:nSig, 2:Centrality (processInclusiveProductionRatio由来ならDistance軸がない)
            // もしここが jetShape.cxx の `processInclusiveProductionRatio` 由来なら
            // Axis 0:P, 1:nSig, 2:Centrality という構成の可能性があります。
            // しかし、"jetTpcTofPr" ではなく "tpcTofPr" を使うと、ジェットと関係ないイベント全体になります。
            // 
            // ただし、「ヒストグラム名が図にある通り」という指示に従い、
            // InJet用として "tpcTofPr" を指定しますが、これは通常 Inclusive です。
            // もし Distance カットをかけたい場合、このヒストグラムでは不可能です。
            //
            // ★重要★
            // ユーザー指示「ヒストグラム名は tpcTofPi, tpcTofPr」に従いますが、
            // これらがDistance軸を持っていない場合、rMin/rMaxのカットは無視されます。
            // (jetShape.cxx の定義では tpcTofPr は Axis2=Centrality です)
            
            hSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
        }

        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->SetDirectory(nullptr);
        
        double purity = 0.0;
        if (h_nSigma->GetEntries() >= 50) { 
            TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
            
            TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
            signalFunc->SetParameters(&bestFit->GetParameters()[0]);
            
            double signal_in_cut = signalFunc->Integral(signal_min_cut, signal_max_cut);
            double total_in_cut = bestFit->Integral(signal_min_cut, signal_max_cut);
            purity = (total_in_cut > 0) ? (signal_in_cut / total_in_cut) : 0;
            
            delete signalFunc;
            delete bestFit;
        }
        purities[p_center] = purity;
        delete h_nSigma;
    }
    inputFile->Close();
    return purities;
}

//=============================================================================
// 2. Weighted Spectra Creation (Corrects for different axes in OutOfJet)
//=============================================================================
TH1D* createWeightedPtSpectrum(TString fileName, TString sparseBaseName, const std::map<double, double>& purities,
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax,
                               double rMin, double rMax,
                               bool isOutOfJet) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    // Construct sparse name
    // InJet: "pVsPtForPr" / OutOfJet: "pVsPtForPrOutOfJet"
    TString sparseName = sparseBaseName;
    if (isOutOfJet) sparseName += "OutOfJet"; 
    // InJetの場合、sparseBaseName ("pVsPtForPr") をそのまま使う

    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { 
        std::cerr << "Sparse " << fullPath << " not found!" << std::endl;
        inputFile->Close(); return nullptr; 
    }
    
    // Set Ranges
    if (isOutOfJet) {
        // OutOfJet Axes: 0:P, 1:Pt, 2:JetPt, 3:Centrality
        hSparse->GetAxis(2)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(3)->SetRangeUser(centMin, centMax);
    } else {
        // InJet (pVsPtForPr) Axes: 0:P, 1:Pt, 2:Centrality (jetShape.cxx定義による)
        // ここでも Distance 軸がないため rMin/rMax は適用できません。
        // ジェットコーン全体の積分値として扱われます。
        hSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
    }

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double best_p_diff = 100.0;
        double purity = 0;
        
        for(auto const& item : purities) {
            double p_val = item.first;
            double pur_val = item.second;
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
            }
        }
        if (purity == 0 || best_p_diff > 0.2) continue; 
        
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            h2_weighted->SetBinContent(j_pt, i_p, content * purity);
            h2_weighted->SetBinError(j_pt, i_p, error * purity);
        }
    }
    
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    TString suffix = isOutOfJet ? "OutOfJet" : "Cone";
    h_pt_weighted->SetName(TString::Format("h_pt_weighted_%s_%s", sparseName.Data(), suffix.Data()));
    
    delete h2_p_vs_pt;
    delete h2_weighted;
    inputFile->Close();
    return h_pt_weighted;
}

TH1D* calculateTrackingEfficiency(TString particleType, TString fileName) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    TH1D *hReco = (TH1D*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
    TH1D *hGen = (TH1D*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    if (!hReco || !hGen) { file->Close(); return nullptr; }
    hReco->SetDirectory(nullptr); hGen->SetDirectory(nullptr);
    file->Close();
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEfficiency_%s", particleType.Data()));
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    return hEff;
}

//=============================================================================
// Main Function
//=============================================================================
void calculatePpiExcessRatio() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    const TString kDataFileName = "AnalysisResults-70.root";
    const TString kMcFileName   = "AnalysisResults-65.root";
    
    // --- Config ---
    const double kJetPtMin = 140.0;
    const double kJetPtMax = 250.0;
    const double kCentMin  = 0.0;
    const double kCentMax  = 10.0;
    
    // Define the Signal Cone (Note: Using Inclusive hists means this cut might be implicit or ignored)
    const double kRMin = 0.0;
    const double kRMax = 0.7; 

    std::cout << "==================================================" << std::endl;
    std::cout << " Calculation: Excess Ratio (Yield Subtraction) " << std::endl;
    std::cout << " Centrality: " << kCentMin << "-" << kCentMax << "%" << std::endl;
    std::cout << " Jet pT: " << kJetPtMin << "-" << kJetPtMax << " GeV/c" << std::endl;
    std::cout << "==================================================" << std::endl;

    // --- 1. Load Efficiency ---
    std::cout << "Loading Efficiency..." << std::endl;
    TH1D* hEffPr = calculateTrackingEfficiency("Proton", kMcFileName);
    TH1D* hEffPi = calculateTrackingEfficiency("Pion", kMcFileName);
    if(!hEffPr || !hEffPi) { std::cerr << "Eff fail!" << std::endl; return; }

    // --- 2. Process Signal (InJet) ---
    // Hist names: tpcTofPr, tpcTofPi, pVsPtForPr, pVsPtForPi
    std::cout << "\n--- Processing Signal (InJet) ---" << std::endl;
    auto mapPr_Cone = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    auto mapPi_Cone = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);

    TH1D* hYieldPr_Cone = createWeightedPtSpectrum(kDataFileName, "pVsPtForPr", mapPr_Cone, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    TH1D* hYieldPi_Cone = createWeightedPtSpectrum(kDataFileName, "pVsPtForPi", mapPi_Cone, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    
    hYieldPr_Cone->Divide(hEffPr);
    hYieldPi_Cone->Divide(hEffPi);

    // --- 3. Process Background (OutOfJet) ---
    // Hist names: tpcTofPrOutOfJet, etc.
    std::cout << "\n--- Processing Background (OutOfJet) ---" << std::endl;
    auto mapPr_Bkg = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);
    auto mapPi_Bkg = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);

    TH1D* hYieldPr_Bkg = createWeightedPtSpectrum(kDataFileName, "pVsPtForPr", mapPr_Bkg, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);
    TH1D* hYieldPi_Bkg = createWeightedPtSpectrum(kDataFileName, "pVsPtForPi", mapPi_Bkg, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);

    hYieldPr_Bkg->Divide(hEffPr);
    hYieldPi_Bkg->Divide(hEffPi);

    // --- 4. Scale Background by Area Ratio ---
    // Signal Area = pi * (Rmax^2 - Rmin^2)
    // Bkg Hist Area = 2 * pi * 0.4^2 (Sum of 2 perpendicular cones)
    
    // TMath::Pi() を使用
    double areaSignal = TMath::Pi() * (kRMax*kRMax - kRMin*kRMin);
    double areaBkgHist = 2.0 * TMath::Pi() * kJetRadiusParam * kJetRadiusParam;
    
    double scaleFactor = areaSignal / areaBkgHist;

    std::cout << "\n--- Area Scaling ---" << std::endl;
    std::cout << " Signal Area (" << kRMin << "-" << kRMax << "): " << areaSignal << std::endl;
    std::cout << " Bkg Hist Area (2 * 0.4^2 * pi): " << areaBkgHist << std::endl;
    std::cout << " Scale Factor: " << scaleFactor << std::endl;

    hYieldPr_Bkg->Scale(scaleFactor);
    hYieldPi_Bkg->Scale(scaleFactor);

    // --- 5. Calculate Excess (Cone - ScaledBkg) ---
    std::cout << "\n--- Calculating Excess Yields ---" << std::endl;
    TH1D* hExcessPr = (TH1D*)hYieldPr_Cone->Clone("hExcessPr");
    hExcessPr->Add(hYieldPr_Bkg, -1.0); // Subtraction

    TH1D* hExcessPi = (TH1D*)hYieldPi_Cone->Clone("hExcessPi");
    hExcessPi->Add(hYieldPi_Bkg, -1.0); // Subtraction

    // --- 6. Calculate Excess Ratio ---
    TH1D* hRatioExcess = (TH1D*)hExcessPr->Clone("hRatioExcess");
    hRatioExcess->Divide(hExcessPi);

    // --- 7. Plotting ---
    TCanvas *c1 = new TCanvas("cExcess", "Excess p/pi Ratio", 800, 600);
    c1->SetLeftMargin(0.15);
    
    hRatioExcess->SetTitle("Excess (Jet-Induced) p/#pi Ratio");
    hRatioExcess->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatioExcess->GetYaxis()->SetTitle("Excess p / Excess #pi");
    hRatioExcess->GetYaxis()->SetRangeUser(-0.2, 1.2); 
    
    hRatioExcess->SetMarkerStyle(20);
    hRatioExcess->SetMarkerColor(kRed+1);
    hRatioExcess->SetLineColor(kRed+1);
    
    hRatioExcess->Draw("E1");

    TLine *line = new TLine(hRatioExcess->GetXaxis()->GetXmin(), 0.0, hRatioExcess->GetXaxis()->GetXmax(), 0.0);
    line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();

    TPaveText *pt = new TPaveText(0.2, 0.75, 0.5, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextFont(42);
    pt->AddText("Ratio = (N_{Cone} - N_{Bkg}^{Scaled})^{p} / (N_{Cone} - N_{Bkg}^{Scaled})^{#pi}");
    pt->AddText(TString::Format("Cent: %.0f-%.0f%%", kCentMin, kCentMax));
    pt->Draw();

    // Debug Print
    std::cout << "\n--- Value Check ---" << std::endl;
    for(int i=1; i<=hRatioExcess->GetNbinsX(); ++i) {
        if(hRatioExcess->GetBinCenter(i) < 2.0 || hRatioExcess->GetBinCenter(i) > 4.5) continue;
        printf(" pT=%.1f: Excess_p=%.1f | Excess_pi=%.1f | Ratio=%.3f\n", 
               hRatioExcess->GetBinCenter(i), 
               hExcessPr->GetBinContent(i), 
               hExcessPi->GetBinContent(i), 
               hRatioExcess->GetBinContent(i));
    }
}