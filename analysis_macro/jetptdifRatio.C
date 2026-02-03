#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TLine.h>
#include <TMath.h>
#include <TError.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility>
#include <iomanip>

// ============================================================================================
// Constants & Settings
// ============================================================================================
const double kJetPtMax = 250.0;
const double kRMin = 0.0;
const double kRMax = 0.1; // R < 0.1 (Core)

// ★ Jet Radius (解析タスクの設定値。UEの面積計算に使用)
const double kJetRadius = 0.4; 

// ★ ビンをまとめる数
const int kRebinFactor = 2; 

// ============================================================================================
// Helper Functions: Fitting (Standard Implementation)
// ============================================================================================

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    int oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning;

    std::vector<double> p4_range; 
    for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 

    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    bool isHighPt = (p_min > 4.1);

    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;

                tempFit->SetParameter(0, pr_init_height);
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                tempFit->SetParameter(3, kaon_init_height); 
                tempFit->SetParameter(4, p4_init); 
                tempFit->SetParameter(5, 1.5);

                tempFit->SetParameter(6, peak_height * 0.7); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 2.0);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 2.5);

                hist->Fit(tempFit, "QNRB0");

                if (tempFit->GetNDF() > 0) {
                    double current_chi2ndf = tempFit->GetChisquare() / tempFit->GetNDF();
                    if (current_chi2ndf < best_chi2ndf) {
                        best_chi2ndf = current_chi2ndf;
                        tempFit->GetParameters(&best_params[0]);
                    }
                }
            }
        }
    }
    delete tempFit;
    
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88);
    } else {
        finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    }
    
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);

    if (isHighPt) {
        finalFit->SetParLimits(3, peak_height * 0.20, peak_height * 1.2); 
    } else {
        finalFit->SetParLimits(3, peak_height * 0.05, peak_height * 1.0);
    }

    finalFit->SetParLimits(4, 0.4, 1.2); 
    finalFit->SetParLimits(5, 0.6, 1.2);

    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.5, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);

    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);

    hist->Fit(finalFit, "RBQ0");

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    int oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning;
    
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    TF1* tempFit = new TF1("tempFitPionSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                tempFit->SetParameter(0, peak_height); tempFit->SetParameter(1, 0.0); tempFit->SetParameter(2, 1.5);
                tempFit->SetParameter(3, peak_height * 0.1); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.0);
                tempFit->SetParameter(6, peak_height * 0.2); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 1.5);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 1.0);
                
                hist->Fit(tempFit, "QNRB0");

                if (tempFit->GetNDF() > 0) {
                    double current_chi2ndf = tempFit->GetChisquare() / tempFit->GetNDF();
                    if (current_chi2ndf < best_chi2ndf) {
                        best_chi2ndf = current_chi2ndf;
                        tempFit->GetParameters(&best_params[0]);
                    }
                }
            }
        }
    }
    delete tempFit;

    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed); 
    finalFit->SetLineWidth(2);

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

    hist->Fit(finalFit, "RBQ0");

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

// ============================================================================================
// Core Logic: Smoothing Purity & PID Efficiency
// ============================================================================================
std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(
    TString fileName, TString particleType, double centMin, double centMax, 
    double jetPtMin, double jetPtMax, double rMin, double rMax) 
{
    std::map<double, double> finalPurityMap;
    std::map<double, double> finalEffMap;
    TString sparseName = (particleType == "proton") ? "jetTpcTofPr" : "jetTpcTofPi"; 
    double sigMin = (particleType == "proton") ? -3.5 : -0.5;
    double sigMax = (particleType == "proton") ? 0.5 : 3.5;

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) return {finalPurityMap, finalEffMap};
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { f->Close(); return {finalPurityMap, finalEffMap}; }

    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);         
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); 
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   

    double p_start = 0.6; double p_end = 8.0; double p_step = 0.2; 
    TGraph *gPurity = new TGraph();
    TGraph *gEff    = new TGraph();
    int pointIdx = 0;

    // std::cout << "--- " << particleType << " (Jet " << jetPtMin << "-" << jetPtMax << ", R " << rMin << "-" << rMax << ") ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        if (p_center <= 3.0) {
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            finalPurityMap[p_center] = 1.0;
            finalEffMap[p_center]    = idealInCut / idealTotal;
            continue;
        }

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        if (h_nSigma->GetEntries() < 50) { delete h_nSigma; continue; }

        TF1* fit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        TF1* fSig = new TF1("fSig", "gaus", -10, 10);
        fSig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));

        double intSigInCut = fSig->Integral(sigMin, sigMax); 
        double intSigTotal = fSig->Integral(-10, 10);        
        double intTotInCut = fit->Integral(sigMin, sigMax); 
        double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
        double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;
        
        if (purity > 1.0) purity = 1.0; if (purity < 0.0) purity = 0.0;
        if (pidEff > 1.0) pidEff = 1.0; if (pidEff < 0.0) pidEff = 0.0;

        gPurity->SetPoint(pointIdx, p_center, purity);
        gEff->SetPoint(pointIdx, p_center, pidEff);
        pointIdx++;
        delete h_nSigma; delete fit; delete fSig;
    }

    if (pointIdx > 2) { 
        TF1 *fPurFit = new TF1("fPurFit", "pol3", 3.0, p_end);
        TF1 *fEffFit = new TF1("fEffFit", "pol3", 3.0, p_end);
        
        gPurity->Fit(fPurFit, "Q R N"); gEff->Fit(fEffFit, "Q R N");

        for (double p_curr = 3.1; p_curr < p_end; p_curr += p_step) {
             double smoothPur = fPurFit->Eval(p_curr);
             double smoothEff = fEffFit->Eval(p_curr);
             if (smoothPur > 1.0) smoothPur = 1.0; if (smoothPur < 0.0) smoothPur = 0.0;
             if (smoothEff > 1.0) smoothEff = 1.0; if (smoothEff < 0.0) smoothEff = 0.0;
             finalPurityMap[p_curr] = smoothPur;
             finalEffMap[p_curr]    = smoothEff;
        }
        delete fPurFit; delete fEffFit;
    } else {
        std::cout << "Warning: Not enough points for smoothing (" << particleType << ")." << std::endl;
        for(int i=0; i<pointIdx; ++i) {
            double x, y; gPurity->GetPoint(i, x, y);
            finalPurityMap[x] = y; gEff->GetPoint(i, x, y); finalEffMap[x] = y;
        }
    }
    delete gPurity; delete gEff;
    f->Close();
    return {finalPurityMap, finalEffMap};
}

// ============================================================================================
// Signal Spectrum Creation (Inside Jet Cone)
// ============================================================================================
TH1D* createCorrectedSpectrum(TString fileName, TString histName, 
                              std::map<double, double>& purityMap, 
                              std::map<double, double>& pidEffMap, 
                              double centMin, double centMax,
                              double jetPtMin, double jetPtMax,
                              double rMin, double rMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", histName.Data()));
    if (!hSparse) {
        std::cout << "Error: Histogram " << histName << " not found." << std::endl;
        if(f) f->Close(); return nullptr;
    }

    // Axis 0: P, 1: Pt, 2: Distance, 3: JetPt, 4: Centrality
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);         
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); 
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   

    TH2D *h2 = (TH2D*)hSparse->Projection(0, 1); // Y=Axis0(P), X=Axis1(Pt)
    
    TH1D *hPt = h2->ProjectionX(TString::Format("hPt_%s_j%.0f", histName.Data(), jetPtMin));
    hPt->Reset();
    hPt->SetDirectory(0);

    for (int by = 1; by <= h2->GetNbinsY(); ++by) {
        double p_val = h2->GetYaxis()->GetBinCenter(by);
        if (p_val > 8.0) continue;

        double purity = 1.0; double pidEff = 1.0; double min_diff = 1.0; 
        for (auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_val) < min_diff) {
                min_diff = std::abs(map_p - p_val);
                purity = map_pur; pidEff = pidEffMap[map_p]; 
            }
        }
        if (min_diff > 0.3 && p_val > 3.0) continue; 
        if (pidEff < 0.01) pidEff = 1.0; 

        double totalCorrectionFactor = purity / pidEff;
        
        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double content = h2->GetBinContent(bx, by);
            double error   = h2->GetBinError(bx, by);
            if (content <= 0) continue;
            double corrContent = content * totalCorrectionFactor;
            double corrError   = error * totalCorrectionFactor;
            hPt->SetBinContent(bx, hPt->GetBinContent(bx) + corrContent);
            hPt->SetBinError(bx, std::sqrt(pow(hPt->GetBinError(bx), 2) + pow(corrError, 2)));
        }
    }
    
    if (kRebinFactor > 1) {
        hPt->Rebin(kRebinFactor);
    }

    delete h2;
    f->Close();
    return hPt;
}

// ============================================================================================
// ★ Background Spectrum Creation (Perp Cone / OutOfJet)
// ============================================================================================
TH1D* createCorrectedBackgroundSpectrum(TString fileName, TString histName, 
                              std::map<double, double>& purityMap, 
                              std::map<double, double>& pidEffMap, 
                              double centMin, double centMax,
                              double jetPtMin, double jetPtMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", histName.Data()));
    if (!hSparse) {
        std::cout << "Error: BG Histogram " << histName << " not found." << std::endl;
        if(f) f->Close(); return nullptr;
    }

    // OutOfJet Structure:
    // Axis 0: P, 1: Pt, 2: JetPt, 3: Centrality
    hSparse->GetAxis(2)->SetRangeUser(jetPtMin, jetPtMax); 
    hSparse->GetAxis(3)->SetRangeUser(centMin, centMax);   

    TH2D *h2 = (TH2D*)hSparse->Projection(0, 1); // Y=P, X=Pt
    
    TH1D *hPt = h2->ProjectionX(TString::Format("hPtBg_%s_j%.0f", histName.Data(), jetPtMin));
    hPt->Reset();
    hPt->SetDirectory(0);

    for (int by = 1; by <= h2->GetNbinsY(); ++by) {
        double p_val = h2->GetYaxis()->GetBinCenter(by);
        if (p_val > 8.0) continue;

        double purity = 1.0; double pidEff = 1.0; double min_diff = 1.0; 
        for (auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_val) < min_diff) {
                min_diff = std::abs(map_p - p_val);
                purity = map_pur; pidEff = pidEffMap[map_p]; 
            }
        }
        if (min_diff > 0.3 && p_val > 3.0) continue; 
        if (pidEff < 0.01) pidEff = 1.0; 

        double totalCorrectionFactor = purity / pidEff;

        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double content = h2->GetBinContent(bx, by);
            double error   = h2->GetBinError(bx, by);
            if (content <= 0) continue;
            double corrContent = content * totalCorrectionFactor;
            double corrError   = error * totalCorrectionFactor;
            hPt->SetBinContent(bx, hPt->GetBinContent(bx) + corrContent);
            hPt->SetBinError(bx, std::sqrt(pow(hPt->GetBinError(bx), 2) + pow(corrError, 2)));
        }
    }
    
    if (kRebinFactor > 1) {
        hPt->Rebin(kRebinFactor);
    }

    delete h2;
    f->Close();
    return hPt;
}

// ============================================================================================
// Tracking Efficiency Calculation (Reco/Gen)
// ============================================================================================
TH1D* calculateEfficiency(TString fileName, TString partName, double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    
    THnSparseD *hRecoSp = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptHistogram%sTof", partName.Data()));
    THnSparseD *hGenSp  = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptGenerated%s", partName.Data()));
    
    if (!hRecoSp || !hGenSp) { if(f) f->Close(); return nullptr; }

    // ★ MCのCentralityは常に0-100% (ユーザー指定)
    hRecoSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    
    TH1D *hReco = hRecoSp->Projection(0);
    TH1D *hGen  = hGenSp->Projection(0);
    
    if (kRebinFactor > 1) {
        hReco->Rebin(kRebinFactor);
        hGen->Rebin(kRebinFactor);
    }

    TH1D *hEff = (TH1D*)hReco->Clone(TString::Format("eff_%s_%.0f_%.0f", partName.Data(), centMin, centMax));
    hEff->Divide(hReco, hGen, 1, 1, "B");
    hEff->SetDirectory(0);
    f->Close();
    return hEff;
}

// ============================================================================================
// Main Function: Jet Analysis with UE Subtraction
// ============================================================================================

void PlotJetPtDependentRatio() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    const TString dataFile = "AnalysisResults-106.root"; 
    const TString mcFile   = "AnalysisResults-91.root";

    const double dataCentMin = 0.0;
    const double dataCentMax = 10.0;

    // ★ Jet Pt Bins
    std::vector<std::pair<double, double>> jetPtBins = {
        {20.0, kJetPtMax},
        {40.0, kJetPtMax},
        {60.0, kJetPtMax},
        {80.0, kJetPtMax},
        {100.0, kJetPtMax}
    };

    std::vector<int> colors  = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1};
    std::vector<int> markers = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond};

    TCanvas *c1 = new TCanvas("cRatio", "Proton/Pion Ratio vs Jet Pt (UE Subtracted)", 900, 700);
    c1->cd();
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetGridy();

    TLegend *leg = new TLegend(0.55, 0.60, 0.88, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);

    // =========================================================
    // STEP 0: Calculate Efficiency 
    // ★ 変更点: MCのCentralityは強制的に 0-100% を使用
    // =========================================================
    std::cout << ">>> Calculating Efficiency (Using MC Centrality 0-100%) <<<" << std::endl;
    TH1D* hEffPrInclusive = calculateEfficiency(mcFile, "Proton", 0.0, 100.0);
    TH1D* hEffPiInclusive = calculateEfficiency(mcFile, "Pion", 0.0, 100.0);

    if (!hEffPrInclusive || !hEffPiInclusive) {
        std::cerr << "Error: Failed to calculate efficiency." << std::endl;
        return;
    }

    bool isFirstDraw = true;

    // =========================================================
    // Area Factor Calculation for UE Subtraction
    // Signal Area = pi * (RMax^2 - RMin^2)
    // BG Area     = 2 * pi * JetRadius^2 (Perp Cones x2)
    // =========================================================
    double areaSignal = TMath::Pi() * (kRMax*kRMax - kRMin*kRMin);
    double areaBg     = 2.0 * TMath::Pi() * kJetRadius * kJetRadius;
    double areaScaleFactor = areaSignal / areaBg;

    std::cout << ">>> Area Scale Factor: " << areaScaleFactor << " <<<" << std::endl;
    std::cout << "    Signal R: " << kRMin << " - " << kRMax << std::endl;
    std::cout << "    BG Jet R: " << kJetRadius << " (x2 cones)" << std::endl;

    for (size_t i = 0; i < jetPtBins.size(); ++i) {
        double jMin = jetPtBins[i].first;
        double jMax = jetPtBins[i].second;

        std::cout << "\n=== Jet pT: " << jMin << "-" << jMax << " GeV/c ===" << std::endl;

        // 1. Get PID Correction Maps (From Signal region data)
        auto mapsPr = getPidCorrectionMaps(dataFile, "proton", dataCentMin, dataCentMax, jMin, jMax, kRMin, kRMax);
        auto mapsPi = getPidCorrectionMaps(dataFile, "pion", dataCentMin, dataCentMax, jMin, jMax, kRMin, kRMax);
        
        // 2. Get Signal Spectra (PID Corrected, 1/PID_Eff applied)
        TH1D* hPrSig = createCorrectedSpectrum(dataFile, "jetpVsPtForPr", mapsPr.first, mapsPr.second, dataCentMin, dataCentMax, jMin, jMax, kRMin, kRMax);
        TH1D* hPiSig = createCorrectedSpectrum(dataFile, "jetpVsPtForPi", mapsPi.first, mapsPi.second, dataCentMin, dataCentMax, jMin, jMax, kRMin, kRMax);

        // 3. Get Background Spectra (PID Corrected)
        TH1D* hPrBg = createCorrectedBackgroundSpectrum(dataFile, "pVsPtForPrOutOfJet", mapsPr.first, mapsPr.second, dataCentMin, dataCentMax, jMin, jMax);
        TH1D* hPiBg = createCorrectedBackgroundSpectrum(dataFile, "pVsPtForPiOutOfJet", mapsPi.first, mapsPi.second, dataCentMin, dataCentMax, jMin, jMax);

        if (!hPrSig || !hPiSig || !hPrBg || !hPiBg) {
            std::cerr << "Skip Jet pT " << jMin << "-" << jMax << " due to missing hists." << std::endl;
            continue;
        }

        // 4. Perform Subtraction: Signal - (Factor * BG)
        TH1D* hPrSub = (TH1D*)hPrSig->Clone(TString::Format("hPrSub_j%.0f", jMin));
        hPrSub->Add(hPrBg, -1.0 * areaScaleFactor);

        TH1D* hPiSub = (TH1D*)hPiSig->Clone(TString::Format("hPiSub_j%.0f", jMin));
        hPiSub->Add(hPiBg, -1.0 * areaScaleFactor);

        // 5. Apply Tracking Efficiency Correction (Divide by Eff)
        TH1D* hPrFinal = (TH1D*)hPrSub->Clone(TString::Format("hPrFinal_j%.0f", jMin));
        hPrFinal->Divide(hEffPrInclusive);
        
        TH1D* hPiFinal = (TH1D*)hPiSub->Clone(TString::Format("hPiFinal_j%.0f", jMin));
        hPiFinal->Divide(hEffPiInclusive);

        // 6. Calculate Ratio
        TH1D* hRatio = (TH1D*)hPrFinal->Clone(TString::Format("hRatio_j%.0f", jMin));
        hRatio->Divide(hPiFinal);

        // Style
        int col = colors[i % colors.size()];
        int mrk = markers[i % markers.size()];
        hRatio->SetMarkerStyle(mrk); hRatio->SetMarkerColor(col); hRatio->SetLineColor(col); hRatio->SetLineWidth(2);
        
        if (isFirstDraw) {
            hRatio->GetXaxis()->SetRangeUser(0.0, 6.0);
            hRatio->GetYaxis()->SetRangeUser(0.0, 1.6); 
            
            hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio (UE Subtracted)");
            hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hRatio->GetYaxis()->SetTitleOffset(1.2);
            hRatio->Draw("PE");
            isFirstDraw = false; 
        } else {
            hRatio->Draw("PE SAME");
        }
        leg->AddEntry(hRatio, TString::Format("#it{p}_{T,jet} > %.0f GeV/#it{c}", jMin), "pe");

        // Clean up temporary hists
        delete hPrSig; delete hPiSig; delete hPrBg; delete hPiBg;
        delete hPrSub; delete hPiSub; delete hPrFinal; delete hPiFinal;
    }

    leg->Draw();

    TPaveText *pt = new TPaveText(0.15, 0.75, 0.5, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextSize(0.035); pt->SetTextAlign(12); pt->SetTextFont(42);
    pt->AddText("#bf{This thesis (UE Subtracted)}");
    pt->AddText(TString::Format("%.0f-%.0f%% Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", dataCentMin, dataCentMax));
    pt->AddText("ch-particle jet, Anti-#it{k}_{T}, #it{R} = 0.4");
    pt->AddText("#it{p}_{T,track} > 0.15 GeV/#it{c}");
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    pt->AddText("|#it{#eta}_{track}| < 0.9");
    pt->Draw();

    //c1->SaveAs("JetPtDependentRatio_UESubtracted_MC0_100.pdf");
}