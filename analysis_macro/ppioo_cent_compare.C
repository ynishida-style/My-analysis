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
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility> 

// ============================================================================================
// Helper Functions: Fitting
// (O-O用にチューニングされたパラメータ設定を使用)
// ============================================================================================
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // --- Grid Search Ranges ---
    std::vector<double> p4_range; 
    for (double p4 = 0.8; p4 <= 1.5; p4 += 0.1) p4_range.push_back(p4); 

    // Pion (Param 6-8)
    std::vector<double> p7_range; 
    double p7_start = (p_min > 5.0) ? 2.5 : 3.5; 
    for (double p7 = p7_start; p7 <= 6.0; p7 += 0.5) p7_range.push_back(p7);

    // Electron (Param 9-11)
    std::vector<double> p10_range; 
    for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // 領域判定
    bool isMidPt  = (p_min >= 3.0 && p_min < 4.8); // Kaon混入が激しい領域
    bool isHighPt = (p_min >= 4.8);               // Pionが近づく領域

    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    // --- Step 1: Grid Search ---
    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                double pr_init_h = peak_height * 0.7;
                double ka_init_h = peak_height * 0.3;
                
                if (isMidPt) {
                     pr_init_h = peak_height * 0.5; 
                     ka_init_h = peak_height * 0.5;
                }

                double params[12] = {
                    pr_init_h, 0.0, 1.0,            
                    ka_init_h, p4_init, 1.0,        
                    peak_height * 0.1, p7_init, 1.5,  
                    peak_height * 0.05, p10_init, 2.0 
                };
                tempFit->SetParameters(params);
                
                tempFit->SetParLimits(0, 0, peak_height * 1.5);
                tempFit->SetParLimits(1, -0.2, 0.2); 
                tempFit->SetParLimits(2, 0.8, 1.1); 

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

    // --- Step 2: Final Fit ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    
    // Limits Tuning
    if (isMidPt) {
        finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 0.90); 
    } else {
        finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 1.5); 
    }
    
    finalFit->SetParLimits(1, -0.1, 0.1); 
    finalFit->SetParLimits(2, 0.95, 1.01); 

    if (isMidPt) {
        finalFit->SetParLimits(3, peak_height * 0.2, peak_height * 0.9);
    } else {
        finalFit->SetParLimits(3, 0.0, peak_height * 1.5);
    }
    
    finalFit->SetParLimits(4, 0.6, 1.8); 
    finalFit->SetParLimits(5, 0.8, 1.3); 

    finalFit->SetParLimits(6, 0.0, peak_height * 1.2);
    if (isHighPt) {
        finalFit->SetParLimits(7, 2.0, 7.0); 
    } else {
        finalFit->SetParLimits(7, 3.0, 7.0); 
    }
    finalFit->SetParLimits(8, 0.5, 2.5);

    finalFit->SetParLimits(9, 0.0, peak_height * 0.5);
    finalFit->SetParLimits(10, 7.0, 15.0);
    finalFit->SetParLimits(11, 0.6, 2.0);

    hist->Fit(finalFit, "QNRB0"); 
    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // std::cout << "--- Starting parameter optimization for PION in p bin " << p_min << "-" << p_max << " ---" << std::endl;
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
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
                delete tempFit;
            }
        }
    }

    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed); 
    finalFit->SetLineWidth(2);

    finalFit->SetParLimits(0, peak_height * 0.7, peak_height * 1.1); finalFit->SetParLimits(1, -0.05, 0.05); finalFit->SetParLimits(2, 0.4, 1.8);
    finalFit->SetParLimits(3, peak_height * 0.01, peak_height * 0.4); finalFit->SetParLimits(4, -3.8, -2.4); finalFit->SetParLimits(5, 0.5, 1.2);
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 0.6); finalFit->SetParLimits(7, -5.0, -3.7); finalFit->SetParLimits(8, 1.2, 1.8);
    finalFit->SetParLimits(9, peak_height * 0.01, peak_height * 0.5); finalFit->SetParLimits(10, 3.0, 8.0); finalFit->SetParLimits(11, 1.0, 1.8);

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

// ============================================================================================
// Core Logic
// ============================================================================================

std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(TString fileName, TString particleType, double centMin, double centMax) {
    std::map<double, double> purityMap;
    std::map<double, double> effMap;
    TString sparseName;
    double sigMin, sigMax; 

    if (particleType == "proton") { sparseName = "tpcTofPr"; sigMin = -3.5; sigMax = 0.5; } 
    else { sparseName = "tpcTofPi"; sigMin = -0.5; sigMax = 3.5; }

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Error: Cannot open file " << fileName << std::endl; return {purityMap, effMap}; }
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { f->Close(); return {purityMap, effMap}; }

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax); 
    double p_start = 0.6; double p_end = 8.0; double p_step = 0.2; 

    std::cout << "--- Calculating Maps for " << particleType << " (" << centMin << "-" << centMax << "%) ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            effMap[p_center] = idealInCut / idealTotal; 
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

        purityMap[p_center] = purity;
        effMap[p_center]    = pidEff;
        delete h_nSigma; delete fit; delete fSig;
    }
    f->Close();
    return {purityMap, effMap};
}

TH1D* createCorrectedSpectrum(TString fileName, TString histName, 
                              std::map<double, double>& purityMap, 
                              std::map<double, double>& pidEffMap, 
                              double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", histName.Data()));
    if (!hSparse) return nullptr;

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
    TH2D *h2 = (TH2D*)hSparse->Projection(0, 1);
    TH1D *hPt = h2->ProjectionX(TString::Format("hPt_%s_%.0f_%.0f", histName.Data(), centMin, centMax));
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
        if (min_diff > 0.3) { if (p_val > 3.0) continue; }
        if (pidEff < 0.01) pidEff = 1.0; 
        double totalCorrectionFactor = purity / pidEff;

        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double content = h2->GetBinContent(bx, by);
            double error   = h2->GetBinError(bx, by);
            if (content <= 0) continue;
            double corrContent = content * totalCorrectionFactor;
            double corrError   = error * totalCorrectionFactor; 
            double curContent = hPt->GetBinContent(bx);
            double curError   = hPt->GetBinError(bx);
            hPt->SetBinContent(bx, curContent + corrContent);
            hPt->SetBinError(bx, std::sqrt(curError*curError + corrError*corrError));
        }
    }
    f->Close();
    return hPt;
}

TH1D* calculateEfficiency(TString fileName, TString partName, double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    THnSparseD *hRecoSp = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptHistogram%sTof", partName.Data()));
    THnSparseD *hGenSp  = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptGenerated%s", partName.Data()));
    if (!hRecoSp || !hGenSp) { f->Close(); return nullptr; }

    hRecoSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    
    TH1D *hReco = hRecoSp->Projection(0);
    TH1D *hGen  = hGenSp->Projection(0);
    TH1D *hEff = (TH1D*)hReco->Clone(TString::Format("eff_%s_%.0f_%.0f", partName.Data(), centMin, centMax));
    hEff->Divide(hReco, hGen, 1, 1, "B");
    hEff->SetDirectory(0);
    f->Close();
    return hEff;
}

// ============================================================================================
// Wrapper Function: Get Ratio for Specific Centrality
// ============================================================================================
TH1D* GetRatioForCentrality(TString dataFile, TString mcFile, double centMin, double centMax, int color, int markerStyle, TString label) {
    
    // 1. Get PID Correction Maps (Data)
    auto mapsPr = getPidCorrectionMaps(dataFile, "proton", centMin, centMax);
    auto mapsPi = getPidCorrectionMaps(dataFile, "pion", centMin, centMax);
    
    if (mapsPr.first.empty() || mapsPi.first.empty()) {
        std::cerr << "Failed to get PID maps for " << label << std::endl;
        return nullptr;
    }

    // 2. Create Corrected Spectra
    TH1D* hPrRaw = createCorrectedSpectrum(dataFile, "pVsPtForPr", mapsPr.first, mapsPr.second, centMin, centMax);
    TH1D* hPiRaw = createCorrectedSpectrum(dataFile, "pVsPtForPi", mapsPi.first, mapsPi.second, centMin, centMax);

    if (!hPrRaw || !hPiRaw) return nullptr;

    // 3. Calculate Efficiency (MC)
    // 統計を確保するため、MCはInclusive (0-100%) を使用
    TH1D* hEffPr = calculateEfficiency(mcFile, "Proton", 0.0, 100.0);
    TH1D* hEffPi = calculateEfficiency(mcFile, "Pion", 0.0, 100.0);

    if (!hEffPr || !hEffPi) return nullptr;

    // 4. Apply Efficiency Correction
    TH1D* hPrCorr = (TH1D*)hPrRaw->Clone(TString::Format("hPrCorr_%.0f", centMin));
    hPrCorr->Divide(hEffPr);
    TH1D* hPiCorr = (TH1D*)hPiRaw->Clone(TString::Format("hPiCorr_%.0f", centMin));
    hPiCorr->Divide(hEffPi);

    // 5. Calculate Ratio (Rebinning included as per original code)
    hPrCorr->Rebin(2);
    hPiCorr->Rebin(2);
    
    TH1D* hRatio = (TH1D*)hPrCorr->Clone(TString::Format("hRatio_%.0f_%.0f", centMin, centMax));
    hRatio->Divide(hPiCorr);

    // 6. Styling
    hRatio->SetMarkerStyle(markerStyle);
    hRatio->SetMarkerColor(color);
    hRatio->SetLineColor(color);
    hRatio->SetLineWidth(2);
    hRatio->SetTitle(label);
    hRatio->SetDirectory(0);

    // Cleanup
    delete hPrRaw; delete hPiRaw; delete hEffPr; delete hEffPi;
    delete hPrCorr; delete hPiCorr;

    return hRatio;
}

// ============================================================================================
// Main
// ============================================================================================

void CompareCentralitiesOO() {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    
    const TString dataFile = "AnalysisResults-100.root";
    const TString mcFile   = "AnalysisResults-91.root";

    std::cout << "=== Processing Centrality 0-30% ===" << std::endl;
    TH1D* hRatio0_30 = GetRatioForCentrality(dataFile, mcFile, 0.0, 30.0, kRed+1, 20, "0-30%");

    std::cout << "\n=== Processing Centrality 60-100% ===" << std::endl;
    TH1D* hRatio60_100 = GetRatioForCentrality(dataFile, mcFile, 60.0, 100.0, kBlue+1, 24, "60-100%");

    if (!hRatio0_30 || !hRatio60_100) {
        std::cerr << "Error: Could not generate one of the ratios." << std::endl;
        return;
    }

    // --- Plotting ---
    TCanvas *c1 = new TCanvas("cCompOO", "Ratio Comparison O-O", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetGridy();

    // Axis settings
    hRatio0_30->SetTitle("");
    hRatio0_30->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hRatio0_30->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio0_30->GetYaxis()->SetTitleOffset(1.3);
    hRatio0_30->GetXaxis()->SetTitleSize(0.045);
    hRatio0_30->GetYaxis()->SetTitleSize(0.045);
    hRatio0_30->GetXaxis()->SetRangeUser(0.0, 6.0);
    
    hRatio0_30->SetMinimum(0.0);
    hRatio0_30->SetMaximum(1.2); 

    hRatio0_30->Draw("E1");
    hRatio60_100->Draw("E1 SAME");

    // Legend
    TLegend *leg = new TLegend(0.6, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Centrality Class");
    leg->AddEntry(hRatio0_30, "0 - 10%", "pe");
    leg->AddEntry(hRatio60_100, "60 - 100%", "pe");
    leg->Draw();

    // Text info
    TPaveText *pt = new TPaveText(0.15, 0.78, 0.5, 0.92, "NDC");
    pt->SetFillStyle(0); 
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42); 
    pt->SetTextSize(0.035);
    
    pt->AddText("#bf{This thesis}");
    pt->AddText("O-O #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText("|#eta_{track}| < 0.9");
    pt->AddText("Inclusive");
    pt->Draw();

    c1->SaveAs("Ratio_Comparison_OO_0-30_vs_60-100.pdf");
}