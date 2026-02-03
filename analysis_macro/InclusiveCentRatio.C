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
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility>
#include <iomanip>

// ============================================================================================
// Constants
// ============================================================================================
const double kJetPtMin = 20.0; 
const double kJetPtMax = 200.0; 

// ============================================================================================
// Helper Functions: Fitting
// ============================================================================================

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // --- Grid Search Ranges (Updated) ---
    // 修正点1: グリッドサーチの範囲を少し現実に即したものに調整
    // KaonはProton(0)に近い位置に寄ってくるため、1.5-3.0は遠すぎます。0.5-1.5付近を重点的に見ます。
    std::vector<double> p4_range; 
    for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 

    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // ★Momentum switch: 4.2 GeV以上で挙動を変えるフラグ
    bool isHighPt = (p_min > 4.1);

    // --- Step 1: Grid Search ---
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                
                // Proton (Signal)
                // 修正点2: 4.2GeV以上なら、Protonの初期値を少し下げておく
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                tempFit->SetParameter(0, pr_init_height); 
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                // Backgrounds
                // 修正点3: 4.2GeV以上なら、Kaonの初期値を高く設定する (0.07 -> 0.30)
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;
                
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
                delete tempFit;
            }
        }
    }
    
    // --- Step 2: Final Fit with Best Parameters ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    // Limits
    // 修正点4: 【最重要】4.2GeV以上の場合、Protonの高さ上限を厳しくする
    // これによりProtonがデータを突き抜けるのを防ぎ、余った部分をKaonに埋めさせます
    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88); // 上限を88%程度に抑える
    } else {
        finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    }
    
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);

    // Kaon Limits
    // 修正点5: 4.2GeV以上の場合、Kaonの下限を引き上げる
    if (isHighPt) {
        finalFit->SetParLimits(3, peak_height * 0.20, peak_height * 1.2); // 最低でも15%はKaonがいると強制
    } else {
        finalFit->SetParLimits(3, peak_height * 0.05, peak_height * 1.0);
    }

    // Kaon Mean Limit: グリッドサーチ範囲と矛盾しないように広げておく
    finalFit->SetParLimits(4, 0.4, 1.2); 
    finalFit->SetParLimits(5, 0.6, 1.2);

    // Pion Limits (変更なし)
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.5, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);

    // Electron Limits (変更なし)
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    int oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;

    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.2) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.4) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.4) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    TF1* tempFit = new TF1("tempFitPionSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                double params[12] = {
                    peak_height, 0.0, 1.5,           // Pion
                    peak_height * 0.1, p4_init, 1.0,  // Kaon
                    peak_height * 0.2, p7_init, 1.5,  // Proton
                    peak_height * 0.1, p10_init, 1.0  // Electron
                };
                tempFit->SetParameters(params);
                tempFit->SetParLimits(0, peak_height * 0.5, peak_height * 1.5);
                tempFit->SetParLimits(1, -0.05, 0.05);
                tempFit->SetParLimits(2, 0.4, 1.8);
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
    finalFit->SetParLimits(0, peak_height * 0.5, peak_height * 1.5);
    finalFit->SetParLimits(1, -0.05, 0.05);
    finalFit->SetParLimits(2, 0.4, 1.8);
    finalFit->SetParLimits(3, 0.0, peak_height * 0.6);
    finalFit->SetParLimits(4, -3.8, -2.4);
    finalFit->SetParLimits(5, 0.5, 1.2);
    finalFit->SetParLimits(6, 0.0, peak_height * 0.8);
    finalFit->SetParLimits(7, -5.0, -3.7);
    finalFit->SetParLimits(8, 1.2, 1.8);
    finalFit->SetParLimits(9, 0.0, peak_height * 0.6);
    finalFit->SetParLimits(10, 3.0, 8.0);
    finalFit->SetParLimits(11, 1.0, 1.8);

    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    hist->Fit(finalFit, "QRB0");
    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

// ============================================================================================
// Core Logic: Purity & PID Efficiency
// ============================================================================================

std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(TString fileName, TString particleType, double centMin, double centMax) {
    std::map<double, double> purityMap;
    std::map<double, double> effMap;
    
    TString sparseName;
    double sigMin, sigMax; 

    if (particleType == "proton") {
        sparseName = "tpcTofPr"; 
        sigMin = -3.5; sigMax = 0.5; 
    } else {
        sparseName = "tpcTofPi"; 
        sigMin = -0.5; sigMax = 3.5; 
    }

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) return {purityMap, effMap};
    
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { f->Close(); return {purityMap, effMap}; }

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax); 

    double p_start = 0.6; 
    double p_end = 8.0; 
    double p_step = 0.2; 
    int totalSteps = (int)((p_end - p_start) / p_step);
    int currentStep = 0;

    std::cout << "--- Processing " << particleType << " (" << centMin << "-" << centMax << "%) ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        currentStep++;
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        // Progress Bar
        int barWidth = 40;
        float progress = (float)currentStep / totalSteps;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            effMap[p_center] = idealInCut / idealTotal; 
            
            // Log Purity for low pT
            printf(" p=%.2f: Purity=%.3f, PID Eff=%.3f (Low pT ideal)\n", p_center, 1.0, effMap[p_center]);
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

        // Log Purity
        printf(" p=%.2f: Purity=%.3f, PID Eff=%.3f\n", p_center, purity, pidEff);

        delete h_nSigma; delete fit; delete fSig;
    }
    std::cout << std::endl;
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

        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0; 
        
        for (auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_val) < min_diff) {
                min_diff = std::abs(map_p - p_val);
                purity = map_pur;
                pidEff = pidEffMap[map_p]; 
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
    f->Close();
    return hPt;
}

TH1D* calculateEfficiency(TString fileName, TString partName, double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    
    THnSparseD *hRecoSp = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptHistogram%sTof", partName.Data()));
    THnSparseD *hGenSp  = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptGenerated%s", partName.Data()));
    
    if (!hRecoSp || !hGenSp) { if(f) f->Close(); return nullptr; }

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
// Main
// ============================================================================================

void PlotInclusiveRatioFinal() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    const TString dataFile = "AnalysisResults-104.root";
    const TString mcFile   = "AnalysisResults-91.root";

    // Centrality Bins Definition
    std::vector<std::pair<double, double>> centBins = {
        {0.0, 10.0},
        {10.0, 30.0},
        {30.0, 50.0},
        {50.0, 100.0} // 4つ目: 50-100%
    };

    // Visualization Settings
    std::vector<int> colors = {kBlack, kRed, kBlue, kGreen+2};
    std::vector<int> markers = {20, 21, 33, 34};

    // Canvas Setup
    TCanvas *c1 = new TCanvas("cSpectra", "Inclusive Spectra", 800, 700);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12); gPad->SetRightMargin(0.05);

    TCanvas *c2 = new TCanvas("cRatio", "Proton/Pion Ratio", 800, 700);
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.12); gPad->SetRightMargin(0.05); gPad->SetGridy();

    TLegend *legSpectra = new TLegend(0.55, 0.65, 0.90, 0.88);
    legSpectra->SetBorderSize(0); legSpectra->SetFillStyle(0); legSpectra->SetHeader("Pb-Pb Cent. (Pion)");

    TLegend *legRatio = new TLegend(0.2, 0.65, 0.5, 0.88);
    legRatio->SetBorderSize(0); legRatio->SetFillStyle(0);

    // =========================================================
    // STEP 0: Calculate Efficiency using 0-100% MC Data ONCE
    // =========================================================
    std::cout << ">>> Calculating Tracking/TOF Efficiency (0-100% Inclusive) <<<" << std::endl;
    TH1D* hEffPrInclusive = calculateEfficiency(mcFile, "Proton", 0.0, 100.0);
    TH1D* hEffPiInclusive = calculateEfficiency(mcFile, "Pion", 0.0, 100.0);

    if (!hEffPrInclusive || !hEffPiInclusive) {
        std::cerr << "Critical Error: Could not calculate inclusive efficiency!" << std::endl;
        return;
    }

    // Storage to keep histograms in scope
    std::vector<TH1D*> storedRatios;
    std::vector<TH1D*> storedPions;

    // --- Loop over Centrality Bins ---
    for (size_t i = 0; i < centBins.size(); ++i) {
        double centMin = centBins[i].first;
        double centMax = centBins[i].second;

        std::cout << "\n=== Processing Centrality: " << centMin << "-" << centMax << "% ===" << std::endl;

        // 1. Get Correction Maps
        auto mapsPr = getPidCorrectionMaps(dataFile, "proton", centMin, centMax);
        auto mapsPi = getPidCorrectionMaps(dataFile, "pion", centMin, centMax);
        
        // 2. Spectra (Corrected for PID Purity/Eff)
        TH1D* hPrRaw = createCorrectedSpectrum(dataFile, "pVsPtForPr", mapsPr.first, mapsPr.second, centMin, centMax);
        TH1D* hPiRaw = createCorrectedSpectrum(dataFile, "pVsPtForPi", mapsPi.first, mapsPi.second, centMin, centMax);

        if (!hPrRaw || !hPiRaw) {
            std::cerr << "Skip Cent " << centMin << "-" << centMax << " due to missing hists." << std::endl;
            continue;
        }

        // 3. Efficiency Correction (Using 0-100% Efficiency)
        TH1D* hPrCorr = (TH1D*)hPrRaw->Clone(TString::Format("hPrCorr_%.0f", centMin));
        hPrCorr->Divide(hEffPrInclusive); // Divide by 0-100% eff
        
        TH1D* hPiCorr = (TH1D*)hPiRaw->Clone(TString::Format("hPiCorr_%.0f", centMin));
        hPiCorr->Divide(hEffPiInclusive); // Divide by 0-100% eff

        // 4. Ratio
        TH1D* hRatio = (TH1D*)hPrCorr->Clone(TString::Format("hRatio_%.0f", centMin));
        hRatio->Divide(hPiCorr);

        // --- Styling ---
        int col = colors[i % colors.size()];
        int mrk = markers[i % markers.size()];

        hPiCorr->SetMarkerStyle(mrk); hPiCorr->SetMarkerColor(col); hPiCorr->SetLineColor(col);
        hRatio->SetMarkerStyle(mrk);  hRatio->SetMarkerColor(col);  hRatio->SetLineColor(col);
        
        // --- Draw Spectra (Pion only to avoid clutter) ---
        c1->cd();
        if (i == 0) {
            hPiCorr->GetXaxis()->SetRangeUser(0.0, 6.0);
            hPiCorr->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}");
            hPiCorr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hPiCorr->Draw("PE");
        } else {
            hPiCorr->Draw("PE SAME");
        }
        legSpectra->AddEntry(hPiCorr, TString::Format("%.0f-%.0f%%", centMin, centMax), "pe");

        // --- Draw Ratio ---
        c2->cd();
        if (i == 0) {
            hRatio->GetXaxis()->SetRangeUser(0.0, 6.0);
            hRatio->GetYaxis()->SetTitle("p / #pi");
            hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hRatio->SetMinimum(0.0);
            hRatio->SetMaximum(1.2); 
            hRatio->Draw("PE");
        } else {
            hRatio->Draw("PE SAME");
        }
        legRatio->AddEntry(hRatio, TString::Format("%.0f-%.0f%%", centMin, centMax), "pe");

        storedRatios.push_back(hRatio);
        storedPions.push_back(hPiCorr);
    }

    c1->cd(); legSpectra->Draw();
    c2->cd(); legRatio->Draw();

    // Thesis Text Box
    TPaveText *pt = new TPaveText(0.55, 0.65, 0.88, 0.88, "NDC");
    pt->SetFillStyle(0); 
    pt->SetBorderSize(0);
    pt->SetTextSize(0.035);
    pt->SetTextAlign(12);

    pt->SetTextFont(42);
    
    pt->AddText("#bf{This thesis}");
    pt->AddText("|#eta_{track}| < 0.9");
    pt->AddText("Pb-Pb inclusive #sqrt{#it{s}_{NN}} = 5.36 TeV");
    
    pt->Draw();

    c1->SaveAs("InclusiveSpectra_CentDep_0-100Eff.pdf");
    c2->SaveAs("InclusiveRatio_CentDep_0-100Eff.pdf");
}