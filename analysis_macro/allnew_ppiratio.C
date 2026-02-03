#include <THnSparse.h>
#include <TH1.h>
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
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility>
#include <Rtypes.h>

//________________________________________________________________________________________________
// Helper Functions: Fitting (Optimized & Cleaned)
// Updated with High-Pt Logic and Adjusted Limits
//________________________________________________________________________________________________

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // --- Grid Search Ranges ---
    // KaonはProton(0)に近い位置に寄るため、0.5-1.2付近を重点的に探索
    std::vector<double> p4_range; 
    for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 

    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // ★Momentum switch: 4.2 GeV以上で挙動を変えるフラグ
    bool isHighPt = (p_min > 4.1);

    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    // --- Step 1: Grid Search ---
    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                // 4.2GeV以上なら、Protonの初期値を少し下げ、Kaonを上げる
                double pr_init_height   = isHighPt ? peak_height * 0.8 : peak_height;
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;

                double params[12] = {
                    pr_init_height, 0.0, 1.5,         // Proton (Signal)
                    kaon_init_height, p4_init, 1.5,   // Kaon
                    peak_height * 0.7, p7_init, 2.0,  // Pion
                    peak_height * 0.1, p10_init, 2.5  // Electron
                };
                tempFit->SetParameters(params);
                
                // Limits for Search
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParLimits(2, 0.7, 2.0);

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

    // --- Step 2: Final Fit with Best Parameters ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    
    // --- Limits applied for final fit (厳格化 & High-Pt対応) ---
    
    // 【最重要】4.2GeV以上の場合、Protonの高さ上限を厳しくする
    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88); // 上限を88%程度に
    } else {
        finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    }

    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);

    // Kaon Limits
    // 4.2GeV以上の場合、Kaonの下限を引き上げる
    if (isHighPt) {
        finalFit->SetParLimits(3, peak_height * 0.20, peak_height * 1.2); // 最低でも20%はKaonがいると強制
    } else {
        finalFit->SetParLimits(3, peak_height * 0.05, peak_height * 1.0);
    }

    finalFit->SetParLimits(4, 0.4, 1.2); // Mean Limit
    finalFit->SetParLimits(5, 0.6, 1.2); // Sigma Limit

    // Pion Limits
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.5, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);

    // Electron Limits
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);

    hist->Fit(finalFit, "QNRB0"); 

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // Pion用は基本的にCode B仕様を使用
    std::cout << "--- Starting parameter optimization for PION in p bin " << p_min << "-" << p_max << " ---" << std::endl;

    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    TF1* tempFit = new TF1("tempFitPionSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                double params[12] = {
                    peak_height, 0.0, 1.5,            // Signal (Pion)
                    peak_height * 0.1, p4_init, 1.0,  // Kaon
                    peak_height * 0.2, p7_init, 1.5,  // Proton
                    peak_height * 0.1, p10_init, 1.0  // Electron
                };
                tempFit->SetParameters(params);
                tempFit->SetParLimits(0, peak_height * 0.7, peak_height * 1.1);
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
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    // Limits (Code B仕様)
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

    hist->Fit(finalFit, "QNRB0");

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

//________________________________________________________________________________________________
// PID Map Calculation (Purity & PID Efficiency)
//________________________________________________________________________________________________

std::pair<std::map<double, double>, std::map<double, double>> calculatePuritiesAndEfficiencies(
    const TString particleType, const TString fileName, 
    double r_min, double r_max, 
    double jetPt_min, double jetPt_max,
    double cent_min, double cent_max) 
{
    // ... 設定部 ...
    TString sparseName;
    double signal_min_cut, signal_max_cut;
    
    if (particleType == "proton") {
        sparseName = "jetTpcTofPr"; 
        signal_min_cut = -3.5; signal_max_cut = 0.5;
    } else if (particleType == "pion") {
        sparseName = "jetTpcTofPi"; 
        signal_min_cut = -0.5; signal_max_cut = 3.5;
    } else {
        std::cerr << "Error: Unknown particle type" << std::endl;
        return {{},{}};
    }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {{},{}}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) {
        std::cerr << "Error: Could not find " << fullPath << std::endl;
        inputFile->Close();
        return {{},{}};
    }

    // Axis Mapping: 0:P, 1:nSigma, 2:Dist, 3:JetPt, 4:Cent
    hSparse->GetAxis(2)->SetRangeUser(r_min, r_max);         
    hSparse->GetAxis(3)->SetRangeUser(jetPt_min, jetPt_max); 
    hSparse->GetAxis(4)->SetRangeUser(cent_min, cent_max);   

    std::map<double, double> purityMap;
    std::map<double, double> effMap;

    const double p_start = 0.6; 
    const double p_end = 8.0;
    const double p_step = 0.2;
    
    std::cout << "--- PID Calc: " << particleType << " (r=" << r_min << "-" << r_max << ") ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        // --- Low pT Region ---
        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            
            // Ideal Efficiency Calculation
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(signal_min_cut, signal_max_cut);
            double idealTotal = gaussian.Integral(-10, 10);
            double idealEff   = idealInCut / idealTotal;
            effMap[p_center] = idealEff;
            continue;
        }

        // --- High pT Region (Fit) ---
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        if (h_nSigma->GetEntries() < 50) {
            delete h_nSigma;
            continue;
        }
        
        TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        TF1* fSig = new TF1("fSig", "gaus", -10, 10);
        fSig->SetParameters(bestFit->GetParameter(0), bestFit->GetParameter(1), bestFit->GetParameter(2));

        double intSigInCut = fSig->Integral(signal_min_cut, signal_max_cut); 
        double intSigTotal = fSig->Integral(-10, 10);
        double intTotInCut = bestFit->Integral(signal_min_cut, signal_max_cut); 

        double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
        double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;

        if (purity > 1.0) purity = 1.0; if (purity < 0.0) purity = 0.0;
        if (pidEff > 1.0) pidEff = 1.0; if (pidEff < 0.0) pidEff = 0.0;

        purityMap[p_center] = purity;
        effMap[p_center] = pidEff;

        printf("  p = %.2f | Purity = %.3f | PID Eff = %.3f\n", p_center, purity, pidEff);
        
        delete h_nSigma; delete fSig; delete bestFit;
    }

    inputFile->Close();
    return {purityMap, effMap};
}

//________________________________________________________________________________________________
// Corrected Spectrum Creation (Fixed Logic)
//________________________________________________________________________________________________

TH1D* createCorrectedPtSpectrum(TString fileName, TString sparseName, 
                                const std::map<double, double>& purityMap, 
                                const std::map<double, double>& effMap,
                                double r_min, double r_max, 
                                double jetPt_min, double jetPt_max,
                                double cent_min, double cent_max) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) {
        inputFile->Close();
        return nullptr;
    }

    // --- 1. まず他軸のカットを適用 ---
    // Axis Mapping: 0:P, 1:nSigma, 2:Dist, 3:JetPt, 4:Cent
    hSparse->GetAxis(2)->SetRangeUser(r_min, r_max);         
    hSparse->GetAxis(3)->SetRangeUser(jetPt_min, jetPt_max); 
    hSparse->GetAxis(4)->SetRangeUser(cent_min, cent_max);   

    // --- 2. 【重要】PIDシグナル軸(Axis 1)にカットを適用する ---
    if (sparseName.Contains("Pr")) { // Proton case
        hSparse->GetAxis(1)->SetRangeUser(-3.5, 0.5); 
    } else { // Pion case
        hSparse->GetAxis(1)->SetRangeUser(-0.5, 3.5);
    }

    // --- 3. 運動量軸(Axis 0)への射影 (Raw Counts in PID Cut) ---
    TH1D *h_p_raw = hSparse->Projection(0); 
    h_p_raw->SetDirectory(nullptr);
    TString histName = TString::Format("%s_corr_r_%.1f_%.1f", sparseName.Data(), r_min, r_max);
    h_p_raw->SetName(histName + "_raw");

    // --- 4. Binごとの補正 ---
    TH1D *h_p_corrected = (TH1D*)h_p_raw->Clone(histName);
    h_p_corrected->Reset();

    for (int i = 1; i <= h_p_raw->GetNbinsX(); ++i) {
        double p_center = h_p_raw->GetBinCenter(i);
        double content  = h_p_raw->GetBinContent(i);
        double error    = h_p_raw->GetBinError(i);

        if (content <= 0) continue;
        if (p_center > 10.0) continue; // 範囲外ガード

        // マップから補正値を取得
        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0;

        for(auto const& [p_val, pur_val] : purityMap) {
            double diff = std::abs(p_val - p_center);
            if (diff < min_diff) { 
                min_diff = diff;
                purity = pur_val;
                if (effMap.count(p_val)) pidEff = effMap.at(p_val);
            }
        }
        
        // マップとのマッチングが遠すぎる場合はスキップ、またはIdeal (低pT)
        if (min_diff > 0.3) {
            if (p_center <= 3.0) { /* use defaults 1.0/ideal */ } 
            else { continue; }
        }

        if (pidEff < 0.01) pidEff = 1.0; 

        // 補正係数: N_corr = N_raw(in cut) * Purity / PID_Eff
        double correctionFactor = purity / pidEff;
        
        h_p_corrected->SetBinContent(i, content * correctionFactor);
        h_p_corrected->SetBinError(i, error * correctionFactor);
    }

    delete h_p_raw;
    inputFile->Close();
    return h_p_corrected;
}

//________________________________________________________________________________________________
// Efficiency Calculation (Tracking + TOF match) - Unchanged
//________________________________________________________________________________________________

TH1D* calculateTotalEfficiency(TString particleType, TString fileName, 
                               double jetPt_min, double jetPt_max,
                               double cent_min, double cent_max) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    
    TString numName = TString::Format("jet-shape-task/ptHistogram%sTof", particleType.Data());
    TString denName = TString::Format("jet-shape-task/ptGenerated%s", particleType.Data());
    
    THnSparseD *hNumSparse = (THnSparseD*)file->Get(numName);
    THnSparseD *hDenSparse = (THnSparseD*)file->Get(denName);
    
    if (!hNumSparse || !hDenSparse) {
        file->Close();
        return nullptr;
    }

    hNumSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hNumSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);

    hDenSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hDenSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);

    TH1D* hNum = hNumSparse->Projection(0);
    TH1D* hDen = hDenSparse->Projection(0);
    hNum->Sumw2(); hDen->Sumw2();
    
    TH1D* hEff = (TH1D*)hNum->Clone(TString::Format("hEffTotal_%s", particleType.Data()));
    hEff->SetDirectory(nullptr);
    hEff->Divide(hNum, hDen, 1.0, 1.0, "B"); 

    delete hNum; delete hDen;
    file->Close();
    return hEff;
}

//________________________________________________________________________________________________
// Main Macro
//________________________________________________________________________________________________

void calculateProtonPionRatioO2() {
    gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTextFont(42);

    const TString kDataFileName = "AnalysisResults-93.root"; 
    const TString kMcFileName   = "AnalysisResults-91.root"; 

    // --- 設定 ---
    const double kJetPtMin = 60.0;
    const double kJetPtMax = 100.0;
    const double kCentMin = 0.0;
    const double kCentMax = 30.0;

    std::cout << "\n=====================================================================" << std::endl;
    std::cout << "=== Processing for Jet pT: " << kJetPtMin << " - " << kJetPtMax << " GeV/c ===" << std::endl;
    std::cout << "=====================================================================\n" << std::endl;

    std::vector<std::pair<double, double>> r_ranges;
    for (double r_min = 0.0; r_min < 0.7; r_min += 0.1) {
        r_ranges.push_back({r_min, r_min + 0.1});
    }

    // --- Efficiency (MC) ---
    std::cout << "\n--- Calculating Tracking+TOF Efficiencies (MC) ---" << std::endl;
    TH1D* hEffTotalProton = calculateTotalEfficiency("Proton", kMcFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    TH1D* hEffTotalPion   = calculateTotalEfficiency("Pion",   kMcFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    
    if (!hEffTotalProton || !hEffTotalPion) {
        std::cerr << "Failed to calculate efficiencies." << std::endl;
        return;
    }

    std::vector<TH1D*> h_ratios;

    for (const auto& range : r_ranges) {
        double r_min = range.first;
        double r_max = range.second;
        std::cout << "\n=== Processing r bin: " << r_min << " - " << r_max << " ===" << std::endl;

        // 1. PID Maps
        auto mapsProton = calculatePuritiesAndEfficiencies("proton", kDataFileName, r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        auto mapsPion   = calculatePuritiesAndEfficiencies("pion",   kDataFileName, r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);

        // 2. Corrected Spectra (Fixed Logic)
        TH1D* hPtProton_pidCorr = createCorrectedPtSpectrum(kDataFileName, "jetpVsPtForPr", mapsProton.first, mapsProton.second, r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        TH1D* hPtPion_pidCorr   = createCorrectedPtSpectrum(kDataFileName, "jetpVsPtForPi", mapsPion.first,   mapsPion.second,   r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        
        if (!hPtProton_pidCorr || !hPtPion_pidCorr) {
            std::cerr << "Skipping r=" << r_min << " due to missing spectra." << std::endl;
            continue;
        }

        // 3. Efficiency Correction
        TH1D* hPtProton_corrected = (TH1D*)hPtProton_pidCorr->Clone(TString::Format("hPtProton_corrected_r_%.1f", r_min));
        hPtProton_corrected->Divide(hEffTotalProton);

        TH1D* hPtPion_corrected = (TH1D*)hPtPion_pidCorr->Clone(TString::Format("hPtPion_corrected_r_%.1f", r_min));
        hPtPion_corrected->Divide(hEffTotalPion);

        // 4. Ratio Calculation
        int rebinFactor = 2; 
        hPtProton_corrected->Rebin(rebinFactor);
        hPtPion_corrected->Rebin(rebinFactor);

        TH1D* hRatio = (TH1D*)hPtProton_corrected->Clone(TString::Format("hRatio_p_pi_r_%.1f", r_min));
        hRatio->Divide(hPtPion_corrected);
        hRatio->SetDirectory(nullptr);
        h_ratios.push_back(hRatio);

        delete hPtProton_pidCorr; delete hPtPion_pidCorr;
        delete hPtProton_corrected; delete hPtPion_corrected;
    }

    delete hEffTotalProton; delete hEffTotalPion;

    // --- 5. Plotting ---
    TCanvas *c_ratio_overlay = new TCanvas("c_ratio_overlay", "Proton/Pion Ratio", 800, 600);
    c_ratio_overlay->cd();
    gPad->SetGridy();
    gPad->SetLeftMargin(0.12);

    std::vector<int> colors = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kCyan, kOrange + 7};
    std::vector<int> markers = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCross, kFullStar, kOpenCircle};
    TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.9);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.03);

    bool first = true;
    for (size_t i = 0; i < h_ratios.size(); ++i) {
        TH1D* hRatio = h_ratios[i];
        int color_idx = i % colors.size();
        hRatio->SetLineColor(colors[color_idx]);
        hRatio->SetMarkerColor(colors[color_idx]);
        hRatio->SetMarkerStyle(markers[i % markers.size()]);
        
        hRatio->GetXaxis()->SetRangeUser(0, 10); 

        if (first) {
            hRatio->SetTitle("p/#pi Ratio within Jets (Corrected)");
            hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
            hRatio->GetYaxis()->SetRangeUser(0.0, 1.5);
            hRatio->Draw("E1");
            first = false;
        } else {
            hRatio->Draw("E1 SAME");
        }
        legend->AddEntry(hRatio, TString::Format("%.1f < #it{r} < %.1f", r_ranges[i].first, r_ranges[i].second), "lpe");
    }
    legend->Draw();

    TPaveText *pt = new TPaveText(0.15, 0.7, 0.5, 0.9, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0);
    pt->SetTextFont(42);
    pt->AddText("ALICE Simulation / Data");
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", kJetPtMin, kJetPtMax));
    pt->Draw();
}