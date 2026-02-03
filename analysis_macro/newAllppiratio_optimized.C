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
#include <TString.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility> // for std::pair

// ============================================================================================
// Helper Functions: Fitting (Logic form Code A)
// ============================================================================================

/**
 * @brief プロトン用のグリッドサーチと最終フィット (Code A仕様)
 * High Pt (4.1 GeV以上) での挙動切り替えやパラメータ制限を含む
 */
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // --- Grid Search Ranges ---
    std::vector<double> p4_range; for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 
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
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                tempFit->SetParameter(0, pr_init_height); 
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                // Backgrounds
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
    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88); // 上限を厳しく
    } else {
        finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    }
    
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);

    // Kaon Limits
    if (isHighPt) {
        finalFit->SetParLimits(3, peak_height * 0.20, peak_height * 1.2);
    } else {
        finalFit->SetParLimits(3, peak_height * 0.05, peak_height * 1.0);
    }

    finalFit->SetParLimits(4, 0.4, 1.2); 
    finalFit->SetParLimits(5, 0.6, 1.2);

    // Pion Limits
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.5, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);

    // Electron Limits
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

/**
 * @brief パイオン用のグリッドサーチと最終フィット (Code A仕様)
 */
TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // Ranges
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.2) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.4) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.4) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // --- Step 1: Grid Search ---
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

    // --- Step 2: Final Fit ---
    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed); 
    finalFit->SetLineWidth(2);

    // Limits
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
    return finalFit;
}

// ============================================================================================
// Core Logic: Purity AND PID Efficiency Calculation (Merged Logic)
// ============================================================================================

/**
 * @brief 運動量ビンごとのPurityとPID Efficiencyを計算する関数
 * Code Aのロジックを使用 (PurityだけでなくEffも計算、低pTの扱いなど)
 * Code BのJet Taggedヒストグラムに対応
 */
std::pair<std::map<double, double>, std::map<double, double>> calculatePidCorrections(
    const TString particleType, const TString fileName, 
    double r_min, double r_max, 
    double jetPt_min, double jetPt_max,
    double cent_min, double cent_max) 
{
    std::map<double, double> purityMap;
    std::map<double, double> effMap;

    TString sparseName;
    double sigMin, sigMax; // カット範囲

    if (particleType == "proton") {
        sparseName = "jetTpcTofPr"; // Jet Tagged
        sigMin = -3.5; sigMax = 0.5;
    } else if (particleType == "pion") {
        sparseName = "jetTpcTofPi"; // Jet Tagged
        sigMin = -0.5; sigMax = 3.5;
    } else {
        std::cerr << "Error: Unknown particle type '" << particleType << "'" << std::endl;
        return {purityMap, effMap};
    }

    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {purityMap, effMap}; }
    
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) {
        std::cerr << "Error: Could not find " << fullPath << " in " << fileName << std::endl;
        inputFile->Close();
        return {purityMap, effMap};
    }

    std::cout << "--- Calculating Correction Maps for " << particleType << " (R=" << r_min << "-" << r_max << ") ---" << std::endl;
    std::cout << "      pT   |  Type   | Purity | PID Eff | Method" << std::endl; // ヘッダー表示

    double p_start = 0.6;
    double p_end = 8.0; 
    double p_step = 0.2;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        // --- Low Pt Logic (Code A) ---
        // 3.0 GeV以下はTOF分離が良いのでPurity=1.0、Efficiencyは理想的なガウス分布から計算
        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); // ideal sigma=1
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            double idealEff   = idealInCut / idealTotal;

            effMap[p_center] = idealEff;

            // ★追加: ログ出力 (Low Pt)
            printf(" %4.2f GeV/c | %-7s | %6.4f | %6.4f  | LowPt(Fixed)\n", 
                   p_center, particleType.Data(), 1.0, idealEff);

            continue;
        }

        // --- High Pt: Fit ---
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);         // Axis 0: p
        hSparse->GetAxis(2)->SetRangeUser(r_min, r_max);         // Axis 2: r
        hSparse->GetAxis(3)->SetRangeUser(jetPt_min, jetPt_max); // Axis 3: jet pT
        hSparse->GetAxis(4)->SetRangeUser(cent_min, cent_max);   // Axis 4: Centrality

        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->Rebin(2);
        
        if (h_nSigma->GetEntries() < 50) { 
            // エントリー不足の場合のログ
            printf(" %4.2f GeV/c | %-7s | ------ | ------  | Skipped(LowStats)\n", p_center, particleType.Data());
            delete h_nSigma; 
            continue; 
        }

        TF1* fit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        // Signal成分抽出
        TF1* fSig = new TF1("fSig", "gaus", -10, 10);
        fSig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));

        // 積分計算
        double intSigInCut = fSig->Integral(sigMin, sigMax); // カット内の信号
        double intSigTotal = fSig->Integral(-10, 10);        // 全信号 (PID効率の分母)
        double intTotInCut = fit->Integral(sigMin, sigMax);  // カット内の全粒子 (Purityの分母)

        // 1. Purity Calculation
        double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
        
        // 2. PID Efficiency Calculation
        double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;
        
        // ガード
        if (purity > 1.0) purity = 1.0;
        if (purity < 0.0) purity = 0.0;
        if (pidEff > 1.0) pidEff = 1.0;
        if (pidEff < 0.0) pidEff = 0.0;

        purityMap[p_center] = purity;
        effMap[p_center]    = pidEff;

        // ★追加: ログ出力 (High Pt / Fitted)
        printf(" %4.2f GeV/c | %-7s | %6.4f | %6.4f  | Fitted\n", 
               p_center, particleType.Data(), purity, pidEff);

        delete h_nSigma; delete fit; delete fSig;
    }

    // Restore ranges
    hSparse->GetAxis(0)->SetRange(1, hSparse->GetAxis(0)->GetNbins());
    hSparse->GetAxis(2)->SetRange(1, hSparse->GetAxis(2)->GetNbins());
    hSparse->GetAxis(3)->SetRange(1, hSparse->GetAxis(3)->GetNbins());
    hSparse->GetAxis(4)->SetRange(1, hSparse->GetAxis(4)->GetNbins());

    inputFile->Close();
    return {purityMap, effMap};
}

/**
 * @brief PurityとEfficiencyの両方を使って補正したスペクトルを作成する関数
 * (Code BのWeightedSpectrumをCode Aのロジックで拡張)
 */
TH1D* createCorrectedJetSpectrum(TString fileName, TString sparseName, 
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
        std::cerr << "Error: Could not find " << fullPath << " in " << fileName << std::endl;
        inputFile->Close();
        return nullptr;
    }

    hSparse->GetAxis(2)->SetRangeUser(r_min, r_max);         // Axis 2: r
    hSparse->GetAxis(3)->SetRangeUser(jetPt_min, jetPt_max); // Axis 3: jet pT
    hSparse->GetAxis(4)->SetRangeUser(cent_min, cent_max);   // Axis 4: Centrality

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1); // Y=p, X=pT
    h2_p_vs_pt->SetDirectory(nullptr);

    // Restore ranges
    hSparse->GetAxis(2)->SetRange(1, hSparse->GetAxis(2)->GetNbins());
    hSparse->GetAxis(3)->SetRange(1, hSparse->GetAxis(3)->GetNbins());
    hSparse->GetAxis(4)->SetRange(1, hSparse->GetAxis(4)->GetNbins());

    TH2D *h2_corrected = (TH2D*)h2_p_vs_pt->Clone();
    TString histName = TString::Format("%s_corr_r_%.1f_%.1f_jet_%.0f_%.0f_cent_%.0f_%.0f", 
                                     sparseName.Data(), r_min, r_max, jetPt_min, jetPt_max, cent_min, cent_max);
    h2_corrected->SetName(histName);
    h2_corrected->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        
        // マップから値を取得
        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0;

        for(auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_center) < min_diff) {
                min_diff = std::abs(map_p - p_center);
                purity = map_pur;
                pidEff = effMap.at(map_p);
            }
        }
        
        // マッチするビンが遠すぎる場合 (マップ範囲外)
        if (min_diff > 0.3) {
             // 3.0GeV以下でマップがないのはおかしい(LowPtLogicがあるはず)が、念のため
             if (p_center > 3.0) continue; 
        }

        if (pidEff < 0.01) pidEff = 1.0; // 安全策

        double totalCorrectionFactor = purity / pidEff;
        
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            
            // 補正を適用
            h2_corrected->SetBinContent(j_pt, i_p, content * totalCorrectionFactor);
            h2_corrected->SetBinError(j_pt, i_p, error * totalCorrectionFactor); // エラーも同様にスケーリング
        }
    }

    TH1D *h_pt_corrected = h2_corrected->ProjectionX();
    h_pt_corrected->SetDirectory(nullptr);
    h_pt_corrected->SetName("h_pt_" + histName);
    
    delete h2_p_vs_pt;
    delete h2_corrected;
    inputFile->Close();
    return h_pt_corrected;
}

// ============================================================================================
// Tracking & TOF Efficiency (Unchanged from Code B)
// ============================================================================================

/**
 * @brief トラッキング効率を計算する関数
 */
TH1D* calculateTrackingEfficiency(TString particleType, TString fileName, 
                                  double jetPt_min, double jetPt_max,
                                  double cent_min, double cent_max) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }
    
    TString recoHistName = TString::Format("jet-shape-task/ptHistogram%s", particleType.Data());
    TString genHistName = TString::Format("jet-shape-task/ptGenerated%s", particleType.Data());
    
    THnSparseD *hRecoSparse = (THnSparseD*) file->Get(recoHistName);
    THnSparseD *hGenSparse = (THnSparseD*) file->Get(genHistName);
    
    if (!hRecoSparse || !hGenSparse) {
        std::cerr << "Error: Could not retrieve THnSparse for tracking efficiency" << std::endl;
        file->Close();
        return nullptr;
    }

    hRecoSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hGenSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hRecoSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);
    hGenSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);

    TH1D* hReco = hRecoSparse->Projection(0);
    TH1D* hGen = hGenSparse->Projection(0);

    hReco->SetDirectory(nullptr);
    hGen->SetDirectory(nullptr);
    file->Close();

    hReco->Sumw2();
    hGen->Sumw2();
    TH1D* hEfficiency = (TH1D*)hReco->Clone(TString::Format("hEfficiency_Track_%s", particleType.Data()));
    hEfficiency->SetDirectory(nullptr);
    hEfficiency->Divide(hReco, hGen, 1.0, 1.0, "B");
    
    delete hReco; delete hGen;
    delete hRecoSparse; delete hGenSparse; 
    
    return hEfficiency;
}

/**
 * @brief TOFマッチング効率を計算する関数
 */
TH1D* calculateTofEfficiency(TString particleType, TString fileName,
                             double jetPt_min, double jetPt_max,
                             double cent_min, double cent_max) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }
    
    TString tofHistName = TString::Format("jet-shape-task/ptHistogram%sTof", particleType.Data());
    TString recoHistName = TString::Format("jet-shape-task/ptHistogram%s", particleType.Data());
    
    THnSparseD *hTofSparse = (THnSparseD*) file->Get(tofHistName);
    THnSparseD *hRecoSparse = (THnSparseD*) file->Get(recoHistName);
    
    if (!hTofSparse || !hRecoSparse) {
        std::cerr << "Error: Could not retrieve THnSparse for TOF efficiency" << std::endl;
        file->Close();
        return nullptr;
    }

    hTofSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hRecoSparse->GetAxis(1)->SetRangeUser(jetPt_min, jetPt_max);
    hTofSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);
    hRecoSparse->GetAxis(2)->SetRangeUser(cent_min, cent_max);

    TH1D* hTof = hTofSparse->Projection(0);
    TH1D* hReco = hRecoSparse->Projection(0);
    
    hTof->SetDirectory(nullptr);
    hReco->SetDirectory(nullptr);
    file->Close();

    hTof->Sumw2();
    hReco->Sumw2();
    TH1D* hEfficiency = (TH1D*)hTof->Clone(TString::Format("hEfficiency_Tof_%s", particleType.Data()));
    hEfficiency->SetDirectory(nullptr);
    hEfficiency->Divide(hTof, hReco, 1.0, 1.0, "B");
    
    delete hTof; delete hReco;
    delete hTofSparse; delete hRecoSparse;

    return hEfficiency;
}

TPaveText* createPaveText(double x1, double y1, double x2, double y2, 
                          double jetPtMin, double jetPtMax,
                          double centMin, double centMax) {
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("#bf{This thesis}");
    pt->AddText(TString::Format("%.0f-%.0f%% Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV",centMin, centMax));
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", jetPtMin, jetPtMax));
    pt->AddText("ch-particle jet, Anti-#it{k}_{T}, #it{R} = 0.3");
    pt->AddText("#it{p}_{T,track} > 0.15 GeV/#it{c}");
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    pt->AddText("|#it{#eta}_{track}| < 0.9");
    return pt;
}

// ============================================================================================
// Main Function
// ============================================================================================

void calculateProtonPionRatio() {
    Int_t oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;

    gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const TString kDataFileName = "AnalysisResults-108.root";
    const TString kMcFileName   = "AnalysisResults-91.root";

    // --- Cuts ---
    const double kJetPtMin = 60.0;
    const double kJetPtMax = 250.0;
    const double kCentMin = 0.0;
    const double kCentMax = 10.0;

    // --- MC Efficiency Cuts ---
    const double kMcJetPtMin = 20.0;
    const double kMcJetPtMax = 250.0;
    const double kMcCentMin = 0.0;
    const double kMcCentMax = 100.0;

    std::vector<std::pair<double, double>> r_ranges;
    for (double r_min = 0.0; r_min < 0.7; r_min += 0.1) {
        r_ranges.push_back({r_min, r_min + 0.1});
    }

    std::cout << "\n--- Calculating Tracking & TOF Efficiencies (using " << kMcFileName << ") ---" << std::endl;
    
    TH1D* hEffTrackProton = calculateTrackingEfficiency("Proton", kMcFileName, kMcJetPtMin, kMcJetPtMax, kMcCentMin, kMcCentMax);
    TH1D* hEffTrackPion   = calculateTrackingEfficiency("Pion",   kMcFileName, kMcJetPtMin, kMcJetPtMax, kMcCentMin, kMcCentMax);
    TH1D* hEffTofProton   = calculateTofEfficiency("Proton", kMcFileName, kMcJetPtMin, kMcJetPtMax, kMcCentMin, kMcCentMax);
    TH1D* hEffTofPion     = calculateTofEfficiency("Pion",   kMcFileName, kMcJetPtMin, kMcJetPtMax, kMcCentMin, kMcCentMax);
    
    if (!hEffTrackProton || !hEffTrackPion || !hEffTofProton || !hEffTofPion) {
        std::cerr << "Efficiency calculation failed." << std::endl;
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        return;
    }

    // 全効率 = Tracking Efficiency * TOF Matching Efficiency
    TH1D* hEffTotalProton = (TH1D*)hEffTrackProton->Clone("hEffTotalProton");
    hEffTotalProton->Multiply(hEffTofProton);
    hEffTotalProton->SetDirectory(nullptr);

    TH1D* hEffTotalPion = (TH1D*)hEffTrackPion->Clone("hEffTotalPion");
    hEffTotalPion->Multiply(hEffTofPion);
    hEffTotalPion->SetDirectory(nullptr);

    std::vector<TH1D*> h_ratios;

    for (const auto& range : r_ranges) {
        double r_min = range.first;
        double r_max = range.second;

        // --- 1. Purity & PID Efficiency Map 計算 (Code Aロジック) ---
        // 戻り値は {PurityMap, PidEffMap} のペア
        auto correctionsPr = calculatePidCorrections("proton", kDataFileName, r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        auto correctionsPi = calculatePidCorrections("pion",   kDataFileName, r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);

        // --- 2. スペクトル作成 (Code Aロジック: Purity/PidEff補正) ---
        // この時点で PID関連の補正は全て完了したスペクトルが返る
        TH1D* hPtProton_pidCorr = createCorrectedJetSpectrum(kDataFileName, "jetpVsPtForPr", 
                                                             correctionsPr.first, correctionsPr.second, 
                                                             r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        TH1D* hPtPion_pidCorr   = createCorrectedJetSpectrum(kDataFileName, "jetpVsPtForPi", 
                                                             correctionsPi.first, correctionsPi.second, 
                                                             r_min, r_max, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
        
        if (!hPtProton_pidCorr || !hPtPion_pidCorr) continue;

        // --- 3. Tracking/TOF 効率補正 (残りの効率) ---
        TH1D* hPtProton_final = (TH1D*)hPtProton_pidCorr->Clone(TString::Format("hPtProton_final_r_%.1f_%.1f", r_min, r_max));
        hPtProton_final->Divide(hEffTotalProton);

        TH1D* hPtPion_final = (TH1D*)hPtPion_pidCorr->Clone(TString::Format("hPtPion_final_r_%.1f_%.1f", r_min, r_max));
        hPtPion_final->Divide(hEffTotalPion);

        // --- 4. Ratio Calculation (Rebin) ---
        Double_t xbins[] = {
            0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
            2.5, 3.0, 3.5, 4.0, 5.0, 6.0
        };
        int nbins = sizeof(xbins) / sizeof(Double_t) - 1;

        TH1D* hPtProton_rebin = (TH1D*)hPtProton_final->Rebin(nbins, TString::Format("hPtProton_rebin_r_%.1f", r_min), xbins);
        TH1D* hPtPion_rebin   = (TH1D*)hPtPion_final->Rebin(nbins, TString::Format("hPtPion_rebin_r_%.1f", r_min), xbins);

        TH1D* hRatio = (TH1D*)hPtProton_rebin->Clone(TString::Format("hRatio_p_pi_r_%.1f_%.1f", r_min, r_max));
        hRatio->Divide(hPtPion_rebin);
        hRatio->SetDirectory(nullptr);
        h_ratios.push_back(hRatio);

        // Cleanup
        delete hPtProton_pidCorr; delete hPtPion_pidCorr;
        delete hPtProton_final;   delete hPtPion_final;
        delete hPtProton_rebin;   delete hPtPion_rebin;
    }

    delete hEffTrackProton; delete hEffTrackPion;
    delete hEffTofProton; delete hEffTofPion;
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

    for (size_t i = 0; i < h_ratios.size(); ++i) {
        TH1D* hRatio = h_ratios[i];
        int color_idx = i % colors.size();
        hRatio->SetLineColor(colors[color_idx]);
        hRatio->SetMarkerColor(colors[color_idx]);
        hRatio->SetMarkerStyle(markers[i % markers.size()]);

        if (i == 0) {
            hRatio->SetTitle("p/#pi Ratio");
            hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
            hRatio->GetYaxis()->SetRangeUser(0.0, 1.1);
            hRatio->Draw("E1");
        } else {
            hRatio->Draw("E1 SAME");
        }
        legend->AddEntry(hRatio, TString::Format("%.1f < #it{r} < %.1f", r_ranges[i].first, r_ranges[i].second), "lpe");
    }
    legend->Draw();

    TPaveText* pave_ratio = createPaveText(0.15, 0.7, 0.5, 0.9, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    pave_ratio->Draw();

    gErrorIgnoreLevel = oldErrorIgnoreLevel;
}