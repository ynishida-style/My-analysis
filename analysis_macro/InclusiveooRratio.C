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
#include <utility> // for std::pair

// ============================================================================================
// Helper Functions: Fitting
// (Parameters and Logic synced with the provided reference code)
// ============================================================================================
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // --- Grid Search Ranges ---
    // Kaon (Param 3-5)
    std::vector<double> p4_range; 
    for (double p4 = 0.8; p4 <= 1.5; p4 += 0.1) p4_range.push_back(p4); 

    // Pion (Param 6-8)
    // 高pT (>5GeV) ではPionがProtonに近づくため、探索範囲を下（3.0付近）まで広げる
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
                
                // 初期値の工夫
                double pr_init_h = peak_height * 0.7;
                double ka_init_h = peak_height * 0.3;
                
                if (isMidPt) {
                     pr_init_h = peak_height * 0.5; // Kaonに譲るため低めに
                     ka_init_h = peak_height * 0.5;
                }

                double params[12] = {
                    pr_init_h, 0.0, 1.0,            // Proton 
                    ka_init_h, p4_init, 1.0,        // Kaon 
                    peak_height * 0.1, p7_init, 1.5,  // Pion
                    peak_height * 0.05, p10_init, 2.0 // Electron
                };
                tempFit->SetParameters(params);
                
                // Grid Search中の制約
                tempFit->SetParLimits(0, 0, peak_height * 1.5);
                tempFit->SetParLimits(1, -0.2, 0.2); 
                tempFit->SetParLimits(2, 0.8, 1.1); // 太りすぎ防止

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
    
    // =========================================================================================
    // ★ Parameter Limits Tuning (Overcounting対策)
    // =========================================================================================

    // --- 1. Proton (Param 0: Height, 1: Mean, 2: Sigma) ---
    // MidPt (3.0-4.8) では Protonがピークの100%になるのを禁止し、Kaonに配分を強制する
    if (isMidPt) {
        finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 0.90); 
    } else {
        finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 1.5); 
    }
    
    finalFit->SetParLimits(1, -0.1, 0.1); 
    // Sigma上限: ここを緩めるとProtonがKaonを食べてPurityが偽増大する
    finalFit->SetParLimits(2, 0.95, 1.01); 

    // --- 2. Kaon (Param 3: Height, 4: Mean, 5: Sigma) ---
    if (isMidPt) {
        // Kaonが消えるのを防ぐため、下限を高く設定
        finalFit->SetParLimits(3, peak_height * 0.2, peak_height * 0.9);
    } else {
        finalFit->SetParLimits(3, 0.0, peak_height * 1.5);
    }
    
    finalFit->SetParLimits(4, 0.6, 1.8); 
    finalFit->SetParLimits(5, 0.8, 1.3); 

    // --- 3. Pion (Param 6-8) ---
    finalFit->SetParLimits(6, 0.0, peak_height * 1.2);
    // 高pTではPionがProtonに近づくため、下限を緩和
    if (isHighPt) {
        finalFit->SetParLimits(7, 2.0, 7.0); 
    } else {
        finalFit->SetParLimits(7, 3.0, 7.0); 
    }
    finalFit->SetParLimits(8, 0.5, 2.5);

    // --- 4. Electron (Param 9-11) ---
    finalFit->SetParLimits(9, 0.0, peak_height * 0.5);
    finalFit->SetParLimits(10, 7.0, 15.0);
    finalFit->SetParLimits(11, 0.6, 2.0);

    hist->Fit(finalFit, "QNRB0"); 

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}



/**
* @brief パイオン用のグリッドサーチと最終フィットを行う関数
* (パラメータ範囲・Limitをご指定のコードに合わせて更新)
*/
TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    std::cout << "--- Starting parameter optimization for PION in p bin " << p_min << "-" << p_max << " ---" << std::endl;
    
    // Ranges (Code B仕様)
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
    std::cout << "Grid search finished. Best preliminary chi2/ndf: " << best_chi2ndf << std::endl;

    // --- Step 2: Final Fit with Best Parameters (Updated Limits) ---
    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed); // Code A logic uses returned fit for calc, color matters less but keeping logic
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

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

// ============================================================================================
// Core Logic: Purity AND PID Efficiency Calculation
// ============================================================================================

// 戻り値: <PurityMap, EfficiencyMap> のペア
std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(TString fileName, TString particleType, double centMin, double centMax) {
    std::map<double, double> purityMap;
    std::map<double, double> effMap;
    
    TString sparseName;
    double sigMin, sigMax; // カット範囲

    if (particleType == "proton") {
        sparseName = "tpcTofPr"; 
        sigMin = -3.5; sigMax = 0.5; // Protonのカット範囲
    } else {
        sparseName = "tpcTofPi"; 
        sigMin = -0.5; sigMax = 3.5; // Pionのカット範囲
    }

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return {purityMap, effMap};
    }
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { 
        std::cerr << "Error: Cannot find histogram " << sparseName << " in " << fileName << std::endl;
        f->Close(); 
        return {purityMap, effMap}; 
    }

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax); 

    double p_start = 0.6; 
    double p_end = 8.0; 
    double p_step = 0.2; 

    std::cout << "--- Calculating Purity & PID Efficiency for " << particleType << " ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        // 低pT (TOFが効いて分離が良い領域) は1.0と仮定
        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); // ideal
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            double idealEff   = idealInCut / idealTotal;

            effMap[p_center] = idealEff; 
            continue;
        }

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        if (h_nSigma->GetEntries() < 50) { delete h_nSigma; continue; }

        TF1* fit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        // Signal成分のみ抽出
        TF1* fSig = new TF1("fSig", "gaus", -10, 10);
        fSig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));

        // 積分計算
        double intSigInCut = fSig->Integral(sigMin, sigMax); // カット内の信号
        double intSigTotal = fSig->Integral(-10, 10);        // 全信号 (PID効率の分母)
        double intTotInCut = fit->Integral(sigMin, sigMax);  // カット内の全粒子 (Purityの分母)

        // 1. Purity Calculation
        double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
        
        // 2. PID Efficiency Calculation (カットで捨てた分の補正)
        double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;
        
        // ガード
        if (purity > 1.0) purity = 1.0;
        if (purity < 0.0) purity = 0.0;
        if (pidEff > 1.0) pidEff = 1.0;
        if (pidEff < 0.0) pidEff = 0.0; // 異常値回避

        purityMap[p_center] = purity;
        effMap[p_center]    = pidEff;

        printf(" p=%.2f: Purity=%.3f, PID Eff=%.3f\n", p_center, purity, pidEff);

        delete h_nSigma; delete fit; delete fSig;
    }
    f->Close();
    return {purityMap, effMap};
}

// PurityとPID効率の両方を適用してスペクトルを補正する関数
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
    
    TH1D *hPt = h2->ProjectionX(TString::Format("hPt_%s", histName.Data()));
    hPt->Reset();
    hPt->SetDirectory(0);

    for (int by = 1; by <= h2->GetNbinsY(); ++by) {
        double p_val = h2->GetYaxis()->GetBinCenter(by);
        
        if (p_val > 8.0) continue;

        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0; 
        
        // PurityとEffMapの検索 (同じキーを持っている前提)
        for (auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_val) < min_diff) {
                min_diff = std::abs(map_p - p_val);
                purity = map_pur;
                pidEff = pidEffMap[map_p]; // 対応する効率を取得
            }
        }

        //purity = 1.0;
        
        // マッチするpTビンが遠すぎる場合の処理
        if (min_diff > 0.3) {
             if (p_val > 3.0) continue; 
        }

        // PID Efficiencyが極端に小さい(0含む)場合の保護
        if (pidEff < 0.01) pidEff = 1.0; 

        // 補正係数: Purity / Efficiency
        double totalCorrectionFactor = purity / pidEff;

        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double content = h2->GetBinContent(bx, by);
            double error   = h2->GetBinError(bx, by);
            
            if (content <= 0) continue;

            // 補正適用
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
    
    // Tof付き、Inclusive
    THnSparseD *hRecoSp = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptHistogram%sTof", partName.Data()));
    THnSparseD *hGenSp  = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptGenerated%s", partName.Data()));
    
    if (!hRecoSp || !hGenSp) {
        std::cout << "Error loading Efficiency Histograms for " << partName << std::endl;
        if(f) f->Close();
        return nullptr;
    }

    hRecoSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    
    TH1D *hReco = hRecoSp->Projection(0);
    TH1D *hGen  = hGenSp->Projection(0);
    
    TH1D *hEff = (TH1D*)hReco->Clone(TString::Format("eff_%s", partName.Data()));
    hEff->Divide(hReco, hGen, 1, 1, "B");
    hEff->SetDirectory(0);
    f->Close();
    return hEff;
}

// ============================================================================================
// Main
// ============================================================================================

void PlotInclusiveRatioFinal() {
    // 統計情報の表示をオフ
    gStyle->SetOptStat(0);
    // フォントのデフォルトを標準（42）に設定して太字化を防止
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    
    const TString dataFile = "AnalysisResults-100.root";
    const TString mcFile   = "AnalysisResults-91.root";

    const double dataCentMin = 0.0;
    const double dataCentMax = 5.0;
    const double mcCentMin = 0.0;
    const double mcCentMax = 100.0;

    auto mapsPr = getPidCorrectionMaps(dataFile, "proton", dataCentMin, dataCentMax);
    auto mapsPi = getPidCorrectionMaps(dataFile, "pion", dataCentMin, dataCentMax);
    
    std::map<double, double> purMapPr = mapsPr.first;
    std::map<double, double> effPidPr = mapsPr.second;
    std::map<double, double> purMapPi = mapsPi.first;
    std::map<double, double> effPidPi = mapsPi.second;

    TH1D* hPrRaw = createCorrectedSpectrum(dataFile, "pVsPtForPr", purMapPr, effPidPr, dataCentMin, dataCentMax);
    TH1D* hPiRaw = createCorrectedSpectrum(dataFile, "pVsPtForPi", purMapPi, effPidPi, dataCentMin, dataCentMax);

    TH1D* hEffPr = calculateEfficiency(mcFile, "Proton", mcCentMin, mcCentMax);
    TH1D* hEffPi = calculateEfficiency(mcFile, "Pion", mcCentMin, mcCentMax);

    if (!hPrRaw || !hPiRaw || !hEffPr || !hEffPi) {
        std::cerr << "Error creating histograms." << std::endl;
        return;
    }

    TH1D* hPrCorr = (TH1D*)hPrRaw->Clone("hPrCorr");
    hPrCorr->Divide(hEffPr);
    TH1D* hPiCorr = (TH1D*)hPiRaw->Clone("hPiCorr");
    hPiCorr->Divide(hEffPi);

    // 比の計算（ご指定のコードに合わせ Rebin(2) を適用）
    hPrCorr->Rebin(2);
    hPiCorr->Rebin(2);
    TH1D* hRatio = (TH1D*)hPrCorr->Clone("hRatio");
    hRatio->Divide(hPiCorr);

    // --- Plotting ---
    hRatio->GetXaxis()->SetRangeUser(0.0, 6.0);

    // Canvas 1: Spectra (見た目の調整)
    TCanvas *c1 = new TCanvas("cSpectra", "Spectra", 800, 600);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    
    // 見た目の設定
    hPiCorr->SetMarkerStyle(20); // kFullCircle
    hPiCorr->SetMarkerColor(kBlue + 1); 
    hPiCorr->SetLineColor(kBlue + 1);
    
    hPrCorr->SetMarkerStyle(21); // kFullSquare
    hPrCorr->SetMarkerColor(kRed + 1);  
    hPrCorr->SetLineColor(kRed + 1);

    hPiCorr->SetTitle("Inclusive p_{T} Spectra");
    hPiCorr->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}");
    hPiCorr->Draw("PE");
    hPrCorr->Draw("PE SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88); 
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hPiCorr, "Pions", "lpe");
    leg->AddEntry(hPrCorr, "Protons", "lpe");
    leg->Draw();

    // Canvas 2: Ratio
    TCanvas *c2 = new TCanvas("cRatio", "Ratio", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetGridy();

    hRatio->SetTitle("p/#pi Ratio");
    hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitleOffset(1.2);
    
    hRatio->SetMarkerStyle(34); // kFullCross
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    
    hRatio->SetMinimum(0.0);
    hRatio->SetMaximum(1.1); 

    hRatio->Draw("E1");

    TPaveText *pt = new TPaveText(0.15, 0.75, 0.5, 0.9, "NDC");
    pt->SetFillStyle(0); 
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42); 
    pt->SetTextSize(0.035);
    
    pt->AddText("#bf{This thesis}");
    pt->AddText("|#eta_{track}| < 0.9");
    pt->AddText(TString::Format("%.0f-%.0f%% O-O inclusive", dataCentMin, dataCentMax));
    pt->AddText("#sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->Draw();

    c1->SaveAs("InclusiveSpectra_FullyCorr.pdf");
    c2->SaveAs("InclusiveRatio_FullyCorr.pdf");
}