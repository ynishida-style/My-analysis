#include <THnSparse.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

// ============================================================================================
// Fitting Function (Grid Search + 4-Gaussian Fit)
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




// ============================================================================================
// Main Function
// ============================================================================================
void CheckProtonFitGridSearch() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // ※ファイル名は環境に合わせて調整してください
    const TString fileName = "AnalysisResults-100.root";
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) { std::cerr << "File Error" << std::endl; return; }

    const TString histName = "jet-shape-task/tpcTofPr";
    THnSparseD *hSparse = (THnSparseD*)file->Get(histName);
    if (!hSparse) { std::cerr << "Hist Error: " << histName << std::endl; return; }

    // Rなどの他軸のカット（必要に応じて）
    hSparse->GetAxis(2)->SetRangeUser(0.0, 10.0);

    // キャンバス作成
    TCanvas *cPr = new TCanvas("cPr", "Proton Fit Grid Search 3.0-5.4 GeV/c", 1600, 1200);
    cPr->Divide(4, 3);

    int pad_idx = 1;
    // ループ範囲: 3.0 GeV/c からスタート
    for (double p_min = 4.0; p_min < 6.0; p_min += 0.2) {
        double p_max = p_min + 0.2;
        
        if (pad_idx > 12) break; // 12枚まで描画
        cPr->cd(pad_idx);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);

        // 射影
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        TString nameH = TString::Format("h_nSigma_Pr_%.1f", p_min);
        h_nSigma->SetName(nameH);
        h_nSigma->SetTitle(TString::Format("%.1f < p < %.1f GeV/c", p_min, p_max));

        if (h_nSigma->GetEntries() < 50) { 
            delete h_nSigma; 
            continue; // エントリが少なすぎる場合はスキップ（padは進めない）
        }

        // --- フィッティング実行 ---
        TF1* fitResult = optimizeAndFitProton(h_nSigma, p_min, p_max);

        // --- Purity Calculation ---
        // Signal: 1番目のガウス関数 [0-2] (Proton)
        TF1* fSigOnly = new TF1("fSigOnly", "gaus", -10, 10);
        fSigOnly->SetParameters(fitResult->GetParameter(0), fitResult->GetParameter(1), fitResult->GetParameter(2));
        fSigOnly->SetLineColor(kRed); 
        fSigOnly->SetLineWidth(2);
        
        // 積分範囲 (-3.5 ~ 0.5 sigma などを想定)
        double signal_min_cut = -3.0; // Proton cut (Typically |nSigma| < 3 or similar)
        double signal_max_cut = 3.0;
        
        // Purity計算 (積分範囲は解析の定義に合わせて調整してください)
        // ここでは全範囲積分での比率を簡易的に計算します
        double sig_integral = fSigOnly->Integral(-10, 10);
        double total_integral = fitResult->Integral(-10, 10); 
        
        // もし特定のCut内でのPurityが見たい場合は積分範囲を変更してください
        // double sig_integral = fSigOnly->Integral(-3, 3);
        // double total_integral = fitResult->Integral(-3, 3);

        double purity = (total_integral > 0) ? sig_integral / total_integral : 0.0;
        purity = std::max(0.0, std::min(1.0, purity));

        // --- Drawing ---
        h_nSigma->GetXaxis()->SetTitle("TPC n#sigma_{p}"); 
        h_nSigma->GetYaxis()->SetTitle("Counts");
        h_nSigma->GetXaxis()->SetRangeUser(-5, 10); // 表示範囲を少しズーム
        h_nSigma->SetMarkerStyle(20); 
        h_nSigma->SetMarkerSize(0.6);
        h_nSigma->Draw("PE");

        fitResult->Draw("SAME");
        fSigOnly->Draw("SAME"); // Proton成分(赤)

        // 背景成分 (Kaon)
        TF1* fBgKaon = new TF1("fBgKaon", "gaus", -10, 10);
        fBgKaon->SetParameters(fitResult->GetParameter(3), fitResult->GetParameter(4), fitResult->GetParameter(5));
        fBgKaon->SetLineColor(kMagenta); 
        fBgKaon->SetLineStyle(2); 
        fBgKaon->SetLineWidth(2);
        fBgKaon->Draw("SAME");

        // 背景成分 (Pion)
        TF1* fBgPion = new TF1("fBgPion", "gaus", -10, 10);
        fBgPion->SetParameters(fitResult->GetParameter(6), fitResult->GetParameter(7), fitResult->GetParameter(8));
        fBgPion->SetLineColor(kBlue); 
        fBgPion->SetLineStyle(2); 
        fBgPion->Draw("SAME");

        // 背景成分 (Electron)
        TF1* fBgEle = new TF1("fBgEle", "gaus", -10, 10);
        fBgEle->SetParameters(fitResult->GetParameter(9), fitResult->GetParameter(10), fitResult->GetParameter(11));
        fBgEle->SetLineColor(kGreen+2); 
        fBgEle->SetLineStyle(2); 
        fBgEle->Draw("SAME");

        // テキスト情報の表示
        TPaveText *pt = new TPaveText(0.55, 0.55, 0.88, 0.88, "NDC");
        pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextAlign(12);
        pt->SetTextSize(0.04);
        pt->AddText(TString::Format("p: %.1f - %.1f", p_min, p_max));
        pt->AddText(TString::Format("Purity: %.3f", purity));
        
        double chi2 = fitResult->GetChisquare();
        int ndf = fitResult->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2/ndf : 0;
        pt->AddText(TString::Format("#chi^{2}/ndf: %.1f", chi2ndf));
        
        // デバッグ用にSigmaを表示 (ここが1.0付近になっているか確認)
        pt->AddText(TString::Format("#sigma_{p}: %.3f", fitResult->GetParameter(2)));
        pt->Draw();

        pad_idx++; 
    }
    cPr->Update();
}