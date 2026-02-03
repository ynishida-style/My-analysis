#include <THnSparse.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TError.h> // エラーレベル制御用
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

// ============================================================================================
// Fitting Function (Provided by User)
// ============================================================================================

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // Grid Search Ranges
    std::vector<double> p4_range; 
    for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 

    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    bool isHighPt = (p_min > 4.1);

    // --- Step 1: Grid Search ---
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                // "Q0" オプションを追加してQuietモードにする
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                
                // Proton (Signal)
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                tempFit->SetParameter(0, pr_init_height); 
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                // Backgrounds (Kaon, Pion, Electron)
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;
                
                tempFit->SetParameter(3, kaon_init_height); 
                tempFit->SetParameter(4, p4_init); 
                tempFit->SetParameter(5, 1.5);

                tempFit->SetParameter(6, peak_height * 0.7); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 2.0);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 2.5);

                // "Q"オプション: Quiet mode (最小限の出力)
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
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    // Limits application
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

    // "Q"オプションを追加
    hist->Fit(finalFit, "QRB0");
    return finalFit;
}

// ============================================================================================
// Main Analysis Function
// ============================================================================================

void AnalyzeJetShapeProtonFit() {
    // --- 警告メッセージの抑制設定 ---
    Int_t originalErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError; 

    // Styling
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    
    // File I/O
    const TString fileName = "AnalysisResults-106.root"; 
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) { 
        std::cerr << "Error: Cannot open file " << fileName << std::endl; 
        gErrorIgnoreLevel = originalErrorLevel; 
        return; 
    }

    const TString histName = "jet-shape-task/jetTpcTofPr";
    THnSparseD *hSparse = (THnSparseD*)file->Get(histName);
    if (!hSparse) { 
        std::cerr << "Error: Histogram " << histName << " not found." << std::endl; 
        gErrorIgnoreLevel = originalErrorLevel;
        return; 
    }

    // ========================================================================================
    // 1. Define Cuts (変数で管理して表示と実際のカットを一致させる)
    // ========================================================================================
    
    // Distance R range
    double r_min = 0.0;
    double r_max = 0.1; // コード内の 0.4001, 0.4999 に対応
    
    // Jet pT range
    double jetPt_min = 60.0;
    double jetPt_max = 199.0;

    // Centrality range
    double cent_min = 0.0;
    double cent_max = 10.0;

    // Apply Cuts
    // わずかな誤差で境界を含める/含めないを制御するため微小値を加減算
    hSparse->GetAxis(2)->SetRangeUser(r_min + 0.0001, r_max - 0.0001); 
    hSparse->GetAxis(3)->SetRangeUser(jetPt_min, jetPt_max);
    hSparse->GetAxis(4)->SetRangeUser(cent_min, cent_max);

    std::cout << ">> Settings Applied:" << std::endl;
    std::cout << "   R (Distance): " << r_min << " - " << r_max << std::endl;
    std::cout << "   Jet pT:       " << jetPt_min << " - " << jetPt_max << " GeV/c" << std::endl;
    std::cout << "   Centrality:   " << cent_min << " - " << cent_max << " %" << std::endl;

    // ========================================================================================
    // Prepare Canvas
    // ========================================================================================

    TCanvas *cOut = new TCanvas("cOut", "Proton Fit", 1600, 1000); // 縦を少し広げました
    // 上部に情報を書くスペースを空けるため、マージンを調整するか、単に上書きします
    cOut->Divide(4, 3, 0.01, 0.01); 

    int pad_idx = 1;
    
    // ... (Purity Calculation Loop - 変更なし) ...
    std::cout << "\n--- Purity Calculation Results ---" << std::endl;
    for (double p_min = 3.0; p_min < 7.0; p_min += 0.4) { 
        if (pad_idx > 12) break;

        double p_max = p_min + 0.4;
        cOut->cd(pad_idx);
        
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->Rebin(2);
        h_nSigma->SetName(TString::Format("h_nSigma_Jet_%.1f_%.1f", p_min, p_max));
        
        if (h_nSigma->GetEntries() < 100) { delete h_nSigma; continue; }

        TF1* fitResult = optimizeAndFitProton(h_nSigma, p_min, p_max);

        // --- Purity Logic ---
        TF1* fSigOnly = new TF1("fSigOnly", "gaus", -10, 10);
        fSigOnly->SetParameters(fitResult->GetParameter(0), fitResult->GetParameter(1), fitResult->GetParameter(2));
        double signal_min_cut = -3.0; 
        double signal_max_cut = 0.5;
        double sig_integral = fSigOnly->Integral(signal_min_cut, signal_max_cut);
        double total_integral = fitResult->Integral(signal_min_cut, signal_max_cut);
        double purity = (total_integral > 0) ? sig_integral / total_integral : 0.0;
        
        printf("  %4.1f   |  %4.1f   |  %6.3f\n", p_min, p_max, purity);

        // --- Draw Plot ---
        h_nSigma->SetTitle(""); // 個別のタイトルは消す
        h_nSigma->GetXaxis()->SetTitle("TPC n#sigma_{p}");
        h_nSigma->GetXaxis()->SetRangeUser(-6, 8); 
        h_nSigma->SetMarkerStyle(20);
        h_nSigma->SetMarkerSize(0.6);
        h_nSigma->Draw("PE");
        fitResult->Draw("SAME");
        fSigOnly->SetLineColor(kRed);
        fSigOnly->SetLineStyle(1);
        fSigOnly->Draw("SAME");

        TF1* fBgKaon = new TF1("fBgKaon", "gaus", -10, 10); 
        fBgKaon->SetParameters(fitResult->GetParameter(3), fitResult->GetParameter(4), fitResult->GetParameter(5));
        fBgKaon->SetLineColor(kMagenta); fBgKaon->SetLineStyle(2); fBgKaon->Draw("SAME");

        // 個別のパッド内の情報（運動量範囲とPurity）
        TPaveText *pt = new TPaveText(0.55, 0.55, 0.88, 0.85, "NDC");
        pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextAlign(12); pt->SetTextSize(0.045);
        pt->AddText(TString::Format("#it{p}: %.1f - %.1f GeV/#it{c}", p_min, p_max));
        pt->AddText(TString::Format("Purity: %.3f", purity));
        pt->Draw();

        pad_idx++;
    }
    
    // ========================================================================================
    // ★ Global Info Display (ここを追加)
    // ========================================================================================
    cOut->cd(0); // メインのキャンバスに戻る
    
    // 上部中央に配置するテキストボックス
    TPaveText *globalInfo = new TPaveText(0.2, 0.96, 0.8, 0.995, "NDC");
    globalInfo->SetFillColor(kWhite); // 背景を白にして下の文字と重ならないようにする
    globalInfo->SetBorderSize(1);
    globalInfo->SetTextAlign(22); // 中央揃え
    globalInfo->SetTextSize(0.025); // 全体が見えるサイズ
    
    // 定義した変数を使って文字列を作成
    TString infoStr = TString::Format("Cut Conditions:  R: %.1f - %.1f  |  Jet #it{p}_{T}: %.0f - %.0f GeV/#it{c}  |  Cent: %.0f - %.0f %%", 
                                      r_min, r_max, jetPt_min, jetPt_max, cent_min, cent_max);
    
    globalInfo->AddText(infoStr);
    globalInfo->Draw();

    cOut->Update();
    gErrorIgnoreLevel = originalErrorLevel;
}