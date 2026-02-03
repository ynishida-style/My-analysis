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
    std::vector<double> p4_range; for (double p4 = 0.8; p4 <= 1.4; p4 += 0.2) p4_range.push_back(p4); 
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();
    bool isLowPt = (p_min < 3.9); 
    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                double ka_init_val = isLowPt ? peak_height * 0.3 : peak_height * 0.4;
                double pr_init_val = isLowPt ? peak_height * 0.7 : peak_height * 0.8;
                double params[12] = {
                    pr_init_val, 0.0, 1.0, ka_init_val, p4_init, 1.0,
                    peak_height * 0.5, p7_init, 1.5, peak_height * 0.05, p10_init, 2.0
                };
                tempFit->SetParameters(params);
                tempFit->SetParLimits(0, 0, peak_height * 1.5);
                tempFit->SetParLimits(1, -0.5, 0.5); 
                tempFit->SetParLimits(3, 0, peak_height * 1.5);
                tempFit->SetParLimits(4, 0.5, 2.0);  
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
    
    if (isLowPt) { finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 0.98); } 
    else { finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.5); }
    finalFit->SetParLimits(1, -0.03, 0.05); 
    finalFit->SetParLimits(2, 0.8, 1.15);

    if (isLowPt) { finalFit->SetParLimits(3, peak_height * 0.01, peak_height * 1.0); } 
    else { finalFit->SetParLimits(3, peak_height * 0.2, peak_height * 1.5); }
    finalFit->SetParLimits(4, 0.8, 1.6); 
    finalFit->SetParLimits(5, 0.5, 1.2); 

    finalFit->SetParLimits(6, 0.0, peak_height * 1.2); 
    finalFit->SetParLimits(7, 3.0, 6.0); 
    finalFit->SetParLimits(8, 0.5, 2.0);

    finalFit->SetParLimits(9, 0.0, peak_height * 0.5); 
    finalFit->SetParLimits(10, 7.0, 12.0); 
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
    const TString fileName = "AnalysisResults-98.root";
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) { std::cerr << "File Error" << std::endl; return; }

    const TString histName = "jet-shape-task/tpcTofPr";
    THnSparseD *hSparse = (THnSparseD*)file->Get(histName);
    if (!hSparse) { std::cerr << "Hist Error" << std::endl; return; }

    // Rなどの他軸のカット（必要に応じて）
    hSparse->GetAxis(2)->SetRangeUser(0.0, 10.0);

    TCanvas *cPr = new TCanvas("cPr", "Proton Fit Grid Search 3.0-5.0 GeV/c", 1600, 900);
    cPr->Divide(4, 3);

    int pad_idx = 1;
    for (double p_min = 3.0; p_min < 5.0; p_min += 0.2) {
        double p_max = p_min + 0.2;
        cPr->cd(pad_idx);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->SetName(TString::Format("h_nSigma_Pr_%.1f", p_min));

        if (h_nSigma->GetEntries() < 50) { 
            delete h_nSigma; 
            pad_idx++; 
            continue; 
        }

        // 新しいフィッティング関数の呼び出し
        TF1* fitResult = optimizeAndFitProton(h_nSigma, p_min, p_max);

        // --- Purity Calculation ---
        // Signal: 1番目のガウス関数 [0-2]
        TF1* fSigOnly = new TF1("fSigOnly", "gaus", -10, 10);
        fSigOnly->SetParameters(fitResult->GetParameter(0), fitResult->GetParameter(1), fitResult->GetParameter(2));
        
        double signal_min_cut = -3.5;
        double signal_max_cut = 0.5;
        double sig_integral = fSigOnly->Integral(signal_min_cut, signal_max_cut);
        double total_integral = fitResult->Integral(signal_min_cut, signal_max_cut);
        double purity = (total_integral > 0) ? sig_integral / total_integral : 0.0;
        purity = std::max(0.0, std::min(1.0, purity));

        // --- Drawing ---
        h_nSigma->GetXaxis()->SetTitle("TPC n#sigma_{p}"); 
        h_nSigma->GetYaxis()->SetTitle("Counts");
        h_nSigma->GetXaxis()->SetRangeUser(-10, 12);
        h_nSigma->SetMarkerStyle(20); 
        h_nSigma->SetMarkerSize(0.6);
        h_nSigma->Draw("PE");

        fitResult->Draw("SAME");

        // 各成分の描画 (背景成分)
        fSigOnly->SetLineColor(kRed); fSigOnly->SetLineWidth(2); fSigOnly->Draw("SAME");

        TF1* fBg1 = new TF1("fBg1", "gaus", -10, 10);
        fBg1->SetParameters(fitResult->GetParameter(3), fitResult->GetParameter(4), fitResult->GetParameter(5));
        fBg1->SetLineColor(kMagenta); fBg1->SetLineStyle(2); fBg1->Draw("SAME");

        TF1* fBg2 = new TF1("fBg2", "gaus", -10, 10);
        fBg2->SetParameters(fitResult->GetParameter(6), fitResult->GetParameter(7), fitResult->GetParameter(8));
        fBg2->SetLineColor(kBlue); fBg2->SetLineStyle(2); fBg2->Draw("SAME");

        TF1* fBg3 = new TF1("fBg3", "gaus", -10, 10);
        fBg3->SetParameters(fitResult->GetParameter(9), fitResult->GetParameter(10), fitResult->GetParameter(11));
        fBg3->SetLineColor(kGreen+2); fBg3->SetLineStyle(2); fBg3->Draw("SAME");

        // テキスト情報の表示
        TPaveText *pt = new TPaveText(0.55, 0.60, 0.88, 0.88, "NDC");
        pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextAlign(12);
        pt->SetTextSize(0.04);
        pt->AddText(TString::Format("p: %.1f - %.1f", p_min, p_max));
        pt->AddText(TString::Format("Purity: %.3f", purity));
        double chi2ndf = (fitResult->GetNDF() > 0) ? fitResult->GetChisquare()/fitResult->GetNDF() : 0;
        pt->AddText(TString::Format("#chi^{2}/ndf: %.1f", chi2ndf));
        pt->Draw();

        pad_idx++; 
        if (pad_idx > 12) break; 
    }
    cPr->Update();
}