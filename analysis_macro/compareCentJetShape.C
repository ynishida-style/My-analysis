#include <iostream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"

// ==========================================
// 共通: THnSparseからProfileを取得する関数 (変更なし)
// ==========================================
TProfile* getProfileFromTHn(TDirectory* dir, TString histName, double minCent, double maxCent, double minPt, double maxPt, TString suffix) {
    THnSparse* hn = (THnSparse*)dir->Get(histName);
    if (!hn) {
        std::cerr << "Error: Cannot find histogram " << histName << std::endl;
        return nullptr;
    }
    hn->GetAxis(2)->SetRangeUser(minPt, maxPt);
    hn->GetAxis(3)->SetRangeUser(minCent, maxCent);
    TH2D* h2 = (TH2D*)hn->Projection(1, 0);
    h2->SetName(Form("h2_%s", suffix.Data()));
    
    return h2->ProfileX(Form("hProf_%s", suffix.Data()));
}

// ==========================================
// 幾何学的補正 (Densityへの変換) (変更なし)
// ==========================================
void applyGeometryNorm(TH1D* h) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double rMin = h->GetBinLowEdge(i);
        double rMax = h->GetXaxis()->GetBinUpEdge(i);
        double dr   = rMax - rMin;
        double area = TMath::Pi() * (rMax * rMax - rMin * rMin);
        
        double scaleFactor = 0.0;
        if (area > 1e-6) scaleFactor = dr / area; 

        h->SetBinContent(i, h->GetBinContent(i) * scaleFactor);
        h->SetBinError(i, h->GetBinError(i) * scaleFactor);
    }
}

// ==========================================
// 実験情報表示 (Centrality引数を削除し汎用化)
// ==========================================
TPaveText* createPaveText(double x1, double y1, double x2, double y2, 
                          double jetPtMin, double jetPtMax) {
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12); // 左揃え・垂直中央
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("#bf{This thesis}");
    //pt->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    //pt->AddText("O-O #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText("ch-particle jet, Anti-k_{T} R=0.4");
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", jetPtMin, jetPtMax));
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    return pt;
}

// ==========================================
// Main Function
// ==========================================
void plotJetShapeOverlay()
{
    // スタイル設定
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // 1. 設定
    TString fileName = "AnalysisResults-106.root";
    TString dirName  = "jet-shape-task";
    
    double minPtCorr = 80.0; 
    double maxPtCorr = 200.0;

    // Centrality 定義
    const int nBins = 3;
    double centBins[nBins+1] = {0.0, 30.0, 60.0, 100.0};
    
    // プロット用設定
    int colors[] = {kBlack, kRed, kBlue};
    int markers[] = {20, 21, 22}; // Circle, Square, Triangle

    // 2. 読み込み
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: File open failed." << std::endl;
        return;
    }
    TDirectory* dir = (TDirectory*)file->Get(dirName);
    if (!dir) {
        std::cerr << "Error: Directory not found." << std::endl;
        return;
    }

    // 結果格納用ベクター
    std::vector<TH1D*> hResults;

    // 3. ループ処理 (0-30, 30-60, 60-100)
    for (int i = 0; i < nBins; ++i) {
        double minCent = centBins[i];
        double maxCent = centBins[i+1];
        TString suffix = TString::Format("C%.0f_%.0f", minCent, maxCent);

        std::cout << "Processing Centrality: " << minCent << " - " << maxCent << "%" << std::endl;

        // Profile取得
        TProfile* pSig = getProfileFromTHn(dir, "ptSum",    minCent, maxCent, minPtCorr, maxPtCorr, "Sig_" + suffix);
        TProfile* pBg1 = getProfileFromTHn(dir, "ptSumBg1", minCent, maxCent, minPtCorr, maxPtCorr, "Bg1_" + suffix);
        TProfile* pBg2 = getProfileFromTHn(dir, "ptSumBg2", minCent, maxCent, minPtCorr, maxPtCorr, "Bg2_" + suffix);

        if (!pSig || !pBg1 || !pBg2) {
            std::cerr << "Skipping bin " << i << " due to missing hist." << std::endl;
            continue;
        }

        // 射影
        TH1D* hSig = pSig->ProjectionX("hSig_" + suffix); 
        TH1D* hBg1 = pBg1->ProjectionX("hBg1_" + suffix);
        TH1D* hBg2 = pBg2->ProjectionX("hBg2_" + suffix);

        // Bg Average計算 (指定コードと同一手法)
        TH1D* hBgAvg = (TH1D*)hBg1->Clone("hBgAvg_" + suffix);
        hBgAvg->Add(hBg2);
        hBgAvg->Scale(0.5);

        // Final Result (引き算)
        TH1D* hFinal = (TH1D*)hSig->Clone("hFinal_" + suffix);
        hFinal->Add(hBgAvg, -1.0);

        // ★ 面積規格化 (指定コードではコメントアウトされていましたが、密度比較のため有効化しました)
        // もしRaw Countの比較が必要な場合はここをコメントアウトしてください。
        //applyGeometryNorm(hFinal);

        // スタイル設定
        hFinal->SetMarkerStyle(markers[i]);
        hFinal->SetMarkerColor(colors[i]);
        hFinal->SetLineColor(colors[i]);
        hFinal->SetMarkerSize(1.2);
        hFinal->SetLineWidth(2);

        hResults.push_back(hFinal);
    }

    // ==========================================
    // 4. 描画: 重ね書き (Overlay)
    // ==========================================
    if (hResults.empty()) return;

    TCanvas* cFinal = new TCanvas("cFinalOverlay", "Jet Shape Overlay", 800, 600); 
    gPad->SetLogy();
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    // 凡例
    TLegend* legend = new TLegend(0.55, 0.60, 0.88, 0.80);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0); 
    legend->SetTextSize(0.04);
    legend->SetHeader("Centrality");

    // 描画ループ
    for (size_t i = 0; i < hResults.size(); ++i) {
        TH1D* h = hResults[i];
        
        if (i == 0) {
            // 最初のヒストグラムで軸設定
            h->SetTitle(""); 
            h->GetXaxis()->SetTitle("#it{r}");
            
            // 指定コードと同じY軸タイトル
            h->GetYaxis()->SetTitle("#rho(#it{r}) = #frac{1}{#it{N}_{jet}} #frac{d#it{p}_{T}}{d#it{r}#kern[0.3]{d}#eta}");
            h->GetYaxis()->SetTitleSize(0.045); 
            h->GetYaxis()->SetTitleOffset(1.5);
            h->GetXaxis()->SetTitleSize(0.045);
            
            h->GetXaxis()->SetRangeUser(0.0, 0.7);
            
            // Y軸レンジ調整 (0番目の最大値を基準に)
            double maxVal = h->GetMaximum();
            if (maxVal <= 0) maxVal = 1.0;
            h->GetYaxis()->SetRangeUser(0.005, maxVal * 10.0);

            h->Draw("PE");
        } else {
            h->Draw("PE SAME");
        }

        // 凡例追加
        TString label = TString::Format("%.0f-%.0f%%", centBins[i], centBins[i+1]);
        legend->AddEntry(h, label, "pe");
    }

    legend->Draw();

    // 実験情報 (Centrality表記は削除したもの)
    TPaveText* paveFinal = createPaveText(0.18, 0.18, 0.5, 0.45, minPtCorr, maxPtCorr);
    paveFinal->Draw();
}

