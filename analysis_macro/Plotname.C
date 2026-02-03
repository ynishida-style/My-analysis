#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

void plotTofBeta() {
    // --- スタイル設定 ---
    gStyle->SetOptStat(0);      // 統計ボックスを非表示
    gStyle->SetOptTitle(0);     // タイトル非表示
    gStyle->SetPalette(kBird);  // 見やすいカラーパレット設定

    // --- 1. ファイルの読み込み ---
    TFile *f = TFile::Open("AnalysisResults-70.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: File 'AnalysisResults-70.root' could not be opened." << std::endl;
        return;
    }

    // --- 2. ヒストグラムの取得 ---
    THnSparse *hSparse = (THnSparse*)f->Get("jet-shape-task/tpcDedx");
    if (!hSparse) {
        std::cerr << "Error: Histogram 'jet-shape-task/tofBeta' not found." << std::endl;
        f->Close();
        return;
    }

// Projection(1, 0): Y軸にAxis 1 (dEdx), X軸にAxis 0 (p) を割り当てて2次元化
    TH2D *h = hSparse->Projection(1, 0);
    h->SetName("hTpcDedxProj"); // 名前を付けておく（推奨）

    // --- 3. 軸ラベルとタイトルの設定 ---
    //h->SetTitle("TOF #beta vs Momentum"); // タイトル
    
    // X軸: p (GeV/c)
    h->GetXaxis()->SetTitle("#it{p}_{track} (GeV/#it{c})");
    h->GetXaxis()->SetTitleOffset(1.0);

    // Y軸: TOF beta
    h->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (arb. units)");
    h->GetYaxis()->SetTitleOffset(1.2); // ラベルが重ならないように調整

    // --- 4. 描画 ---
    TCanvas *c1 = new TCanvas("c1", "TOF Beta vs Momentum", 800, 600);
    
    // 余白調整（右側にカラーパレットを表示するスペースを確保）
    gPad->SetRightMargin(0.12);
    gPad->SetLeftMargin(0.12);
    
    // Z軸（頻度）を対数表示にする（PIDのラインが見やすくなります）
    gPad->SetLogz(); 

    // COLZオプションで2次元カラープロット描画
    h->Draw("COLZ");

    // 必要に応じて保存
    // c1->SaveAs("TofBetaPlot.pdf");

    // =========================================================
    // その他の説明書き（図の中）
    // =========================================================
    TPaveText *ptInfo = new TPaveText(0.15, 0.75, 0.45, 0.88, "NDC");
    ptInfo->SetBorderSize(0);
    
    // ★半透明設定
    ptInfo->SetFillStyle(1001);             
    ptInfo->SetFillColorAlpha(kWhite, 0.8); 

    ptInfo->SetTextAlign(12);
    ptInfo->SetTextFont(42);
    ptInfo->SetTextSize(0.028);
    ptInfo->AddText("#bf{This thesis}");
    ptInfo->AddText("0-10% Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    ptInfo->Draw();
}