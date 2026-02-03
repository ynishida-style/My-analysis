#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TString.h>

// __________________________________________________________________________________
TPaveText* createPaveText(double x1, double y1, double x2, double y2, 
                          double jetPtMin, double jetPtMax,
                          double centMin, double centMax) {
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12); // 左揃え・垂直中央
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    
    pt->AddText("#bf{This thesis}");
    pt->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    //pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    //pt->AddText(TString::Format("Centrality %.0f-%.0f%%", centMin, centMax));
    // #{k} を #it{k} に修正しました (kをイタリックにする標準的な記法)
    pt->AddText("ch-particle jet, Anti-#it{k}_{T}, R = 0.4"); 
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    
    return pt;
}

// __________________________________________________________________________________
void CompareJetPt_NoMarker() {
    // --- スタイル設定 ---
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // --- 1. ファイル読み込み ---
    TFile *f = TFile::Open("AnalysisResults-83.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: File 'AnalysisResults-82.root' could not be opened." << std::endl;
        return;
    }

    // --- 2. ヒストグラム取得 ---
    TH1 *hRaw = (TH1*)f->Get("jet-shape-task/jetPt");
    TH1 *hSub = (TH1*)f->Get("jet-shape-task/ptCorr");

    if (!hRaw || !hSub) {
        std::cerr << "Error: Histograms not found." << std::endl;
        f->Close();
        return;
    }

    // --- 3. 見た目の設定（マーカーなし・線のみ） ---
    
    // jetPt (Raw) -> 青線
    hRaw->SetLineColor(kBlue+1);
    hRaw->SetLineWidth(2);
    hRaw->SetMarkerStyle(0); 
    hRaw->SetMarkerSize(0);

    // ptCorr (UE subtracted) -> 赤線
    hSub->SetLineColor(kRed+1);
    hSub->SetLineWidth(2);
    hSub->SetMarkerStyle(0); 
    hSub->SetMarkerSize(0);

    // --- 4. 軸の設定 ---
    hRaw->GetXaxis()->SetRangeUser(0.0, 200.0);
    
    // Log表示にするため、最小値は0より大きい値にする必要があります
    // ジェットスペクトルは見やすさを考慮して 0.1, 0.5, 1.0, 10.0 あたりを設定します
    hRaw->SetMinimum(1.0); 

    hRaw->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");
    hRaw->GetYaxis()->SetTitle("Counts");
    hRaw->GetXaxis()->SetTitleOffset(1.2);
    hRaw->GetYaxis()->SetTitleOffset(1.3);

    // --- 5. 描画 ---
    TCanvas *c1 = new TCanvas("c1", "Jet Pt Comparison", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    
    // ★ Logスケール設定
    gPad->SetLogy();

    hRaw->Draw("HIST");       
    hSub->Draw("HIST SAME"); 

    // --- 6. 凡例 (Legend) ---
    TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hRaw, "Raw #it{p}_{T}", "l");      
    leg->AddEntry(hSub, "UE subtracted", "l");
    leg->Draw();

    // --- 7. テキストボックスの追加 ---
    TPaveText* ptInfo = createPaveText(0.55, 0.45, 0.88, 0.68, 0.0, 200.0, 0.0, 10.0);
    ptInfo->Draw();

    // c1->SaveAs("JetPt_LineOnly_Log.pdf");
}