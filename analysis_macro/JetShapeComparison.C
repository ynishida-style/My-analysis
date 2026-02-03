#include <iostream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"

// 共通: THnSparseから引き算済みのJet Shapeを取得する関数
TH1D* getJetShape(TString fileName, TString dirName, double minCent, double maxCent, double minPt, double maxPt, TString suffix) {
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << fileName << std::endl;
        return nullptr;
    }
    TDirectory* dir = (TDirectory*)file->Get(dirName);
    if (!dir) return nullptr;

    auto getProf = [&](TString hName, TString s) {
        THnSparse* hn = (THnSparse*)dir->Get(hName);
        if (!hn) return (TProfile*)nullptr;
        hn->GetAxis(2)->SetRangeUser(minPt, maxPt);
        hn->GetAxis(3)->SetRangeUser(minCent, maxCent);
        TH2D* h2 = (TH2D*)hn->Projection(1, 0);
        return h2->ProfileX(Form("hProf_%s_%s", hName.Data(), s.Data()));
    };

    TProfile* pSig = getProf("ptSum", suffix);
    TProfile* pBg1 = getProf("ptSumBg1", suffix);
    TProfile* pBg2 = getProf("ptSumBg2", suffix);

    if (!pSig || !pBg1 || !pBg2) return nullptr;

    TH1D* hFinal = pSig->ProjectionX(Form("hFinal_%s", suffix.Data()));
    TH1D* hBgAvg = pBg1->ProjectionX(Form("hBgAvg_%s", suffix.Data()));
    hBgAvg->Add(pBg2->ProjectionX("tmp"));
    hBgAvg->Scale(0.5);

    hFinal->Add(hBgAvg, -1.0);
    
    // 統計誤差のみ考慮（必要に応じてシステム誤差を追加）
    hFinal->SetDirectory(0); // ファイルを閉じてもデータが残るようにする
    file->Close();
    return hFinal;
}

void compareJetShape() {
    gStyle->SetOptStat(0);
    
    // --- 1. 設定 ---
    double minPt = 60.0; 
    double maxPt = 200.0;
    TString dirName = "jet-shape-task";

    // --- 2. データの取得 ---
    // pp データ
    TH1D* hPP = getJetShape("AnalysisResults-83.root", dirName, 0, 100, minPt, maxPt, "pp");
    // Pb-Pb 0-10% データ
    TH1D* hPbPb = getJetShape("AnalysisResults-106.root", dirName, 0, 30, minPt, maxPt, "PbPb");

    if (!hPP || !hPbPb) return;

    // --- 3. 描画設定 ---
    TCanvas* c = new TCanvas("cCompare", "Compare pp and PbPb", 800, 700);
    c->SetLogy();
    c->SetTicks(1, 1);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.13);

    // スタイル適用
    hPP->SetMarkerStyle(20);
    hPP->SetMarkerColor(kBlack);
    hPP->SetLineColor(kBlack);
    hPP->SetMarkerSize(1.2);

    hPbPb->SetMarkerStyle(21);
    hPbPb->SetMarkerColor(kRed+1);
    hPbPb->SetLineColor(kRed+1);
    hPbPb->SetMarkerSize(1.2);

    // 軸ラベル
    hPP->SetTitle("");
    hPP->GetXaxis()->SetTitle("#it{r}");
    hPP->GetYaxis()->SetTitle("#rho(#it{r}) = #frac{1}{#it{N}_{jet}} #frac{d#it{p}_{T}}{d#it{r}d#eta}");
    hPP->GetXaxis()->SetTitleSize(0.04);
    hPP->GetYaxis()->SetTitleSize(0.04);
    hPP->GetYaxis()->SetTitleOffset(1.4);
    hPP->GetXaxis()->SetRangeUser(0.0, 0.7);
    hPP->GetYaxis()->SetRangeUser(0.01, hPbPb->GetMaximum() * 10.0);

    hPP->Draw("PE");
    hPbPb->Draw("PE SAME");

    // --- 4. 凡例とテキスト ---
    TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.82);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hPP, "pp #sqrt{#it{s}} = 13.6 TeV", "lp");
    leg->AddEntry(hPbPb, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", "lp");
    leg->Draw();

    TPaveText *pt = new TPaveText(0.18, 0.15, 0.45, 0.35, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("#bf{This thesis}");
    pt->AddText("ch-particle jet, Anti-#it{k}_{T}, #it{R} = 0.4");
    pt->AddText(Form("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", minPt, maxPt));
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    pt->Draw();
}