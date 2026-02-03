#include <iostream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"

// 共通: ファイルを指定してJet Shapeを取得する関数
TH1D* getJetShapeFromFile(TString fileName, TString dirName, double minCent, double maxCent, double minPt, double maxPt, TString suffix) {
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << fileName << std::endl;
        return nullptr;
    }
    TDirectory* dir = (TDirectory*)file->Get(dirName);
    if (!dir) {
        file->Close();
        return nullptr;
    }

    auto getProf = [&](TString hName, TString s) {
        THnSparse* hn = (THnSparse*)dir->Get(hName);
        if (!hn) return (TProfile*)nullptr;
        // 軸の範囲設定
        hn->GetAxis(2)->SetRangeUser(minPt, maxPt);
        hn->GetAxis(3)->SetRangeUser(minCent, maxCent);
        
        // 射影。名前が被らないように suffix を付与
        TH2D* h2 = (TH2D*)hn->Projection(1, 0);
        h2->SetName(Form("h2_%s_%s", hName.Data(), s.Data()));
        
        TProfile* prof = h2->ProfileX(Form("hProf_%s_%s", hName.Data(), s.Data()));
        delete h2; // 不要になった中間オブジェクトを削除
        return prof;
    };

    TProfile* pSig = getProf("ptSum", suffix);
    TProfile* pBg1 = getProf("ptSumBg1", suffix);
    TProfile* pBg2 = getProf("ptSumBg2", suffix);

    if (!pSig || !pBg1 || !pBg2) {
        std::cerr << "Error: Histograms not found in " << fileName << std::endl;
        file->Close();
        return nullptr;
    }

    // 各 Profile から TH1D を作成（名前をユニークに）
    TH1D* hFinal = pSig->ProjectionX(Form("hFinal_%s", suffix.Data()));
    TH1D* hBg1   = pBg1->ProjectionX(Form("hBg1_%s", suffix.Data()));
    TH1D* hBg2   = pBg2->ProjectionX(Form("hBg2_%s", suffix.Data()));

    // 背景の平均計算
    TH1D* hBgAvg = (TH1D*)hBg1->Clone(Form("hBgAvg_%s", suffix.Data()));
    hBgAvg->Add(hBg2);
    hBgAvg->Scale(0.5);

    // 引き算
    hFinal->Add(hBgAvg, -1.0);
    
    // 後片付け
    hFinal->SetDirectory(0); 
    delete hBg1;
    delete hBg2;
    delete hBgAvg;
    
    file->Close();
    return hFinal;
}

void compareJetShapeFiles() {
    gStyle->SetOptStat(0);
    
    double minPt = 60.0; 
    double maxPt = 200.0;
    double minCent = 0.0;
    double maxCent = 30.0;
    TString dirName = "jet-shape-task";

    std::vector<TString> fileNames = {"AnalysisResults-97.root", "AnalysisResults-93.root", "AnalysisResults-82.root"};
    std::vector<TString> labels    = {"R = 0.2", "R = 0.3", "R = 0.4"};
    std::vector<int> colors        = {kBlue+1, kRed+1, kGreen+2};
    std::vector<int> markers       = {20, 21, 22};

    TCanvas* c = new TCanvas("cCompR", "Compare Jet Radii", 800, 700);
    c->SetLogy();
    c->SetTicks(1, 1);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.13);

    TLegend* leg = new TLegend(0.60, 0.65, 0.88, 0.82);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    for (size_t i = 0; i < fileNames.size(); ++i) {
        // suffix を変えることで、内部のヒストグラム名が衝突しないようにする
        TH1D* h = getJetShapeFromFile(fileNames[i], dirName, minCent, maxCent, minPt, maxPt, Form("Ridx%zu", i));
        if (!h) continue;

        h->SetMarkerStyle(markers[i]);
        h->SetMarkerColor(colors[i]);
        h->SetLineColor(colors[i]);
        h->SetMarkerSize(1.0);
        h->SetLineWidth(2);

        if (i == 0) {
            h->SetTitle("");
            h->GetXaxis()->SetTitle("#it{r}");
            h->GetYaxis()->SetTitle("#rho(#it{r}) = #frac{1}{#it{N}_{jet}} #frac{d#it{p}_{T}}{d#it{r}d#eta}");
            h->GetXaxis()->SetTitleSize(0.03);
            h->GetYaxis()->SetTitleSize(0.03);
            h->GetYaxis()->SetTitleOffset(1.4);
            h->GetXaxis()->SetRangeUser(0.0, 0.7);
            h->GetYaxis()->SetRangeUser(0.01, h->GetMaximum() * 20.0);
            h->Draw("PE");
        } else {
            h->Draw("PE SAME");
        }
        leg->AddEntry(h, labels[i], "lp");
    }

    leg->Draw();

    TPaveText *pt = new TPaveText(0.18, 0.15, 0.50, 0.35, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText(Form("Centrality %.0f-%.0f%%", minCent, maxCent));
    pt->AddText(Form("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", minPt, maxPt));
    pt->Draw();
}