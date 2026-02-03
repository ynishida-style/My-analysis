#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"

// __________________________________________________________________________________
// ★修正: THnSparseからの射影ではなく、TH1Dを直接取得するように変更
TH1D* getProjection(TFile* f, TString listName, TString histName) 
{
    // リスト（ディレクトリ）を取得
    TDirectory* dir = (TDirectory*)f->Get(listName); 
    if (!dir) {
        std::cerr << "Error: Directory '" << listName << "' not found." << std::endl;
        return nullptr;
    }

    // ヒストグラムを取得 (THnSparseD ではなく TH1D として取得)
    TH1D* hRaw = (TH1D*)dir->Get(histName);
    if (!hRaw) {
        std::cerr << "Error: Histogram '" << histName << "' not found inside '" << listName << "'." << std::endl;
        return nullptr;
    }

    // ヒストグラムをコピーして名前を変更 (メモリ管理のため)
    TH1D* hRet = (TH1D*)hRaw->Clone(Form("%s_clone", histName.Data()));
    
    // 必要であればここでRebinを行う (例: 2ビンを1つにまとめる)
    // hRet->Rebin(2); 

    hRet->Sumw2();
    hRet->SetDirectory(0); // ファイルを閉じても消えないようにする
    
    // 元のポインタはGetで取得したものなのでdeleteしてはいけない（Fileが管理しているため）
    // Cloneした方(hRet)を返す
    return hRet;
}

// __________________________________________________________________________________
void setStyle(TH1D* h, int color, int markerStyle, TString yTitle) {
    if (!h) return;
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle); 
    h->SetMarkerSize(0.8);          

    h->SetTitle("");
    h->GetYaxis()->SetTitle(yTitle);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetRangeUser(0.0, 1.3); // 範囲指定
    
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.1);
    
    // 描画範囲 (必要に応じて調整)
    h->GetXaxis()->SetRangeUser(0.15, 10.0); 
}

// __________________________________________________________________________________
TPaveText* createPaveText(double x1, double y1, double x2, double y2) {
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12); 
    pt->SetTextFont(42); 
    pt->SetTextSize(0.04); 
    
    pt->AddText("#bf{This thesis}");
    pt->AddText("pp inclusive #sqrt{#it{s}} = 13.6 TeV");
    pt->AddText("Monte Carlo data");
    pt->AddText("#it{p}_{T,track} > 0.15 GeV/#it{c}");
    pt->AddText("|#it{#eta}_{track}| < 0.9");
    
    return pt;
}

// __________________________________________________________________________________
void DrawEfficiency_Separate() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    // ★注意: ファイル名が正しいか確認してください
    TString fileName = "AnalysisResults-58.root"; 
    TString listName = "jet-shape-task";
    
    TFile *file = new TFile(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    TH1D *hGenPion    = getProjection(file, listName, "ptGeneratedPion");
    TH1D *hGenKaon    = getProjection(file, listName, "ptGeneratedKaon");
    TH1D *hGenProton  = getProjection(file, listName, "ptGeneratedProton");

    TH1D *hRecoPion   = getProjection(file, listName, "ptHistogramPion");
    TH1D *hRecoKaon   = getProjection(file, listName, "ptHistogramKaon");
    TH1D *hRecoProton = getProjection(file, listName, "ptHistogramProton");

    TH1D *hTofPion    = getProjection(file, listName, "ptHistogramPionTof");
    TH1D *hTofKaon    = getProjection(file, listName, "ptHistogramKaonTof");
    TH1D *hTofProton  = getProjection(file, listName, "ptHistogramProtonTof");

    // エラーチェック
    if (!hGenPion || !hRecoPion || !hTofPion) {
        std::cerr << "Critical Error: One or more histograms could not be loaded." << std::endl;
        file->Close();
        return;
    }

    // --- Efficiency Calculation ---
    // Reco / Gen
    TH1D* hEffTrackPi = (TH1D*)hRecoPion->Clone("hEffTrackPi"); hEffTrackPi->Divide(hRecoPion, hGenPion, 1.0, 1.0, "B");
    TH1D* hEffTrackKa = (TH1D*)hRecoKaon->Clone("hEffTrackKa"); hEffTrackKa->Divide(hRecoKaon, hGenKaon, 1.0, 1.0, "B");
    TH1D* hEffTrackPr = (TH1D*)hRecoProton->Clone("hEffTrackPr"); hEffTrackPr->Divide(hRecoProton, hGenProton, 1.0, 1.0, "B");

    // TOF / Reco
    TH1D* hEffMatchPi = (TH1D*)hTofPion->Clone("hEffMatchPi"); hEffMatchPi->Divide(hTofPion, hRecoPion, 1.0, 1.0, "B");
    TH1D* hEffMatchKa = (TH1D*)hTofKaon->Clone("hEffMatchKa"); hEffMatchKa->Divide(hTofKaon, hRecoKaon, 1.0, 1.0, "B");
    TH1D* hEffMatchPr = (TH1D*)hTofProton->Clone("hEffMatchPr"); hEffMatchPr->Divide(hTofProton, hRecoProton, 1.0, 1.0, "B");

    // TOF / Gen
    TH1D* hEffTotalPi = (TH1D*)hTofPion->Clone("hEffTotalPi"); hEffTotalPi->Divide(hTofPion, hGenPion, 1.0, 1.0, "B");
    TH1D* hEffTotalKa = (TH1D*)hTofKaon->Clone("hEffTotalKa"); hEffTotalKa->Divide(hTofKaon, hGenKaon, 1.0, 1.0, "B");
    TH1D* hEffTotalPr = (TH1D*)hTofProton->Clone("hEffTotalPr"); hEffTotalPr->Divide(hTofProton, hGenProton, 1.0, 1.0, "B");

    // --- Style Settings ---
    int colPi = kRed+1;
    int colKa = kGreen+2;
    int colPr = kBlue+1;
    int mkPi = 20; 
    int mkKa = 21; 
    int mkPr = 22; 

    setStyle(hEffTrackPi, colPi, mkPi, "Tracking Efficiency");
    setStyle(hEffTrackKa, colKa, mkKa, "Tracking Efficiency");
    setStyle(hEffTrackPr, colPr, mkPr, "Tracking Efficiency");

    setStyle(hEffMatchPi, colPi, mkPi, "TOF Matching Efficiency");
    setStyle(hEffMatchKa, colKa, mkKa, "TOF Matching Efficiency");
    setStyle(hEffMatchPr, colPr, mkPr, "TOF Matching Efficiency");

    setStyle(hEffTotalPi, colPi, mkPi, "Total Efficiency");
    setStyle(hEffTotalKa, colKa, mkKa, "Total Efficiency");
    setStyle(hEffTotalPr, colPr, mkPr, "Total Efficiency");

    TPaveText* pt = createPaveText(0.45, 0.20, 0.88, 0.55);

    int cW = 800; int cH = 600;

    // --- Draw ---
    TCanvas *cTrack = new TCanvas("cTrack", "Tracking Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffTrackPi->Draw("PE"); hEffTrackKa->Draw("PE same"); hEffTrackPr->Draw("PE same");
    
    TLegend *leg1 = new TLegend(0.2, 0.7, 0.4, 0.85); 
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.04);
    leg1->AddEntry(hEffTrackPi, "#pi", "pe"); 
    leg1->AddEntry(hEffTrackKa, "K", "pe"); 
    leg1->AddEntry(hEffTrackPr, "p", "pe");
    leg1->Draw();
    pt->Draw("SAME");
    //cTrack->SaveAs("Efficiency_Tracking.pdf");

    TCanvas *cMatch = new TCanvas("cMatch", "Matching Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffMatchPi->Draw("PE"); hEffMatchKa->Draw("PE same"); hEffMatchPr->Draw("PE same");
    leg1->Draw(); 
    pt->Draw("SAME");
    //cMatch->SaveAs("Efficiency_Matching.pdf");

    TCanvas *cTotal = new TCanvas("cTotal", "Total Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffTotalPi->Draw("PE"); hEffTotalKa->Draw("PE same"); hEffTotalPr->Draw("PE same");
    leg1->Draw();
    pt->Draw("SAME");
    //cTotal->SaveAs("Efficiency_Total.pdf");
    
    // ファイルは閉じる（ヒストグラムはCloneしているので安全）
    // file->Close(); // new TFileしているので、プログラム終了まで開けておいても良いが、閉じるならここ
}