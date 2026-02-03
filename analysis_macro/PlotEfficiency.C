#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"

// __________________________________________________________________________________
TH1D* getProjection(TFile* f, TString listName, TString histName) 
{
    TDirectory* dir = (TDirectory*)f->Get(listName); 
    if (!dir) return nullptr;

    THnSparseD* hn = (THnSparseD*)dir->Get(histName);
    if (!hn) return nullptr;

    TH1D* hRaw = (TH1D*)hn->Projection(0); 
    
    int nBins = 100;
    double ptMin = 0.0;
    double ptMax = 20.0;
    TH1D* hRebinned = new TH1D(Form("%s_proj", histName.Data()), "", nBins, ptMin, ptMax);
    
    // Sumw2を事前に作成
    hRebinned->Sumw2();

    for (int i = 1; i <= hRaw->GetNbinsX(); ++i) {
        double center  = hRaw->GetBinCenter(i);
        if (center < ptMin || center > ptMax) continue;

        double content = hRaw->GetBinContent(i);
        double error   = hRaw->GetBinError(i);
        
        int newBin = hRebinned->FindBin(center);
        double currentContent = hRebinned->GetBinContent(newBin);
        
        hRebinned->SetBinContent(newBin, currentContent + content);
        hRebinned->SetBinError(newBin, sqrt(pow(hRebinned->GetBinError(newBin), 2) + error*error));
    }
    
    hRebinned->SetDirectory(0);
    delete hRaw;
    return hRebinned;
}

// __________________________________________________________________________________
// ★修正: 引数に markerStyle を追加しました
void setStyle(TH1D* h, int color, int markerStyle, TString yTitle) {
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle); // ★ここで適用
    h->SetMarkerSize(1.0);          // 形を変えるときは少し大きめが見やすいです

    h->SetTitle("");
    h->GetYaxis()->SetTitle(yTitle);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetRangeUser(0.0, 1.2);
    
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetRangeUser(0.0, 6.0); 
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
    pt->AddText("Pb-Pb inclusive #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText("O-O inclusive #sqrt{#it{s}_{NN}} = 5.36 TeV");
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
    
    TString fileName = "AnalysisResults_mc.root"; 
    TString listName = "jet-shape-task";
    
    TFile *file = new TFile(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
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

    if (!hGenPion || !hRecoPion || !hTofPion) {
        std::cerr << "Error: Histograms not found." << std::endl;
        return;
    }

    // --- Efficiency Calculation ---
    TH1D* hEffTrackPi = (TH1D*)hRecoPion->Clone("hEffTrackPi"); hEffTrackPi->Divide(hRecoPion, hGenPion, 1.0, 1.0, "B");
    TH1D* hEffTrackKa = (TH1D*)hRecoKaon->Clone("hEffTrackKa"); hEffTrackKa->Divide(hRecoKaon, hGenKaon, 1.0, 1.0, "B");
    TH1D* hEffTrackPr = (TH1D*)hRecoProton->Clone("hEffTrackPr"); hEffTrackPr->Divide(hRecoProton, hGenProton, 1.0, 1.0, "B");

    TH1D* hEffMatchPi = (TH1D*)hTofPion->Clone("hEffMatchPi"); hEffMatchPi->Divide(hTofPion, hRecoPion, 1.0, 1.0, "B");
    TH1D* hEffMatchKa = (TH1D*)hTofKaon->Clone("hEffMatchKa"); hEffMatchKa->Divide(hTofKaon, hRecoKaon, 1.0, 1.0, "B");
    TH1D* hEffMatchPr = (TH1D*)hTofProton->Clone("hEffMatchPr"); hEffMatchPr->Divide(hTofProton, hRecoProton, 1.0, 1.0, "B");

    TH1D* hEffTotalPi = (TH1D*)hTofPion->Clone("hEffTotalPi"); hEffTotalPi->Divide(hTofPion, hGenPion, 1.0, 1.0, "B");
    TH1D* hEffTotalKa = (TH1D*)hTofKaon->Clone("hEffTotalKa"); hEffTotalKa->Divide(hTofKaon, hGenKaon, 1.0, 1.0, "B");
    TH1D* hEffTotalPr = (TH1D*)hTofProton->Clone("hEffTotalPr"); hEffTotalPr->Divide(hTofProton, hGenProton, 1.0, 1.0, "B");

    // --- Style Settings ---
    int colPi = kGreen+1;
    int colKa = kRed+1;
    int colPr = kBlue+1;

    // ★ マーカーの形状を指定
    int mkPi = 20; // Full Circle (●)
    int mkKa = 21; // Full Square (■)
    int mkPr = 22; // Full Diamond (◆) 

    setStyle(hEffTrackPi, colPi, mkPi, "Tracking efficiency");
    setStyle(hEffTrackKa, colKa, mkKa, "Tracking efficiency");
    setStyle(hEffTrackPr, colPr, mkPr, "Tracking efficiency");

    setStyle(hEffMatchPi, colPi, mkPi, "TOF Matching efficiency");
    setStyle(hEffMatchKa, colKa, mkKa, "TOF Matching efficiency");
    setStyle(hEffMatchPr, colPr, mkPr, "TOF Matching efficiency");

    setStyle(hEffTotalPi, colPi, mkPi, "Total efficiency (Track #times TOF)");
    setStyle(hEffTotalKa, colKa, mkKa, "Total efficiency (Track #times TOF)");
    setStyle(hEffTotalPr, colPr, mkPr, "Total efficiency (Track #times TOF)");

    TPaveText* pt = createPaveText(0.45, 0.20, 0.88, 0.50);

    // --- Draw ---
    int cW = 600; int cH = 500;

    // 凡例のオプションも "lpe" (Line + Point + Error) に変更してマーカーを表示させます

    TCanvas *cTrack = new TCanvas("cTrack", "Tracking Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffTrackPi->Draw("PE"); hEffTrackKa->Draw("PE same"); hEffTrackPr->Draw("PE same");
    
    TLegend *leg1 = new TLegend(0.55, 0.55, 0.75, 0.75);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.05);
    leg1->AddEntry(hEffTrackPi, "#pi", "pe"); 
    leg1->AddEntry(hEffTrackKa, "K", "pe"); 
    leg1->AddEntry(hEffTrackPr, "p", "pe");
    leg1->Draw();
    pt->Draw("SAME");
    cTrack->SaveAs("Efficiency_Tracking.pdf");

    TCanvas *cMatch = new TCanvas("cMatch", "Matching Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffMatchPi->Draw("PE"); hEffMatchKa->Draw("PE same"); hEffMatchPr->Draw("PE same");
    
    TLegend *leg2 = new TLegend(0.2, 0.6, 0.4, 0.8);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.05);
    leg2->AddEntry(hEffMatchPi, "#pi", "pe"); 
    leg2->AddEntry(hEffMatchKa, "K", "pe"); 
    leg2->AddEntry(hEffMatchPr, "p", "pe");
    leg2->Draw();
    pt->Draw("SAME");
    cMatch->SaveAs("Efficiency_Matching.pdf");

    TCanvas *cTotal = new TCanvas("cTotal", "Total Efficiency", cW, cH);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.05);
    hEffTotalPi->Draw("PE"); hEffTotalKa->Draw("PE same"); hEffTotalPr->Draw("PE same");
    
    TLegend *leg3 = new TLegend(0.55, 0.55, 0.75, 0.75);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0); leg3->SetTextSize(0.05);
    leg3->AddEntry(hEffTotalPi, "#pi", "pe"); 
    leg3->AddEntry(hEffTotalKa, "K", "pe"); 
    leg3->AddEntry(hEffTotalPr, "p", "pe");
    leg3->Draw();
    pt->Draw("SAME");
    cTotal->SaveAs("Efficiency_Total.pdf");
}