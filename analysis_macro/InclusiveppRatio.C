#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h> // 追加
#include <TPaveText.h>
#include <TString.h>
#include <TLine.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility> // for std::pair

// ============================================================================================
// Helper Functions: Fitting
// (Collision System Specific Logic - PRESERVED)
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

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // (元のコードと同じPion用関数)
    std::cout << "--- Starting parameter optimization for PION in p bin " << p_min << "-" << p_max << " ---" << std::endl;
    
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.2) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.4) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.4) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                
                tempFit->SetParameter(0, peak_height); tempFit->SetParameter(1, 0.0); tempFit->SetParameter(2, 1.5);
                tempFit->SetParameter(3, peak_height * 0.1); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.0);
                tempFit->SetParameter(6, peak_height * 0.2); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 1.5);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 1.0);
                
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

    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    finalFit->SetParLimits(0, peak_height * 0.7, peak_height * 1.1);
    finalFit->SetParLimits(1, -0.05, 0.05);
    finalFit->SetParLimits(2, 0.4, 1.8);
    
    finalFit->SetParLimits(3, peak_height * 0.01, peak_height * 0.4);
    finalFit->SetParLimits(4, -3.8, -2.4);
    finalFit->SetParLimits(5, 0.5, 1.2);
    
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 0.6);
    finalFit->SetParLimits(7, -5.0, -3.7);
    finalFit->SetParLimits(8, 1.2, 1.8);
    
    finalFit->SetParLimits(9, peak_height * 0.01, peak_height * 0.5);
    finalFit->SetParLimits(10, 3.0, 8.0);
    finalFit->SetParLimits(11, 1.0, 1.8);

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

// ============================================================================================
// Core Logic: Purity AND PID Efficiency Calculation
// ============================================================================================

std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(TString fileName, TString particleType, double centMin, double centMax) {
    std::map<double, double> purityMap;
    std::map<double, double> effMap;
    
    TString sparseName;
    double sigMin, sigMax; 

    if (particleType == "proton") {
        sparseName = "tpcTofPr"; 
        sigMin = -3.5; sigMax = 0.5; 
    } else {
        sparseName = "tpcTofPi"; 
        sigMin = -0.5; sigMax = 3.5; 
    }

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return {purityMap, effMap};
    }
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { 
        std::cerr << "Error: Cannot find histogram " << sparseName << " in " << fileName << std::endl;
        f->Close(); 
        return {purityMap, effMap}; 
    }

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax); 

    double p_start = 0.6; 
    double p_end = 8.0; 
    double p_step = 0.2; 

    std::cout << "--- Calculating Purity & PID Efficiency for " << particleType << " ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;

        if (p_center <= 3.0) {
            purityMap[p_center] = 1.0;
            
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            double idealEff   = idealInCut / idealTotal;

            effMap[p_center] = idealEff; 
            continue;
        }

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        if (h_nSigma->GetEntries() < 50) { delete h_nSigma; continue; }

        TF1* fit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        TF1* fSig = new TF1("fSig", "gaus", -10, 10);
        fSig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));

        double intSigInCut = fSig->Integral(sigMin, sigMax); 
        double intSigTotal = fSig->Integral(-10, 10);        
        double intTotInCut = fit->Integral(sigMin, sigMax);  

        double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
        double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;
        
        if (purity > 1.0) purity = 1.0;
        if (purity < 0.0) purity = 0.0;
        if (pidEff > 1.0) pidEff = 1.0;
        if (pidEff < 0.0) pidEff = 0.0; 

        purityMap[p_center] = purity;
        effMap[p_center]    = pidEff;

        printf(" p=%.2f: Purity=%.3f, PID Eff=%.3f\n", p_center, purity, pidEff);

        delete h_nSigma; delete fit; delete fSig;
    }
    f->Close();
    return {purityMap, effMap};
}

TH1D* createCorrectedSpectrum(TString fileName, TString histName, 
                              std::map<double, double>& purityMap, 
                              std::map<double, double>& pidEffMap, 
                              double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    THnSparseD *hSparse = (THnSparseD*)f->Get(TString::Format("jet-shape-task/%s", histName.Data()));
    if (!hSparse) return nullptr;

    hSparse->GetAxis(2)->SetRangeUser(centMin, centMax);

    TH2D *h2 = (TH2D*)hSparse->Projection(0, 1);
    
    TH1D *hPt = h2->ProjectionX(TString::Format("hPt_%s", histName.Data()));
    hPt->Reset();
    hPt->SetDirectory(0);

    for (int by = 1; by <= h2->GetNbinsY(); ++by) {
        double p_val = h2->GetYaxis()->GetBinCenter(by);
        
        if (p_val > 8.0) continue;

        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0; 
        
        for (auto const& [map_p, map_pur] : purityMap) {
            if (std::abs(map_p - p_val) < min_diff) {
                min_diff = std::abs(map_p - p_val);
                purity = map_pur;
                pidEff = pidEffMap[map_p]; 
            }
        }
        
        if (min_diff > 0.3) {
             if (p_val > 3.0) continue; 
        }

        if (pidEff < 0.01) pidEff = 1.0; 

        double totalCorrectionFactor = purity / pidEff;

        for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
            double content = h2->GetBinContent(bx, by);
            double error   = h2->GetBinError(bx, by);
            
            if (content <= 0) continue;

            double corrContent = content * totalCorrectionFactor;
            double corrError   = error * totalCorrectionFactor; 

            double curContent = hPt->GetBinContent(bx);
            double curError   = hPt->GetBinError(bx);

            hPt->SetBinContent(bx, curContent + corrContent);
            hPt->SetBinError(bx, std::sqrt(curError*curError + corrError*corrError));
        }
    }
    f->Close();
    return hPt;
}

TH1D* calculateEfficiency(TString fileName, TString partName, double centMin, double centMax) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f) return nullptr;
    
    THnSparseD *hRecoSp = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptHistogram%sTof", partName.Data()));
    THnSparseD *hGenSp  = (THnSparseD*)f->Get(TString::Format("jet-shape-task/ptGenerated%s", partName.Data()));
    
    if (!hRecoSp || !hGenSp) {
        std::cout << "Error loading Efficiency Histograms for " << partName << std::endl;
        if(f) f->Close();
        return nullptr;
    }

    hRecoSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSp->GetAxis(2)->SetRangeUser(centMin, centMax);
    
    TH1D *hReco = hRecoSp->Projection(0);
    TH1D *hGen  = hGenSp->Projection(0);
    
    TH1D *hEff = (TH1D*)hReco->Clone(TString::Format("eff_%s", partName.Data()));
    hEff->Divide(hReco, hGen, 1, 1, "B");
    hEff->SetDirectory(0);
    f->Close();
    return hEff;
}

// ============================================================================================
// Main
// ============================================================================================
void PlotInclusiveRatioFinal() {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    
    // パラメータ保持 (AnalysisResults-98.root, 60-100%)
    const TString dataFile = "AnalysisResults-99.root";
    const TString mcFile   = "AnalysisResults-91.root";

    const double dataCentMin = 0.0;
    const double dataCentMax = 100.0;
    const double mcCentMin = 0.0;
    const double mcCentMax = 100.0;

    auto mapsPr = getPidCorrectionMaps(dataFile, "proton", dataCentMin, dataCentMax);
    auto mapsPi = getPidCorrectionMaps(dataFile, "pion", dataCentMin, dataCentMax);
    
    std::map<double, double> purMapPr = mapsPr.first;
    std::map<double, double> effPidPr = mapsPr.second;
    std::map<double, double> purMapPi = mapsPi.first;
    std::map<double, double> effPidPi = mapsPi.second;

    TH1D* hPrRaw = createCorrectedSpectrum(dataFile, "pVsPtForPr", purMapPr, effPidPr, dataCentMin, dataCentMax);
    TH1D* hPiRaw = createCorrectedSpectrum(dataFile, "pVsPtForPi", purMapPi, effPidPi, dataCentMin, dataCentMax);

    TH1D* hEffPr = calculateEfficiency(mcFile, "Proton", mcCentMin, mcCentMax);
    TH1D* hEffPi = calculateEfficiency(mcFile, "Pion", mcCentMin, mcCentMax);

    if (!hPrRaw || !hPiRaw || !hEffPr || !hEffPi) {
        std::cerr << "Error creating histograms." << std::endl;
        return;
    }

    TH1D* hPrCorr = (TH1D*)hPrRaw->Clone("hPrCorr");
    hPrCorr->Divide(hEffPr);
    TH1D* hPiCorr = (TH1D*)hPiRaw->Clone("hPiCorr");
    hPiCorr->Divide(hEffPi);

    // 比の計算
    hPrCorr->Rebin(2);
    hPiCorr->Rebin(2);
    TH1D* hRatio = (TH1D*)hPrCorr->Clone("hRatio");
    hRatio->Divide(hPiCorr);

    // --- Plotting ---
    hRatio->GetXaxis()->SetRangeUser(0.0, 6.0);
    hPiCorr->GetXaxis()->SetRangeUser(0.0, 6.0); 

    // =========================================================================
    // Canvas 1: Spectra
    // =========================================================================
    TCanvas *c1 = new TCanvas("cSpectra", "Spectra", 800, 600);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    
    hPiCorr->SetMarkerStyle(20); 
    hPiCorr->SetMarkerColor(kBlue + 1); 
    hPiCorr->SetLineColor(kBlue + 1);
    
    hPrCorr->SetMarkerStyle(21); 
    hPrCorr->SetMarkerColor(kRed + 1);  
    hPrCorr->SetLineColor(kRed + 1);

    hPiCorr->SetTitle("Inclusive p_{T} Spectra");
    hPiCorr->GetYaxis()->SetTitle("dN/dp_{T}");
    hPiCorr->Draw("PE");
    hPrCorr->Draw("PE SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88); 
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hPiCorr, "Pions", "lpe");
    leg->AddEntry(hPrCorr, "Protons", "lpe");
    leg->Draw();

    TPaveText *pt1 = new TPaveText(0.20, 0.20, 0.55, 0.40, "NDC");
    pt1->SetFillStyle(0); pt1->SetBorderSize(0); pt1->SetTextAlign(12); pt1->SetTextFont(42); pt1->SetTextSize(0.035);
    pt1->AddText("#bf{This thesis}");
    pt1->AddText("|#eta_{track}| < 0.9");
    pt1->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    pt1->AddText(TString::Format("centrality : %.0f-%.0f%%", dataCentMin, dataCentMax));
    pt1->AddText("inclusive");
    pt1->Draw();

    // =========================================================================
    // Canvas 2: Ratio
    // =========================================================================
    TCanvas *c2 = new TCanvas("cRatio", "Ratio", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetGridy();

    hRatio->SetTitle("p/#pi Ratio");
    hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitleOffset(1.2);
    hRatio->SetMarkerStyle(34); hRatio->SetMarkerColor(kBlack); hRatio->SetLineColor(kBlack);
    hRatio->SetMinimum(0.0); hRatio->SetMaximum(1.1); 
    hRatio->Draw("E1");


    TPaveText *pt2 = new TPaveText(0.15, 0.75, 0.5, 0.9, "NDC");
    pt2->SetFillStyle(0); pt2->SetBorderSize(0); pt2->SetTextAlign(12); pt2->SetTextFont(42); pt2->SetTextSize(0.035);
    pt2->AddText("#bf{This thesis}");
    pt2->AddText("|#eta_{track}| < 0.9");
    pt2->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    pt2->AddText(TString::Format("centrality : %.0f-%.0f%%", dataCentMin, dataCentMax));
    pt2->AddText("inclusive");
    pt2->Draw();

    // =========================================================================
    // Canvas 3: Purity Check (New Addition)
    // =========================================================================
    TCanvas *c3 = new TCanvas("cPurity", "Purity Check", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetGridy();

    // TGraphErrors を使用 (横方向のエラーバーを描画するため)
    TGraphErrors *gPurPr = new TGraphErrors();
    TGraphErrors *gPurPi = new TGraphErrors();

    double p_step = 0.2; // ビン幅
    double x_err = p_step / 2.0; // 誤差は幅の半分

    int nPr = 0;
    for (auto const& [p, val] : purMapPr) {
        gPurPr->SetPoint(nPr, p, val);
        gPurPr->SetPointError(nPr, x_err, 0.0); // 縦誤差は0
        nPr++;
    }

    int nPi = 0;
    for (auto const& [p, val] : purMapPi) {
        gPurPi->SetPoint(nPi, p, val);
        gPurPi->SetPointError(nPi, x_err, 0.0);
        nPi++;
    }

    // 描画範囲の設定 (3 GeV/c - 6 GeV/c)
    TH1F *frame = c3->DrawFrame(3.0, 0.0, 6.0, 1.2); 
    frame->SetTitle("PID Purity");
    frame->GetXaxis()->SetTitle("#it{p}_{track} (GeV/#it{c})");
    frame->GetYaxis()->SetTitle("Purity");
    frame->GetYaxis()->SetTitleOffset(1.2);

    gPurPr->SetMarkerStyle(21); gPurPr->SetMarkerColor(kRed + 1); gPurPr->SetLineColor(kRed + 1); gPurPr->SetMarkerSize(1.2);
    gPurPi->SetMarkerStyle(20); gPurPi->SetMarkerColor(kBlue + 1); gPurPi->SetLineColor(kBlue + 1); gPurPi->SetMarkerSize(1.2);

    // 描画 (P E1 オプション等でエラーバーを描画可能ですが、TGraphErrorsはデフォルトで描画します)
    gPurPr->Draw("P SAME"); 
    gPurPi->Draw("P SAME");

    TLegend *legPur = new TLegend(0.6, 0.2, 0.88, 0.35); 
    legPur->SetBorderSize(0); 
    legPur->SetFillStyle(0);
    legPur->AddEntry(gPurPi, "Pion Purity", "p");
    legPur->AddEntry(gPurPr, "Proton Purity", "p");
    legPur->Draw();

    TPaveText *pt3 = new TPaveText(0.20, 0.20, 0.55, 0.35, "NDC");
    pt3->SetFillStyle(0); pt3->SetBorderSize(0); pt3->SetTextAlign(12); pt3->SetTextFont(42); pt3->SetTextSize(0.035);

    pt3->AddText("#bf{This thesis}");
    pt3->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    pt3->AddText("|#eta_{track}| < 0.9");
    pt3->Draw();

    //c1->SaveAs("InclusiveSpectra_FullyCorr.pdf");
    //c2->SaveAs("InclusiveRatio_FullyCorr.pdf");
    //c3->SaveAs("InclusivePurity_Check.pdf");
}