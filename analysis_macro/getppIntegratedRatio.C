#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility>

//________________________________________________________________________________________________
// Fitting Functions (Ported from Code B - Robust Version)
//________________________________________________________________________________________________

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max) {
    // --- Grid Search Ranges ---
    std::vector<double> p4_range; for (double p4 = 0.8; p4 <= 1.4; p4 += 0.2) p4_range.push_back(p4); 
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // 領域判定 (Code B Logic)
    bool isLowPt = (p_min < 3.9); 

    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError; // 一時的にエラー出力を抑制

    // --- Step 1: Grid Search ---
    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                double ka_init_val = isLowPt ? peak_height * 0.3 : peak_height * 0.4;
                double pr_init_val = isLowPt ? peak_height * 0.7 : peak_height * 0.8;

                double params[12] = {
                    pr_init_val, 0.0, 1.0,            // Proton 
                    ka_init_val, p4_init, 1.0,        // Kaon 
                    peak_height * 0.5, p7_init, 1.5,  // Pion
                    peak_height * 0.05, p10_init, 2.0 // Electron
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

    // --- Step 2: Final Fit ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    
    // --- Parameter Limits Tuning (Key for Stability) ---
    if (isLowPt) {
        // Protonが全てを食うのを防ぐ (上限92%)
        finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 0.92); 
        // Kaonの下限を保証 (下限30%)
        finalFit->SetParLimits(3, peak_height * 0.30, peak_height * 1.0);
    } else {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 1.5); 
        finalFit->SetParLimits(3, peak_height * 0.1, peak_height * 1.5);
    }
    finalFit->SetParLimits(1, -0.1, 0.1);
    finalFit->SetParLimits(2, 0.8, 1.15); // Sigma Constraint

    finalFit->SetParLimits(4, 0.7, 1.6); // Kaon Mean
    finalFit->SetParLimits(5, 0.8, 1.2); 

    finalFit->SetParLimits(6, 0.0, peak_height * 1.2); // Pion
    finalFit->SetParLimits(7, 3.0, 6.0); 
    finalFit->SetParLimits(8, 0.5, 2.0);

    finalFit->SetParLimits(9, 0.0, peak_height * 0.5); // Electron
    finalFit->SetParLimits(10, 7.0, 12.0);
    finalFit->SetParLimits(11, 0.6, 2.0);

    hist->Fit(finalFit, "QNRB0"); 

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max) {
    // Code B Logic
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    
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

//________________________________________________________________________________________________
// PID Map Calculation (Combination of Code A and B)
// Code AのJetSparse構造と、Code BのPID効率計算ロジックを統合
//________________________________________________________________________________________________

std::pair<std::map<double, double>, std::map<double, double>> getPidCorrectionMaps(
    const TString particleType, const TString fileName,
    double jetPtMin, double jetPtMax,
    double centMin, double centMax,
    double rMin, double rMax) 
{
    std::map<double, double> purityMap;
    std::map<double, double> pidEffMap; // Code B: PID Efficiency Map

    TString sparseName;
    double sigMin, sigMax;
    if (particleType == "proton") {
        sparseName = "jetTpcTofPr"; // 0:p, 1:nSigma, 2:distance, 3:jetPt, 4:cent
        sigMin = -3.5; sigMax = 0.5;
    } else {
        sparseName = "jetTpcTofPi"; 
        sigMin = -0.5; sigMax = 3.5;
    }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {purityMap, pidEffMap}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { inputFile->Close(); return {purityMap, pidEffMap}; }
    
    // 軸設定 (Code A: r-dependence)
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);         
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); 
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    const double p_start = 2.0;
    const double p_end = 6.5;
    const double p_step = 0.2;
    
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->SetDirectory(nullptr);
        
        double purity = 0.0;
        double pidEff = 1.0;

        if (h_nSigma->GetEntries() >= 50) {
            TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
            
            // --- Code B Logic: Purity + PID Efficiency Calculation ---
            TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
            signalFunc->SetParameters(bestFit->GetParameter(0), bestFit->GetParameter(1), bestFit->GetParameter(2));
            
            double intSigInCut = signalFunc->Integral(sigMin, sigMax); // カット内信号
            double intSigTotal = signalFunc->Integral(-10, 10);        // 全信号
            double intTotInCut = bestFit->Integral(sigMin, sigMax);    // カット内全粒子

            if (intTotInCut > 0) purity = intSigInCut / intTotInCut;
            if (intSigTotal > 0) pidEff = intSigInCut / intSigTotal;

            // 安全策
            purity = std::max(0.0, std::min(1.0, purity));
            pidEff = std::max(0.0, std::min(1.0, pidEff));

            delete signalFunc;
            delete bestFit;
        }
        
        purityMap[p_center] = purity;
        pidEffMap[p_center] = pidEff; // Store PID Eff
        delete h_nSigma;
    }
    inputFile->Close();
    return {purityMap, pidEffMap};
}

//________________________________________________________________________________________________
// Spectrum Weighting (Updated to use PID Efficiency)
//________________________________________________________________________________________________

TH1D* createWeightedPtSpectrum(TString fileName, TString particleType, 
                               const std::map<double, double>& purityMap,
                               const std::map<double, double>& pidEffMap, // Added
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax,
                               double rMin, double rMax) {
    TString sparseName;
    if (particleType == "proton") sparseName = "jetpVsPtForPr";
    else sparseName = "jetpVsPtForPi";

    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { inputFile->Close(); return nullptr; }
    
    // Axes: 0:P, 1:Pt, 2:Distance, 3:JetPt, 4:Centrality
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    // Projection(1, 0) -> X:Pt, Y:P
    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(1, 0); 
    h2_p_vs_pt->SetDirectory(nullptr);
    
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        
        double best_p_diff = 100.0;
        double purity = 0;
        double pidEff = 1.0;

        for(auto const& [p_val, pur_val] : purityMap) {
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
                pidEff = pidEffMap.at(p_val); // 対応するPID効率
            }
        }
        
        if (purity == 0 || best_p_diff > 0.2) continue; 
        if (pidEff < 0.01) pidEff = 1.0; // ゼロ除算防止

        // 補正係数: Purity / PID Efficiency
        double correctionFactor = purity / pidEff;

        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            h2_weighted->SetBinContent(j_pt, i_p, content * correctionFactor);
            h2_weighted->SetBinError(j_pt, i_p, error * correctionFactor);
        }
    }
    
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    h_pt_weighted->SetName(TString::Format("h_pt_weighted_%s_r%.1f_%.1f", sparseName.Data(), rMin, rMax));
    
    delete h2_p_vs_pt;
    delete h2_weighted;
    inputFile->Close();
    return h_pt_weighted;
}

TH1D* calculateTrackingEfficiency(TString particleType, TString fileName, 
                                  double jetPtMin, double jetPtMax, 
                                  double centMin, double centMax) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;

    THnSparseD *hRecoSparse = (THnSparseD*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
    THnSparseD *hGenSparse = (THnSparseD*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    
    if (!hRecoSparse || !hGenSparse) { file->Close(); return nullptr; }

    hRecoSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax);
    hRecoSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax);
    hGenSparse->GetAxis(2)->SetRangeUser(centMin, centMax);

    TH1D *hReco = hRecoSparse->Projection(0);
    hReco->SetDirectory(nullptr);
    TH1D *hGen = hGenSparse->Projection(0);
    hGen->SetDirectory(nullptr);
    file->Close();
    
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEfficiency_%s", particleType.Data()));
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    delete hReco; delete hGen;
    return hEff;
}

//________________________________________________________________________________________________
// Main Macro
//________________________________________________________________________________________________

void calculateRatioVsR_Scan() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    const TString kDataFileName = "AnalysisResults-83.root"; // Code Aのデータ
    const TString kMcFileName   = "AnalysisResults-91.root";

    const double kJetPtMin = 20.0; // Code Aの設定にあわせる
    const double kJetPtMax = 250.0;
    const double kCentMin  = 0.0;
    const double kCentMax  = 100.0;
    
    const double kTrackPtMinForRatio = 2.0; 
    const double kTrackPtMaxForRatio = 5.0;

    const int nRBins = 7;
    double rBins[nRBins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};

    TH1D* hRatioVsR = new TH1D("hRatioVsR", "", nRBins, rBins);
    hRatioVsR->Sumw2();

    std::cout << "Loading Efficiency maps..." << std::endl;
    TH1D* hEffPr = calculateTrackingEfficiency("Proton", kMcFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    TH1D* hEffPi = calculateTrackingEfficiency("Pion", kMcFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);

    if (!hEffPr || !hEffPi) { std::cerr << "Eff Error" << std::endl; return; }

    std::cout << "=== Starting r-scan Loop (" << nRBins << " bins) ===" << std::endl;

    for (int i = 0; i < nRBins; ++i) {
        double rMin = rBins[i];
        double rMax = rBins[i+1];
        
        std::cout << "\nProcessing bin " << i+1 << "/" << nRBins << " : r = [" << rMin << ", " << rMax << "]" << std::endl;

        // 1. Calculate Maps (Purity & PID Efficiency)
        auto mapsPr = getPidCorrectionMaps("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, rMin, rMax);
        auto mapsPi = getPidCorrectionMaps("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, rMin, rMax);

        // 2. Create Weighted Spectra (Correction = Purity / PID_Efficiency)
        TH1D* hSpecPr = createWeightedPtSpectrum(kDataFileName, "proton", mapsPr.first, mapsPr.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, rMin, rMax);
        TH1D* hSpecPi = createWeightedPtSpectrum(kDataFileName, "pion",   mapsPi.first, mapsPi.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, rMin, rMax);

        if (!hSpecPr || !hSpecPi) {
            std::cerr << "  Warning: Spectra creation failed." << std::endl;
            if (hSpecPr) delete hSpecPr;
            if (hSpecPi) delete hSpecPi;
            continue;
        }

        // 3. Efficiency Correction (Manual Division)
        auto divideManually = [](TH1D* hSpec, TH1D* hEff) {
            for (int b = 1; b <= hSpec->GetNbinsX(); ++b) {
                double pt = hSpec->GetBinCenter(b);
                double content = hSpec->GetBinContent(b);
                double error = hSpec->GetBinError(b);
                int binEff = hEff->FindBin(pt);
                if (hEff->IsBinUnderflow(binEff) || hEff->IsBinOverflow(binEff)) continue;
                double eff = hEff->GetBinContent(binEff);
                if (eff > 1e-4) { 
                    hSpec->SetBinContent(b, content / eff);
                    hSpec->SetBinError(b, error / eff); 
                } else {
                    hSpec->SetBinContent(b, 0.0); hSpec->SetBinError(b, 0.0);
                }
            }
        };
        divideManually(hSpecPr, hEffPr);
        divideManually(hSpecPi, hEffPi);

        // 4. Integrate Yields
        int binStartPr = hSpecPr->FindBin(kTrackPtMinForRatio + 0.001);
        int binEndPr   = hSpecPr->FindBin(kTrackPtMaxForRatio - 0.001);
        double errPr;
        double yieldPr = hSpecPr->IntegralAndError(binStartPr, binEndPr, errPr);

        int binStartPi = hSpecPi->FindBin(kTrackPtMinForRatio + 0.001);
        int binEndPi   = hSpecPi->FindBin(kTrackPtMaxForRatio - 0.001);
        double errPi;
        double yieldPi = hSpecPi->IntegralAndError(binStartPi, binEndPi, errPi);

        // 5. Calculate Ratio
        if (yieldPi > 0 && yieldPr > 0) {
            double ratio = yieldPr / yieldPi;
            double ratioErr = ratio * std::sqrt(std::pow(errPr/yieldPr, 2) + std::pow(errPi/yieldPi, 2));
            hRatioVsR->SetBinContent(i + 1, ratio);
            hRatioVsR->SetBinError(i + 1, ratioErr);
            printf("  Result: p/pi = %.4f +/- %.4f (Pr:%.0f, Pi:%.0f)\n", ratio, ratioErr, yieldPr, yieldPi);
        } else {
            std::cout << "  Warning: Yield is zero (Pr:" << yieldPr << ", Pi:" << yieldPi << ")" << std::endl;
        }
        delete hSpecPr; delete hSpecPi;
    }

    // =========================================================
    // 【追加部分】 結果をROOTファイルに保存
    // =========================================================
    TString outputFileName = "pp_reference.root";
    TFile *fOut = new TFile(outputFileName, "RECREATE");
    
    // ヒストグラムに分かりやすい名前をつけて保存
    hRatioVsR->SetName("hRatio_pp"); 
    hRatioVsR->Write();
    
    fOut->Close();
    
    std::cout << "Analysis saved to " << outputFileName << std::endl;

    // ... (描画や後処理) .

    // --- Drawing ---
    TCanvas *c1 = new TCanvas("cRatioScan", "Ratio vs r Scan", 800, 600);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.12);
    c1->SetGridy();

    hRatioVsR->SetTitle("");
    hRatioVsR->GetXaxis()->SetTitle("Distance from Jet Axis #it{r}");
    hRatioVsR->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hRatioVsR->GetYaxis()->SetRangeUser(0.0, 0.5); // ppなので比率は低めと予想されます
    
    hRatioVsR->SetMarkerStyle(kFullCircle);
    hRatioVsR->SetMarkerSize(1.2);
    hRatioVsR->SetMarkerColor(kBlue+1);
    hRatioVsR->SetLineColor(kBlue+1);
    hRatioVsR->SetLineWidth(2);
    hRatioVsR->Draw("E1");

    TPaveText *pt = new TPaveText(0.2, 0.75, 0.6, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextFont(42);
    pt->AddText("#bf{This thesis}");
    pt->AddText(TString::Format("pp #sqrt{s} = 13.6 TeV, %.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", kJetPtMin, kJetPtMax));
    pt->AddText("ch-particle jet, Anti-k_{T} R=0.4");
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    pt->AddText(TString::Format("%.1f < #it{p}_{T,track} < %.1f GeV/#it{c}", kTrackPtMinForRatio, kTrackPtMaxForRatio));
    pt->AddText("|#it{#eta}_{track}| < 0.9");
    pt->Draw();
    
    // c1->SaveAs("RatioVsR_pp.pdf");
    delete hEffPr; delete hEffPi;
}