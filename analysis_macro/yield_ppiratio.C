#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TPad.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <Rtypes.h> // For gErrorIgnoreLevel

//________________________________________________________________________________________________
// The following helper functions must be defined BEFORE they are called by the main function.
//________________________________________________________________________________________________

/**
* @brief プロトン用のグリッドサーチと最終フィットを行う関数
*/
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // ... (内容は変更なし)
    std::cout << "--- Starting grid search for PROTON in p bin " << p_min << "-" << p_max << " ---" << std::endl;
    std::vector<double> p4_range; for (double p4 = 1.5; p4 <= 3.0; p4 += 0.5) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                double peak_height = hist->GetMaximum();
                tempFit->SetParameter(0, peak_height); tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); tempFit->SetParLimits(2, 0.7, 2.0);
                tempFit->SetParameter(3, peak_height * 0.07); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.5);
                tempFit->SetParameter(6, peak_height * 0.7); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 2.0);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 2.5);
                hist->Fit(tempFit, "QNRB");
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
    std::cout << "Grid search finished for Proton. Best preliminary chi2/ndf: " << best_chi2ndf << std::endl;
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);
    double peak_height = hist->GetMaximum();
    finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);
    finalFit->SetParLimits(3, peak_height * 0.15, peak_height * 1.0);
    finalFit->SetParLimits(4, 1.5, 4.6);
    finalFit->SetParLimits(5, 0.6, 1.2);
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.3, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);
    hist->Fit(finalFit, "RBQ");
    return finalFit;
}

//________________________________________________________________________________________________
/**
* @brief パイオン用のグリッドサーチと最終フィットを行う関数
*/
TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // ... (内容は変更なし)
    std::cout << "--- Starting parameter optimization for PION in p bin " << p_min << "-" << p_max << " ---" << std::endl;
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                double peak_height = hist->GetMaximum();
                tempFit->SetParameter(0, peak_height); tempFit->SetParameter(1, 0.0); tempFit->SetParameter(2, 1.5);
                tempFit->SetParameter(3, peak_height * 0.1); tempFit->SetParameter(4, p4_init); tempFit->SetParameter(5, 1.0);
                tempFit->SetParameter(6, peak_height * 0.2); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 1.5);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 1.0);
                hist->Fit(tempFit, "QNRB");
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
    std::cout << "Grid search finished. Best preliminary chi2/ndf: " << best_chi2ndf << std::endl;
    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    double peak_height = hist->GetMaximum();
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
    hist->Fit(finalFit, "BQ");
    return finalFit;
}

//________________________________________________________________________________________________
/**
* @brief 運動量ビンごとの純度を計算し、Chi2/NDFも表示する
* ★変更: Jet pT, Centrality 引数を追加
*/
std::map<double, double> calculatePurities(const TString particleType, const TString fileName,
                                           double jetPtMin, double jetPtMax,
                                           double centMin, double centMax) {
    TString sparseName;
    double signal_min_cut, signal_max_cut;
    if (particleType == "proton") {
        sparseName = "tpcTofPr";
        signal_min_cut = -3.5;
        signal_max_cut = 0.5;
    } else if (particleType == "pion") {
        sparseName = "tpcTofPi";
        signal_min_cut = -0.5;
        signal_max_cut = 3.5;
    } else {
        std::cerr << "Error: Unknown particle type '" << particleType << "'" << std::endl;
        return {};
    }
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) {
        std::cerr << "Error: Could not find " << fullPath << " in " << fileName << std::endl;
        inputFile->Close();
        return {};
    }
    
    // Axis 0: p, Axis 1: nSigma, Axis 2: r
    // Axis 3: Jet pT, Axis 4: Centrality
    
    std::map<double, double> purities;
    const double p_start = 3.0;
    const double p_end = 7.0;
    const double p_step = 0.2;
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        // カットを適用
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        hSparse->GetAxis(2)->SetRangeUser(0.6, 0.7);
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); // ★ Jet pT
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   // ★ Centrality

        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->SetDirectory(nullptr);
        if (h_nSigma->GetEntries() < 50) {
            delete h_nSigma;
            continue;
        }
        TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
        signalFunc->SetParameters(&bestFit->GetParameters()[0]);
        double signal_in_cut = signalFunc->Integral(signal_min_cut, signal_max_cut);
        double total_in_cut = bestFit->Integral(signal_min_cut, signal_max_cut);
        double purity = (total_in_cut > 0) ? (signal_in_cut / total_in_cut) : 0;
        double chi2 = bestFit->GetChisquare();
        double ndf = bestFit->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
        purities[p_center] = purity;
        printf(" p = %.2f GeV/c | Purity (%-6s) = %.4f | Chi2/NDF = %.2f\n", p_center, particleType.Data(), purity, chi2ndf);
        delete h_nSigma;
        delete signalFunc;
        delete bestFit;
    }
    
    // 軸のリセット
    hSparse->GetAxis(0)->SetRange(1, hSparse->GetAxis(0)->GetNbins());
    hSparse->GetAxis(2)->SetRange(1, hSparse->GetAxis(2)->GetNbins());
    hSparse->GetAxis(3)->SetRange(1, hSparse->GetAxis(3)->GetNbins());
    hSparse->GetAxis(4)->SetRange(1, hSparse->GetAxis(4)->GetNbins());
    
    inputFile->Close();
    return purities;
}

//________________________________________________________________________________________________
/**
* @brief 純度で重み付けしたpTスペクトルを作成する関数
* ★変更: Jet pT, Centrality 引数を追加
*/
TH1D* createWeightedPtSpectrum(TString fileName, TString sparseName, const std::map<double, double>& purities,
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) {
        std::cerr << "Error: Could not find " << fullPath << " in " << fileName << std::endl;
        inputFile->Close();
        return nullptr;
    }
    
    // Axis 0: p, Axis 1: pT, Axis 2: r
    // Axis 3: Jet pT, Axis 4: Centrality
    
    // カットを適用
    hSparse->GetAxis(2)->SetRangeUser(0.6, 0.7);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); // ★ Jet pT
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   // ★ Centrality

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    
    // 軸のリセット
    hSparse->GetAxis(2)->SetRange(1, hSparse->GetAxis(2)->GetNbins());
    hSparse->GetAxis(3)->SetRange(1, hSparse->GetAxis(3)->GetNbins());
    hSparse->GetAxis(4)->SetRange(1, hSparse->GetAxis(4)->GetNbins());

    TH2D *h2_p_vs_pt_weighted = (TH2D*)h2_p_vs_pt->Clone();
    // 名前が重複しないように情報を付加
    TString histName = TString::Format("%s_weighted_jet%.0f_cent%.0f", sparseName.Data(), jetPtMin, centMin);
    h2_p_vs_pt_weighted->SetName(histName);
    h2_p_vs_pt_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double best_p_diff = std::numeric_limits<double>::max();
        double purity = 0;
        for(auto const& [p_val, pur_val] : purities) {
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
            }
        }
        if (purity == 0) continue;
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            h2_p_vs_pt_weighted->SetBinContent(j_pt, i_p, content * purity);
            h2_p_vs_pt_weighted->SetBinError(j_pt, i_p, error * purity);
        }
    }
    TH1D *h_pt_weighted = h2_p_vs_pt_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    h_pt_weighted->SetName(TString::Format("h_pt_%s", histName.Data()));
    delete h2_p_vs_pt;
    delete h2_p_vs_pt_weighted;
    inputFile->Close();
    return h_pt_weighted;
}

//________________________________________________________________________________________________
/**
* @brief トラッキング効率を計算する関数
*/
TH1D* calculateTrackingEfficiency(TString particleType, TString fileName) {
    // ... (内容は変更なし)
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }
    TString recoHistName = TString::Format("jet-shape-task/ptHistogram%s", particleType.Data());
    TString genHistName = TString::Format("jet-shape-task/ptGenerated%s", particleType.Data());
    TH1D *hReco = (TH1D*) file->Get(recoHistName);
    TH1D *hGen = (TH1D*) file->Get(genHistName);
    if (!hReco || !hGen) {
        std::cerr << "Error: Could not retrieve histograms for " << particleType
        << " (Tried: " << recoHistName << " and " << genHistName << ")" << std::endl;
        file->Close();
        return nullptr;
    }
    hReco->SetDirectory(nullptr);
    hGen->SetDirectory(nullptr);
    file->Close();
    hReco->Sumw2();
    hGen->Sumw2();
    TH1D* hEfficiency = (TH1D*)hReco->Clone(TString::Format("hEfficiency_%s", particleType.Data()));
    hEfficiency->SetDirectory(nullptr);
    hEfficiency->Divide(hReco, hGen, 1.0, 1.0, "B");
    return hEfficiency;
}

//________________________________________________________________________________________________
/**
* @brief 指定されたテキストを持つTPaveTextオブジェクトを作成するヘルパー関数
* ★変更: Jet pT, Centrality 引数を追加
*/
TPaveText* createPaveText(double x1, double y1, double x2, double y2,
                          double jetPtMin, double jetPtMax,
                          double centMin, double centMax) {
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("ALICE work in progress");
    //pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    //pt->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    pt->AddText("O-O #sqrt{#it{s}_{NN}} = 5.36 TeV");
    // ★ Centrality と Jet pT の情報を追加 (%%でエスケープ)
    pt->AddText(TString::Format("Centrality %.0f-%.0f%%", centMin, centMax));
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", jetPtMin, jetPtMax));
    pt->AddText("#it{p}_{T,track} > 0.15 GeV/#it{c}");
    return pt;
}

//________________________________________________________________________________________________
/**
* @brief メイン関数: p/pi 比を計算しプロットする
*/
void calculateProtonPionRatio() {
    Int_t oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;

    gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const TString kDataFileName = "AnalysisResults-70.root";
    const TString kMcFileName = "AnalysisResults-65.root";

    //const TString kDataFileName = "AnalysisResults-96.root";
    //const TString kMcFileName = "AnalysisResults-58.root";

    // --- ★ 新規: カットの設定 ---
    const double kJetPtMin = 140.0;
    const double kJetPtMax = 250.0; // 上限を十分に大きく設定
    const double kCentMin  = 0.0;
    const double kCentMax  = 10.0;

    std::cout << "\n=====================================================================" << std::endl;
    std::cout << "=== Jet pT cut: " << kJetPtMin << " - " << kJetPtMax << " GeV/c ===" << std::endl;
    std::cout << "=== Centrality cut: " << kCentMin << " - " << kCentMax << " % ===" << std::endl;
    std::cout << "=====================================================================\n" << std::endl;

    // --- 1. 純度計算 ---
    std::cout << "\n--- Calculating Proton & Pion Purities (using " << kDataFileName << ") ---" << std::endl;
    // ★引数追加
    std::map<double, double> protonPurities = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    std::map<double, double> pionPurities   = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax);

    // --- 純度のサマリー表を表示 ---
    std::cout << "\n\n--- Purity Calculation Summary ---" << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << " p Range (GeV/c) | Proton Purity | Pion Purity " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    const double p_step = 0.2;
    for (auto const& [p_center, proton_purity] : protonPurities) {
        double p_min = p_center - p_step / 2.0;
        double p_max = p_center + p_step / 2.0;
        double pion_purity = pionPurities.count(p_center) ? pionPurities.at(p_center) : 0.0;
        printf(" %.2f - %.2f | %.4f | %.4f\n", p_min, p_max, proton_purity, pion_purity);
    }
    std::cout << "===================================================\n" << std::endl;

    // --- 2. 純度で重み付けしたスペクトル作成 ---
    std::cout << "--- Creating Purity-Weighted pT Spectra (using " << kDataFileName << ") ---" << std::endl;
    // ★引数追加
    TH1D* hPtProton_purity = createWeightedPtSpectrum(kDataFileName, "pVsPtForPr", protonPurities, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    TH1D* hPtPion_purity   = createWeightedPtSpectrum(kDataFileName, "pVsPtForPi",   pionPurities, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    
    if (!hPtProton_purity || !hPtPion_purity) {
        std::cerr << "Error: Failed to create spectra. Check if data exists for these cuts." << std::endl;
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        return;
    }

    // --- 3. 効率計算 ---
    std::cout << "\n--- Calculating Tracking Efficiencies (using " << kMcFileName << ") ---" << std::endl;
    // Note: Efficiency histograms in jetShape.cxx are typically 1D (Inclusive). 
    // If exact centrality-dependent efficiency is needed, MC file must contain such histograms.
    TH1D* hEffProton = calculateTrackingEfficiency("Proton", kMcFileName);
    TH1D* hEffPion = calculateTrackingEfficiency("Pion", kMcFileName);
    if (!hEffProton || !hEffPion) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        return;
    }

    // --- 4. 効率で収量を補正 ---
    std::cout << "--- Applying Efficiency Correction to Yields ---" << std::endl;
    TH1D* hPtProton_corrected = (TH1D*)hPtProton_purity->Clone("hPtProton_corrected");
    hPtProton_corrected->Divide(hEffProton);
    TH1D* hPtPion_corrected = (TH1D*)hPtPion_purity->Clone("hPtPion_corrected");
    hPtPion_corrected->Divide(hEffPion);

    // --- 5. 効率補正後のグラフを描画 ---
    // ★ createPaveText にもカット情報を渡す
    
    TCanvas *c_pi = new TCanvas("c_pion_pt_final", "Pion pT Spectrum (Corrected)", 800, 600);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15);
    hPtPion_corrected->SetTitle("Pion #it{p}_{T} Spectrum");
    hPtPion_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtPion_corrected->GetYaxis()->SetTitle("Corrected Counts");
    hPtPion_corrected->SetLineColor(kBlue + 1); hPtPion_corrected->SetMarkerColor(kBlue + 1); hPtPion_corrected->SetMarkerStyle(kFullCircle);
    hPtPion_corrected->Draw("PE");
    TPaveText* pave_pi = createPaveText(0.5, 0.65, 0.88, 0.88, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    pave_pi->Draw();

    TCanvas *c_pr = new TCanvas("c_proton_pt_final", "Proton pT Spectrum (Corrected)", 800, 600);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15);
    hPtProton_corrected->SetTitle("Proton #it{p}_{T} Spectrum");
    hPtProton_corrected->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPtProton_corrected->GetYaxis()->SetTitle("Corrected Counts");
    hPtProton_corrected->SetLineColor(kRed + 1); hPtProton_corrected->SetMarkerColor(kRed + 1); hPtProton_corrected->SetMarkerStyle(kFullSquare);
    hPtProton_corrected->Draw("PE");
    TPaveText* pave_pr = createPaveText(0.5, 0.65, 0.88, 0.88, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    pave_pr->Draw();

    TCanvas *c_overlaid = new TCanvas("c_overlaid_pt_final", "Overlaid pT Spectra (Corrected)", 800, 600);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15);
    double max_y = std::max(hPtPion_corrected->GetMaximum(), hPtProton_corrected->GetMaximum());
    hPtPion_corrected->GetYaxis()->SetRangeUser(hPtPion_corrected->GetMinimum(0.1), max_y * 1.5);
    hPtPion_corrected->SetTitle("Purity & Efficiency Corrected #it{p}_{T} Spectra");
    hPtPion_corrected->Draw("PE");
    hPtProton_corrected->Draw("PE SAME");
    TLegend* leg_spectra = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg_spectra->SetBorderSize(0); leg_spectra->SetFillStyle(0);
    leg_spectra->AddEntry(hPtPion_corrected, "Pions", "lpe");
    leg_spectra->AddEntry(hPtProton_corrected, "Protons", "lpe");
    leg_spectra->Draw();
    TPaveText* pave_overlaid = createPaveText(0.18, 0.7, 0.55, 0.88, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    pave_overlaid->Draw();

    // --- 6. 比を計算してプロット (補正後) ---
    hPtProton_corrected->Rebin(2);
    hPtPion_corrected->Rebin(2);
    TH1D* hRatio = (TH1D*)hPtProton_corrected->Clone("hRatio_p_pi_final");
    hRatio->Divide(hPtPion_corrected);
    
    TCanvas *c_ratio = new TCanvas("c_ratio_final", "Proton/Pion Ratio (Corrected)", 800, 600);
    c_ratio->cd();
    gPad->SetGridy();
    gPad->SetLeftMargin(0.12);
    hRatio->SetTitle("p/#pi Ratio");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(kFullCross);
    hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hRatio->GetYaxis()->SetTitleOffset(1.2);
    hRatio->GetYaxis()->SetRangeUser(0.0, 1.1);
    hRatio->Draw("E1");
    TPaveText* pave_ratio = createPaveText(0.15, 0.75, 0.5, 0.9, kJetPtMin, kJetPtMax, kCentMin, kCentMax);
    pave_ratio->Draw();

    // 抑制していたメッセージレベルを元に戻す
    gErrorIgnoreLevel = oldErrorIgnoreLevel;
}