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
#include <TLine.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <Rtypes.h>

//________________________________________________________________________________________________
// ヘルパー関数群
//________________________________________________________________________________________________

/**
* @brief プロトン用のグリッドサーチと最終フィットを行う関数
*/
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // ★警告抑制: Infoレベルのメッセージを消す
    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError; 

    std::vector<double> p4_range; for (double p4 = 1.5; p4 <= 3.0; p4 += 0.5) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);

    // グリッドサーチ (Quietモード)
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
                
                // "Q": Quiet (最低限), "N": No draw
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

    // 最終フィット
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
    hist->Fit(finalFit, "QRB"); // Qオプション追加

    // 警告レベルを戻す
    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

/**
* @brief パイオン用のグリッドサーチと最終フィットを行う関数
*/
TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max)
{
    // ★警告抑制
    Int_t oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

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
    hist->Fit(finalFit, "QB"); // Qオプション追加
    
    // 警告レベルを戻す
    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

/**
* @brief 運動量ビンごとの純度を計算 (進捗表示付き)
*/
std::map<double, double> calculatePurities(const TString particleType, const TString fileName,
                                           double jetPtMin, double jetPtMax,
                                           double centMin, double centMax,
                                           double rMin, double rMax) {
    TString sparseName;
    double signal_min_cut, signal_max_cut;
    if (particleType == "proton") {
        sparseName = "tpcTofPr";
        signal_min_cut = -3.5; signal_max_cut = 0.5;
    } else if (particleType == "pion") {
        sparseName = "tpcTofPi";
        signal_min_cut = -0.5; signal_max_cut = 3.5;
    } else {
        return {};
    }
    
    std::cout << ">>> Start Purity Calc for " << particleType << "..." << std::endl;

    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { std::cerr << "File Error: " << fileName << std::endl; return {}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { inputFile->Close(); return {}; }
    
    std::map<double, double> purities;
    const double p_start = 2.0; 
    const double p_end = 8.0; // ★変更: 7 GeV/c を超える範囲まで計算 (2.0 ~ 8.0)
    const double p_step = 0.2;
    
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        // 進捗と純度計算の様子を表示
        // std::cout << "\r Processing p = " << p_center << " GeV/c ... " << std::flush;

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        hSparse->GetAxis(2)->SetRangeUser(rMin, rMax); 
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

        TH1D *h_nSigma = hSparse->Projection(1);
        h_nSigma->SetName(TString::Format("h_sigma_%s_p%.1f", particleType.Data(), p_center));
        h_nSigma->SetDirectory(nullptr);
        
        if (h_nSigma->GetEntries() < 50) { 
            delete h_nSigma;
            // エントリ不足の場合は0にするが、ログは残す
            // printf("\n [Info] Low stats at p=%.1f\n", p_center);
            purities[p_center] = 0.0;
            continue;
        }
        
        TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
        
        TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
        signalFunc->SetParameters(&bestFit->GetParameters()[0]);
        double signal_in_cut = signalFunc->Integral(signal_min_cut, signal_max_cut);
        double total_in_cut = bestFit->Integral(signal_min_cut, signal_max_cut);
        double purity = (total_in_cut > 0) ? (signal_in_cut / total_in_cut) : 0;
        
        purities[p_center] = purity;
        
        // ★変更: 純度の計算結果を表示
        printf(" p = %.1f | Purity = %.3f | Chi2/NDF = %.1f\n", p_center, purity, (bestFit->GetNDF()>0 ? bestFit->GetChisquare()/bestFit->GetNDF() : 0));

        delete h_nSigma; delete signalFunc; delete bestFit;
    }
    
    std::cout << ">>> Finished Purity Calc for " << particleType << ".\n" << std::endl;

    for(int i=0; i<5; i++) hSparse->GetAxis(i)->SetRange(0, 0);
    inputFile->Close();
    return purities;
}

/**
* @brief Purity重み付けなしのRawスペクトルを作成する関数
*/
TH1D* createRawPtSpectrum(TString fileName, TString sparseName,
                          double jetPtMin, double jetPtMax,
                          double centMin, double centMax,
                          double rMin, double rMax) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { inputFile->Close(); return nullptr; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    TH1D *h_pt_raw = hSparse->Projection(1);
    h_pt_raw->SetDirectory(nullptr);
    TString histName = TString::Format("h_pt_raw_%s_c%.0f_r%.1f", sparseName.Data(), centMin, rMin);
    h_pt_raw->SetName(histName);
    h_pt_raw->SetTitle(histName);
    
    delete hSparse;
    inputFile->Close();
    return h_pt_raw;
}

/**
* @brief 純度で重み付けしたpTスペクトルを作成
*/
TH1D* createWeightedPtSpectrum(TString fileName, TString sparseName, const std::map<double, double>& purities,
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax,
                               double rMin, double rMax) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { inputFile->Close(); return nullptr; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted_temp");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double purity = 0;
        double min_dist = 1.0;
        for(auto const& [p_val, pur_val] : purities) {
            if (std::abs(p_val - p_center) < min_dist) {
                min_dist = std::abs(p_val - p_center);
                purity = pur_val;
            }
        }
        if (min_dist > 0.3) purity = 0;
        if (purity <= 0) continue;
        
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            h2_weighted->SetBinContent(j_pt, i_p, content * purity);
            h2_weighted->SetBinError(j_pt, i_p, error * purity);
        }
    }
    
    TH1D *h_pt = h2_weighted->ProjectionX();
    h_pt->SetDirectory(nullptr);
    h_pt->SetName(TString::Format("h_pt_weighted_%s", sparseName.Data()));
    
    delete h2_p_vs_pt; delete h2_weighted; delete hSparse;
    inputFile->Close();
    return h_pt;
}

/**
* @brief トラッキング効率を計算する関数
*/
TH1D* calculateTrackingEfficiency(TString particleType, TString fileName) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    
    TH1D *hReco = (TH1D*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
    TH1D *hGen = (TH1D*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    
    if (!hReco || !hGen) { file->Close(); return nullptr; }
    
    hReco->SetDirectory(nullptr); hGen->SetDirectory(nullptr);
    file->Close();
    
    hReco->Sumw2(); hGen->Sumw2();
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEff_%s", particleType.Data()));
    hEff->SetDirectory(nullptr);
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    delete hReco; delete hGen;
    return hEff;
}


//________________________________________________________________________________________________
// メイン関数: 補正の各段階を比較表示する (3分割レイアウト)
//________________________________________________________________________________________________

void compareCorrectionSteps() {
    // --- Settings ---
    const TString kDataFileName = "AnalysisResults-70.root";
    const TString kMcFileName   = "AnalysisResults-65.root";
    
    const double kJetPtMin = 30.0;
    const double kJetPtMax = 250.0;
    const double kCentMin  = 0.0;
    const double kCentMax  = 10.0;
    const double kRMin     = 0.0;
    const double kRMax     = 0.7; 

    // グローバル設定
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kFALSE); // ウィンドウを表示する
    
    std::cout << "=== Comparing Correction Steps ===" << std::endl;

    // 1. Calculate Purities
    // ログ抑制のため、calculatePurities内部で出力制御済み
    auto pPurities = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);
    auto piPurities = calculatePurities("pion", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);

    // 2. Load Efficiencies
    std::cout << "Loading Efficiencies..." << std::endl;
    TH1D* effPr = calculateTrackingEfficiency("Proton", kMcFileName);
    TH1D* effPi = calculateTrackingEfficiency("Pion", kMcFileName);

    if (!effPr || !effPi) { std::cerr << "Efficiency load failed." << std::endl; return; }

    // --- Prepare Histograms ---
    // Proton
    TH1D* hPr_Raw = createRawPtSpectrum(kDataFileName, "tpcTofPr", kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);
    TH1D* hPr_Pur = createWeightedPtSpectrum(kDataFileName, "pVsPtForPr", pPurities, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);
    TH1D* hPr_Final = (TH1D*)hPr_Pur->Clone("hPr_Final");
    hPr_Final->Divide(effPr);

    // Pion
    TH1D* hPi_Raw = createRawPtSpectrum(kDataFileName, "tpcTofPi", kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);
    TH1D* hPi_Pur = createWeightedPtSpectrum(kDataFileName, "pVsPtForPi", piPurities, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax);
    TH1D* hPi_Final = (TH1D*)hPi_Pur->Clone("hPi_Final");
    hPi_Final->Divide(effPi);

    // --- Drawing : 3分割レイアウト ---
    // キャンバス作成 (横長)
    TCanvas *cCheck = new TCanvas("cCheck", "Correction Steps Comparison", 1500, 500);
    cCheck->Divide(3, 1); // 3列1行

    // --- Common Styles ---
    // Proton: Red Square, Pion: Blue Circle
    auto setStyle = [](TH1D* h, int color, int marker) {
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(marker);
        h->SetLineWidth(2);
        h->SetMarkerSize(1.0);
    };
    
    setStyle(hPr_Raw, kRed+1, kFullSquare);
    setStyle(hPi_Raw, kBlue+1, kFullCircle);
    
    setStyle(hPr_Pur, kRed+1, kFullSquare);
    setStyle(hPi_Pur, kBlue+1, kFullCircle);
    
    setStyle(hPr_Final, kRed+1, kFullSquare);
    setStyle(hPi_Final, kBlue+1, kFullCircle);

    double xMin = 2.0;
    double xMax = 8.0; // ★ 8GeVまで表示

    // ==========================================
    // Pad 1: Raw Spectra (重ね書き)
    // ==========================================
    cCheck->cd(1);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
    
    hPi_Raw->SetTitle("Step 1: Raw Spectra (Counts)");
    hPi_Raw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPi_Raw->GetYaxis()->SetTitle("Raw Counts");
    hPi_Raw->GetXaxis()->SetRangeUser(xMin, xMax);
    
    // Y軸調整 (Logで見えるように)
    double maxRaw = std::max(hPi_Raw->GetMaximum(), hPr_Raw->GetMaximum());
    hPi_Raw->SetMaximum(maxRaw * 5.0);
    hPi_Raw->SetMinimum(1.0); // Logの最小値

    hPi_Raw->Draw("PE");
    hPr_Raw->Draw("PE SAME"); // 重ね書き

    TLegend *leg1 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->AddEntry(hPi_Raw, "Pion (Raw)", "lp");
    leg1->AddEntry(hPr_Raw, "Proton (Raw)", "lp");
    leg1->Draw();

    // ==========================================
    // Pad 2: Purity Corrected (重ね書き)
    // ==========================================
    cCheck->cd(2);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);

    hPi_Pur->SetTitle("Step 2: Purity Corrected");
    hPi_Pur->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPi_Pur->GetYaxis()->SetTitle("Counts #times Purity");
    hPi_Pur->GetXaxis()->SetRangeUser(xMin, xMax);
    
    double maxPur = std::max(hPi_Pur->GetMaximum(), hPr_Pur->GetMaximum());
    hPi_Pur->SetMaximum(maxPur * 5.0);
    hPi_Pur->SetMinimum(0.1);

    hPi_Pur->Draw("PE");
    hPr_Pur->Draw("PE SAME");

    TLegend *leg2 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->AddEntry(hPi_Pur, "Pion (Pur)", "lp");
    leg2->AddEntry(hPr_Pur, "Proton (Pur)", "lp");
    leg2->Draw();

    // ==========================================
    // Pad 3: Fully Corrected (重ね書き)
    // ==========================================
    cCheck->cd(3);
    gPad->SetLogy(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);

    hPi_Final->SetTitle("Step 3: Fully Corrected (Yield)");
    hPi_Final->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPi_Final->GetYaxis()->SetTitle("Corrected Yield");
    hPi_Final->GetXaxis()->SetRangeUser(xMin, xMax);

    double maxFin = std::max(hPi_Final->GetMaximum(), hPr_Final->GetMaximum());
    hPi_Final->SetMaximum(maxFin * 5.0);
    hPi_Final->SetMinimum(maxFin * 1e-4);

    hPi_Final->Draw("PE");
    hPr_Final->Draw("PE SAME");

    TLegend *leg3 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0);
    leg3->AddEntry(hPi_Final, "Pion (Final)", "lp");
    leg3->AddEntry(hPr_Final, "Proton (Final)", "lp");
    leg3->Draw();

    // 全体タイトル的な情報を描画
    cCheck->cd(0);
    TPaveText *pt = new TPaveText(0.3, 0.94, 0.7, 0.99, "NDC");
    pt->SetBorderSize(0); pt->SetFillStyle(0);
    pt->AddText(TString::Format("Comparison Steps: Cent %.0f-%.0f%%, %.1f < r < %.1f", kCentMin, kCentMax, kRMin, kRMax));
    pt->Draw();
}