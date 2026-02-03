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
#include <Rtypes.h>

//________________________________________________________________________________________________
// Helper Functions
//________________________________________________________________________________________________

/**
* @brief プロトン用のグリッドサーチと最終フィットを行う関数
*/
TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max)
{
    // ... (前回の内容と同様)
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
    // std::cout << "Grid search finished for Proton. Best preliminary chi2/ndf: " << best_chi2ndf << std::endl;
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
    // ... (前回の内容と同様)
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
    // std::cout << "Grid search finished. Best preliminary chi2/ndf: " << best_chi2ndf << std::endl;
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
    
    std::map<double, double> purities;
    const double p_start = 3.0;
    const double p_end = 7.0;
    const double p_step = 0.2;
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        hSparse->GetAxis(2)->SetRangeUser(0.6, 0.7);
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); // Jet pT
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   // Centrality

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
        
        purities[p_center] = purity;
        
        delete h_nSigma;
        delete signalFunc;
        delete bestFit;
    }
    inputFile->Close();
    return purities;
}

//________________________________________________________________________________________________
/**
* @brief 純度で重み付けしたpTスペクトルを作成する関数
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
    
    hSparse->GetAxis(2)->SetRangeUser(0.6, 0.7);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); // Jet pT
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);   // Centrality

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    
    // 重み付け用ヒストグラム
    TH2D *h2_p_vs_pt_weighted = (TH2D*)h2_p_vs_pt->Clone();
    TString histName = TString::Format("%s_weighted_jet%.0f_cent%.0f_%.0f", sparseName.Data(), jetPtMin, centMin, centMax);
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
* @brief トラッキング効率を計算する関数 (Inclusiveと仮定)
*/
TH1D* calculateTrackingEfficiency(TString particleType, TString fileName) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return nullptr;
    }
    
    // ※ MCファイル内のヒストグラム名は環境に合わせて確認・修正してください
    // 一般的な名前: "jet-shape-task/mcTpcTofPr" / "jet-shape-task/mcGenPr" など
    // もしくは TH1D が直接保存されている場合:
    TString recoHistName, genHistName;
    if (particleType == "Proton") {
        recoHistName = "jet-shape-task/mcTpcTofPr"; // または "jet-shape-task/ptHistogramProton"
        genHistName  = "jet-shape-task/mcGenPr";    // または "jet-shape-task/ptGeneratedProton"
    } else {
        recoHistName = "jet-shape-task/mcTpcTofPi";
        genHistName  = "jet-shape-task/mcGenPi";
    }

    // THnSparseD の場合を想定して Projection する例
    THnSparseD *hSparseRec = (THnSparseD*)file->Get(recoHistName);
    THnSparseD *hSparseGen = (THnSparseD*)file->Get(genHistName);
    
    TH1D *hReco = nullptr;
    TH1D *hGen = nullptr;

    if (hSparseRec && hSparseGen) {
        hReco = hSparseRec->Projection(1); // 1: pT
        hGen  = hSparseGen->Projection(1);
    } else {
        // 見つからない場合は TH1D として再試行 (以前のコードとの互換性)
        hReco = (TH1D*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
        hGen  = (TH1D*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    }

    if (!hReco || !hGen) {
        std::cerr << "Error: Could not retrieve efficiency histograms for " << particleType << std::endl;
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
    hEfficiency->Divide(hReco, hGen, 1.0, 1.0, "B"); // B: Binomial errors
    
    delete hSparseRec; delete hSparseGen;
    delete hReco; delete hGen;
    
    return hEfficiency;
}

//________________________________________________________________________________________________
/**
* @brief メイン関数: 複数のCentralityビンでp/pi比を計算し重ねてプロットする
*/
void calculateProtonPionRatio() {
    Int_t oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;

    gROOT->SetBatch(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const TString kDataFileName = "AnalysisResults-70.root";
    const TString kMcFileName = "AnalysisResults-65.root";
    const double kJetPtMin = 30.0;
    const double kJetPtMax = 250.0;

    // --- ★ 新規: Centrality ビンの定義 ---
    struct CentBin {
        double min;
        double max;
        int color;
        int marker;
    };

    std::vector<CentBin> centBins = {
        {0.0,  10.0, kBlack,      kFullCircle},
        {10.0, 30.0, kRed + 1,    kFullSquare},
        {30.0, 50.0, kGreen + 2,  kFullTriangleUp},
        {50.0, 100.0, kBlue + 1,  kFullTriangleDown}
    };

    // 結果のヒストグラムを保存するコンテナ
    std::vector<TH1D*> ratioHists;

    // --- 効率計算 (MCファイル) ---
    // ※効率がCentralityに依存しない(Inclusive)場合の処理です。
    // もしMCもCentralityごとに分かれている場合はループ内で計算してください。
    std::cout << "\n--- Calculating Tracking Efficiencies (Inclusive) ---" << std::endl;
    TH1D* hEffProton = calculateTrackingEfficiency("Proton", kMcFileName);
    TH1D* hEffPion   = calculateTrackingEfficiency("Pion", kMcFileName);
    
    if (!hEffProton || !hEffPion) {
        std::cerr << "Error: Failed to calculate efficiency." << std::endl;
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        return;
    }

    // --- ループ処理開始 ---
    for (const auto& bin : centBins) {
        std::cout << "\n=====================================================================" << std::endl;
        std::cout << "=== Processing Centrality: " << bin.min << " - " << bin.max << " % ===" << std::endl;
        std::cout << "=====================================================================\n" << std::endl;

        // 1. 純度計算
        std::map<double, double> protonPurities = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, bin.min, bin.max);
        std::map<double, double> pionPurities   = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, bin.min, bin.max);

        // 2. 純度で重み付けしたスペクトル作成
        TH1D* hPtProton = createWeightedPtSpectrum(kDataFileName, "pVsPtForPr", protonPurities, kJetPtMin, kJetPtMax, bin.min, bin.max);
        TH1D* hPtPion   = createWeightedPtSpectrum(kDataFileName, "pVsPtForPi",   pionPurities, kJetPtMin, kJetPtMax, bin.min, bin.max);

        if (!hPtProton || !hPtPion) {
            std::cerr << "Skipping this bin due to error." << std::endl;
            continue;
        }

        // 3. 効率補正
        // Cloneして補正用ヒストグラムを作成 (名前が被らないようにする)
        TH1D* hPtProton_corr = (TH1D*)hPtProton->Clone(TString::Format("hPtProton_corr_%.0f_%.0f", bin.min, bin.max));
        TH1D* hPtPion_corr   = (TH1D*)hPtPion->Clone(TString::Format("hPtPion_corr_%.0f_%.0f", bin.min, bin.max));
        
        hPtProton_corr->Divide(hEffProton);
        hPtPion_corr->Divide(hEffPion);

        // 4. 比の計算 (Rebin処理含む)
        hPtProton_corr->Rebin(2);
        hPtPion_corr->Rebin(2);

        TH1D* hRatio = (TH1D*)hPtProton_corr->Clone(TString::Format("hRatio_%.0f_%.0f", bin.min, bin.max));
        hRatio->Divide(hPtPion_corr);

        // 5. スタイル設定と保存
        hRatio->SetLineColor(bin.color);
        hRatio->SetMarkerColor(bin.color);
        hRatio->SetMarkerStyle(bin.marker);
        hRatio->SetLineWidth(2);
        hRatio->SetStats(0);
        
        // 凡例用にタイトルにCentrality情報を入れておく（後で凡例作成時に使用）
        hRatio->SetTitle(TString::Format("Cent. %.0f-%.0f%%", bin.min, bin.max));

        ratioHists.push_back(hRatio);
    }

    // --- 重ね書きプロット ---
    if (ratioHists.empty()) {
        std::cerr << "No histograms to plot." << std::endl;
        return;
    }

    TCanvas *c_overlay = new TCanvas("c_overlay", "Proton/Pion Ratio vs Centrality", 900, 700);
    c_overlay->cd();
    gPad->SetGridy();
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);

    // 凡例の作成
    TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    // 最初のヒストグラムを描画して枠を作る
    TH1D* firstHist = ratioHists[0];
    firstHist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    firstHist->GetYaxis()->SetTitle("Proton / Pion Ratio");
    firstHist->GetYaxis()->SetTitleOffset(1.2);
    // Y軸の範囲を適宜調整
    firstHist->GetYaxis()->SetRangeUser(0.0, 1.2); 
    
    // タイトルを消して、きれいなグラフにする
    firstHist->SetTitle(""); 

    firstHist->Draw("E1");
    leg->AddEntry(firstHist, TString::Format("%.0f-%.0f%%", centBins[0].min, centBins[0].max), "lp");

    // 残りのヒストグラムを重ねる
    for (size_t i = 1; i < ratioHists.size(); ++i) {
        ratioHists[i]->Draw("SAME E1");
        leg->AddEntry(ratioHists[i], TString::Format("%.0f-%.0f%%", centBins[i].min, centBins[i].max), "lp");
    }

    leg->Draw();

    // テキストボックスの追加 (共通情報)
    TPaveText *pt = new TPaveText(0.15, 0.75, 0.5, 0.88, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("ALICE work in progress");
    pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV"); // システム情報は必要に応じて修正してください
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", kJetPtMin, kJetPtMax));
    pt->Draw();

    gErrorIgnoreLevel = oldErrorIgnoreLevel;
}