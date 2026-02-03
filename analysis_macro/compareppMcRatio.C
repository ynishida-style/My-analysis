// File Name: CompareRatios.C
// Run command: root -l CompareRatios.C

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <string>

// ヘルパー関数: ファイル名とスタイルを指定して比率(p/pi)ヒストグラムを返す
TH1D* GetRatioHistogram(std::string fileName, std::string label, int color, int markerStyle) {
    TFile *f = TFile::Open(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cout << "Error: File not found -> " << fileName << std::endl;
        return nullptr;
    }

    TH1D *hp = (TH1D*)f->Get("h_num_p");
    TH1D *hpi = (TH1D*)f->Get("h_num_pi");

    if (!hp || !hpi) {
        std::cout << "Error: Histograms not found in " << fileName << std::endl;
        f->Close();
        return nullptr;
    }

    // ヒストグラムのクローン作成（名前が被らないようにLabelを付加）
    TH1D *hRatio = (TH1D*)hp->Clone(Form("hRatio_%s", label.c_str()));
    
    // ディレクトリから切り離す（ファイルを閉じてもヒストグラムがメモリに残るようにする）
    hRatio->SetDirectory(0); 

    // 比率計算: p / pi
    hRatio->Divide(hpi);

    // スタイル設定
    hRatio->SetLineColor(color);
    hRatio->SetMarkerColor(color);
    hRatio->SetMarkerStyle(markerStyle);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetLineWidth(2);
    hRatio->SetTitle(label.c_str()); // 凡例用タイトル

    f->Close(); // ファイルを閉じる
    return hRatio;
}

void CompareRatios() {
    // スタイル設定
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // ---------------------------------------------------------
    // 1. 各パターンのヒストグラムを取得
    // ---------------------------------------------------------
    // 色コード: kBlack(黒), kRed(赤), kBlue(青)
    // マーカー: 24(白抜き丸), 20(塗りつぶし丸), 21(四角) など
    
    // CR Off (Baseline) -> 黒
    TH1D *hCRoff = GetRatioHistogram("Results_CRoff.root", "CR Off", kBlack, 24);
    
    // CR Mode 0 (Default MPI-based) -> 赤
    TH1D *hMode0 = GetRatioHistogram("Results_Mode0.root", "CR Mode 0 (Default)", kRed, 20);

    // CR Mode 2 (Gluon-move / Baryon Enhancement) -> 青
    TH1D *hMode2 = GetRatioHistogram("Results_Mode2.root", "CR Mode 2 (Gluon Move)", kBlue, 21);

    if (!hCRoff || !hMode0 || !hMode2) return; // 読み込み失敗時は終了

    // ---------------------------------------------------------
    // 2. 描画設定
    // ---------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "p/pi Ratio Comparison", 800, 600);
    gPad->SetGridy(); // Y軸グリッドだけあると比率が見やすい

    // 軸のタイトルなどを設定（代表してhCRoffに行わせる）
    hCRoff->SetTitle("");
    hCRoff->GetXaxis()->SetTitle("Distance from Jet Axis #it{r}");
    hCRoff->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hCRoff->GetXaxis()->SetRangeUser(0.0, 0.7);
    hCRoff->GetYaxis()->SetRangeUser(0.0, 0.6); // データに合わせて調整してください

    // 重ね書き (Draw "SAME")
    // 描画順序: エラーバー(E)と点(P)または曲線(C)など
    hCRoff->Draw("EP");      // 最初に座標軸ごと描画
    hMode0->Draw("EP SAME"); // 重ねて描画
    hMode2->Draw("EP SAME"); 

    // ---------------------------------------------------------
    // 3. 凡例 (Legend)
    // ---------------------------------------------------------
    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.85); // 右上に配置
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->SetHeader("PYTHIA 8 pp #sqrt{s} = 13.6 TeV");
    
    leg->AddEntry(hCRoff, "CR Off", "lep");
    leg->AddEntry(hMode0, "CR Mode 0 (MPI-based)", "lep");
    leg->AddEntry(hMode2, "CR Mode 2 (Gluon-move)", "lep");
    
    leg->Draw();

    // 保存
    // c1->SaveAs("Comparison_Result.pdf");
}