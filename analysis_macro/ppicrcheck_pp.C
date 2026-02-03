void PlotRatio() {
    // スタイル設定（余白やフォントサイズなどを見やすくする）
    gStyle->SetOptStat(0); // 統計ボックスを非表示
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // 1. ファイルを開く
    TFile *f = TFile::Open("Results_Pythia_pp.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error: File not found!" << std::endl;
        return;
    }

    // 2. ヒストグラムを取得
    TH1D *hp = (TH1D*)f->Get("h_num_p");
    TH1D *hpi = (TH1D*)f->Get("h_num_pi");

    if (!hp || !hpi) {
        std::cout << "Error: Histograms not found!" << std::endl;
        return;
    }

    // 3. 比率の計算 (p / pi)
    // Cloneして新しいヒストグラムを作成（名前を "hRatio" に）
    TH1D *hRatio = (TH1D*)hp->Clone("hRatio");
    hRatio->Divide(hpi); // hp / hpi を計算

    // 4. 見た目の設定
    hRatio->SetTitle(""); // タイトルは空（必要なら入れる）
    
    // 軸ラベル
    hRatio->GetXaxis()->SetTitle("Distance from Jet Axis r");
    hRatio->GetYaxis()->SetTitle("Proton / Pion Ratio");
    
    // 範囲設定
    hRatio->GetXaxis()->SetRangeUser(0.0, 0.7);
    hRatio->GetYaxis()->SetRangeUser(0.0, 0.6); // 論文の図に合わせて調整

    // マーカー設定
    hRatio->SetMarkerStyle(20); // 塗りつぶし丸
    hRatio->SetMarkerSize(1.2);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);

    // 5. 描画
    TCanvas *c1 = new TCanvas("c1", "p/pi Ratio", 800, 600);
    
    // グリッド線を入れると見やすい
    gPad->SetGridy();
    
    hRatio->Draw("EP"); // E:エラーバー, P:マーカー

    // 凡例 (Legend) の追加
    TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85); // 位置 (x1, y1, x2, y2)
    leg->SetBorderSize(0);
    leg->AddEntry(hRatio, "PYTHIA 8 (pp #sqrt{s}=13.6 TeV)", "lep");
    // leg->AddEntry((TObject*)0, "CR Mode: Gluon Move", ""); // 設定をメモしたい場合
    leg->Draw();

    // 6. 画像として保存（必要であればコメントアウトを外す）
    // c1->SaveAs("Ratio_Result.pdf");
}