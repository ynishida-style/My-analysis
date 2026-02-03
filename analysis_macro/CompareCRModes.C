//Filiname:CompareCRModes.C

void CompareModes() {
    gStyle->SetOptStat(0);
    
    // 設定: ファイル名とラベル、色
    const int nModes = 3;
    const char* files[]  = {"Results_CRoff.root", "Results_Mode0.root", "Results_Mode2.root"};
    const char* labels[] = {"CR Off", "CR Mode 0 (Default)", "CR Mode 2 (Gluon-move)"};
    int colors[] = {kBlue, kBlack, kRed};
    int markers[] = {25, 20, 21}; // Open square, Full circle, Full square

    TCanvas *c1 = new TCanvas("c1", "CR Mode Comparison", 800, 600);
    TLegend *leg = new TLegend(0.5, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);

    TH1D* hFrame = new TH1D("hFrame", "", 10, 0, 0.7);
    hFrame->GetYaxis()->SetRangeUser(0.0, 0.6); // 縦軸範囲
    hFrame->GetXaxis()->SetTitle("Distance from Jet Axis r");
    hFrame->GetYaxis()->SetTitle("Proton / Pion Ratio");
    hFrame->Draw();

    for (int i = 0; i < nModes; ++i) {
        TFile *f = TFile::Open(files[i]);
        if(!f) { std::cout << "File not found: " << files[i] << std::endl; continue; }

        TH1D *hp  = (TH1D*)f->Get("h_num_p");
        TH1D *hpi = (TH1D*)f->Get("h_num_pi");
        
        // Ratio計算
        TH1D *hRatio = (TH1D*)hp->Clone(Form("hRatio_%d", i));
        hRatio->Divide(hpi);
        hRatio->SetDirectory(0); // ファイルを閉じてもメモリに残す

        hRatio->SetLineColor(colors[i]);
        hRatio->SetMarkerColor(colors[i]);
        hRatio->SetMarkerStyle(markers[i]);
        hRatio->SetMarkerSize(1.2);
        
        hRatio->Draw("SAME EP");
        leg->AddEntry(hRatio, labels[i], "lp");
        
        f->Close();
    }
    
    leg->Draw();
    c1->SaveAs("Comparison_CR_Modes.png");
}