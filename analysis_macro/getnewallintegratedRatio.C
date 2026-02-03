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
#include <TMath.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <Rtypes.h>

//=============================================================================
// Configuration
//=============================================================================
const double kJetRadiusParam = 0.4; 

//=============================================================================
// Fitting Helper Functions (Updated to match the provided source code)
//=============================================================================

// Updated: Added p_min argument for High Pt logic
TF1* optimizeAndFitProton(TH1D* hist, double p_min) {
    // --- Grid Search Ranges (From Source Code) ---
    std::vector<double> p4_range; for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // ★Momentum switch: 4.1 GeV以上で挙動を変えるフラグ (Source Code logic)
    bool isHighPt = (p_min > 4.1);

    // --- Step 1: Grid Search ---
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                
                // Proton (Signal)
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                tempFit->SetParameter(0, pr_init_height); 
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                // Backgrounds
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;
                tempFit->SetParameter(3, kaon_init_height); 
                tempFit->SetParameter(4, p4_init); 
                tempFit->SetParameter(5, 1.5);

                tempFit->SetParameter(6, peak_height * 0.7); tempFit->SetParameter(7, p7_init); tempFit->SetParameter(8, 2.0);
                tempFit->SetParameter(9, peak_height * 0.1); tempFit->SetParameter(10, p10_init); tempFit->SetParameter(11, 2.5);

                hist->Fit(tempFit, "QNRB0"); // Q:Quiet, N:NoDraw, R:Range, B:Bound, 0:NoStore

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
    
    // --- Step 2: Final Fit with Best Parameters ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    // Limits
    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88); // 上限を厳しく
    } else {
        finalFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
    }
    
    finalFit->SetParLimits(1, -0.1, 0.05);
    finalFit->SetParLimits(2, 0.5, 1.8);

    // Kaon Limits
    if (isHighPt) {
        finalFit->SetParLimits(3, peak_height * 0.20, peak_height * 1.2);
    } else {
        finalFit->SetParLimits(3, peak_height * 0.05, peak_height * 1.0);
    }

    finalFit->SetParLimits(4, 0.4, 1.2); 
    finalFit->SetParLimits(5, 0.6, 1.2);

    // Pion Limits
    finalFit->SetParLimits(6, peak_height * 0.01, peak_height * 1.0);
    finalFit->SetParLimits(7, 2.5, 5.0);
    finalFit->SetParLimits(8, 0.3, 1.7);

    // Electron Limits
    finalFit->SetParLimits(9, peak_height * 0.002, peak_height * 0.3);
    finalFit->SetParLimits(10, 7.0, 10.5);
    finalFit->SetParLimits(11, 0.6, 1.5);

    hist->Fit(finalFit, "RBQ0");
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist) {
    // --- Grid Search Ranges (From Source Code) ---
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.2) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.4) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.4) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // --- Step 1: Grid Search ---
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

    // --- Step 2: Final Fit ---
    TF1* finalFit = new TF1("finalFitPion", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed); 
    finalFit->SetLineWidth(2);

    // Limits
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

//=============================================================================
// Calculation Functions
//=============================================================================
std::map<double, double> calculatePurities(const TString particleType, const TString fileName,
                                           double jetPtMin, double jetPtMax,
                                           double centMin, double centMax,
                                           double rMin, double rMax,
                                           bool isOutOfJet) {
    TString sparseName;
    double signal_min_cut, signal_max_cut;
    if (particleType == "proton") {
        sparseName = isOutOfJet ? "tpcTofPrOutOfJet" : "jetTpcTofPr";
        signal_min_cut = -3.5; signal_max_cut = 0.5;
    } else {
        sparseName = isOutOfJet ? "tpcTofPiOutOfJet" : "jetTpcTofPi";
        signal_min_cut = -0.5; signal_max_cut = 3.5;
    }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { std::cout << "Error: " << fullPath << " not found!" << std::endl; inputFile->Close(); return {}; }
    
    std::map<double, double> purities;
    const double p_start = 2.0;
    const double p_end = 6.5; 
    const double p_step = 0.2;
    
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        // --- Low Pt Logic (From Source Code) ---
        // 3.0 GeV以下はTOF分離が良いのでPurity=1.0と仮定
        if (p_center <= 3.0) {
            purities[p_center] = 1.0;
            continue;
        }

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        
        if (isOutOfJet) {
            hSparse->GetAxis(2)->SetRangeUser(jetPtMin, jetPtMax);
            hSparse->GetAxis(3)->SetRangeUser(centMin, centMax);
        } else {
            hSparse->GetAxis(2)->SetRangeUser(rMin, rMax); 
            hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
            hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);
        }



        TH1D *h_nSigma_raw = hSparse->Projection(1); // rawとして定義する
TH1D *h_nSigma = nullptr;                    // 最終的に使う変数を初期化
int nBins = h_nSigma_raw->GetNbinsX();
        
        // ビン数が奇数の場合 (例: 101)、偶数 (例: 100) にトリミングする
        if (nBins % 2 != 0) {
            double xmin = h_nSigma_raw->GetXaxis()->GetXmin();
            double binWidth = h_nSigma_raw->GetBinWidth(1);
            // 新しい最大値 = 最小値 + (今のビン数 - 1) * ビン幅
            double xmax_new = xmin + (nBins - 1) * binWidth;

            h_nSigma = new TH1D(TString::Format("%s_even", h_nSigma_raw->GetName()), 
                                h_nSigma_raw->GetTitle(), 
                                nBins - 1, xmin, xmax_new);
            
            // データコピー (最後の1ビンは除外される)
            for (int b = 1; b < nBins; ++b) { 
                h_nSigma->SetBinContent(b, h_nSigma_raw->GetBinContent(b));
                h_nSigma->SetBinError(b, h_nSigma_raw->GetBinError(b));
            }
            delete h_nSigma_raw; // 生データは削除
        } else {
            h_nSigma = h_nSigma_raw; // 偶数ならそのまま使う
        }

        // ここで安全にRebin(2)を実行
        h_nSigma->Rebin(2);
        h_nSigma->SetDirectory(nullptr);
        
        double purity = 0.0;
        if (h_nSigma->GetEntries() >= 50) { 
            // Pass p_min for Proton High Pt Logic
            TF1* bestFit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min) : optimizeAndFitPion(h_nSigma);
            
            TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
            signalFunc->SetParameters(bestFit->GetParameter(0), bestFit->GetParameter(1), bestFit->GetParameter(2));
            
            double signal_in_cut = signalFunc->Integral(signal_min_cut, signal_max_cut);
            double total_in_cut = bestFit->Integral(signal_min_cut, signal_max_cut);
            purity = (total_in_cut > 0) ? (signal_in_cut / total_in_cut) : 0;
            
            delete signalFunc; delete bestFit;
        }
        
        if (purity > 1.0) purity = 1.0;
        if (purity < 0.0) purity = 0.0;
        
        purities[p_center] = purity;
        delete h_nSigma;
    }
    inputFile->Close();
    return purities;
}

TH1D* createWeightedPtSpectrum(TString fileName, TString sparseName, const std::map<double, double>& purities,
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax,
                               double rMin, double rMax,
                               bool isOutOfJet) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { std::cout << "Error: " << fullPath << " not found!" << std::endl; inputFile->Close(); return nullptr; }
    
    if (isOutOfJet) {
        hSparse->GetAxis(2)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(3)->SetRangeUser(centMin, centMax);
    } else {
        hSparse->GetAxis(2)->SetRangeUser(rMin, rMax); 
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);
    }

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double best_p_diff = 100.0;
        double purity = 1.0; // Default 1.0 if not found (though code below filters)
        
        // Find matching purity
        double found = false;
        for(auto const& item : purities) {
            double p_val = item.first;
            double pur_val = item.second;
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
                found = true;
            }
        }
        
        // Map Matching Constraint
        if (!found || best_p_diff > 0.2) {
             // Exception: if low pt logic applies and no map entry exists, we might want to keep it?
             // But calculatePurities populates map for low pT now.
             continue; 
        }
        
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            h2_weighted->SetBinContent(j_pt, i_p, content * purity);
            h2_weighted->SetBinError(j_pt, i_p, error * purity);
        }
    }
    
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    TString suffix = isOutOfJet ? "OutOfJet" : TString::Format("Cone_%.1f_%.1f", rMin, rMax);
    h_pt_weighted->SetName(TString::Format("h_pt_weighted_%s_%s_C%.0f", sparseName.Data(), suffix.Data(), centMin));
    
    delete h2_p_vs_pt; delete h2_weighted; inputFile->Close();
    return h_pt_weighted;
}

TH1D* calculateTrackingEfficiency(TString particleType, TString fileName, double jetPtMin, double jetPtMax) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    
    TString recoName = TString::Format("jet-shape-task/ptHistogram%s", particleType.Data());
    TString genName  = TString::Format("jet-shape-task/ptGenerated%s", particleType.Data());
    
    THnSparseD *hRecoSparse = (THnSparseD*) file->Get(recoName);
    THnSparseD *hGenSparse = (THnSparseD*) file->Get(genName);
    
    if (!hRecoSparse || !hGenSparse) { 
        std::cout << "Eff Hists not found: " << recoName << std::endl;
        file->Close(); return nullptr; 
    }

    hRecoSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax);
    hGenSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax);
    hRecoSparse->GetAxis(2)->SetRangeUser(0.0, 100.0);
    hGenSparse->GetAxis(2)->SetRangeUser(0.0, 100.0);

    TH1D *hReco = hRecoSparse->Projection(0);
    TH1D *hGen = hGenSparse->Projection(0);
    
    hReco->SetDirectory(nullptr); 
    hGen->SetDirectory(nullptr);
    file->Close();
    
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEfficiency_%s", particleType.Data()));
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    
    delete hReco; delete hGen;
    return hEff;
}

void getIntegratedYield(TH1D* hSpec, double ptMin, double ptMax, double& yield, double& error) {
    if (!hSpec) { yield = 0; error = 0; return; }
    int binStart = hSpec->FindBin(ptMin + 0.001);
    int binEnd   = hSpec->FindBin(ptMax - 0.001);
    yield = hSpec->IntegralAndError(binStart, binEnd, error);
}

// Function to calculate ratio for one centrality
TH1D* getExcessRatioMinusPP(double centMin, double centMax, 
                            double jetPtMin, double jetPtMax,
                            double trackPtMin, double trackPtMax,
                            TH1D* hEffPr, TH1D* hEffPi, TH1D* hRatioPP,
                            const TString& dataFileName) {
    
    const int nRBins = 7;
    double rBins[nRBins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    TH1D* hFinal = new TH1D(TString::Format("hFinal_C%.0f_%.0f", centMin, centMax), "", nRBins, rBins);
    hFinal->SetDirectory(0);

    // 1. Calculate Background Yields
    auto mapPr_Bkg = calculatePurities("proton", dataFileName, jetPtMin, jetPtMax, centMin, centMax, 0, 0, true);
    auto mapPi_Bkg = calculatePurities("pion",   dataFileName, jetPtMin, jetPtMax, centMin, centMax, 0, 0, true);

    TH1D* hSpecPr_Bkg = createWeightedPtSpectrum(dataFileName, "pVsPtForPrOutOfJet", mapPr_Bkg, jetPtMin, jetPtMax, centMin, centMax, 0, 0, true);
    TH1D* hSpecPi_Bkg = createWeightedPtSpectrum(dataFileName, "pVsPtForPiOutOfJet", mapPi_Bkg, jetPtMin, jetPtMax, centMin, centMax, 0, 0, true);
    
    if(!hSpecPr_Bkg || !hSpecPi_Bkg) return hFinal;

    hSpecPr_Bkg->Divide(hEffPr);
    hSpecPi_Bkg->Divide(hEffPi);

    double yieldPr_Bkg, errPr_Bkg, yieldPi_Bkg, errPi_Bkg;
    getIntegratedYield(hSpecPr_Bkg, trackPtMin, trackPtMax, yieldPr_Bkg, errPr_Bkg);
    getIntegratedYield(hSpecPi_Bkg, trackPtMin, trackPtMax, yieldPi_Bkg, errPi_Bkg);
    
    delete hSpecPr_Bkg; delete hSpecPi_Bkg;

    // 2. Loop r bins
    for (int i = 0; i < nRBins; ++i) {
        double rMin = rBins[i];
        double rMax = rBins[i+1];
        
        auto mapPr = calculatePurities("proton", dataFileName, jetPtMin, jetPtMax, centMin, centMax, rMin, rMax, false);
        auto mapPi = calculatePurities("pion",   dataFileName, jetPtMin, jetPtMax, centMin, centMax, rMin, rMax, false);

        TH1D* hSpecPr = createWeightedPtSpectrum(dataFileName, "jetpVsPtForPr", mapPr, jetPtMin, jetPtMax, centMin, centMax, rMin, rMax, false);
        TH1D* hSpecPi = createWeightedPtSpectrum(dataFileName, "jetpVsPtForPi", mapPi, jetPtMin, jetPtMax, centMin, centMax, rMin, rMax, false);

        if (!hSpecPr || !hSpecPi) { delete hSpecPr; delete hSpecPi; continue; }

        hSpecPr->Divide(hEffPr);
        hSpecPi->Divide(hEffPi);

        double yieldPr_Sig, errPr_Sig, yieldPi_Sig, errPi_Sig;
        getIntegratedYield(hSpecPr, trackPtMin, trackPtMax, yieldPr_Sig, errPr_Sig);
        getIntegratedYield(hSpecPi, trackPtMin, trackPtMax, yieldPi_Sig, errPi_Sig);

        double areaSig = TMath::Pi() * (rMax*rMax - rMin*rMin);
        double areaBkg = 2.0 * TMath::Pi() * kJetRadiusParam * kJetRadiusParam;
        double scale = areaSig / areaBkg;
        
        double excessPr = yieldPr_Sig - (yieldPr_Bkg * scale);
        double excessPi = yieldPi_Sig - (yieldPi_Bkg * scale);
        double errExcessPr = std::sqrt(errPr_Sig*errPr_Sig + std::pow(errPr_Bkg * scale, 2));
        double errExcessPi = std::sqrt(errPi_Sig*errPi_Sig + std::pow(errPi_Bkg * scale, 2));

        if (excessPi > 0 && excessPr > 0) {
            double ratioPbPb = excessPr / excessPi;
            double ratioPbPbErr = ratioPbPb * std::sqrt(std::pow(errExcessPr/excessPr, 2) + std::pow(errExcessPi/excessPi, 2));
            
            double ratioPP = hRatioPP->GetBinContent(i+1);
            double errPP   = hRatioPP->GetBinError(i+1);
            
            hFinal->SetBinContent(i+1, ratioPbPb - ratioPP);
            hFinal->SetBinError(i+1, std::sqrt(ratioPbPbErr*ratioPbPbErr + errPP*errPP));
        }
        delete hSpecPr; delete hSpecPi;
    }
    return hFinal;
}

//=============================================================================
// Main Function
//=============================================================================
void calculateExcessRatioMultiCent() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    const TString kDataFileName = "AnalysisResults-108.root";
    const TString kMcFileName   = "AnalysisResults-91.root";
    const TString kPPRefFileName = "pp_reference.root"; 

    const double kJetPtMin = 80.0;
    const double kJetPtMax = 250.0;
    const double kTrackPtMin = 2.0;
    const double kTrackPtMax = 5.0;
    
    // Load PP
    TFile *fPP = TFile::Open(kPPRefFileName, "READ");
    if (!fPP || fPP->IsZombie()) { std::cout << "PP file missing" << std::endl; return; }
    TH1D *hRatioPP = (TH1D*)fPP->Get("hRatio_pp");
    hRatioPP->SetDirectory(0); fPP->Close();

    // Load Efficiency
    TH1D* hEffPr = calculateTrackingEfficiency("Proton", kMcFileName, kJetPtMin, kJetPtMax);
    TH1D* hEffPi = calculateTrackingEfficiency("Pion", kMcFileName, kJetPtMin, kJetPtMax);

    if (!hEffPr || !hEffPi) {
        std::cout << "Failed to calculate efficiency." << std::endl;
        return;
    }

    // Centrality Bins (0-10, 10-30, 30-50)
    std::vector<std::pair<double, double>> centBins = {{0, 10}, {10, 30}, {30, 50}};
    std::vector<int> colors = {kRed+1, kGreen+2, kBlue+1};
    std::vector<int> markers = {20, 21, 22};

    TCanvas *c1 = new TCanvas("cMulti", "Excess Ratio - pp", 800, 600);
    c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.12);
    
    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);

    std::vector<TH1D*> hists;

    for(size_t i=0; i<centBins.size(); ++i) {
        double cMin = centBins[i].first;
        double cMax = centBins[i].second;
        
        std::cout << "\n=== Processing Centrality " << cMin << "-" << cMax << "% ===" << std::endl;
        
        TH1D* h = getExcessRatioMinusPP(cMin, cMax, kJetPtMin, kJetPtMax, kTrackPtMin, kTrackPtMax, hEffPr, hEffPi, hRatioPP, kDataFileName);
        
        h->SetLineColor(colors[i]);
        h->SetMarkerColor(colors[i]);
        h->SetMarkerStyle(markers[i]);
        h->SetMarkerSize(1.5);
        h->SetLineWidth(2);
        
        h->GetXaxis()->SetTitle("Distance from Jet Axis #it{r}");
        h->GetYaxis()->SetTitle("(p/#pi)^{Excess}_{Pb-Pb} - (p/#pi)_{pp}");
        h->GetXaxis()->SetLabelSize(0.045); h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetLabelSize(0.045); h->GetYaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->GetYaxis()->SetRangeUser(-0.5, 1.0);

        if(i==0) h->Draw("E1");
        else h->Draw("E1 SAME");
        
        leg->AddEntry(h, TString::Format("%.0f-%.0f%%", cMin, cMax), "p");
        hists.push_back(h);
    }

    TLine *line = new TLine(0.0, 0.0, 0.7, 0.0);
    line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
    leg->Draw();

    TPaveText *pt = new TPaveText(0.18, 0.65, 0.55, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); 
    pt->SetTextFont(42); pt->SetTextSize(0.035); pt->SetTextAlign(12);
    pt->AddText("#bf{This thesis}");
    pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText("pp #sqrt{#it{s}} = 13.6 TeV");
    pt->AddText("ch-particle jet, Anti-k_{T} R=0.4");
    pt->AddText("|#it{#eta}_{jet}| < 0.4");
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", kJetPtMin, kJetPtMax));
    pt->AddText(TString::Format("%.1f < #it{p}_{T,track} < %.1f GeV/#it{c}", kTrackPtMin, kTrackPtMax));
    pt->AddText("|#it{#eta}_{track}| < 0.9");
    pt->Draw();
}