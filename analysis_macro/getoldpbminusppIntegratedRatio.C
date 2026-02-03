#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility>
#include <array>

//===============================================================================================
// PART 1: Pb-Pb Specific Functions (Logic from Code 1, updated with PID Eff)
//===============================================================================================

TF1* PbPb_FitProton(TH1D* hist, double p_min, double p_max) {
    // Code 1 Fit Logic
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
    TF1* finalFit = new TF1("finalFitProtonPbPb", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
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
    hist->Fit(finalFit, "QRN0"); 
    return finalFit;
}

TF1* PbPb_FitPion(TH1D* hist, double p_min, double p_max) {
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
    TF1* finalFit = new TF1("finalFitPionPbPb", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
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
    hist->Fit(finalFit, "QRN0");
    return finalFit;
}

// 変更点: 戻り値を pair<purityMap, pidEffMap> に変更し、PID効率も計算する
std::pair<std::map<double, double>, std::map<double, double>> PbPb_GetCorrectionMaps(
    const TString particleType, const TString fileName,
    double jetPtMin, double jetPtMax,
    double centMin, double centMax,
    double rMin, double rMax) 
{
    std::map<double, double> purityMap;
    std::map<double, double> pidEffMap;

    TString sparseName;
    double signal_min_cut, signal_max_cut;
    if (particleType == "proton") { sparseName = "tpcTofPr"; signal_min_cut = -3.5; signal_max_cut = 0.5; } 
    else { sparseName = "tpcTofPi"; signal_min_cut = -0.5; signal_max_cut = 3.5; }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {purityMap, pidEffMap}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { inputFile->Close(); return {purityMap, pidEffMap}; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    const double p_start = 2.0; const double p_end = 6.5; const double p_step = 0.2;
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        double purity = 0.0;
        double pidEff = 1.0;

        if (h_nSigma->GetEntries() >= 100) {
            TF1* bestFit = (particleType == "proton") ? PbPb_FitProton(h_nSigma, p_min, p_max) : PbPb_FitPion(h_nSigma, p_min, p_max);
            
            // PID Efficiency Calculation (Same logic as pp)
            TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
            signalFunc->SetParameters(&bestFit->GetParameters()[0]); // 0,1,2 params are signal gaus
            
            double intSigInCut = signalFunc->Integral(signal_min_cut, signal_max_cut); // Signal inside cut
            double intSigTotal = signalFunc->Integral(-10, 10);                        // Total Signal
            double intTotInCut = bestFit->Integral(signal_min_cut, signal_max_cut);    // All particles inside cut

            if (intTotInCut > 0) purity = intSigInCut / intTotInCut;
            if (intSigTotal > 0) pidEff = intSigInCut / intSigTotal;

            purity = std::max(0.0, std::min(1.0, purity));
            pidEff = std::max(0.0, std::min(1.0, pidEff));

            delete signalFunc; delete bestFit;
        }
        purityMap[p_center] = purity;
        pidEffMap[p_center] = pidEff;
        delete h_nSigma;
    }
    inputFile->Close();
    return {purityMap, pidEffMap};
}

// 変更点: pidEffMapを受け取り、補正を適用する
TH1D* PbPb_GetSpectrum(TString fileName, TString sparseName, 
                       const std::map<double, double>& purityMap,
                       const std::map<double, double>& pidEffMap,
                       double jetPtMin, double jetPtMax,
                       double centMin, double centMax,
                       double rMin, double rMax) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { inputFile->Close(); return nullptr; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted_pbpb");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double best_p_diff = 100.0; double purity = 0; double pidEff = 1.0;
        
        for(auto const& [p_val, pur_val] : purityMap) {
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
                pidEff = pidEffMap.at(p_val);
            }
        }
        if (purity == 0 || best_p_diff > 0.2) continue; 
        if (pidEff < 0.01) pidEff = 1.0; // Zero division guard

        double correctionFactor = purity / pidEff; // Apply both Purity and PID Efficiency

        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            h2_weighted->SetBinContent(j_pt, i_p, h2_p_vs_pt->GetBinContent(j_pt, i_p) * correctionFactor);
            h2_weighted->SetBinError(j_pt, i_p, h2_p_vs_pt->GetBinError(j_pt, i_p) * correctionFactor);
        }
    }
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    delete h2_p_vs_pt; delete h2_weighted; inputFile->Close();
    return h_pt_weighted;
}

TH1D* PbPb_LoadEff(TString particleType, TString fileName) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    TH1D *hReco = (TH1D*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
    TH1D *hGen = (TH1D*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    if (!hReco || !hGen) { file->Close(); return nullptr; }
    hReco->SetDirectory(nullptr); hGen->SetDirectory(nullptr);
    file->Close();
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEff_PbPb_%s", particleType.Data()));
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    return hEff;
}

//===============================================================================================
// PART 2: pp Specific Functions (Logic from Code 2 - Robust)
//===============================================================================================

TF1* PP_FitProton(TH1D* hist, double p_min, double p_max) {
    std::vector<double> p4_range; for (double p4 = 0.8; p4 <= 1.4; p4 += 0.2) p4_range.push_back(p4); 
    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();
    bool isLowPt = (p_min < 3.9); 

    Int_t oldLevel = gErrorIgnoreLevel; gErrorIgnoreLevel = kError; 
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                TF1* tempFit = new TF1("tempFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
                double ka_init_val = isLowPt ? peak_height * 0.3 : peak_height * 0.4;
                double pr_init_val = isLowPt ? peak_height * 0.7 : peak_height * 0.8;
                double params[12] = { pr_init_val, 0.0, 1.0, ka_init_val, p4_init, 1.0, peak_height * 0.5, p7_init, 1.5, peak_height * 0.05, p10_init, 2.0 };
                tempFit->SetParameters(params);
                tempFit->SetParLimits(0, 0, peak_height * 1.5);
                tempFit->SetParLimits(1, -0.5, 0.5); 
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
    TF1* finalFit = new TF1("finalFitProtonPP", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    if (isLowPt) { finalFit->SetParLimits(0, peak_height * 0.3, peak_height * 0.92); finalFit->SetParLimits(3, peak_height * 0.30, peak_height * 1.0); }
    else { finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 1.5); finalFit->SetParLimits(3, peak_height * 0.1, peak_height * 1.5); }
    finalFit->SetParLimits(1, -0.1, 0.1); finalFit->SetParLimits(2, 0.8, 1.15); 
    hist->Fit(finalFit, "QNRB0"); 
    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* PP_FitPion(TH1D* hist, double p_min, double p_max) {
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
    TF1* finalFit = new TF1("finalFitPionPP", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    hist->Fit(finalFit, "QNRB0");
    return finalFit;
}

std::pair<std::map<double, double>, std::map<double, double>> PP_GetCorrectionMaps(
    const TString particleType, const TString fileName,
    double jetPtMin, double jetPtMax,
    double centMin, double centMax,
    double rMin, double rMax) 
{
    std::map<double, double> purityMap, pidEffMap;
    TString sparseName;
    double sigMin, sigMax;
    if (particleType == "proton") { sparseName = "jetTpcTofPr"; sigMin = -3.5; sigMax = 0.5; } 
    else { sparseName = "jetTpcTofPi"; sigMin = -0.5; sigMax = 3.5; }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return {purityMap, pidEffMap};
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { inputFile->Close(); return {purityMap, pidEffMap}; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);         
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax); 
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    const double p_start = 2.0; const double p_end = 6.5; const double p_step = 0.2;
    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        double purity = 0.0, pidEff = 1.0;
        if (h_nSigma->GetEntries() >= 50) {
            TF1* bestFit = (particleType == "proton") ? PP_FitProton(h_nSigma, p_min, p_max) : PP_FitPion(h_nSigma, p_min, p_max);
            TF1* signalFunc = new TF1("signalFunc", "gaus", -10, 10);
            signalFunc->SetParameters(bestFit->GetParameter(0), bestFit->GetParameter(1), bestFit->GetParameter(2));
            
            double intSigInCut = signalFunc->Integral(sigMin, sigMax);
            double intSigTotal = signalFunc->Integral(-10, 10);
            double intTotInCut = bestFit->Integral(sigMin, sigMax);

            if (intTotInCut > 0) purity = intSigInCut / intTotInCut;
            if (intSigTotal > 0) pidEff = intSigInCut / intSigTotal;
            purity = std::max(0.0, std::min(1.0, purity));
            pidEff = std::max(0.0, std::min(1.0, pidEff));
            delete signalFunc; delete bestFit;
        }
        purityMap[p_center] = purity; pidEffMap[p_center] = pidEff;
        delete h_nSigma;
    }
    inputFile->Close();
    return {purityMap, pidEffMap};
}

TH1D* PP_GetSpectrum(TString fileName, TString particleType, 
                     const std::map<double, double>& purityMap,
                     const std::map<double, double>& pidEffMap,
                     double jetPtMin, double jetPtMax,
                     double centMin, double centMax,
                     double rMin, double rMax) {
    TString sparseName = (particleType == "proton") ? "jetpVsPtForPr" : "jetpVsPtForPi";
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(TString::Format("jet-shape-task/%s", sparseName.Data()));
    if (!hSparse) { inputFile->Close(); return nullptr; }
    
    hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(1, 0); 
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted_pp");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        double best_p_diff = 100.0; double purity = 0; double pidEff = 1.0;
        for(auto const& [p_val, pur_val] : purityMap) {
            if (std::abs(p_val - p_center) < best_p_diff) {
                best_p_diff = std::abs(p_val - p_center);
                purity = pur_val;
                pidEff = pidEffMap.at(p_val);
            }
        }
        if (purity == 0 || best_p_diff > 0.2) continue; 
        if (pidEff < 0.01) pidEff = 1.0;
        double correctionFactor = purity / pidEff;

        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            h2_weighted->SetBinContent(j_pt, i_p, h2_p_vs_pt->GetBinContent(j_pt, i_p) * correctionFactor);
            h2_weighted->SetBinError(j_pt, i_p, h2_p_vs_pt->GetBinError(j_pt, i_p) * correctionFactor);
        }
    }
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    delete h2_p_vs_pt; delete h2_weighted; inputFile->Close();
    return h_pt_weighted;
}

TH1D* PP_LoadEff(TString particleType, TString fileName, double jetPtMin, double jetPtMax, double centMin, double centMax) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) return nullptr;
    THnSparseD *hRecoSparse = (THnSparseD*) file->Get(TString::Format("jet-shape-task/ptHistogram%s", particleType.Data()));
    THnSparseD *hGenSparse = (THnSparseD*) file->Get(TString::Format("jet-shape-task/ptGenerated%s", particleType.Data()));
    if (!hRecoSparse || !hGenSparse) { file->Close(); return nullptr; }
    
    hRecoSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax); hRecoSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
    hGenSparse->GetAxis(1)->SetRangeUser(jetPtMin, jetPtMax); hGenSparse->GetAxis(2)->SetRangeUser(centMin, centMax);
    TH1D *hReco = hRecoSparse->Projection(0); TH1D *hGen = hGenSparse->Projection(0);
    hReco->SetDirectory(nullptr); hGen->SetDirectory(nullptr);
    file->Close();
    TH1D* hEff = (TH1D*)hReco->Clone(TString::Format("hEff_PP_%s", particleType.Data()));
    hEff->Divide(hReco, hGen, 1.0, 1.0, "B");
    delete hReco; delete hGen;
    return hEff;
}

//===============================================================================================
// MAIN: Calculate Ratio PbPb - Ratio pp (Centrality Dependent + PID Eff Corrected)
//===============================================================================================

void calculateRatioExcessCentralityDep() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    const TString kPbPbData = "AnalysisResults-70.root";
    const TString kPbPbMC   = "AnalysisResults-65.root";
    const TString kPPData   = "AnalysisResults-83.root";
    const TString kPPMC     = "AnalysisResults-91.root";

    const double kJetPtMin = 40.0;
    const double kJetPtMax = 60.0;
    const double kTrackPtMin = 2.5; 
    const double kTrackPtMax = 4.5;
    
    std::vector<std::pair<double, double>> centBins = {
        {0.0, 10.0},
        {10.0, 30.0},
        {30.0, 50.0}
    };
    std::vector<int> colors = {kBlack, kRed+1, kBlue+1};
    std::vector<int> markers = {kFullCircle, kFullSquare, kFullTriangleUp};

    const int nRBins = 7;
    double rBins[nRBins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};

    std::cout << "Loading Efficiencies..." << std::endl;
    TH1D* hEffPr_PbPb = PbPb_LoadEff("Proton", kPbPbMC);
    TH1D* hEffPi_PbPb = PbPb_LoadEff("Pion", kPbPbMC);
    TH1D* hEffPr_PP   = PP_LoadEff("Proton", kPPMC, kJetPtMin, kJetPtMax, 0, 100);
    TH1D* hEffPi_PP   = PP_LoadEff("Pion", kPPMC, kJetPtMin, kJetPtMax, 0, 100);

    if (!hEffPr_PbPb || !hEffPi_PbPb || !hEffPr_PP || !hEffPi_PP) {
        std::cerr << "Error loading efficiencies." << std::endl; return;
    }

    // =========================================================
    // STEP 1: Calculate pp Baseline
    // =========================================================
    std::cout << "=== Calculating pp Baseline ===" << std::endl;
    std::vector<double> ppRatioVec(nRBins, 0.0);
    std::vector<double> ppErrorVec(nRBins, 0.0);

    for (int i = 0; i < nRBins; ++i) {
        double rMin = rBins[i], rMax = rBins[i+1];
        auto mapsPr = PP_GetCorrectionMaps("proton", kPPData, kJetPtMin, kJetPtMax, 0, 100, rMin, rMax);
        auto mapsPi = PP_GetCorrectionMaps("pion",   kPPData, kJetPtMin, kJetPtMax, 0, 100, rMin, rMax);
        
        TH1D* hSpecPr = PP_GetSpectrum(kPPData, "proton", mapsPr.first, mapsPr.second, kJetPtMin, kJetPtMax, 0, 100, rMin, rMax);
        TH1D* hSpecPi = PP_GetSpectrum(kPPData, "pion",   mapsPi.first, mapsPi.second, kJetPtMin, kJetPtMax, 0, 100, rMin, rMax);
        
        if (hSpecPr && hSpecPi) {
            auto divideManually = [](TH1D* hS, TH1D* hE) {
                for (int b=1; b<=hS->GetNbinsX(); ++b) {
                    double eff = hE->GetBinContent(hE->FindBin(hS->GetBinCenter(b)));
                    if(eff > 1e-4) { hS->SetBinContent(b, hS->GetBinContent(b)/eff); hS->SetBinError(b, hS->GetBinError(b)/eff); }
                    else { hS->SetBinContent(b, 0); hS->SetBinError(b, 0); }
                }
            };
            divideManually(hSpecPr, hEffPr_PP);
            divideManually(hSpecPi, hEffPi_PP);

            double ePr, ePi;
            double yPr = hSpecPr->IntegralAndError(hSpecPr->FindBin(kTrackPtMin+0.001), hSpecPr->FindBin(kTrackPtMax-0.001), ePr);
            double yPi = hSpecPi->IntegralAndError(hSpecPi->FindBin(kTrackPtMin+0.001), hSpecPi->FindBin(kTrackPtMax-0.001), ePi);
            
            if (yPi > 0) {
                ppRatioVec[i] = yPr / yPi;
                ppErrorVec[i] = ppRatioVec[i] * std::sqrt(std::pow(ePr/yPr, 2) + std::pow(ePi/yPi, 2));
            }
            delete hSpecPr; delete hSpecPi;
        }
    }

    // =========================================================
    // STEP 2: Loop Centrality & Calculate Excess
    // =========================================================
    TCanvas *c1 = new TCanvas("cExcessCent", "Excess vs R Centrality Dep", 700, 600);
    c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.12);

    TLegend *leg = new TLegend(0.6, 0.65, 0.88, 0.85);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);

    std::vector<TH1D*> hExcessList;

    for (size_t ic = 0; ic < centBins.size(); ++ic) {
        double centMin = centBins[ic].first;
        double centMax = centBins[ic].second;
        std::cout << "\n=== Processing Centrality: " << centMin << "-" << centMax << "% ===" << std::endl;
        
        TH1D* hExcess = new TH1D(TString::Format("hExcess_%.0f_%.0f", centMin, centMax), "", nRBins, rBins);
        hExcess->Sumw2();
        hExcess->SetMarkerStyle(markers[ic]);
        hExcess->SetMarkerColor(colors[ic]);
        hExcess->SetLineColor(colors[ic]);
        hExcess->SetMarkerSize(1.5);
        hExcess->SetLineWidth(2);

        for (int i = 0; i < nRBins; ++i) {
            double rMin = rBins[i]; double rMax = rBins[i+1];
            
            // PbPb Calculation (Now using Updated Map with PID Eff)
            double ratioPbPb = 0, errPbPb = 0;
            auto mapsPr = PbPb_GetCorrectionMaps("proton", kPbPbData, kJetPtMin, kJetPtMax, centMin, centMax, rMin, rMax);
            auto mapsPi = PbPb_GetCorrectionMaps("pion",   kPbPbData, kJetPtMin, kJetPtMax, centMin, centMax, rMin, rMax);
            
            // Pass the PID Eff Map to Spectrum calculator
            TH1D* hSpecPr = PbPb_GetSpectrum(kPbPbData, "pVsPtForPr", mapsPr.first, mapsPr.second, kJetPtMin, kJetPtMax, centMin, centMax, rMin, rMax);
            TH1D* hSpecPi = PbPb_GetSpectrum(kPbPbData, "pVsPtForPi", mapsPi.first, mapsPi.second, kJetPtMin, kJetPtMax, centMin, centMax, rMin, rMax);

            if (hSpecPr && hSpecPi) {
                hSpecPr->Divide(hEffPr_PbPb);
                hSpecPi->Divide(hEffPi_PbPb);
                double ePr, ePi;
                double yPr = hSpecPr->IntegralAndError(hSpecPr->FindBin(kTrackPtMin+0.001), hSpecPr->FindBin(kTrackPtMax-0.001), ePr);
                double yPi = hSpecPi->IntegralAndError(hSpecPi->FindBin(kTrackPtMin+0.001), hSpecPi->FindBin(kTrackPtMax-0.001), ePi);
                if (yPi > 0) {
                    ratioPbPb = yPr / yPi;
                    errPbPb = ratioPbPb * std::sqrt(std::pow(ePr/yPr, 2) + std::pow(ePi/yPi, 2));
                }
                delete hSpecPr; delete hSpecPi;
            }

            if (ratioPbPb > 0 && ppRatioVec[i] > 0) {
                double excess = ratioPbPb - ppRatioVec[i];
                double excessErr = std::sqrt(errPbPb*errPbPb + ppErrorVec[i]*ppErrorVec[i]);
                hExcess->SetBinContent(i+1, excess);
                hExcess->SetBinError(i+1, excessErr);
            }
        }
        hExcessList.push_back(hExcess);
        
        if (ic == 0) {
            hExcess->GetXaxis()->SetTitle("Distance from Jet Axis #it{r}");
            hExcess->GetYaxis()->SetTitle("(p/#pi)_{Pb-Pb} - (p/#pi)_{pp}");
            hExcess->GetYaxis()->SetRangeUser(-0.2, 0.8);
            hExcess->Draw("E1");
        } else {
            hExcess->Draw("E1 SAME");
        }
        leg->AddEntry(hExcess, TString::Format("%.0f-%.0f%%", centMin, centMax), "p");
    }

    TLine *line = new TLine(0.0, 0.0, 0.7, 0.0);
    line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
    leg->Draw();

    TPaveText *pt = new TPaveText(0.2, 0.75, 0.5, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextFont(42);
    pt->AddText("#bf{This thesis}");
    pt->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    pt->AddText(TString::Format("%.0f < #it{p}_{T,jet} < %.0f GeV/#it{c}", kJetPtMin, kJetPtMax));
    pt->AddText(TString::Format("%.1f < #it{p}_{T,track} < %.1f GeV/#it{c}", kTrackPtMin, kTrackPtMax));
    pt->Draw();
    
    delete hEffPr_PbPb; delete hEffPi_PbPb; delete hEffPr_PP; delete hEffPi_PP;
}