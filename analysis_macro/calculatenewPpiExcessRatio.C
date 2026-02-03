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
#include <TError.h>  // 追加
#include <TGraph.h>  // 追加
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <utility> // for std::pair
#include <iomanip> // 追加

//=============================================================================
// Configuration
//=============================================================================
const double kComplexFitThreshold = 500.0;
const double kJetRadiusParam = 0.4; 

//=============================================================================
// Fitting Helper Functions (Updated to Code B Logic / Smoothing Ready)
//=============================================================================

TF1* optimizeAndFitProton(TH1D* hist, double p_min, double p_max) {
    // コンソールスパム防止
    int oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning;

    // --- Grid Search Ranges ---
    std::vector<double> p4_range; 
    for (double p4 = 0.5; p4 <= 1.2; p4 += 0.3) p4_range.push_back(p4); 

    std::vector<double> p7_range; for (double p7 = 3.5; p7 <= 5.5; p7 += 0.5) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 8.0; p10 <= 12.0; p10 += 1.0) p10_range.push_back(p10);

    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    // ★Momentum switch: 4.2 GeV以上で挙動を変えるフラグ
    bool isHighPt = (p_min > 4.1);

    // --- Step 1: Grid Search ---
    TF1* tempFit = new TF1("tempFitProtonSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
                double pr_init_height = isHighPt ? peak_height * 0.8 : peak_height;
                double kaon_init_height = isHighPt ? peak_height * 0.30 : peak_height * 0.07;

                // Set Parameters
                tempFit->SetParameter(0, pr_init_height);
                tempFit->SetParLimits(0, peak_height * 0.2, peak_height * 1.1);
                
                tempFit->SetParameter(1, 0.0); 
                tempFit->SetParLimits(1, -0.2, 0.5);
                
                tempFit->SetParameter(2, 1.5); 
                tempFit->SetParLimits(2, 0.7, 2.0);

                // Backgrounds
                tempFit->SetParameter(3, kaon_init_height); 
                tempFit->SetParameter(4, p4_init); 
                tempFit->SetParameter(5, 1.5);

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
            }
        }
    }
    delete tempFit;
    
    // --- Step 2: Final Fit with Best Parameters ---
    TF1* finalFit = new TF1("finalFitProton", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);
    finalFit->SetParameters(&best_params[0]);
    finalFit->SetLineColor(kRed);
    finalFit->SetLineWidth(2);

    // Limits
    if (isHighPt) {
        finalFit->SetParLimits(0, peak_height * 0.1, peak_height * 0.88); // 上限を88%程度に抑える
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

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

TF1* optimizeAndFitPion(TH1D* hist, double p_min, double p_max) {
    int oldLevel = gErrorIgnoreLevel; 
    gErrorIgnoreLevel = kWarning;
    
    std::vector<double> p4_range; for (double p4 = -2.5; p4 <= -1.5; p4 += 0.1) p4_range.push_back(p4);
    std::vector<double> p7_range; for (double p7 = -5.0; p7 <= -3.5; p7 += 0.2) p7_range.push_back(p7);
    std::vector<double> p10_range; for (double p10 = 1.5; p10 <= 2.5; p10 += 0.2) p10_range.push_back(p10);
    
    double best_chi2ndf = std::numeric_limits<double>::max();
    std::vector<double> best_params(12);
    double peak_height = hist->GetMaximum();

    TF1* tempFit = new TF1("tempFitPionSearch", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -10, 10);

    // --- Step 1: Grid Search ---
    for (double p4_init : p4_range) {
        for (double p7_init : p7_range) {
            for (double p10_init : p10_range) {
                
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
            }
        }
    }
    delete tempFit;

    // --- Step 2: Final Fit ---
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

    gErrorIgnoreLevel = oldLevel;
    return finalFit;
}

//=============================================================================
// 1. Purity Calculation (Using Smoothing Logic)
//=============================================================================
std::pair<std::map<double, double>, std::map<double, double>> calculatePurities(
    const TString particleType, const TString fileName,
    double jetPtMin, double jetPtMax,
    double centMin, double centMax,
    double rMin, double rMax,
    bool isOutOfJet) 
{
    std::map<double, double> finalPurityMap;
    std::map<double, double> finalEffMap;

    TString sparseName;
    double sigMin, sigMax;
    if (particleType == "proton") {
        sparseName = isOutOfJet ? "jetTpcTofPr" : "jetTpcTofPr"; 
        sigMin = -3.5; sigMax = 0.5;
    } else {
        sparseName = isOutOfJet ? "jetTpcTofPi" : "jetTpcTofPi"; 
        sigMin = -0.5; sigMax = 3.5;
    }
    
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) { return {finalPurityMap, finalEffMap}; }
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { 
        std::cerr << "Hist " << fullPath << " not found!" << std::endl;
        inputFile->Close(); return {finalPurityMap, finalEffMap}; 
    }
    
    // Axis Settings
    if (isOutOfJet) {
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);
        hSparse->GetAxis(2)->SetRangeUser(0.0, 0.0); 
    } else {
        hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
        hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
        hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);
    }

    double p_start = 0.6; double p_end = 8.0; double p_step = 0.2; 
    TGraph *gPurity = new TGraph();
    TGraph *gEff    = new TGraph();
    int pointIdx = 0;

    std::cout << "--- " << particleType << " (Jet " << jetPtMin << "-" << jetPtMax << ", " << (isOutOfJet ? "Bkg" : "Cone") << ") ---" << std::endl;

    for (double p_min = p_start; p_min < p_end; p_min += p_step) {
        double p_max = p_min + p_step;
        double p_center = (p_min + p_max) / 2.0;
        
        // Low Pt Logic (Fixed)
        if (p_center <= 3.0) {
            TF1 gaussian("ideal", "gaus", -10, 10);
            gaussian.SetParameters(1, 0, 1); 
            double idealInCut = gaussian.Integral(sigMin, sigMax);
            double idealTotal = gaussian.Integral(-10, 10);
            finalPurityMap[p_center] = 1.0;
            finalEffMap[p_center]    = idealInCut / idealTotal;
            continue;
        }

        hSparse->GetAxis(0)->SetRangeUser(p_min, p_max);
        TH1D *h_nSigma = hSparse->Projection(1);
        
        if (h_nSigma->GetEntries() >= 50) { 
            TF1* fit = (particleType == "proton") ? optimizeAndFitProton(h_nSigma, p_min, p_max) : optimizeAndFitPion(h_nSigma, p_min, p_max);
            
            TF1* fSig = new TF1("fSig", "gaus", -10, 10);
            fSig->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
            
            double intSigInCut = fSig->Integral(sigMin, sigMax);
            double intSigTotal = fSig->Integral(-10, 10);
            double intTotInCut = fit->Integral(sigMin, sigMax);
            
            double purity = (intTotInCut > 0) ? intSigInCut / intTotInCut : 0.0;
            double pidEff = (intSigTotal > 0) ? intSigInCut / intSigTotal : 0.0;
            
            if (purity > 1.0) purity = 1.0; if (purity < 0.0) purity = 0.0;
            if (pidEff > 1.0) pidEff = 1.0; if (pidEff < 0.0) pidEff = 0.0;

            gPurity->SetPoint(pointIdx, p_center, purity);
            gEff->SetPoint(pointIdx, p_center, pidEff);
            pointIdx++;
            
            delete fSig;
            delete fit;
        }
        delete h_nSigma;
    }

    // Smoothing (pol3 fit for > 3.0 GeV)
    if (pointIdx > 2) { 
        TF1 *fPurFit = new TF1("fPurFit", "pol3", 3.0, p_end);
        TF1 *fEffFit = new TF1("fEffFit", "pol3", 3.0, p_end);
        
        // std::cout << ">>> Smoothing..." << std::endl;
        gPurity->Fit(fPurFit, "Q R N"); 
        gEff->Fit(fEffFit, "Q R N");

        for (double p_curr = 3.1; p_curr < p_end; p_curr += p_step) {
             double smoothPur = fPurFit->Eval(p_curr);
             double smoothEff = fEffFit->Eval(p_curr);
             
             if (smoothPur > 1.0) smoothPur = 1.0; if (smoothPur < 0.0) smoothPur = 0.0;
             if (smoothEff > 1.0) smoothEff = 1.0; if (smoothEff < 0.0) smoothEff = 0.0;
             
             finalPurityMap[p_curr] = smoothPur;
             finalEffMap[p_curr]    = smoothEff;
        }
        delete fPurFit; delete fEffFit;
    } else {
        std::cout << "Warning: Not enough points for smoothing. Using raw values." << std::endl;
        for(int i=0; i<pointIdx; ++i) {
            double x, y; 
            gPurity->GetPoint(i, x, y); finalPurityMap[x] = y; 
            gEff->GetPoint(i, x, y);    finalEffMap[x] = y;
        }
    }

    delete gPurity; delete gEff;
    inputFile->Close();
    return {finalPurityMap, finalEffMap};
}

//=============================================================================
// 2. Weighted Spectra Creation (Corrects for different axes in OutOfJet)
//=============================================================================
TH1D* createWeightedPtSpectrum(TString fileName, TString sparseBaseName, 
                               const std::map<double, double>& purityMap,
                               const std::map<double, double>& pidEffMap,
                               double jetPtMin, double jetPtMax,
                               double centMin, double centMax,
                               double rMin, double rMax,
                               bool isOutOfJet) {
    TFile *inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) return nullptr;
    
    // Construct sparse name (Always using InJet name as requested)
    TString sparseName = sparseBaseName; 
    TString fullPath = TString::Format("jet-shape-task/%s", sparseName.Data());
    THnSparseD *hSparse = (THnSparseD*)inputFile->Get(fullPath);
    if (!hSparse) { 
        std::cerr << "Sparse " << fullPath << " not found!" << std::endl;
        inputFile->Close(); return nullptr; 
    }
    
    // Set Ranges
    hSparse->GetAxis(3)->SetRangeUser(jetPtMin, jetPtMax);
    hSparse->GetAxis(4)->SetRangeUser(centMin, centMax);
    if (!isOutOfJet) {
        hSparse->GetAxis(2)->SetRangeUser(rMin, rMax);
    }

    TH2D *h2_p_vs_pt = (TH2D*)hSparse->Projection(0, 1);
    h2_p_vs_pt->SetDirectory(nullptr);
    
    TH2D *h2_weighted = (TH2D*)h2_p_vs_pt->Clone("h2_weighted");
    h2_weighted->Reset();
    
    for (int i_p = 1; i_p <= h2_p_vs_pt->GetNbinsY(); ++i_p) {
        double p_center = h2_p_vs_pt->GetYaxis()->GetBinCenter(i_p);
        
        if (p_center > 8.0) continue;

        double purity = 1.0;
        double pidEff = 1.0;
        double min_diff = 1.0;
        
        for(auto const& item : purityMap) {
            double p_val = item.first;
            if (std::abs(p_val - p_center) < min_diff) {
                min_diff = std::abs(p_val - p_center);
                purity = item.second;
                pidEff = pidEffMap.at(p_val);
            }
        }
        
        // Skip if map match is too far (except low pt where we hardcode)
        if (min_diff > 0.3 && p_center > 3.0) continue; 
        if (pidEff < 0.01) pidEff = 1.0;

        double totalCorrectionFactor = purity / pidEff;
        
        for (int j_pt = 1; j_pt <= h2_p_vs_pt->GetNbinsX(); ++j_pt) {
            double content = h2_p_vs_pt->GetBinContent(j_pt, i_p);
            double error = h2_p_vs_pt->GetBinError(j_pt, i_p);
            
            if (content <= 0) continue;

            h2_weighted->SetBinContent(j_pt, i_p, h2_weighted->GetBinContent(j_pt, i_p) + content * totalCorrectionFactor);
            // Error propagation: sqrt( (sigma*factor)^2 )
            double newErr = std::sqrt( pow(h2_weighted->GetBinError(j_pt, i_p), 2) + pow(error * totalCorrectionFactor, 2) );
            h2_weighted->SetBinError(j_pt, i_p, newErr);
        }
    }
    
    TH1D *h_pt_weighted = h2_weighted->ProjectionX();
    h_pt_weighted->SetDirectory(nullptr);
    TString suffix = isOutOfJet ? "OutOfJet" : "Cone";
    h_pt_weighted->SetName(TString::Format("h_pt_weighted_%s_%s", sparseName.Data(), suffix.Data()));
    
    delete h2_p_vs_pt;
    delete h2_weighted;
    inputFile->Close();
    return h_pt_weighted;
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

//=============================================================================
// Main Function
//=============================================================================
void calculatePpiExcessRatio() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    const TString kDataFileName = "AnalysisResults-106.root";
    const TString kMcFileName   = "AnalysisResults-91.root";
    
    // --- Config ---
    const double kJetPtMin = 60.0;
    const double kJetPtMax = 250.0;
    const double kCentMin  = 0.0;
    const double kCentMax  = 10.0;
    
    const double kRMin = 0.0;
    const double kRMax = 0.7; 

    std::cout << "==================================================" << std::endl;
    std::cout << " Calculation: Excess Ratio (Yield Subtraction) " << std::endl;
    std::cout << " Centrality: " << kCentMin << "-" << kCentMax << "%" << std::endl;
    std::cout << " Jet pT: " << kJetPtMin << "-" << kJetPtMax << " GeV/c" << std::endl;
    std::cout << "==================================================" << std::endl;

    // --- 1. Load Efficiency ---
    std::cout << "Loading Efficiency..." << std::endl;
    TH1D* hEffPr = calculateEfficiency(kMcFileName, "Proton", 0.0, 100.0);
    TH1D* hEffPi = calculateEfficiency(kMcFileName, "Pion", 0.0, 100.0);
    if(!hEffPr || !hEffPi) { std::cerr << "Eff fail!" << std::endl; return; }

    // --- 2. Process Signal (InJet) ---
    // Note: getPidCorrectionMaps returns pair<PurityMap, PidEffMap>
    std::cout << "\n--- Processing Signal (InJet) ---" << std::endl;
    auto mapsPr_Cone = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    auto mapsPi_Cone = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);

    TH1D* hYieldPr_Cone = createWeightedPtSpectrum(kDataFileName, "jetpVsPtForPr", mapsPr_Cone.first, mapsPr_Cone.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    TH1D* hYieldPi_Cone = createWeightedPtSpectrum(kDataFileName, "jetpVsPtForPi", mapsPi_Cone.first, mapsPi_Cone.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, kRMin, kRMax, false);
    
    hYieldPr_Cone->Divide(hEffPr);
    hYieldPi_Cone->Divide(hEffPi);

    // --- 3. Process Background (OutOfJet) ---
    std::cout << "\n--- Processing Background (OutOfJet) ---" << std::endl;
    auto mapsPr_Bkg = calculatePurities("proton", kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);
    auto mapsPi_Bkg = calculatePurities("pion",   kDataFileName, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);

    TH1D* hYieldPr_Bkg = createWeightedPtSpectrum(kDataFileName, "jetpVsPtForPr", mapsPr_Bkg.first, mapsPr_Bkg.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);
    TH1D* hYieldPi_Bkg = createWeightedPtSpectrum(kDataFileName, "jetpVsPtForPi", mapsPi_Bkg.first, mapsPi_Bkg.second, kJetPtMin, kJetPtMax, kCentMin, kCentMax, 0, 0, true);

    hYieldPr_Bkg->Divide(hEffPr);
    hYieldPi_Bkg->Divide(hEffPi);

    // --- 4. Scale Background by Area Ratio ---
    double areaSignal = TMath::Pi() * (kRMax*kRMax - kRMin*kRMin);
    double areaBkgHist = 2.0 * TMath::Pi() * kJetRadiusParam * kJetRadiusParam;
    
    double scaleFactor = areaSignal / areaBkgHist;

    std::cout << "\n--- Area Scaling ---" << std::endl;
    std::cout << " Scale Factor: " << scaleFactor << std::endl;

    hYieldPr_Bkg->Scale(scaleFactor);
    hYieldPi_Bkg->Scale(scaleFactor);

    // --- 5. Calculate Excess (Cone - ScaledBkg) ---
    std::cout << "\n--- Calculating Excess Yields ---" << std::endl;
    TH1D* hExcessPr = (TH1D*)hYieldPr_Cone->Clone("hExcessPr");
    hExcessPr->Add(hYieldPr_Bkg, -1.0); // Subtraction

    TH1D* hExcessPi = (TH1D*)hYieldPi_Cone->Clone("hExcessPi");
    hExcessPi->Add(hYieldPi_Bkg, -1.0); // Subtraction

    // --- 6. Calculate Excess Ratio ---
    TH1D* hRatioExcess = (TH1D*)hExcessPr->Clone("hRatioExcess");
    hRatioExcess->Divide(hExcessPi);

    // --- 7. Plotting ---
    TCanvas *c1 = new TCanvas("cExcess", "Excess p/pi Ratio", 800, 600);
    c1->SetLeftMargin(0.15);
    
    hRatioExcess->SetTitle("Excess (Jet-Induced) p/#pi Ratio");
    hRatioExcess->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatioExcess->GetYaxis()->SetTitle("Excess p / Excess #pi");
    hRatioExcess->GetYaxis()->SetRangeUser(-0.2, 1.2); 
    
    hRatioExcess->SetMarkerStyle(20);
    hRatioExcess->SetMarkerColor(kRed+1);
    hRatioExcess->SetLineColor(kRed+1);
    
    hRatioExcess->Draw("E1");

    TLine *line = new TLine(hRatioExcess->GetXaxis()->GetXmin(), 0.0, hRatioExcess->GetXaxis()->GetXmax(), 0.0);
    line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();

    TPaveText *pt = new TPaveText(0.2, 0.75, 0.5, 0.88, "NDC");
    pt->SetFillStyle(0); pt->SetBorderSize(0); pt->SetTextFont(42);
    pt->AddText("Ratio = (N_{Cone} - N_{Bkg}^{Scaled})^{p} / (N_{Cone} - N_{Bkg}^{Scaled})^{#pi}");
    pt->AddText(TString::Format("Cent: %.0f-%.0f%%", kCentMin, kCentMax));
    pt->Draw();

    // Debug Print
    std::cout << "\n--- Value Check ---" << std::endl;
    for(int i=1; i<=hRatioExcess->GetNbinsX(); ++i) {
        if(hRatioExcess->GetBinCenter(i) < 2.0 || hRatioExcess->GetBinCenter(i) > 4.5) continue;
        printf(" pT=%.1f: Excess_p=%.1f | Excess_pi=%.1f | Ratio=%.3f\n", 
               hRatioExcess->GetBinCenter(i), 
               hExcessPr->GetBinContent(i), 
               hExcessPi->GetBinContent(i), 
               hRatioExcess->GetBinContent(i));
    }
}