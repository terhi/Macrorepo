/*
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooRealVar.h"#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "TPaveText.h"
#include "RooAbsReal.h"
*/
void SetStyle()
{

 gStyle->SetCanvasBorderMode(0);
 gStyle->SetCanvasColor(kWhite);
 gStyle->SetCanvasDefH(600); 
 gStyle->SetCanvasDefW(800);
 gStyle->SetCanvasDefX(0);   
 gStyle->SetCanvasDefY(0);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadColor(kWhite);
 gStyle->SetFrameBorderMode(0);
 gStyle->SetFrameBorderSize(0.1);
 gStyle->SetFrameFillColor(0);
 gStyle->SetFrameFillStyle(0);
 gStyle->SetFrameLineColor(1);
 gStyle->SetFrameLineStyle(1);
 gStyle->SetFrameLineWidth(0.1);
 gStyle->SetMarkerStyle(20);
// gStyle->SetMarkerColor(2);
 gStyle->SetMarkerSize(1.3);
 gStyle->SetEndErrorSize(2);
 gStyle->SetOptStat(111);
 gStyle->SetOptTitle(0);
 gStyle->SetTitleColor(1, "XYZ");
 gStyle->SetTitleFont(42, "XYZ");
// gStyle->SetTitleSize(0.05, "XYZ");

gStyle->SetTitleSize(0.07, "XYZ");
 gStyle->SetTitleXOffset(0.9);
 gStyle->SetTitleYOffset(1.2);

// gStyle->SetTitleXOffset(1.05);
// gStyle->SetTitleYOffset(1.3);
 gStyle->SetLabelFont(42, "XYZ");
 gStyle->SetLabelOffset(0.007, "XYZ");
 gStyle->SetLabelSize(0.05, "XYZ");
 gROOT->ForceStyle();
}

void FitDoubleGaus(){
using namespace RooFit;
using namespace std;
SetStyle();
//TFile *McBplusFile = TFile::Open("BpMass.root");
//TH1F* BpMassHisto = (TH1F*)McBplusFile->Get("canv5");

/*
TRandom *rand = new TRandom();

TH1F* BpMassHisto = new TH1F("Bpm","Bpm",100,4.8,6);
	for(int i=0; i < 10000; i++){
		BpMassHisto->Fill(rand->Gaus(5.3,0.15));
		BpMassHisto->Fill(rand->Uniform(4.8,6)); 
	}
BpMassHisto->Draw();
*/

//TFile *McBplusFile = TFile::Open("/tmp/terhi/ReRecoTotal.root");
//TFile *McBplusFile = TFile::Open("/tmp/terhi/BplusMCFiles/BplusMCCleanedNew.root");
//TTree* myBsTree = (TTree*)McBplusFile->Get("BsTree");
//TFile *McBplusFile = TFile::Open("SlimmedBpMC.root");
//TFile *McBplusFile = TFile::Open("SlimmedBpMC.root");
//TFile *McBplusFile = TFile::Open("CutOptFileMCBpNewTrig.root");
//TFile *McBplusFile = TFile::Open("CutOptFileMCBsNewTrig2.root");
//TFile *McBplusFile = TFile::Open("CutOptFiles/CutOptFileMCBdNewTrig.root");

//TFile *McBplusFile= TFile::Open("CutOptFiles/CutOptFileBdDataFinal.root");
//TFile *McBplusFile= TFile::Open("CutOptFiles/CutOptFileBdMCFinal.root");
//TFile *McBplusFile= TFile::Open("CutOptFiles/CutOptFileBsMCFinal2.root");
//TFile *McBplusFile= TFile::Open("CutOptFiles/CutOptFileBpMCFinal.root");
TFile *McBplusFile= TFile::Open("CutOptFiles/CutOptFileBpDataFinal.root");
TTree* myBsTree = (TTree*)McBplusFile->Get("cutTree");


cout << myBsTree->GetEntries() << endl;	
/*
TCut cut1 = "DeltaR > 0.3";
	TCut cut2 = "MuPt > 1.9";
	TCut cut3 = "MuIP < 0.1";
	TCut cut4 = "PFMuon == 1";
	TCut Precut = cut1+cut2+cut3+cut4;	

 RooRealVar MuCharge("MuCharge","MuCharge",-1,1);
 RooRealVar MuPt("MuPt","MuPt",0,20);
 RooRealVar MuIP("MuIP","MuIP",0,0.1);
 RooRealVar DeltaR("DeltaR","DeltaR",0,5); 
 RooRealVar PFMuon("PFMuon","PFMuon",0,1); 

*/
// Build Gaussian PDF //
	RooRealVar BpMass("BpMass","B_{s} mass [GeV/c^{2}]",5.,5.5);
	RooRealVar BdCt2D("BdCt2D","Bd ct",0,0.3);
//	RooArgSet ntuples(BpMass);   
//	RooDataSet dataset("mass ntuple","mass ntuple",myBsTree,ntuples);
	
//	RooDataSet dataset("MassDataSet","mass set", RooArgSet(BpMass,MuIP,MuPt,DeltaR,PFMuon),Import(*myBsTree));
//	TString BdCut = "BdCt2D>0.02";
RooDataSet dataset("MassDataSet","mass set", RooArgSet(BpMass),Import(*myBsTree));

// Build Gaussian PDF //
//	RooRealVar BpMass("BpMass","B^{+} mass [GeV/c^{2}]",5.,6);
//	RooArgSet ntuples(BpMass);   
//	RooDataSet dataset("mass ntuple","mass ntuple",myBsTree,ntuples);
	RooRealVar mean("mean","mean of gaussian",5.25,5.2,5.4);
	//RooRealVar mean("mean","mean of gaussian",5.25,5.2,5.35);
	RooRealVar sigma("sigma","width of gaussian",0.01,0.,0.05);
	RooGaussian gauss("gauss","gaussian PDF",BpMass,mean,sigma);
 	RooRealVar nsig("nsig","signal events",0,myBsTree->GetEntries());
	RooRealVar fGG("fGG","fraction of signal events",0.5,0.,1.);
 
 	
//	RooExtendPdf egauss("gauss","gaussian PDF",gauss,nsig);

	RooRealVar sigma2("sigma2","width of 2nd gaussian",0.04,0.,0.09); 
	RooGaussian gauss2("gauss2","gaussian PDF",BpMass,mean,sigma2);	
//	RooRealVar nsig2("nsig2","signal 2 events",11500,0,23544); 

//	RooExtendPdf egauss2("egauss2","extended gaussian",gauss2,nsig2));

	RooAddPdf doubleGauss("two gaussians","gauss1+gauss2", RooArgList(gauss,gauss2),fGG);

	RooRealVar coeff1("coeff1","coeff1",-1.,-10.,10);
	RooExponential expo("expo","expo",BpMass,coeff1);
	RooRealVar nbkg("nbkg","background events",0,myBsTree->GetEntries());

	RooRealVar fEER("fEER","fEER",0.,1);
	RooRealVar div("div","div",0.035,0.03,0.06);
	RooGenericPdf erfPdf("erf","(TMath::Erf((-BpMass+5.15)/div)+1)",RooArgSet(div,BpMass)); 
	RooAddPdf bkgtot("BkgtotR","erf+expo R", RooArgList(erfPdf,expo),fEER);


//	RooExtendPdf ebkg("ebkg","extended expo",expo,nbkg);
	// RooPoly = 1 + coeff1*x +coeff2*x^2 + ..	

	RooAddPdf gaussians("gaussians"," double gaus+poly",RooArgList(bkgtot,doubleGauss),RooArgSet(nbkg,nsig)); 


	//RooDataHist dataset("data","data",BpMass,BpMassHisto);
	
// Fit the model to the data
gaussians.fitTo(dataset);

// Plot PDF and toy data overlaid
TCanvas *c2 = new TCanvas("c2","c2");
c2->cd();
//c2->SetLeftMargin(0.1457286);
//c2->SetBottomMargin(0.1741463);
   c2->SetLeftMargin(0.1796482);
   c2->SetBottomMargin(0.1748252);
RooPlot* xframe = BpMass.frame();
dataset.plotOn(xframe);
xframe->GetYaxis()->SetTitle("Events / (5 MeV/c^{2})");
//gaussians.paramOn(xframe);
//gaussians_paramBox.SetFillStyle(0);
//dataset.statOn(xframe);
/*
   TPaveText *pt = new TPaveText(0.1771357,0.6045296,0.4560302,0.8867596,"br");
   pt->SetLineColor(0);	
   pt->SetFillStyle(0);	
   pt->AddText("B^{+} MC");
   pt->AddText("#sqrt{s} = 8 TeV");
   pt->AddText("#int L dt = ~10 fb^{-1}");
   pt->AddText("");   
   
*/
 


 //  TLegend *leg = new TLegend(0.1771357,0.6045296,0.4560302,0.8867596,NULL,"brNDC");

//TLegend *leg = new TLegend(0.1721106,0.6045296,0.4246231,0.8867596,NULL,"brNDC");
//TLegend *leg = new TLegend(0.1771357,0.6045296,0.4560302,0.7867596,NULL,"brNDC");
//TLegend *leg = new TLegend(0.1746231,0.7692308,0.410804,0.8863636,NULL,"brNDC");
//TLegend *leg = new TLegend(0.1771357,0.7902098,0.4133166,0.8828671,NULL,"brNDC");
//TLegend *leg = new TLegend(0.1871859,0.729021,0.3994975,0.8653846,NULL,"brNDC");

TLegend *leg = new TLegend(0.2311558,0.7307692,0.4434673,0.8671329,NULL,"brNDC");
leg->SetLineColor(0);
leg->SetFillStyle(0);
leg->AddEntry((TObject*)0, "MC","");
//leg->AddEntry((TObject*)0, "","");
//leg->AddEntry((TObject*)0, "#sqrt{s} = 8 TeV","");
//leg->AddEntry((TObject*)0, "#int L dt = ~10 fb^{-1}",""); 

gaussians.plotOn(xframe,Components(expo),LineStyle(kDashed),LineColor(3)); //green
//gaussians.plotOn(xframe,Components(doubleGauss),LineStyle(kDashed),LineColor(6)); //magenta
gaussians.plotOn(xframe,Components(gauss),LineStyle(kDashed),LineColor(6)); //skyblue
gaussians.plotOn(xframe,Components(gauss2),LineStyle(kDashed),LineColor(7)); //skyblue
//gaussians.plotOn(xframe);
//gaussians.plotOn(xframe,Components("erf"),LineStyle(kDashed),LineColor(2)); 
gaussians.plotOn(xframe); 
xframe->Draw();
//leg->Draw("same");

// pt->Draw("same");
// Print final value of parameters

mean.Print();
sigma.Print();
sigma2.Print();

nsig.Print();
//nsig2.Print();
nbkg.Print();
coeff1.Print();

double sigmaTot = (sigma.getVal()*sigma.getError() + sigma2.getVal()*sigma2.getError() )/(sigma.getError() + sigma2.getError());

cout << "sigma tot " <<sigmaTot << endl;
cout << "limits low: " << mean.getVal()-2*sigmaTot << " high: " << mean.getVal()+3*sigmaTot << endl;

BpMass.setRange("Signal",mean.getVal()-2*sigmaTot,mean.getVal()+2*sigmaTot);
RooAbsReal* integralSig = doubleGauss.createIntegral(BpMass,NormSet(BpMass),Range("Signal")) ;
cout  << "sig " << integralSig->getVal() << endl;

RooAbsReal* integralBkg = bkgtot.createIntegral(BpMass,NormSet(BpMass),Range("Signal")) ;

cout  << "bkg  asfds dsf gssfg "<< integralBkg->getVal() << endl;

} 
