#include "TMath.h"
#include <vector>
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include "TEventList.h"
#include "TEfficiency.h"
#include "TCut.h"
#include "TString.h"
#include "TLegend.h"
#include <string>
#include "FitDoubleGausModule3.C"

using namespace std;

void SetStyle();
void DrawEfficiencies(TH1F *wrongTagsBp,TH1F *totalTagsBp,TH1F* wrongTagsBs,TH1F* totalTagsBs,TH1F* wrongTagsBd,TH1F* totalTagsBd,TH1F* wrongTagsBdata,TH1F* totalTagsBdata,Double_t min,Double_t max, TString XaxisName);
void DrawLegend(TEfficiency* pEffBp, TEfficiency* pEffBs,TEfficiency* pEffBData,TEfficiency* pEffBd);
void FillTagHistos(TTree* tree, TH1F* WrongTagsH,TH1F* TotalTagsH,int histoflag);
void SetCanvMargin(TCanvas *c1);
/*
Double_t GetTaggingPower(TH1F* wrong,TH1F* total, Int_t allEvents){

	Int_t RightandWrongTags = total->GetEntries();
	Double_t TaggingPower=0;
	

	if(RightandWrongTags!=0){
		Double_t efficiency = Double_t(RightandWrongTags)/Double_t(allEvents);
	
		cout << "efficiency " << total->GetName() << ": " << efficiency << endl;
		Int_t wrongTags = wrong->GetEntries();

		Double_t WtagFrac = Double_t(wrongTags)/Double_t(RightandWrongTags);

		TaggingPower = efficiency*(1.0-2*WtagFrac)*(1-2*WtagFrac);
	}
	return TaggingPower; 
}

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
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetTitleColor(1, "XYZ");
 gStyle->SetTitleFont(42, "XYZ");
// gStyle->SetTitleSize(0.05, "XYZ");
// gStyle->SetTitleXOffset(1.05);
// gStyle->SetTitleYOffset(1.2);
 gStyle->SetTitleSize(0.07, "XYZ");
 gStyle->SetTitleXOffset(0.9);
//for bin fit plots
// gStyle->SetTitleYOffset(0.8);
//for bin wtag plots
 gStyle->SetTitleYOffset(0.9);
 gStyle->SetTitleSize(0.07,"t"); 
  gStyle->SetTextSize(0.07); 
 gStyle->SetLabelFont(42, "XYZ");
 gStyle->SetLabelOffset(0.007, "XYZ");
 gStyle->SetLabelSize(0.05, "XYZ");
 gROOT->ForceStyle();



}
/*
void SelectEvents(TTree *tree, TH1F* WtagsMass, TH1F* RtagsMass, TH1F* WrongIP, TH1F *WrongPt, TH1F* WrongdR, TH1F* RightIP,TH1F* RightPt, TH1F* RightdR){

const int Arraysize = 25;
 Double_t MuPt[Arraysize];
 Double_t MuIP3D[Arraysize];

 Int_t isTrackerTagMuon[Arraysize];
 Int_t isPFTagMuon[Arraysize];
 Int_t isGlobalTagMuon[Arraysize];
 Double_t MuPhi[Arraysize];
 Double_t MuEta[Arraysize];
 Double_t MuCharge[Arraysize]; 
 Double_t BplusMass;
 Double_t DeltaR=-1;

 Int_t MuArraySize;

 Double_t BplusEta;
 Double_t BplusPhi;

 Int_t nentries = tree->GetEntries();

 tree->SetBranchAddress("TagMuRecoPt", MuPt);
 tree->SetBranchAddress("TagMuRecoIP3D", MuIP3D);

 tree->SetBranchAddress("TagMuListSize", &MuArraySize);
 tree->SetBranchAddress("TagMuRecoChg", MuCharge);
 tree->SetBranchAddress("GlobalTagMuon", isGlobalTagMuon);
 tree->SetBranchAddress("TrackerTagMuon", isTrackerTagMuon);
 tree->SetBranchAddress("PFTagMuon", isPFTagMuon);
 
 tree->SetBranchAddress("TagMuRecoEta",MuEta);
 tree->SetBranchAddress("TagMuRecoPhi",MuPhi); 

 tree->SetBranchAddress("BplusEta",&BplusEta);
 tree->SetBranchAddress("BplusPhi",&BplusPhi);   
 tree->SetBranchAddress("BplusM_fit",&BplusMass); 

	Int_t noMuons =0;
	 		 //nentries
	for(Int_t i = 0; i<nentries; i++){
			
		tree->GetEntry(i);
			 
		if(MuArraySize <= 0){
		noMuons++;
		continue;
		}					 
		//muonArraySize->Fill(MuArraySize);	
	
		Int_t maxPtMuonIndex = maximumIndex(MuPt,MuArraySize);
		if(maxPtMuonIndex > 24) continue;
		
	  	if( isPFTagMuon[maxPtMuonIndex] != 1  ) continue;
		
		if(MuIP3D[maxPtMuonIndex] > 0.1 || MuIP3D[maxPtMuonIndex] < 0) continue;
		
		if(MuPt[maxPtMuonIndex] < 1.9) continue; // 1.9/1.7 (G&&PF) /2.1 PF
		
		DeltaR = deltaR(BplusEta,BplusPhi,MuEta[maxPtMuonIndex],MuPhi[maxPtMuonIndex]);
		if(DeltaR < 0.3) continue;

		if(MuCharge[maxPtMuonIndex] !=-1 && MuCharge[maxPtMuonIndex] !=1){
			cout << MuCharge[maxPtMuonIndex] << endl;
			continue;
		}
	
		if(MuCharge[maxPtMuonIndex] == 1){
			WtagsMass->Fill(BplusMass);
			WrongIP->Fill(MuIP3D[maxPtMuonIndex]);
			WrongPt->Fill(MuPt[maxPtMuonIndex]);
			WrongdR->Fill(DeltaR);
		}		

		if(MuCharge[maxPtMuonIndex] == -1){
			RtagsMass->Fill(BplusMass);
			RightIP->Fill(MuIP3D[maxPtMuonIndex]);
			RightPt->Fill(MuPt[maxPtMuonIndex]);
			RightdR->Fill(DeltaR);			
		}
  }

	cout << "Number of events with no tag muons " << noMuons << endl; 

	
}*/
/*
void DrawEfficiencies(TH1F *wrongTagsBp,TH1F *totalTagsBp,Double_t min,Double_t max, TString XaxisName){

TEfficiency* pEffBp = 0;
TEfficiency* pEffBpdata = 0;

TH1F* frame = new TH1F("frame","frame",2,min,max);

   frame->GetXaxis()->SetTitle(XaxisName);
   frame->GetYaxis()->SetTitle("Wrong tag fraction #omega");
   frame->SetMinimum(0.15);
   frame->SetMaximum(1);
   frame->Draw();

	if(TEfficiency::CheckConsistency(*wrongTagsBp,*totalTagsBp)){
			
	 	pEffBp = new TEfficiency(*wrongTagsBp,*totalTagsBp);
		pEffBp->SetLineColor(2);
   		pEffBp->SetLineWidth(2);
   		pEffBp->SetMarkerColor(2);
		pEffBp->Draw("samep");
		
		//DrawLegend(pEffBp,pEffBs,pEffBpdata);
	}
	cout << "JEE" << endl;
}*/


void DrawEfficiencies(TH1F *wrongTagsBp,TH1F *totalTagsBp,TH1F* wrongTagsBs,TH1F* totalTagsBs,TH1F* wrongTagsBd,TH1F* totalTagsBd,TH1F* wrongTagsBdata,TH1F* totalTagsBdata,Double_t min,Double_t max, TString XaxisName){

TEfficiency* pEffBd = 0;
TEfficiency* pEffBp = 0;
TEfficiency* pEffBs = 0;
TEfficiency* pEffBpdata = 0;

TH1F* frame = new TH1F("frame","",int(max-min),min,max);

   frame->GetXaxis()->SetTitle(XaxisName);
   frame->GetYaxis()->SetTitle("Wrong tag fraction #omega");
   frame->SetMinimum(0.15);
   frame->SetMaximum(0.55);
   frame->Draw();

	if(TEfficiency::CheckConsistency(*wrongTagsBp,*totalTagsBp) && TEfficiency::CheckConsistency(*wrongTagsBs,*totalTagsBs) && TEfficiency::CheckConsistency(*wrongTagsBdata,*totalTagsBdata) && TEfficiency::CheckConsistency(*wrongTagsBd,*totalTagsBd)){
	
		pEffBd = new TEfficiency(*wrongTagsBd,*totalTagsBd);		
	 	pEffBp = new TEfficiency(*wrongTagsBp,*totalTagsBp);
		pEffBs = new TEfficiency(*wrongTagsBs,*totalTagsBs);
		pEffBpdata = new TEfficiency(*wrongTagsBdata,*totalTagsBdata);

//	pEffBp->SetStatisticOption(TEfficiency::kFCP);
//	pEffBs->SetStatisticOption(TEfficiency::kFCP);
//	pEffBpdata->SetStatisticOption(TEfficiency::kFCP);

		pEffBd->SetLineColor(3);
		pEffBd->SetLineWidth(2);
		pEffBd->SetMarkerColor(3);
		pEffBd->SetMarkerStyle(33);


		pEffBp->SetLineColor(2);
   		pEffBp->SetLineWidth(2);
   		pEffBp->SetMarkerColor(2);

		pEffBs->SetLineColor(1);
   		pEffBs->SetLineWidth(1);
   		pEffBs->SetMarkerColor(1);
		pEffBs->SetMarkerStyle(21);

		pEffBpdata->SetLineColor(4);
		pEffBpdata->SetLineWidth(1);
   		pEffBpdata->SetMarkerColor(4);

		pEffBp->Draw("samep");
		pEffBpdata->Draw("samep");
//		pEffBd->Draw("samep");
		pEffBs->Draw("samep");
		DrawLegend(pEffBp,pEffBs,pEffBpdata,pEffBd);
	}
}

void DrawLegend(TEfficiency* pEffBp, TEfficiency* pEffBs,TEfficiency* pEffBData,TEfficiency* pEffBd){

//TLegend *leg = new TLegend(0.7173367,0.7108014,0.8969849,0.8954704,NULL,"brNDC");
TLegend *leg = new TLegend(0.6984925,0.6153846,0.8781407,0.8741259,NULL,"brNDC"); 
  leg->SetBorderSize(1);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->SetBorderSize(0);  
   leg->AddEntry(pEffBs,"B_{s} MC","lep");
//   leg->AddEntry(pEffBd,"B^{0} MC","lep");
   leg->AddEntry(pEffBp,"B^{+} MC","ple");
   leg->AddEntry(pEffBData,"B^{+} Data","ple");
   leg->Draw("same");
}


void FillTagHistos(TTree* tree, TH1F* WrongTagsH,TH1F* TotalTagsH,int histoflag){
	TCut cut1 = "DeltaR > 0.3";
	TCut cut2 = "MuPt > 2.2";
	TCut cut3 = "MuIP < 0.1";
	TCut cut4 = "PFMuon == 1";	

	TCut PreCuts = cut1+cut2+cut3+cut4;
	
	TString BincutR;
	Double_t nWrongTags,nRightTags;
	Double_t BinLow;
	Double_t BinHigh;
	Double_t upperlim_sigma2R=0.04; //sigma2R upperlimit for fit
	Double_t lowerlim_sigmaR=0.009;
//sigma2  limit 0.009-0.0.06 works for IP first bin
//upper limit 0.009-0.07 works for second,third and 4th IP bin
		
	vector<Int_t> tags(2,-1);
	
			//WrongTagsH->GetNbinsX()+1
	for(int i = 1; i<WrongTagsH->GetNbinsX()+1; i++){ //Pt bin 4th problem with nWtags  
		cout << endl;
		cout << "BIN " << i << endl;	
		cout << endl;
	
	BinLow = WrongTagsH->GetBinLowEdge(i);
	BinHigh = BinLow + WrongTagsH->GetBinWidth(i);

	if(histoflag == 1){
	PreCuts = "MuPt > 2.9 && DeltaR > 0.3 && PFMuon == 1";
	BincutR = TString(PreCuts + Form("MuIP>%f && MuIP< %f ",BinLow,BinHigh) );
	cout << BincutR << endl;

	if(i==2) {upperlim_sigma2R=0.045;} //for Bs cuts case 
	if(i==3) {upperlim_sigma2R=0.043;} // for Bs cuts case

	
	}

	else if(histoflag == 2){
	PreCuts = "MuIP < 0.3 && DeltaR > 0.3 && PFMuon == 1";
	BincutR = TCut(PreCuts + Form("MuPt>%f && MuPt< %f ",BinLow,BinHigh) );
	//if(i == 1) upperlim_sigma2R=0.043; //for Bs cuts case 
        if(i == 1) upperlim_sigma2R=0.045; //for Bd cuts case 
	}
	
	else { //histoflag == 3
	PreCuts = "MuIP < 0.3 && MuPt > 2.9 && PFMuon == 1";
        BincutR = TCut( PreCuts + Form("DeltaR>%f && DeltaR< %f ",BinLow,BinHigh));

	//if(i == 2) lowerlim_sigmaR = 0.015;	//for Bs cuts case 
        if(i == 2) lowerlim_sigmaR = 0.007;	//for Bd cuts case 
	}

	cout << BincutR << endl;
	cout << "histoflag " << histoflag << endl;

	tags = FitDoubleGausModule3(tree,BincutR,upperlim_sigma2R,lowerlim_sigmaR,histoflag,BinLow,BinHigh,i);
	//nRightTags = tags[0]			
	//nWrongTags = tags[1]

	WrongTagsH->SetBinContent(i,tags[1]);
	TotalTagsH->SetBinContent(i,tags[1]+tags[0]);
	cout << nRightTags << tags[0]<< endl;
	cout << nWrongTags << tags[1] << endl;	
	
 	upperlim_sigma2R=0.04; 
	lowerlim_sigmaR=0.009;
	}
}

void SetCanvMargin(TCanvas *c1){
c1->SetLeftMargin(0.1457286);
c1->SetBottomMargin(0.1541463);
}
void WrongtagfracBData(){
using namespace std;

SetStyle();
//Old samples
 
//TFile *BpFile = TFile::Open("CutOptFileMCBpNewTrig.root");
//TFile *BpFile = TFile::Open("CutOptFiles/SlimmedBpMC.root");
//TFile *BsFile = TFile::Open("CutOptFiles/CutOptFileMCBsNew.root");
//TFile *BsFile = TFile::Open("CutOptFileMCBsNewTrig2.root"); //trigger + matching ok
//TFile *BdataFile = TFile::Open("SlimmedBdataNew.root");

//TFile *BdataFile = TFile::Open("CutOptFiles/SlimmedBdataNew.root");
//TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileMCBdNewTrig.root");

//new samples
TFile *BpFile = TFile::Open("CutOptFiles/CutOptFileBpMCFinal.root");
TFile *BsFile = TFile::Open("CutOptFiles/CutOptFileBsMCFinal2.root");
TFile *BdataFile = TFile::Open("CutOptFiles/CutOptFileBpDataFinal.root");
TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileBdMCFinal.root");

TTree* myBdTree = (TTree*)BdFile->Get("cutTree");
TTree* myBpTree = (TTree*)BpFile->Get("cutTree");
TTree* myBsTree = (TTree*)BsFile->Get("cutTree");
TTree* myBdataTree = (TTree*)BdataFile->Get("cutTree"); 

cout << "Bp "<< myBpTree->GetEntries() << endl;
cout << "Bs "<<myBsTree->GetEntries() << endl;
cout << "B data "<<myBdataTree->GetEntries() << endl;

   //Double_t IPBinArray[] = {0,0.0055,0.015,0.035,0.1}; //for Bs cuts
   Double_t IPBinArray[] = {0,0.0055,0.015,0.035,0.1,0.3}; //for Bd
   int BinNumIP = sizeof(IPBinArray)/sizeof(Double_t)-1;

//   Double_t DeltaRBinArray[] = {0.3,1.1,1.68,2.07,2.35,2.6,2.79,2.96,3.095,3.24,3.53,5.0};
//   Double_t DeltaRBinArray[] = {0.3,1.8,2.4,2.75,2.96,3.15,3.53,5.0};
     Double_t DeltaRBinArray[] = {0.3,2.2,2.81,3.15,5.0};	
   int BinNumDeltaR = sizeof(DeltaRBinArray)/sizeof(Double_t)-1;
    
 //  Double_t PtBinArray[] = {2,4,4.5,5.2,6.05,7.35,9.7,20.0};
   //Double_t PtBinArray[] = {2,4.5,6.05,9.7,20.0};
   Double_t PtBinArray[] = {2.9,3.3,4.45,6.8,20.0};		
   int BinNumPt = sizeof(PtBinArray)/sizeof(Double_t)-1; 
	 
   TH1F *WrongTagsIPBs = new TH1F("WrongTagsIPBs","WrongTagsIPBs",BinNumIP,IPBinArray);
   TH1F *TotalTagsIPBs = new TH1F("TotalTagsIPBs","TotalTagsIPBs",BinNumIP,IPBinArray);

   TH1F *WrongTagsIPBp = new TH1F("WrongTagsIPBp","WrongTagsIPBp",BinNumIP,IPBinArray);
   TH1F *TotalTagsIPBp = new TH1F("TotalTagsIPBp","TotalTagsIPBp",BinNumIP,IPBinArray);	

   TH1F *WrongTagsIPBd = new TH1F("WrongTagsIPBd","WrongTagsIPBd",BinNumIP,IPBinArray);
   TH1F *TotalTagsIPBd = new TH1F("TotalTagsIPBd","TotalTagsIPBd",BinNumIP,IPBinArray);

   TH1F *WrongTagsIPBpdata = new TH1F("WrongTagsIPBpdata","WrongTagsIPBpdata",BinNumIP,IPBinArray);
   TH1F *TotalTagsIPBpdata = new TH1F("TotalTagsIPBpdata","TotalTagsIPBpdata",BinNumIP,IPBinArray);

  TH1F *WrongTagsdRBs = new TH1F("WrongTagsDeltaRBs","WrongTagsDeltaRBs",BinNumDeltaR,DeltaRBinArray);
   TH1F *TotalTagsdRBs = new TH1F("TotalTagsDeltaRBs","TotalTagsDeltaRBs",BinNumDeltaR,DeltaRBinArray);

  TH1F *WrongTagsdRBd = new TH1F("WrongTagsDeltaRBd","WrongTagsDeltaRBd",BinNumDeltaR,DeltaRBinArray);
  TH1F *TotalTagsdRBd = new TH1F("TotalTagsDeltaRBd","TotalTagsDeltaRBd",BinNumDeltaR,DeltaRBinArray);

   TH1F *WrongTagsdRBp = new TH1F("WrongTagsDeltaRBp","WrongTagsDeltaRBp",BinNumDeltaR,DeltaRBinArray);
   TH1F *TotalTagsdRBp = new TH1F("TotalTagsDeltaRBp","TotalTagsDeltaRBp",BinNumDeltaR,DeltaRBinArray); 

   TH1F *WrongTagsdRBpdata = new TH1F("WrongTagsDeltaRBpdata","WrongTagsDeltaRBpdata",BinNumDeltaR,DeltaRBinArray);
   TH1F *TotalTagsdRBpdata = new TH1F("TotalTagsDeltaRBpdata","TotalTagsDeltaRBpdata",BinNumDeltaR,DeltaRBinArray);

   
   TH1F *WrongTagsPtBs = new TH1F("WrongTagsPtBs","WrongTagsPtBs",BinNumPt,PtBinArray);
   TH1F *TotalTagsPtBs = new TH1F("TotalTagsPtBs","TotalTagsPtBs",BinNumPt,PtBinArray);

   TH1F *WrongTagsPtBp = new TH1F("WrongTagsPtBp","WrongTagsPtBp",BinNumPt,PtBinArray);
   TH1F *TotalTagsPtBp = new TH1F("TotalTagsPtBp","TotalTagsPtBp",BinNumPt,PtBinArray);

   TH1F *WrongTagsPtBpdata = new TH1F("WrongTagsPtBpdata","WrongTagsPtBpdata",BinNumPt,PtBinArray);
   TH1F *TotalTagsPtBpdata = new TH1F("TotalTagsPtBpdata","TotalTagsPtBpdata",BinNumPt,PtBinArray);	

   TH1F *WrongTagsPtBd = new TH1F("WrongTagsPtBd","WrongTagsPtBd",BinNumPt,PtBinArray);
   TH1F *TotalTagsPtBd = new TH1F("TotalTagsPtBd","TotalTagsPtBd",BinNumPt,PtBinArray);	   

   TCut cut1 = "DeltaR > 0.3 && DeltaR < 5";
   TCut cut2 = "MuPt > 2.9 && MuPt<20";
   TCut cut3 = "MuIP < 0.3 && MuIP > 0";
   TCut cut4 = "PFMuon == 1";
   TCut WrongTagsBp = cut1+cut2+cut3+cut4 + "MuCharge == 1";
   //TCut RightTagsBp = cut1+cut2+cut3+cut4 + "MuCharge == -1"; 
   TCut WrongTagsBs = cut1+cut2+cut3+cut4 + "( (BsFlavour == 1) && (MuCharge == 1) ) || ( (MuCharge == -1) && (BsFlavour == -1 ) )";

  TCut WrongTagsBd = cut1+cut2+cut3+cut4 + "( (BdFlavour == 1) && (MuCharge == 1) ) || ( (MuCharge == -1) && (BdFlavour == -1 ) )";
 
//    TCut cut5 = "BpMass >5.18 && BpMass < 5.5";
//cout << "data good PF muons "<<myBdataTree->GetEntries(cut1+cut2+cut3+cut4) << endl;
//cout << "data all PF muons "<<myBdataTree->GetEntries(cut4) << endl;

   myBsTree->Draw("MuPt>>TotalTagsPtBs", cut1+cut2+cut3+cut4); 	
   myBpTree->Draw("MuPt>>WrongTagsPtBp", WrongTagsBp);
   myBsTree->Draw("MuPt>>WrongTagsPtBs", WrongTagsBs);   
   myBpTree->Draw("MuPt>>TotalTagsPtBp", cut1+cut2+cut3+cut4);

   myBdTree->Draw("MuPt>>TotalTagsPtBd", cut1+cut2+cut3+cut4);
   myBdTree->Draw("MuPt>>WrongTagsPtBd", WrongTagsBd);

FillTagHistos(myBdataTree,WrongTagsPtBpdata,TotalTagsPtBpdata,2);
TCanvas *Ptcanv = new TCanvas("Ptcanv","Ptcanv");
SetCanvMargin(Ptcanv);  
Ptcanv->cd();  
//WrongTagsPtBd->Divide(TotalTagsPtBd);
//WrongTagsPtBd->Draw();
DrawEfficiencies(WrongTagsPtBp,TotalTagsPtBp,WrongTagsPtBs,TotalTagsPtBs,WrongTagsPtBd,TotalTagsPtBd,WrongTagsPtBpdata,TotalTagsPtBpdata,1,21,"Tag muon p_{T} [GeV/c]");

/*


   myBsTree->Draw("MuIP>>TotalTagsIPBs", cut1+cut2+cut3+cut4); 	
   myBpTree->Draw("MuIP>>WrongTagsIPBp", WrongTagsBp);
   myBsTree->Draw("MuIP>>WrongTagsIPBs", WrongTagsBs);   
   myBpTree->Draw("MuIP>>TotalTagsIPBp", cut1+cut2+cut3+cut4);

   myBdTree->Draw("MuIP>>TotalTagsIPBd", cut1+cut2+cut3+cut4);
   myBdTree->Draw("MuIP>>WrongTagsIPBd", WrongTagsBd);	

FillTagHistos(myBdataTree,WrongTagsIPBpdata,TotalTagsIPBpdata,1);
TCanvas *IPcanv = new TCanvas("IPcanv","IPcanv");
SetCanvMargin(IPcanv);
DrawEfficiencies(WrongTagsIPBp,TotalTagsIPBp,WrongTagsIPBs,TotalTagsIPBs,WrongTagsIPBd,TotalTagsIPBd, WrongTagsIPBpdata,TotalTagsIPBpdata,-0.0050,0.305, "Tag muon 3D IP [cm]");



myBsTree->Draw("DeltaR>>TotalTagsDeltaRBs", cut1+cut2+cut3+cut4); 	
myBpTree->Draw("DeltaR>>TotalTagsDeltaRBp",cut1+cut2+cut3+cut4);
myBsTree->Draw("DeltaR>>WrongTagsDeltaRBs", WrongTagsBs);   
myBpTree->Draw("DeltaR>>WrongTagsDeltaRBp",WrongTagsBp);

myBdTree->Draw("DeltaR>>TotalTagsDeltaRBd", cut1+cut2+cut3+cut4);
myBdTree->Draw("DeltaR>>WrongTagsDeltaRBd", WrongTagsBd);	

FillTagHistos(myBdataTree,WrongTagsdRBpdata,TotalTagsdRBpdata,3);
TCanvas *dRcanv = new TCanvas("dRcanv","dRcanv");
SetCanvMargin(dRcanv);


DrawEfficiencies(WrongTagsdRBp,TotalTagsdRBp,WrongTagsdRBs,TotalTagsdRBs,WrongTagsdRBd,TotalTagsdRBd,WrongTagsdRBpdata,TotalTagsdRBpdata,0.2,5.2,"Tag muon isolation #DeltaR");
*/

/*
TCanvas *c2 = new TCanvas("c2","c2");
c2->Divide(3,1);
c2->cd(1);
//myBdataTree->Draw("MuIP>>IP(100,0,20)",cut1+cut2+cut3+cut4);
TotalTagsIPBs->Draw();
/*c2->cd(2);
TotalTagsPt->Draw();
c2->cd(3); 
TotalTagsdRBp->Draw();
*/
}
