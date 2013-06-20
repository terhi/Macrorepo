#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>

 
void DrawLeg(TH1F* BpMCHisto, TH1F* BdataHisto,TH1F* BsHisto=0, TH1F* BdHisto =0, TH1F* BdDataHisto =0){
//TLegend *leg = new TLegend(0.7173367,0.7108014,0.8969849,0.8954704,NULL,"brNDC");
TLegend *leg = new TLegend(0.7223618,0.6341463,0.8932161,0.8885017,NULL,"brNDC");
leg->SetFillStyle(0);
leg->SetBorderSize(0);
//leg->SetBorderMode(0);
 if(BsHisto!=0) {
  leg->AddEntry(BsHisto,"B_{s} MC");
 }
if(BdHisto!=0) leg->AddEntry(BdHisto,"B^{0} MC");

leg->AddEntry(BpMCHisto,"B^{+} MC");
leg->AddEntry(BdataHisto,"B^{+} Data");
if(BdDataHisto!=0) leg->AddEntry(BdDataHisto,"B^{0} Data");
leg->Draw("same");
}


void DrawAxis(TH1F*histo, TString xAxis, TString yAxis){
histo->GetXaxis()->SetTitle(xAxis);
histo->GetYaxis()->SetTitle(yAxis);

}


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
 //gStyle->SetMarkerColor(2);
 gStyle->SetMarkerSize(1.3);
 gStyle->SetEndErrorSize(2);
 //gStyle->SetLineColor(2);
 //gStyle->SetLineWidth(2);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
 gStyle->SetTitleColor(1, "XYZ");
 gStyle->SetTitleFont(42, "XYZ");

// gStyle->SetTitleSize(0.05, "XYZ");
// gStyle->SetTitleXOffset(1.);
// gStyle->SetTitleYOffset(1.35);
 gStyle->SetTitleSize(0.07, "XYZ");
 gStyle->SetTitleXOffset(0.9);
 gStyle->SetTitleYOffset(0.9);
 gStyle->SetTitleSize(0.07,"t");  

gStyle->SetLabelFont(42, "XYZ");
 gStyle->SetLabelOffset(0.007, "XYZ");
 gStyle->SetLabelSize(0.05, "XYZ");
 gROOT->ForceStyle();
}

void Errors(TH1F* histo1, TH1F* histo2, TH1F* histo3=0,TH1F* histo4=0, TH1F* histo5 =0){
histo1->Sumw2();
histo2->Sumw2();
   histo3->Sumw2();

 if(histo4!=0){

   histo4->Sumw2();
  
 }
 if(histo5!=0){histo5->Sumw2();}

}

void MarkerStyles(TH1F* BpMC,TH1F* Bdata=0, TH1F* BsMC=0,TH1F* BdMC=0,TH1F* BdData=0){
 if(Bdata!=0){
Bdata->SetMarkerColor(4);
Bdata->SetLineColor(4);
 }
BpMC->SetMarkerColor(2);
BpMC->SetLineColor(2);
 if(BsMC!=0){

  BsMC->SetMarkerColor(1);
  BsMC->SetLineColor(1);
  BsMC->SetMarkerStyle(21);
 }
if(BdMC!=0){

  BdMC->SetMarkerColor(3);
  BdMC->SetLineColor(3);
  BdMC->SetMarkerStyle(33);
 }

if(BdData!=0){

  BdData->SetMarkerColor(6);
  BdData->SetLineColor(6);
  BdData->SetMarkerStyle(33);
 }


}

void Scale2(TH1F* histoSmaller){

 if(histoSmaller->GetEntries() != 0 ){

histoSmaller->Scale( 1.0/double(histoSmaller->GetEntries()) );
 }
}

void Scale(TH1F* histoSmaller, TH1F* histoBigger){

 if(histoSmaller->GetEntries() != 0 ){
histoSmaller->Scale( double(histoBigger->GetEntries())/double(histoSmaller->GetEntries()) );
 }
}

void SetCanvMargin(TCanvas *c1){
c1->SetLeftMargin(0.1457286);
c1->SetBottomMargin(0.1541463);
}

void BsBpComparison(){

using namespace std;




SetStyle();
/*
TFile *McBplusFile = TFile::Open("/tmp/terhi/BplusMCFiles/BplusMCCleanedNew.root");
TTree* myBpTree = (TTree*)McBplusFile->Get("BsTree");

TFile *McBsFile = TFile::Open("/tmp/terhi/BsMCFiles/BsMCCleanedNew.root");
TTree* myBsTree = (TTree*)McBsFile->Get("BsTree");

TFile *BdataFile = TFile::Open("/tmp/terhi/ReRecoTotal.root");
TTree* myBdataTree = (TTree*)BdataFile->Get("BsTree");
*/

/*
TCanvas *MassCanv = new TCanvas("MassCanv","MassCanv");
MassCanv->cd();
TH1F* BpMass = new TH1F("BpMass","BpMass",100,5.,5.5);
TH1F* BpMassMC = new TH1F("BpMassMC","BpMassMC",100,5.1,5.5);
myBpTree->Draw("BplusM_fit>>BpMassMC");
myBdataTree->Draw("BplusM_fit>>BpMass","","same");
cout << BpMass->GetEntries() << endl;
Errors(BpMass,BpMassMC);
MarkerStyles(BpMassMC,BpMass);
Scale(BpMassMC,BpMass);
SetCanvMargin(MassCanv);
DrawLeg(BpMassMC,BpMass);
DrawAxis(BpMassMC,"B^{+} mass [GeV/c^{2}]","Events");
MassCanv->Update();
 */
/*
TCanvas* Coscanv = new TCanvas("Coscanv","Coscanv"); 
Coscanv->cd();
TH1F* CosBp = new TH1F("CosBp","CosBp",100,0.85,1.1);
myBpTree->Draw("BplusCosTheta>>CosBp");
CosBp->Scale( 1.0/double(CosBp->GetEntries()) );

TH1F* CosBs = new TH1F("CosBs","CosBs",100,0.85,1.1);
CosBs->SetLineColor(2);
myBsTree->Draw("CosDeltaAlpha>>CosBs","","same");
CosBs->Scale( 1.0/double(CosBs->GetEntries()) );

DrawLeg(CosBp,CosBs);
DrawAxis(CosBp,"Cos(#alpha) w.r.t signal B meson p_{xyz} and L_{xyz} ","Events");
//IPcanv->Modified();
Coscanv->Update();
*/

/*
TCanvas* Kaoncanv = new TCanvas("Kaoncanv","Kaoncanv"); 
Kaoncanv->cd();
TH1F* BpDataK1pt = new TH1F("BpDataK1pt","BpDataK1pt",100,0,20);
myBdataTree->Draw("KplusPt>>BpDataK1pt");

TH1F* BsK1pt = new TH1F("BsK1pt","BsK1pt",100,0,20);
myBsTree->Draw("K1Pt_fit>>BsK1pt","","same");

TH1F* BsK2pt = new TH1F("BsK2pt","BsK2pt",100,0,20);
myBsTree->Draw("K2Pt_fit>>BsK2pt","","same");
BsK2pt->SetMarkerColor(1);
BsK2pt->SetLineColor(1);

TH1F* BpKpt = new TH1F("BpKpt","BpKpt",100,0,20);
myBpTree->Draw("KplusPt>>BpKpt","","same");

//Errors(BpDataK1pt,BsK2pt);
//Errors(BpKpt,BsK1pt);
Scale(BpDataK1pt,BsK2pt);
Scale(BpKpt,BsK2pt);
Scale(BsK1pt,BsK2pt);
MarkerStyles(BpKpt,BpDataK1pt,BsK1pt);
DrawLeg(BpKpt,BpDataK1pt,BsK1pt);
SetCanvMargin(KaonCanv);
DrawAxis(BsK2pt,"Kaon p_{T} ","Events");
Kaoncanv->Update();
*/

/*
TCanvas* JpsiMucanv = new TCanvas("JpsiMucanv","JpsiMucanv"); 
JpsiMucanv->cd();

TH1F* BsJpsiMu2pt = new TH1F("BsJpsiMu2pt","BsJpsiMu2pt",100,0,20);
myBsTree->Draw("JpsiMu2Pt_alone>>BsJpsiMu2pt","");
BsJpsiMu2pt->SetMaximum(0.08);

TH1F* BpJpsiMu2pt = new TH1F("BpJpsiMu2pt","BpJpsiMu2pt",100,0,20);
myBpTree->Draw("BplusMu2Pt>>BpJpsiMu2pt","","same");


TH1F* BsJpsiMu1pt = new TH1F("BsJpsiMu1pt","BsJpsiMu1pt",100,0,20);
BsJpsiMu1pt->SetLineColor(2);
myBsTree->Draw("JpsiMu1Pt_alone>>BsJpsiMu1pt","","same");

TH1F* BpJpsiMu1pt = new TH1F("BpJpsiMu1pt","BpJpsiMu1pt",100,0,20);
myBpTree->Draw("BplusMu1Pt>>BpJpsiMu1pt","","same"); //sini

TH1F* BDataJpsiMu1pt = new TH1F("BDataJpsiMu1pt","BDataJpsiMu1pt",100,0,20);
myBdataTree->Draw("BplusMu1Pt>>BDataJpsiMu1pt","","same"); 

TH1F* BDataJpsiMu2pt = new TH1F("BDataJpsiMu2pt","BDataJpsiMu2pt",100,0,20);
myBdataTree->Draw("BplusMu2Pt>>BDataJpsiMu2pt","","same"); 

Scale2(BpJpsiMu1pt);
Scale2(BsJpsiMu2pt);
Scale2(BpJpsiMu2pt);
Scale2(BsJpsiMu1pt);
Scale2(BDataJpsiMu2pt);
Scale2(BDataJpsiMu1pt);

MarkerStyles(BpJpsiMu1pt,BDataJpsiMu1pt,BsJpsiMu1pt);
MarkerStyles(BpJpsiMu2pt,BDataJpsiMu2pt,BsJpsiMu2pt);
DrawLeg(BpJpsiMu1pt,BDataJpsiMu1pt,BsJpsiMu1pt);
SetCanvMargin(JpsiMucanv);
DrawAxis(BsJpsiMu2pt,"Jpsi Mu p_{T} [GeV/c]","Events");
JpsiMucanv->Update();
*/

/*
TCanvas* JpsiMucanv2 = new TCanvas("JpsiMucanv2","JpsiMucanv2"); 
JpsiMucanv2->cd();

TH1F* BpJpsiMu1p = new TH1F("BpJpsiMu1p","BpJpsiMu1p",100,0,60);
myBpTree->Draw("BplusMu1Ptot>>BpJpsiMu1p");
BpJpsiMu1p->Scale( 1.0/double(BpJpsiMu1p->GetEntries()) );
BpJpsiMu1p->SetMinimum(0.);
BpJpsiMu1p->SetMaximum(0.1);

TH1F* BpJpsiMu2p = new TH1F("BpJpsiMu2p","BpJpsiMu2p",100,0,60);
myBpTree->Draw("BplusMu2Ptot>>BpJpsiMu2p","","same");
BpJpsiMu2p->Scale( 1.0/double(BpJpsiMu2p->GetEntries()) );
BpJpsiMu2p->SetLineColor(7);

TH1F* BsJpsiMu2p = new TH1F("BsJpsiMu2p","BsJpsiMu2p",100,0,60);
BsJpsiMu2p->SetLineColor(6);
myBsTree->Draw("TMath::Sqrt(JpsiMu2Pt_alone*JpsiMu2Pt_alone + Mu2Pz_beffit*Mu2Pz_beffit) >>BsJpsiMu2p","","same");
BsJpsiMu2p->Scale( 1.0/double(BsJpsiMu2p->GetEntries()) );

TH1F* BsJpsiMu1p = new TH1F("BsJpsiMu1p","BsJpsiMu1p",100,0,60);
BsJpsiMu1p->SetLineColor(2);
myBsTree->Draw("TMath::Sqrt(JpsiMu1Pt_alone*JpsiMu1Pt_alone + Mu1Pz_beffit*Mu1Pz_beffit) >>BsJpsiMu1p","","same");
BsJpsiMu1p->Scale( 1.0/double(BsJpsiMu1p->GetEntries()) );

//DrawLeg(BpJpsiMu1p,BsJpsiMu1p);
//DrawLeg(BpJpsiMu2p,BsJpsiMu2p);
//DrawAxis(BpJpsiMu1p,"Jpsi Mu P ","Events");

JpsiMucanv2->Update(); 
*/

/*
TCanvas *Jpsicanv = new TCanvas("Jpsicanv","Jpsicanv");
Jpsicanv->cd();

TH1F* BdataJpsiPt = new TH1F("BdataJpsiPt","BdataJpsiPt",100,0.,70);
myBdataTree->Draw("JpsiPt_bplus>>BdataJpsiPt","");

TH1F* BsJpsiPt = new TH1F("BsJpsiPt","BsJpsiPt",100,0.,70);
myBsTree->Draw("JpsiPt_alone>>BsJpsiPt","","same");

TH1F* BpJpsiPt = new TH1F("BpJpsiPt","BpJpsiPt",100,0.,70);
myBpTree->Draw("JpsiPt_bplus>>BpJpsiPt","","same");

MarkerStyles(BpJpsiPt,BdataJpsiPt,BsJpsiPt);
Scale2(BdataJpsiPt);
Scale2(BsJpsiPt);
Scale2(BpJpsiPt);
SetCanvMargin(Jpsicanv);
DrawLeg(BpJpsiPt,BdataJpsiPt,BsJpsiPt);
DrawAxis(BdataJpsiPt,"Jpsi p_{T} [GeV/c]","Events");
Jpsicanv->Update();
*/

/*
TCanvas *KandJpsiVtxcanv = new TCanvas("KandJpsiVtxcanv","KandJpsiVtxcanv");
KandJpsiVtxcanv->cd();
TH1F* KandJpsiVtx = new TH1F("KandJpsiVtx","KandJpsiVtx",100,0.,0.08);
myBpTree->Draw("IP3DKandJpsiVtx>>KandJpsiVtx");
KandJpsiVtx->Scale( 1.0/double(KandJpsiVtx->GetEntries()) );
//DrawLeg(BpJpsiPt,BsJpsiPt);
DrawAxis(KandJpsiVtx,"K track Jpsi vertex distance [cm]","Events");
KandJpsiVtxcanv->Update();

TCanvas *KandBpVtxcanv = new TCanvas("KandBpVtxcanv","KandBpVtxcanv");
KandBpVtxcanv->cd();
TH1F* KandBpVtx = new TH1F("KandBpVtx","KandBpVtx",100,0.,0.9);
myBpTree->Draw("BplusKIP3D>>KandBpVtx");
KandBpVtx->Scale( 1.0/double(KandBpVtx->GetEntries()) );

TH1F* BpLxyz = new TH1F("BpLxyz","BpLxyz",100,0.,0.9);
myBpTree->Draw("BplusLxyz>>BpLxyz","","same");
BpLxyz->Scale( 1.0/double(BpLxyz->GetEntries()) );
BpLxyz->SetLineColor(7);
//DrawLeg(BpJpsiPt,BsJpsiPt);
//DrawAxis(KandBpVtx,"K track Bp PV distance [cm]","Events");
KandBpVtxcanv->Update();
*/



//Bd optimised cuts
/*
TCut IPcut = "MuIP<0.4";
TCut Ptcut = IPcut && "MuPt>3.2";
TCut dRcut = Ptcut && "DeltaR>0.3";
*/

//Bd optimised cuts

TCut IPcut = "MuIP<0.3";
TCut Ptcut = IPcut && "MuPt>2.9 && MuPt < 20";
TCut dRcut = Ptcut && "DeltaR>0.3 && DeltaR <5";

TCut PFMuCut = dRcut && "PFMuon == 1 && MuCharge !=0";

// HUOM KAYTETTY LEIKATTUJA TIEDOSTOJA:  
// CutOptFileMCBpNew_PFMuons.root ja CutOptFileMCBsNew_PFMuons.root

//TFile *BpcutFile = TFile::Open("CutOptFiles/CutOptFileMCBpNew_PFMuons.root");
//TFile *BpcutFile = TFile::Open("CutOptFileMCBpNewTrig.root"); (old)
TFile *BpcutFile = TFile::Open("CutOptFiles/CutOptFileBpMCFinal.root");
TTree* BpcutTree = (TTree*)BpcutFile->Get("cutTree");

//TFile* BsPFMuFile = TFile::Open("CutOptFiles/CutOptFileMCBsNewTrig.root");
//TFile* BsPFMuFile = TFile::Open("CutOptFileMCBsNewTrig2.root"); (old)
TFile* BsPFMuFile = TFile::Open("CutOptFiles/CutOptFileBsMCFinal2.root");

TTree* BscutTree = (TTree*)BsPFMuFile->Get("cutTree");

//TFile* BdPFMuFile = TFile::Open("CutOptFiles/CutOptFileMCBdNewTrig.root");
//TFile* BdPFMuFile = TFile::Open("/tmp/terhi/MCBd2/CutOptFileMCBdNewTrig.root"); (old)
TFile* BdPFMuFile = TFile::Open("CutOptFiles/CutOptFileBdMCFinal.root");
TTree* BdcutTree = (TTree*)BdPFMuFile->Get("cutTree");

//TFile* BdDataMuFile = TFile::Open("/tmp/terhi/BdData/CutOptFileBdDataNewTrig.root"); (old)
TFile* BdDataMuFile = TFile::Open("CutOptFiles/CutOptFileBdDataFinal.root");
TTree* BdDatacutTree = (TTree*)BdDataMuFile->Get("cutTree");

//TFile *BdataCuttedFile = TFile::Open("SlimmedBdataNew.root"); (old)
TFile *BdataCuttedFile = TFile::Open("CutOptFiles/CutOptFileBpDataFinal.root");
TTree* BDatacutTree = (TTree*)BdataCuttedFile->Get("cutTree");

cout<< BDatacutTree->GetEntries() << endl; 
//BscutTree->Draw("MuPt>>pt(100,0,20)","(BsFlavour == -1) || (BsFlavour == 1)");


// IP Canvas
//impact parameter w.r.t signal B meson PV and tag muon

TCanvas* IPcanv = new TCanvas("IPcanv","IPcanv"); 
IPcanv->cd();

TH1F* IPBs = new TH1F("IPBs","IPBs",32,-0.01,0.31); // bin width 0.005 with Bs cuts 
BscutTree->Draw("MuIP>>IPBs",PFMuCut);

TH1F* IPBp = new TH1F("IPBp","IPBp",32,-0.01,0.31);
BpcutTree->Draw("MuIP>>IPBp",PFMuCut,"same");

TH1F* IPBpdata = new TH1F("IPBpdata","IPBpdata",32,-0.01,0.31);
BDatacutTree->Draw("MuIP>>IPBpdata",PFMuCut,"same");

TH1F* IPBd = new TH1F("IPBd","IPBd",32,-0.01,0.31);
BdcutTree->Draw("MuIP>>IPBd",PFMuCut,"same");

TH1F* IPBdData = new TH1F("IPBdData","IPBdData",32,-0.01,0.31);
//BdDatacutTree->Draw("MuIP>>IPBdData",PFMuCut,"same");

//Errors(IPBp,IPBs,IPBpdata,IPBd,IPBdData);
Errors(IPBp,IPBs,IPBpdata,IPBd);

Scale2(IPBs);
Scale2(IPBd);
Scale2(IPBp);
Scale2(IPBpdata);
Scale2(IPBdData);

MarkerStyles(IPBp,IPBpdata,IPBs,IPBd,IPBdData);
//DrawLeg(IPBp,IPBpdata,IPBs,IPBd,IPBdData);
DrawLeg(IPBp,IPBpdata,IPBs,IPBd);
SetCanvMargin(IPcanv);
DrawAxis(IPBs,"Tag muon 3D IP [cm]","Events / (0.01 cm)");
IPBs->GetYaxis()->SetRangeUser(0,0.2);
IPcanv->Update();
IPcanv->SaveAs("BpBsPFTagMuIP_Bdcuts2.pdf");


// Pt Canvas

TCanvas* PtPFMuCanv = new TCanvas("PtPFMuCanv","PtPFMuCanv");
PtPFMuCanv->cd();

TH1F* BsTagMuPt = new TH1F("BsTagMuPt","BsTagMuPt",22,1,21);
BscutTree->Draw("MuPt>>BsTagMuPt",PFMuCut);

TH1F* BplusTagMuPt = new TH1F("BplusTagMuPt","BplusTagMuPt",22,1.,21);
BpcutTree->Draw("MuPt>>BplusTagMuPt",PFMuCut,"same");

TH1F* BpdataTagMuPt = new TH1F("BpdataTagMuPt","BpdataTagMuPt",22,1.,21);
BDatacutTree->Draw("MuPt>>BpdataTagMuPt",PFMuCut,"same");

TH1F* BdTagMuPt = new TH1F("BdTagMuPt","BdTagMuPt",22,1,21);
BdcutTree->Draw("MuPt>>BdTagMuPt",PFMuCut,"same");

TH1F* BddataTagMuPt = new TH1F("BddataTagMuPt","BddataTagMuPt",22,1.,21);
//BdDatacutTree->Draw("MuPt>>BddataTagMuPt",PFMuCut,"same");

Errors(BplusTagMuPt,BsTagMuPt,BpdataTagMuPt,BdTagMuPt,BddataTagMuPt);
Scale2(BplusTagMuPt);
Scale2(BsTagMuPt);
Scale2(BdTagMuPt);
Scale2(BpdataTagMuPt);
Scale2(BddataTagMuPt);
//DrawLeg(BplusTagMuPt,BpdataTagMuPt,BsTagMuPt,BdTagMuPt,BddataTagMuPt);
MarkerStyles(BplusTagMuPt,BpdataTagMuPt,BsTagMuPt,BdTagMuPt,BddataTagMuPt);
DrawLeg(BplusTagMuPt,BpdataTagMuPt,BsTagMuPt,BdTagMuPt);
SetCanvMargin(PtPFMuCanv);
DrawAxis(BsTagMuPt,"Tag muon p_{T} [GeV/c]","Events / (GeV/c)");
//BplusTagMuPt->GetYaxis()->SetRangeUser(0,0.055);
PtPFMuCanv->Update();
PtPFMuCanv->SaveAs("BpBsPFTagMuPt_Bdcuts2.pdf");



//dR Canvas
// isolation between tag muon and signal B meson 
TCanvas* dRcanv = new TCanvas("dRcanv","dRcanv");
//dRcanv->cd();


TH1F* dRBs = new TH1F("dRBs","dRBs",20,0.1,5.1);
BscutTree->Draw("DeltaR>>dRBs",PFMuCut);

TH1F* dRBp = new TH1F("dRBp","dRBp",20,0.1,5.1);
BpcutTree->Draw("DeltaR>>dRBp",PFMuCut,"same");

TH1F* dRBpdata = new TH1F("dRBpdata","dRBpdata",20,0.1,5.1);
BDatacutTree->Draw("DeltaR>>dRBpdata",PFMuCut,"same");

TH1F* dRBd = new TH1F("dRBd","dRBd",20,0.1,5.1);
BdcutTree->Draw("DeltaR>>dRBd",PFMuCut,"same");

TH1F* dRBddata = new TH1F("dRBddata","dRBddata",20,0.1,5.1);
//BdDatacutTree->Draw("DeltaR>>dRBddata",PFMuCut,"same");

Errors(dRBs,dRBp,dRBpdata,dRBd,dRBddata);
Scale2(dRBs);
Scale2(dRBd);
Scale2(dRBp);
Scale2(dRBpdata);
Scale2(dRBddata);
MarkerStyles(dRBp,dRBpdata,dRBs,dRBd,dRBddata);
DrawLeg(dRBp,dRBpdata,dRBs,dRBd);
//DrawLeg(dRBp,dRBpdata,dRBs,dRBd,dRBddata);
SetCanvMargin(dRcanv);
dRBs->GetYaxis()->SetRangeUser(0,0.2);
DrawAxis(dRBs,"Tag muon isolation #DeltaR","Events / (0.25)");
dRcanv->Update();
dRcanv->SaveAs("BpBsPFTagMudR_Bdcuts2.pdf");

}
