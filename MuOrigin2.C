#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include <vector>
#include <iostream>
#include <set>

using namespace std;

Int_t maximumIndex(Double_t array[],Int_t isPF[],Int_t size){
 bool PFMuExists =false;
   
vector<int> PFMuIndices;
vector<double> PFPt;

     Int_t maxIndex = 0;       // start with max = first element
     if(size>25){size =25;}	    	    		

     for(Int_t i = 0; i<size; i++)
     {	
	 if(isPF[i] == 1 && array[i] >= 0){		
 		PFMuExists = true;
		PFMuIndices.push_back(i);
		PFPt.push_back(array[i]);           
	}	
     }

	
	for(Int_t k = 0; k<PFPt.size(); k++){
		
		if(PFPt[maxIndex] <= PFPt[k]){ 
			maxIndex = PFMuIndices[k];
	  	}
	}	
	
	if(PFMuExists == true){		
	return maxIndex;                // return index of max value  in array
	}
	
	else return -1;

}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
	
	
	Double_t deltaEta = eta1-eta2;
	Double_t deltaPhi = phi1-phi2;

//	cout << "deltaPhi " << deltaPhi << endl;
//	cout << "deltaEta " << deltaEta << endl;
 	while (deltaPhi > TMath::Pi()) deltaPhi -= TMath::TwoPi();
	while (deltaPhi <= -TMath::Pi()) deltaPhi += TMath::TwoPi();

	Double_t DeltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi );

	return DeltaR;
}


double FillHistos(TTree* myBsTree, bool BsStudy,bool BpStudy, bool BdStudy, TH1F* MuonID,TH1F* TrueMuonMom ){

// this is a list of all the PDG ids of B mesons
  std::set<int> listOfBmesonIds;
  listOfBmesonIds.insert(511 );   // Bd
  listOfBmesonIds.insert(521 );   // B+
  listOfBmesonIds.insert(10511 );    // B_0*0
  listOfBmesonIds.insert(10521 );    // B_0*+
  listOfBmesonIds.insert(513 );   // B*d
  listOfBmesonIds.insert(523 );   // B*d+
  listOfBmesonIds.insert(10513 );   // B1(L)0
  listOfBmesonIds.insert(10523 );   // B1(L)+
  listOfBmesonIds.insert(20513 );   // B1(H)0
  listOfBmesonIds.insert(20523 );   // B1(H)+
  listOfBmesonIds.insert(515 );    // B2*_0
  listOfBmesonIds.insert(525 );    // B2*_+
  listOfBmesonIds.insert(531 );   // Bs
  listOfBmesonIds.insert(10531 );    // B_s0*_0
  listOfBmesonIds.insert(533 );   // B*s
  listOfBmesonIds.insert(10533 );   // Bs1(L)0
  listOfBmesonIds.insert(20533 );   // Bs1(H)0
  listOfBmesonIds.insert(535 );    // Bs2*_0
  listOfBmesonIds.insert(541 );   // Bc+
  listOfBmesonIds.insert(10541 );   // B*c0+
  listOfBmesonIds.insert(543 );   // B*c+
  listOfBmesonIds.insert(10543 );   // Bc1(L)+
  listOfBmesonIds.insert(20543 );   // Bc1(H)+
  listOfBmesonIds.insert(545 );    // Bc2*_0
  
  listOfBmesonIds.insert(551 );   // etab(1S)
  listOfBmesonIds.insert(10551 );   // chib(1P)
  listOfBmesonIds.insert(100551 );   // etab(2S)
  listOfBmesonIds.insert(110551 );   // chib(2P)
  listOfBmesonIds.insert(200551 );   // etab(3S)
  listOfBmesonIds.insert(210551 );   // chib(3P)
  listOfBmesonIds.insert(553 );   // upsilon(1S)
  listOfBmesonIds.insert(10553 );   // hb(1P)
  listOfBmesonIds.insert(20553 );   // chib1(1P)
  listOfBmesonIds.insert(30553 );   // upsilon1(1D)
  listOfBmesonIds.insert(100553 );   // upsilon(2S)
  listOfBmesonIds.insert(110553 );   // hb(2P)
  listOfBmesonIds.insert(120553 );   // chib1(2P)
  listOfBmesonIds.insert(130553 );   // upsilon1(2D)
  listOfBmesonIds.insert(200553 );   // upsilon(3S)
  listOfBmesonIds.insert(210553 );   // hb(3P)
  listOfBmesonIds.insert(220553 );   // chib1(3P)
  listOfBmesonIds.insert(300553 );   // upsilon(4S)
  listOfBmesonIds.insert(9000553 );   // upsilon(10860)
  listOfBmesonIds.insert(9010553 );   // upsilon(11020)
  listOfBmesonIds.insert(555 );   // chib2(1P)
  listOfBmesonIds.insert(10555 );   // etab2(1D)
  listOfBmesonIds.insert(20555 );   // upsilon2(1D)
  listOfBmesonIds.insert(100555 );   // chib2(2P)
  listOfBmesonIds.insert(110555 );   // etab2(2D)
  listOfBmesonIds.insert(120555 );   // upsilon2(2D)
  listOfBmesonIds.insert(200555 );   // chib2(3P)
  listOfBmesonIds.insert(557 );   // upsilon3(1D)
  listOfBmesonIds.insert(100557 );   // upsilon3(2D)	

 const int Arraysize = 25;
 Double_t TagMuPt[Arraysize];
 Double_t MuP[Arraysize];
 Double_t MuIP3D[Arraysize];
 Double_t MuSimu[Arraysize];
 Int_t MuMomID[Arraysize];
 Int_t MuGrandMomID[Arraysize];
 Int_t isGlobalTagMuon[Arraysize];
 Int_t isTrackerTagMuon[Arraysize];
 Int_t isPFTagMuon[Arraysize];
 Double_t MuPhi[Arraysize];
 Double_t MuEta[Arraysize]; 
 Double_t MuCharge[Arraysize];
 Double_t DeltaR=-1;

 Int_t MuArraySize;
 
 Double_t BEta,BPhi;
 Int_t BsFlavour=0;
 Int_t BdFlavour=0;
 
 Int_t nentries = myBsTree->GetEntries();


 myBsTree->SetBranchAddress("TagMuRecoPt", TagMuPt);
 myBsTree->SetBranchAddress("TagMuRecoP", MuP); 	
 myBsTree->SetBranchAddress("TagMuRecoIP3D", MuIP3D);

 myBsTree->SetBranchAddress("TagMuListSize", &MuArraySize);
 myBsTree->SetBranchAddress("TagMuRecoChg", MuCharge);
 myBsTree->SetBranchAddress("GlobalTagMuon", isGlobalTagMuon);
 myBsTree->SetBranchAddress("TrackerTagMuon",isTrackerTagMuon);
 myBsTree->SetBranchAddress("PFTagMuon", isPFTagMuon);
 
 myBsTree->SetBranchAddress("TagMuRecoEta",MuEta);
 myBsTree->SetBranchAddress("TagMuRecoPhi",MuPhi); 


 if(BpStudy == true){
 myBsTree->SetBranchAddress("BplusEta",&BEta);
 myBsTree->SetBranchAddress("BplusPhi",&BPhi); 
}

 if(BsStudy == true){
 myBsTree->SetBranchAddress("BsIniFlavour",&BsFlavour);
 myBsTree->SetBranchAddress("BsFitEta",&BEta);
 myBsTree->SetBranchAddress("BsFitPhi",&BPhi);  
   
}


Double_t BdEta_Hyp1,BdEta_Hyp2, KstarMass_Hyp1, KstarMass_Hyp2,BdFitVtxProb_Hyp1;
Double_t BdPhi_Hyp1,BdPhi_Hyp2, BdFitM_Hyp1,BdFitM_Hyp2,BdFitVtxProb_Hyp2;

 if(BdStudy == true){
//Bd variables
 
 myBsTree->SetBranchAddress("BdIniFlavour",&BdFlavour); 
 myBsTree->SetBranchAddress("KstarMass_nofit_Hyp1",&KstarMass_Hyp1); 
 myBsTree->SetBranchAddress("KstarMass_nofit_Hyp2",&KstarMass_Hyp2);

 myBsTree->SetBranchAddress("BdFitEta_Hyp1",&BdEta_Hyp1);
 myBsTree->SetBranchAddress("BdFitEta_Hyp2",&BdEta_Hyp2);
 myBsTree->SetBranchAddress("BdFitPhi_Hyp1",&BdPhi_Hyp1);
 myBsTree->SetBranchAddress("BdFitPhi_Hyp2",&BdPhi_Hyp2);
 myBsTree->SetBranchAddress("BdFitM_Hyp1",&BdFitM_Hyp1);
 myBsTree->SetBranchAddress("BdFitM_Hyp2",&BdFitM_Hyp2);

 myBsTree->SetBranchAddress("BdFitVtxProb_Hyp1",&BdFitVtxProb_Hyp1);
 myBsTree->SetBranchAddress("BdFitVtxProb_Hyp2",&BdFitVtxProb_Hyp2);
}

//Bd variables

 //MC Matching
 myBsTree->SetBranchAddress("TagMuSimu",MuSimu);
 myBsTree->SetBranchAddress("MuRecoMCmother",MuMomID);
 myBsTree->SetBranchAddress("MuRecoMCMotherMother",MuGrandMomID);

	int Picounter = 0;
	int Kcounter = 0;
	int Dmomcounter =0;
 	int trueMuoncounter=0;  	
	int allcutpassedtagmuons =0;   

	int Bpmomcounter = 0;
	int Bsmomcounter = 0;
	int Bdmomcounter = 0;
	int Bmomcounter = 0;

	int mysterycounter =0;
	int noMCMatchCounter=0;
	int alltagmuons =0; 		
	int JPsicounter =0; 
	int BgrandmomCounter =0;


	Double_t BdFitM,BdFitVtxProb,KstarMass_nofit;
	Double_t nominalKstarMass = 0.89166; 

 		 //nentries
	for(Int_t i = 0; i<nentries; i++){
		
		myBsTree->GetEntry(i);
		
		//Bd study start
   if(BdStudy == true){
		BdFitM= -999;
		BdFitVtxProb = -999;
		KstarMass_nofit = -999;

	 if(TMath::Abs(KstarMass_Hyp1 - nominalKstarMass) < TMath::Abs(KstarMass_Hyp2 - nominalKstarMass) ) {
	       KstarMass_nofit = KstarMass_Hyp1;
	       BdFitVtxProb = BdFitVtxProb_Hyp1;
	       BdFitM = BdFitM_Hyp1; 
	       BEta = BdEta_Hyp1;
	       BPhi = BdPhi_Hyp1;
         }
   
         else {
      	 KstarMass_nofit = KstarMass_Hyp2;
      	 BdFitVtxProb = BdFitVtxProb_Hyp2;
         BdFitM = BdFitM_Hyp2; 
	 BEta = BdEta_Hyp2;
	 BPhi = BdPhi_Hyp2;
   
     
     }
		if(BdFitVtxProb < 0.02) continue;
		if(BdFitM > 5.55) continue;
		if(BdFitM < 5.0) continue;
		if(KstarMass_nofit < 0.796) continue;	
		if(KstarMass_nofit > 0.996) continue;		


		//Bd study end	
  }		


		if(MuArraySize <= 0) continue;		 
		
		Int_t maxPtMuonIndex = maximumIndex(TagMuPt,isPFTagMuon,MuArraySize);
		if(maxPtMuonIndex == -1) {
			continue;
		}
			
		alltagmuons++;
			

		if(BEta < -999 || MuEta[maxPtMuonIndex] == -999 ) continue; 
		if(BPhi < -999 || MuPhi[maxPtMuonIndex] == -999 ) continue;

		DeltaR = deltaR(BEta,BPhi,MuEta[maxPtMuonIndex],MuPhi[maxPtMuonIndex]);

		if(DeltaR < 0.3) continue;
//		if(MuIP3D[maxPtMuonIndex] < 0) continue;		
		//if(MuPt[maxPtMuonIndex] < 0 ) continue; 
		if(MuIP3D[maxPtMuonIndex] > 0.3 || MuIP3D[maxPtMuonIndex] < 0) continue;
		
		if(TagMuPt[maxPtMuonIndex] < 2.9 ) continue; //1.7/1.9 

		allcutpassedtagmuons++;
		
		if(MuSimu[maxPtMuonIndex] == -999){ 
			noMCMatchCounter++;
			continue;
		}
		
		MuonID->Fill(TMath::Abs(MuSimu[maxPtMuonIndex])); 
//		MuonMomID->Fill(TMath::Abs(MuMomID[maxPtMuonIndex]) ); 
 
		//Muon ID
		if( TMath::Abs(MuSimu[maxPtMuonIndex]) == 211 || TMath::Abs(MuSimu[maxPtMuonIndex]) == 213 || TMath::Abs(MuSimu[maxPtMuonIndex]) == 331 ){
			Picounter++;
/*			PionP->Fill(MuP[maxPtMuonIndex]);
			PionPt->Fill(TagMuPt[maxPtMuonIndex]);
			PiondR->Fill( DeltaR );
			PionIP->Fill(MuIP3D[maxPtMuonIndex]);
*/			
		}

		if( TMath::Abs(MuSimu[maxPtMuonIndex]) == 321 || TMath::Abs(MuSimu[maxPtMuonIndex]) == 323  ){
		
/*			KaonP->Fill(MuP[maxPtMuonIndex]);
			KaonPt->Fill(TagMuPt[maxPtMuonIndex]);
			KaondR->Fill( DeltaR );
			KaonIP->Fill(MuIP3D[maxPtMuonIndex]);
*/
			Kcounter++;
		}

//		if(TMath::Abs(MuSimu[maxPtMuonIndex]) != 211 && TMath::Abs(MuSimu[maxPtMuonIndex]) != 321 && TMath::Abs(MuSimu[maxPtMuonIndex]) != 13){ cout << " Simu MC ID " << MuSimu[maxPtMuonIndex] << endl; }


		if(TMath::Abs(MuSimu[maxPtMuonIndex]) == 13 ){
//		if( (MuSimu[maxPtMuonIndex]  == -13 && BdFlavour == -1) || (MuSimu[maxPtMuonIndex]  == 13 && BdFlavour == 1) ){
			trueMuoncounter++;

			if(MuMomID[maxPtMuonIndex] != -999){
				TrueMuonMom->Fill(TMath::Abs(MuMomID[maxPtMuonIndex]) );
				
//				MuonGMomID->Fill(TMath::Abs(MuGrandMomID[maxPtMuonIndex]));

			if(TMath::Abs(MuMomID[maxPtMuonIndex] ) == 443 ) {
				JPsicounter++ ;	
			}

			   if( TMath::Abs(MuMomID[maxPtMuonIndex] ) == 411 || TMath::Abs(MuMomID[maxPtMuonIndex] ) == 421 ||  TMath::Abs(MuMomID[maxPtMuonIndex] ) == 431){


				Dmomcounter++;
/*				DmuonPt->Fill(TagMuPt[maxPtMuonIndex]);
				DmuonP->Fill(MuP[maxPtMuonIndex]);
				DmuondR->Fill( DeltaR );
				DmuonIP->Fill(MuIP3D[maxPtMuonIndex]);
*/			
				

				if(listOfBmesonIds.find( TMath::Abs(MuGrandMomID[maxPtMuonIndex]) ) != listOfBmesonIds.end()) {BgrandmomCounter++;}

			   }

			   if( listOfBmesonIds.find( TMath::Abs(MuMomID[maxPtMuonIndex]) ) != listOfBmesonIds.end() ){
/*				BmuonP->Fill(MuP[maxPtMuonIndex]);
				BmuonPt->Fill(TagMuPt[maxPtMuonIndex]);
				BmuondR->Fill( DeltaR );
				BmuonIP->Fill(MuIP3D[maxPtMuonIndex]);
*/
				Bmomcounter++;
			   }
		
			if(TMath::Abs(MuMomID[maxPtMuonIndex]) == 531) {Bsmomcounter++;}
			if(TMath::Abs(MuMomID[maxPtMuonIndex]) == 511) {Bdmomcounter++;}
			if(TMath::Abs(MuMomID[maxPtMuonIndex]) == 521) {Bpmomcounter++;}
			
			}	
		}	
		//Muon mother ID
		//should all the neutral B mesons be excluded?
		//cout << "ad" << endl;
 		if(TMath::Abs( MuMomID[maxPtMuonIndex] )== 92){
			mysterycounter++;
		}
		
	}

  if(BdStudy == true){ 
	cout << "Bd Study " << endl;
	cout << "" <<endl;
  } 

  if(BpStudy == true){ 
	cout << "Bp Study " << endl;
	cout << "" <<endl;
  } 

  if(BsStudy == true){ 
	cout << "Bs Study " << endl;
	cout << "" <<endl;
  }


  Double_t xAxis[7] = {-0.05, 0, 0.0055, 0.015, 0.035, 0.1, 0.3}; 
   TEfficiency *muonEff = new TEfficiency("muonEff","",6,xAxis);
     
  muonEff->SetTotalEvents(0,allcutpassedtagmuons); 
  muonEff->SetPassedEvents(0,noMCMatchCounter);   
  double dnomatched = (muonEff->GetEfficiencyErrorLow(0) + muonEff->GetEfficiencyErrorUp(0) )/2.0;

  muonEff->SetTotalEvents(1,allcutpassedtagmuons);  
  muonEff->SetPassedEvents(1,trueMuoncounter);   

  double dtruemuons = (muonEff->GetEfficiencyErrorLow(1) + muonEff->GetEfficiencyErrorUp(1) )/2.0;

  muonEff->SetTotalEvents(2,allcutpassedtagmuons);  
  muonEff->SetPassedEvents(2,Kcounter);   

  double dkaons = (muonEff->GetEfficiencyErrorLow(2) + muonEff->GetEfficiencyErrorUp(2) )/2.0;

  muonEff->SetTotalEvents(4,allcutpassedtagmuons);  
  muonEff->SetPassedEvents(4,Picounter);   

  double dpions = (muonEff->GetEfficiencyErrorLow(4) + muonEff->GetEfficiencyErrorUp(4) )/2.0;

  muonEff->SetTotalEvents(5,allcutpassedtagmuons);  
  muonEff->SetPassedEvents(5,Bmomcounter+BgrandmomCounter);   

  double dBmom = (muonEff->GetEfficiencyErrorLow(5) + muonEff->GetEfficiencyErrorUp(5) )/2.0;

  cout << "no matched tag muons " << muonEff->GetEfficiency(0) << " +- " << dnomatched << endl;
  cout << "true tag muons " << muonEff->GetEfficiency(1) << " +- " << dtruemuons << endl;

  cout << "misid' kaons " << muonEff->GetEfficiency(2) << " +- " << dkaons << endl;

 cout << "misid' pions " << muonEff->GetEfficiency(4) << " +- " << dpions << endl;

 cout << "muons from B mom " << muonEff->GetEfficiency(5) << "+- " << dBmom  << endl; 

cout << "MC matched muons among cut passed muons " << allcutpassedtagmuons-noMCMatchCounter << endl;
//cout << "no MC matched muons among cut passed muons " << noMCMatchCounter << endl;
cout << "true tag muons " << trueMuoncounter << endl;
cout << "muons from Bplus mom " << double(Bpmomcounter)  << endl; 
cout << "muons from Bs mom " << double(Bsmomcounter)  << endl; 
cout << "muons from Bd mom " << double(Bdmomcounter) << endl; 

//cout << "Kaons in tag muons " << Kcounter << endl; 
//cout << "Pions in tag muons " << Picounter << endl; 
cout << "all cut passed tag muons " << allcutpassedtagmuons << endl; 	

cout << "tag muon B moms" << Bmomcounter << endl;
//cout << "tag muon D moms" << Dmomcounter << " with B grandmoms " << BgrandmomCounter << endl;
//cout << "tag muon Jpsi moms" << JPsicounter << endl;
//cout << "all tag muons " << alltagmuons << endl;


return allcutpassedtagmuons-noMCMatchCounter;
}

void SetCanvMargin(TCanvas *c1){
c1->SetLeftMargin(0.1457286);
c1->SetBottomMargin(0.1541463);
}

void Errors(TH1I* histo1, TH1I* histo2, TH1I* histo3=0){
histo1->Sumw2();
histo2->Sumw2();
histo3->Sumw2();
}

void Scale2(TH1I* histoSmaller){

 if(histoSmaller->GetEntries() != 0 ){

histoSmaller->Scale( 1.0/double(histoSmaller->GetEntries()) );
 }
}

void LineStyles(TH1* BpMC, TH1* BsMC=0,TH1* BdMC=0){

BpMC->SetLineColor(2);
BpMC->SetLineWidth(2);
BpMC->SetFillColor(2);
 if(BsMC!=0){
  BsMC->SetFillColor(1); 
  BsMC->SetLineColor(1);
  BsMC->SetLineWidth(2); 
 }
if(BdMC!=0){
  BdMC->SetFillColor(3); 
  BdMC->SetLineColor(3);
  BdMC->SetLineWidth(2);
 }

}

void DrawLeg(TH1* BpMCHisto,TH1* BsHisto=0, TH1* BdHisto =0){
//TLegend *leg = new TLegend(0.7173367,0.7108014,0.8969849,0.8954704,NULL,"brNDC");
TLegend *leg = new TLegend(0.7223618,0.6341463,0.8932161,0.8885017,NULL,"brNDC");
leg->SetFillStyle(0);
leg->SetBorderSize(0);
 if(BsHisto!=0) {
  leg->AddEntry(BsHisto,"B_{s} MC");
 }
leg->AddEntry(BdHisto,"B^{0} MC");
leg->AddEntry(BpMCHisto,"B^{+} MC");
leg->Draw("same");
}

void DrawAxis(THStack *histo, TString xAxis, TString yAxis){
histo->GetXaxis()->SetTitle(xAxis);
histo->GetYaxis()->SetTitle(yAxis);
}

void MuOrigin2(){

  gStyle->SetOptStat(11111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetTitleSize(0.07,"t");  

 gStyle->SetLabelFont(42, "XYZ");
 gStyle->SetLabelOffset(0.007, "XYZ");
 gStyle->SetLabelSize(0.05, "XYZ");

  
//  TFile *BpFile = TFile::Open("/tmp/terhi/JetVectorTestCleaned.root");
  TFile *BpFile = TFile::Open("/tmp/terhi/BplusMCCleanedNew.root");
  TTree* BpTree = (TTree*)BpFile->Get("BsTree");
 
  TFile *BsFile = TFile::Open("/tmp/terhi/BsMCCleanedNew.root");
 // TFile *BsFile = TFile::Open("/tmp/terhi/BsMC/BsMCRun2Cleaned.root");
  TTree* BsTree = (TTree*)BsFile->Get("BsTree"); 
/*
  TFile *BdFile = TFile::Open("/tmp/terhi/BdMCPreCleaned.root");
  TTree* BdTree = (TTree*)BdFile->Get("BsTree"); 
*/
	
	TH1F* MuonIDBs = new TH1F("MuonIDBs","Muon ID Bs",351,-0.5,350.5);
	TH1F* TrueMuonMomBs = new TH1F("TrueMuonMomBs","mother of true muons Bs", 1001,-0.5,1000.5);
//	TH1I* MuonMomIDBs = new TH1I("MuonMomIDBs","Muon mother ID Bs", 601,-0.5,600.5);
//	TH1I* MuonGMomIDBs = new TH1I("MuonGMomIDBs","Muon grand mother ID Bs", 1000,-1000,1000);

	TH1F* MuonIDBp = new TH1F("MuonIDBp","Muon ID Bp",351,-0.5,350.5);
	TH1F* TrueMuonMomBp = new TH1F("TrueMuonMomBp","mother of true muons Bp", 1001,-0.5,1000.5);

	TH1F* MuonIDBd = new TH1F("MuonIDBd","Muon ID Bd",351,-0.5,350.5);
	TH1F* TrueMuonMomBd = new TH1F("TrueMuonMomBd","mother of true muons Bd", 1001,-0.5,1000.5);

double MCmathchedBs = FillHistos(BsTree, true, false, false, MuonIDBs, TrueMuonMomBs);
double MCmathchedBp = FillHistos(BpTree, false, true, false, MuonIDBp, TrueMuonMomBp);
//double MCmathchedBd = FillHistos(BdTree, false, false, true, MuonIDBd, TrueMuonMomBd);

TCanvas* MuID = new TCanvas("MuonMCID","MuonMCID");
  
THStack *hs = new THStack("MuIDstack","MuIDstack"); 
/*
cout << MCmathchedBs << endl;

cout << "jee 1" << endl;
//Errors(MuonIDBp,MuonIDBs,MuonIDBd);
MuonIDBs->Scale(1.0/double(MCmathchedBs));
MuonIDBd->Scale(1.0/double(MCmathchedBd));
MuonIDBp->Scale(1.0/double(MCmathchedBp));

LineStyles(MuonIDBp,MuonIDBs,MuonIDBd);

cout << "jee 2" << endl;

hs->Add(MuonIDBs);
hs->Add(MuonIDBd);
hs->Add(MuonIDBp);

cout << "jee 3" << endl;

SetCanvMargin(MuID);

hs->Draw();
DrawAxis(hs,"Tag muon MC ID","Events");
DrawLeg(MuonIDBp,MuonIDBs,MuonIDBd);
cout << "jee 4" << endl;
MuID->Update();

TCanvas* MuMomCanv = new TCanvas("MuonMomMCID","MuonMomMCID");
THStack *hs2 = new THStack("MuMomIDstack","MuMomIDstack"); 

//cout << MCmathchedBs << endl;

cout << "jee 1" << endl;
//Errors(MuonIDBp,MuonIDBs,MuonIDBd);
TrueMuonMomBs->Scale(1.0/double(4324));
TrueMuonMomBd->Scale(1.0/double(1996));
TrueMuonMomBp->Scale(1.0/double(3580));

LineStyles(TrueMuonMomBp,TrueMuonMomBs,TrueMuonMomBd);

cout << "jee 2" << endl;

hs2->Add(TrueMuonMomBs);
hs2->Add(TrueMuonMomBd);
hs2->Add(TrueMuonMomBp);

cout << "jee 3" << endl;
//DrawLeg(MuonIDBp,MuonIDBs,MuonIDBd);
SetCanvMargin(MuMomCanv);

hs2->Draw();
DrawAxis(hs2,"Tag muon mother MC ID","% from number of true muons");
DrawLeg(TrueMuonMomBp,TrueMuonMomBs,TrueMuonMomBd);
cout << "jee 4" << endl;
MuMomCanv->Update();
*/
/*
	TH1F* KaonPt = new TH1F("KaonPt","Pt of Kaons",100,0,40);
	TH1F* KaonP = new TH1F("KaonP","P of Kaons",120,0,60);	
	TH1F* KaondR = new TH1F("KaondR","dR of Kaons",100,0,5);
	TH1F* KaonIP = new TH1F("KaonIP","IP of Kaons",100,0,0.5);

	TH1F* PionPt = new TH1F("PionPt","Pt of Pions",100,0,40);
	TH1F* PionP = new TH1F("PionP","P of Pions",120,0,60);	
	TH1F* PiondR = new TH1F("PiondR","dR of Pions",100,0,5);
	TH1F* PionIP = new TH1F("PionIP","IP of Pions",100,0,0.5);
	

	TH1F* DmuonPt = new TH1F("DmuonPt","Pt of D mesons muons",100,0,40);
	TH1F* DmuonP = new TH1F("DmuonP","P of D mesons muons",120,0,60);	
	TH1F* DmuondR = new TH1F("DmuondR","dR of D mesons muons",100,0,5);
	TH1F* DmuonIP = new TH1F("DmuonIP","IP of D mesons muons",100,0,0.5);

	TH1F* BmuonPt = new TH1F("BmuonPt","Pt of B mesons muons",100,0,40);
	TH1F* BmuonP = new TH1F("BmuonP","P of B mesons muons",120,0,60);	
	TH1F* BmuondR = new TH1F("BmuondR","dR of B mesons muons",100,0,5);
	TH1F* BmuonIP = new TH1F("BmuonIP","IP of B mesons muons",100,0,0.5);
*/
		
 

/*
 TCanvas* c1 = new TCanvas("Muon ID","Muon ID",10,10,700,700);	
 c1->Divide(3,2);

 c1->cd(1); 
// totalTagsIP->Draw();
MuonID->Draw();
 c1->cd(2); 
TrueMuonMom->Draw();
 c1->cd(3); 
MuonGMomID->Draw();
*/
/*
DmuonP->Scale(1.0/double(DmuonP->GetEntries()));
DmuonP->Draw();
BmuonP->Scale(1.0/double(BmuonP->GetEntries()));
BmuonP->SetLineColor(2); 
BmuonP->Draw("same");
*/

/*
 c1->cd(4); 
DmuonIP->Scale(1.0/double(DmuonIP->GetEntries()));
DmuonIP->Draw();
BmuonIP->Scale(1.0/double(BmuonIP->GetEntries()));
BmuonIP->SetLineColor(2); 
BmuonIP->Draw("same");

 c1->cd(5);
DmuonPt->Scale(1.0/double(DmuonPt->GetEntries()));
DmuonPt->Draw();
BmuonPt->Scale(1.0/double(BmuonPt->GetEntries()));
BmuonPt->SetLineColor(2); 
BmuonPt->Draw("same");
 c1->cd(6);
DmuondR->Scale(1.0/double(DmuondR->GetEntries()));
DmuondR->Draw();
BmuondR->Scale(1.0/double(BmuondR->GetEntries()));
BmuondR->SetLineColor(2); 
BmuondR->Draw("same");

 TCanvas* c2 = new TCanvas("Muon misID'ed","Muon misID'ed");	
 c2->Divide(3,2);

 c2->cd(1);
PionIP->Scale( 1.0/double(PionIP->GetEntries()) );
PionIP->Draw();
BmuonIP->Draw("same");

 c2->cd(2);
PionPt->Scale( 1.0/double(PionPt->GetEntries()) );
PionPt->Draw();
BmuonPt->Draw("same");

 c2->cd(3);
PiondR->Scale( 1.0/double(PiondR->GetEntries()) );
PiondR->Draw();
BmuondR->Draw("same");

 c2->cd(4);
KaonIP->Scale( 1.0/double(KaonIP->GetEntries()) );
KaonIP->Draw();
BmuonIP->Draw("same");

 c2->cd(5);
KaonPt->Scale( 1.0/double(KaonPt->GetEntries()) );
KaonPt->Draw(); 
BmuonPt->Draw("same");

 c2->cd(6);
KaondR->Scale( 1.0/double(KaondR->GetEntries()) );
KaondR->Draw();
BmuondR->Draw("same");

*/
}
