#include "TMath.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TROOT.h"
#include "TStyle.h"

Int_t DauCounter = 0;

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

	
	for(size_t k = 0; k<PFPt.size(); k++){
		
		if(PFPt[maxIndex] <= PFPt[k]){ 
			maxIndex = PFMuIndices[k];
	  	}
	}	
	
	if(PFMuExists == true){		
	return maxIndex;                // return index of max value  in array
	}
	
	else return -1;

}


Int_t isOppoBdMuon(Int_t GenNumberOfBdecays_, Int_t GenNumberOfDaughters_[],Int_t BDauIdMC_[10][15], Int_t BmesonsId_[]){


 bool JpsiFound = false;
 bool KstarFound = false;

 int signalflag = 0;
 int BdSigCounter =0;
 int BdOppoCounter=0;


//	cout << " inside oppo method 1 " << endl; 

	for(int x=0; x < GenNumberOfBdecays_; x++){
	    int B_id = abs(BmesonsId_[x]);
            int BpId = BmesonsId_[x]; 
//	cout << "B Id " << BpId << endl;
	
	    JpsiFound = false;	
	    KstarFound =false;	

		if(GenNumberOfDaughters_[x] > 5) {DauCounter++;}
	   		
		for(int j=0; j < GenNumberOfDaughters_[x]; j++){
			int Dau_id = BDauIdMC_[x][j];
     
//			cout << Dau_id << endl;		
			if(Dau_id == 443) JpsiFound = true;
			if(abs(Dau_id) == 313) KstarFound = true;
		} 

	  if( B_id == 511 && JpsiFound  == true && KstarFound == true && GenNumberOfDaughters_[x] == 2 && signalflag ==0){ 
	BdSigCounter++; signalflag =1; // cout << "B+ channel!" << endl;
	  }
	
	 else if(B_id == 511){ 
	   BdOppoCounter++;	   
	  }			
	}

//	cout << " BdOppoCounter " <<BdOppoCounter << endl; 

	if(BdOppoCounter!=0){return 1;}
	else return 0;
}



Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
	

	Double_t deltaEta = eta1-eta2;
	Double_t deltaPhi = phi1-phi2;

 	while (deltaPhi > TMath::Pi()) deltaPhi -= TMath::TwoPi();
	while (deltaPhi <= -TMath::Pi()) deltaPhi += TMath::TwoPi();

	Double_t DeltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi );

	return DeltaR;
}


void RootFileProducerBdMCTrueMuons(){
 
gROOT->SetStyle("Plain");
//gStyle->SetFillColor(0);  //white fill color
gStyle->SetFrameBorderMode(0);  //no frame border
gStyle->SetCanvasBorderMode(0);  //no canvas border
gStyle->SetOptTitle(0);  //no title
gStyle->SetOptStat(1111111);  //no statistics _or_
gStyle->SetCanvasColor(0);
gStyle->SetStatBorderSize(0);

TFile *MCPUFile = TFile::Open("/tmp/terhi/BdMCPreCleaned.root");
TTree* myBsTree = (TTree*)MCPUFile->Get("BsTree");
//myBsTree->Print(); 
 Double_t nominalKstarMass = 0.89166;  
 const int Arraysize = 25;

 Double_t TagMuPt[Arraysize];
 Double_t TagMuPtRel[Arraysize];
 Double_t MuIP3D[Arraysize];
 Double_t TagMuCharge[Arraysize];

 Int_t MuSimuMom[Arraysize];
 Int_t MuSimuGMom[Arraysize];
 Int_t isPFTagMuon[Arraysize];
 Double_t MuPhi[Arraysize];
 Double_t MuEta[Arraysize]; 
 Double_t MuSimu[Arraysize];

 Int_t MuArraySize;
 Double_t BdCt2DMC; 
 Double_t BdFitVtxProb_Hyp1,BdFitVtxProb_Hyp2;
 Double_t BdEta_Hyp1,BdEta_Hyp2, KstarMass_Hyp1, KstarMass_Hyp2;
 Double_t BdPhi_Hyp1,BdPhi_Hyp2, BdFitM_Hyp1,BdFitM_Hyp2;
 Int_t BdIniFlavour,track1Charge, triggerv11,triggerv12; 
 Int_t nentries = myBsTree->GetEntries();
 
 myBsTree->SetBranchAddress("TagMuRecoPt", TagMuPt);
 myBsTree->SetBranchAddress("TagMuRecoPtRel", TagMuPtRel); 
 myBsTree->SetBranchAddress("TagMuRecoIP3D", MuIP3D); 
 myBsTree->SetBranchAddress("TagMuListSize", &MuArraySize);
 myBsTree->SetBranchAddress("TagMuRecoChg", TagMuCharge);

 myBsTree->SetBranchAddress("PFTagMuon", isPFTagMuon);
 myBsTree->SetBranchAddress("TagMuRecoEta",MuEta);
 myBsTree->SetBranchAddress("TagMuRecoPhi",MuPhi);
 myBsTree->SetBranchAddress("BdIniFlavour", &BdIniFlavour); 
 myBsTree->SetBranchAddress("BdCt2DMC", &BdCt2DMC); 
 myBsTree->SetBranchAddress("triggerbit_4Jpsi_Displaced_v11", &triggerv11); 
 myBsTree->SetBranchAddress("triggerbit_4Jpsi_Displaced_v12", &triggerv12);  

 myBsTree->SetBranchAddress("BdTrack1Charge", &track1Charge); 
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
 myBsTree->SetBranchAddress("MuRecoMCmother",MuSimuMom);
 myBsTree->SetBranchAddress("MuRecoMCMotherMother",MuSimuGMom); 
 myBsTree->SetBranchAddress("TagMuSimu",MuSimu);

 Int_t BDauIdMC_[10][15];
 Int_t BmesonsId_[10];
 Int_t GenNumberOfBdecays_;
 Int_t GenNumberOfDaughters_[10];
 myBsTree->SetBranchAddress("BDauIdMC", BDauIdMC_ );
 myBsTree->SetBranchAddress("GenNumberOfBdecays", &GenNumberOfBdecays_ );
 myBsTree->SetBranchAddress("BmesonsId", BmesonsId_  );
 myBsTree->SetBranchAddress("GenNumberOfDaughters", GenNumberOfDaughters_  );


	Double_t DeltaR, MuPt, MuIP, MuPtRel,BdFitM,BdecayFlavour=-999,SimuMuon;
	Int_t MuCharge,PFMuon,MCMuonGMom,MCMuonMom,OppoBdMeson;

 	TFile f("MuMcOppoBd.root","RECREATE");
        TTree tree("cutTree","cutTree");
        tree.Branch("DeltaR",&DeltaR,"DeltaR/D");  
	tree.Branch("MuPt",&MuPt,"MuPt/D");
	tree.Branch("MuPtRel",&MuPtRel,"MuPtRel/D");    
	tree.Branch("MuIP",&MuIP,"MuIP/D");  
	tree.Branch("MuCharge",&MuCharge,"MuCharge/I");
	tree.Branch("BdecayFlavour",&BdecayFlavour,"BdecayFlavour/I");
	tree.Branch("BdFlavour",&BdIniFlavour,"BdFlavour/I");
	tree.Branch("PFMuon", &PFMuon,"PFMuon/I");
//	tree.Branch("GlobalMuon", &GlobalMuon,"GlobalMuon/I");
//	tree.Branch("TrackerMuon", &TrackerMuon,"TrackerMuon/I");
	tree.Branch("BpMass", &BdFitM,"BpMass/D");
	tree.Branch("BdCt2DMC",&BdCt2DMC,"BdCt2DMC/D");	
	tree.Branch("SimuMuon", &SimuMuon,"SimuMuon/D");
	tree.Branch("MCMuonMom", &MCMuonMom,"MCMuonMom/I");
	tree.Branch("MCMuonGMom", &MCMuonGMom,"MCMuonGMom/I");
	tree.Branch("OppoBdMeson",&OppoBdMeson,"OppoBdMeson/I");  

	Double_t KstarMass_nofit=-999, BdFitVtxProb=-999,BdEta=-999, BdPhi = -999;

	Int_t couterMuArraySize=0, CounterCharge=0, CounterMudR =0, CounterMuIP=0, CounterMaxInd =0 ,CounterFlavour=0 ;

			//nentries
	for(Int_t i = 0; i<nentries; i++){

		BdFitM = -999;
		DeltaR= -999;
		MuPt= -999;
		PFMuon = 0;
		MuCharge = 0;
		MuIP = -999;
		BdecayFlavour = -999;
		MuPtRel = -999;
		SimuMuon = -999;
		MCMuonGMom = -999;
		MCMuonMom = -999;
		OppoBdMeson =0;
		myBsTree->GetEntry(i);
		
		
	 if(TMath::Abs(KstarMass_Hyp1 - nominalKstarMass) < TMath::Abs(KstarMass_Hyp2 - nominalKstarMass) ) {
	       KstarMass_nofit = KstarMass_Hyp1;
	       BdFitVtxProb = BdFitVtxProb_Hyp1;
	       BdFitM = BdFitM_Hyp1; 
	       BdEta = BdEta_Hyp1;
	       BdPhi = BdPhi_Hyp1;

		if(track1Charge==1){BdecayFlavour=1;}
		else BdecayFlavour=-1;     
         }
   
         else {
      	 KstarMass_nofit = KstarMass_Hyp2;
      	 BdFitVtxProb = BdFitVtxProb_Hyp2;
         BdFitM = BdFitM_Hyp2; 
	 BdEta = BdEta_Hyp2;
	 BdPhi = BdPhi_Hyp2;

	if(track1Charge==1){BdecayFlavour=-1;}
	else BdecayFlavour=1;     
     
     }
		if(BdFitVtxProb < 0.02) continue;
		if(BdFitM > 5.55) continue;
		if(BdFitM < 5.0) continue;
		if(KstarMass_nofit < 0.796) continue;	
		if(KstarMass_nofit > 0.996) continue;	

		OppoBdMeson = isOppoBdMuon(GenNumberOfBdecays_, GenNumberOfDaughters_,BDauIdMC_,BmesonsId_);
		
		if(MuArraySize <= 0) {
			couterMuArraySize++;
			tree.Fill(); 			
			continue;
		} 
		
					
		Int_t maxPtMuonIndex = maximumIndex(TagMuPt,isPFTagMuon,MuArraySize);
		if(maxPtMuonIndex == -1) {
			CounterMaxInd++;
			tree.Fill();  
			continue;
		}
 
		SimuMuon = MuSimu[maxPtMuonIndex];
		MCMuonMom = MuSimuMom[maxPtMuonIndex];
		MCMuonGMom = MuSimuGMom[maxPtMuonIndex];

		MuCharge = TagMuCharge[maxPtMuonIndex];
		if(MuCharge!=-1 && MuCharge!=1){
			MuCharge = 0;
			CounterCharge++;
			tree.Fill();  
			continue;
		}		

		if(BdIniFlavour!=1 && BdIniFlavour!=-1){
			MuCharge = 0;
			CounterFlavour++;
			tree.Fill();  
			continue;
		}
		 
		if(MuIP3D[maxPtMuonIndex] < 0 || TagMuPt[maxPtMuonIndex] < 0){ 
			MuCharge = 0; 
			CounterMuIP++;
			tree.Fill();  
			continue;	
		}
			
		if(MuEta[maxPtMuonIndex] == -999 || MuPhi[maxPtMuonIndex] == -999){
			MuCharge = 0;
			CounterMudR++;
			tree.Fill();  
			continue;		
		}

		if(BdEta == -999 || BdPhi == -999) {
			CounterMudR++;
			MuCharge =0;
			tree.Fill();  
			continue;
		}

 
		DeltaR = deltaR(BdEta,BdPhi,MuEta[maxPtMuonIndex],MuPhi[maxPtMuonIndex]);
//		TrackerMuon = isTrackerTagMuon[maxPtMuonIndex];
//		GlobalMuon = isGlobalTagMuon[maxPtMuonIndex]; 
		PFMuon = isPFTagMuon[maxPtMuonIndex];
		MuIP = MuIP3D[maxPtMuonIndex];		
		MuPt = TagMuPt[maxPtMuonIndex];
		MuPtRel = TagMuPtRel[maxPtMuonIndex];		
	
		tree.Fill();
		
		
	}
	tree.Print();
	tree.Write(); 

	cout << "stopped by Muarraysize " << couterMuArraySize << endl;
	cout << "stopped by dR " << CounterMudR << endl;
	cout << "stopped by IP and Pt " << CounterMuIP << endl;
	cout << "stopped by noPF " << CounterMaxInd << endl;
	cout << "stopped by MuCharge " << CounterCharge << endl;
	cout << "stopped by BdIniFlavour " << CounterFlavour << endl;
	cout << "entries " << nentries << endl;
	cout << "B Events with more than 5 daughter particles" << DauCounter << endl;
}
