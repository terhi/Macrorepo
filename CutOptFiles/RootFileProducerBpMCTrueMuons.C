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

	
	for(Int_t k = 0; k<PFPt.size(); k++){
		
		if(PFPt[maxIndex] <= PFPt[k]){ 
			maxIndex = PFMuIndices[k];
	  	}
	}	
	
	if(maxIndex > 24){cout <<"max index > 24!" << endl;}	

	if(PFMuExists == true){		
	return maxIndex;                // return index of max value  in array
	}
	
	else return -1;

}

Int_t isOppoBpMuon(Int_t GenNumberOfBdecays_, Int_t GenNumberOfDaughters_[],Int_t BDauIdMC_[10][15], Int_t BmesonsId_[]){


 bool JpsiFound = false;
 bool KplusFound = false;

 int signalflag = 0;
 int BpSigCounter =0;
 int BplusOppoCounter=0;


//	cout << " inside oppo method 1 " << endl; 

	for(int x=0; x < GenNumberOfBdecays_; x++){
	    int B_id = abs(BmesonsId_[x]);
            int BpId = BmesonsId_[x]; 

//	if(GenNumberOfDaughters_[x] > 5) { DauCounter++; cout << " inside oppo method 1 " << endl;  cout <<  BpId << endl;}
	    JpsiFound = false;	
	    KplusFound =false;	
	   		
		for(int j=0; j < GenNumberOfDaughters_[x]; j++){
			int Dau_id = BDauIdMC_[x][j];
      			//int absDau_id = abs(Dau_id);
			
//			if(GenNumberOfDaughters_[x] > 5) {cout <<  Dau_id << endl;}

			if(Dau_id == 443) JpsiFound = true;
			if(Dau_id == 321) KplusFound = true;
		} 

	  if( BpId == 521 && JpsiFound  == true && KplusFound == true && GenNumberOfDaughters_[x] == 2 && signalflag ==0){ 
	BpSigCounter++; signalflag =1; // cout << "B+ channel!" << endl;
	  }
	
	 else if(B_id == 521){ 
	   BplusOppoCounter++;	   
	  }			
	}

//	cout << " BplusOppoCounter " <<BplusOppoCounter << endl; 

	if(BplusOppoCounter!=0){return 1;}
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


void RootFileProducerBpMCTrueMuons(){
 
gROOT->SetStyle("Plain");
//gStyle->SetFillColor(0);  //white fill color
gStyle->SetFrameBorderMode(0);  //no frame border
gStyle->SetCanvasBorderMode(0);  //no canvas border
gStyle->SetOptTitle(0);  //no title
gStyle->SetOptStat(1111111);  //no statistics _or_
gStyle->SetCanvasColor(0);
gStyle->SetStatBorderSize(0);

//TFile *MCPUFile = TFile::Open("/tmp/terhi/JetVectorTestCleaned.root");
TFile *MCPUFile = TFile::Open("/tmp/terhi/BplusMCCleanedNew.root");

TTree* myBsTree = (TTree*)MCPUFile->Get("BsTree");

 const int Arraysize = 25;
 Double_t TagMuPt[Arraysize];
 Double_t MuIP3D[Arraysize];
 Double_t TagMuCharge[Arraysize];
 Int_t isGlobalTagMuon[Arraysize];
 Int_t isTrackerTagMuon[Arraysize];
 Int_t isPFTagMuon[Arraysize];
 Int_t MuSimuMom[Arraysize];
 Int_t MuSimuGMom[Arraysize];
 Double_t MuPhi[Arraysize];
 Double_t MuEta[Arraysize]; 
 Double_t MuSimu[Arraysize];
 Double_t TagMuPtRel[Arraysize];

 vector<int> *JetTrackCharge = new vector<int>();
 Int_t MuArraySize;
 Double_t BplusEta;
 Double_t BplusPhi, BplusMass;
 Int_t nentries = myBsTree->GetEntries();
 
 myBsTree->SetBranchAddress("TagMuRecoPt", TagMuPt);
 myBsTree->SetBranchAddress("TagMuRecoIP3D", MuIP3D); 
 myBsTree->SetBranchAddress("TagMuListSize", &MuArraySize);
 myBsTree->SetBranchAddress("TagMuRecoChg", TagMuCharge);
 myBsTree->SetBranchAddress("GlobalTagMuon", isGlobalTagMuon);
 myBsTree->SetBranchAddress("TrackerTagMuon", isTrackerTagMuon);
 myBsTree->SetBranchAddress("PFTagMuon", isPFTagMuon);
 myBsTree->SetBranchAddress("TagMuRecoEta",MuEta);
 myBsTree->SetBranchAddress("TagMuRecoPhi",MuPhi);
 myBsTree->SetBranchAddress("BplusEta",&BplusEta);
 myBsTree->SetBranchAddress("BplusPhi",&BplusPhi);
 myBsTree->SetBranchAddress("BplusM_fit",&BplusMass);
 myBsTree->SetBranchAddress("TagMuSimu",MuSimu);
 myBsTree->SetBranchAddress("MuRecoMCmother",MuSimuMom);
 myBsTree->SetBranchAddress("MuRecoMCMotherMother",MuSimuGMom);
// myBsTree->SetBranchAddress("JetTrkCharge", &JetTrackCharge);
//myBsTree->SetBranchAddress("TagMuRecoPtRel",TagMuPtRel);

 Int_t BDauIdMC_[10][15];
 Int_t BmesonsId_[10];
 Int_t GenNumberOfBdecays_;
 Int_t GenNumberOfDaughters_[10];
 myBsTree->SetBranchAddress("BDauIdMC", BDauIdMC_ );
 myBsTree->SetBranchAddress("GenNumberOfBdecays", &GenNumberOfBdecays_ );
 myBsTree->SetBranchAddress("BmesonsId", BmesonsId_  );
 myBsTree->SetBranchAddress("GenNumberOfDaughters", GenNumberOfDaughters_  );



	Double_t DeltaR, MuPt, MuIP,SimuMuon,MuPtRel;
	Int_t MuCharge,PFMuon,MCMuonGMom,MCMuonMom,OppoBpMuon,OppoBpMeson;

 	TFile f("MuMcOppo.root","RECREATE");
        TTree tree("cutTree","cutTree");
        tree.Branch("DeltaR",&DeltaR,"DeltaR/D");  
	tree.Branch("MuPt",&MuPt,"MuPt/D");  
	tree.Branch("MuIP",&MuIP,"MuIP/D");  
	tree.Branch("MuCharge",&MuCharge,"MuCharge/I");
	tree.Branch("BpMass", &BplusMass,"BpMass/D");
	tree.Branch("PFMuon", &PFMuon,"PFMuon/I");
	tree.Branch("SimuMuon", &SimuMuon,"SimuMuon/D");
	tree.Branch("MCMuonMom", &MCMuonMom,"MCMuonMom/I");
	tree.Branch("MCMuonGMom", &MCMuonGMom,"MCMuonGMom/I");
	tree.Branch("MuPtRel",&MuPtRel,"MuPtRel/D");
	tree.Branch("OppoBpMuon",&OppoBpMuon,"OppoBpMuon/I");
	tree.Branch("OppoBpMeson",&OppoBpMeson,"OppoBpMeson/I");  
//	tree.Branch("JetTrackCharge", JetTrackCharge);

	Int_t couterMuArraySize=0, CounterCharge=0, CounterMudR =0, CounterMuIP=0, CounterMaxInd =0,RtagBdcounter=0;

			//nentries
	for(Int_t i = 0; i<nentries; i++){
		
		DeltaR= -999;
		MuPt= -999;
		PFMuon = 0;
		MuCharge = 0;
		MuIP = -999;
		SimuMuon = -999;
		MCMuonGMom = -999;
		MCMuonMom = -999;
		MuPtRel = -999;
		OppoBpMuon = 0;
		OppoBpMeson = 0;

		myBsTree->GetEntry(i);

		OppoBpMeson = isOppoBpMuon(GenNumberOfBdecays_, GenNumberOfDaughters_,BDauIdMC_,BmesonsId_);
 
							
		if(MuArraySize <= 0) {
			couterMuArraySize++;
			tree.Fill(); 			
			continue;
		} 
		
		if(MuArraySize > 25) cout <<"Mu array size exceeded! " <<i << endl; 			
		Int_t maxPtMuonIndex = maximumIndex(TagMuPt,isPFTagMuon,MuArraySize);

		

		if(maxPtMuonIndex == -1){
			CounterMaxInd++;
			tree.Fill();  
			continue;
		}

		SimuMuon = MuSimu[maxPtMuonIndex];
		MCMuonMom = MuSimuMom[maxPtMuonIndex];
		MCMuonGMom = MuSimuGMom[maxPtMuonIndex];	
//		MuPtRel = TagMuPtRel[maxPtMuonIndex];	

		MuCharge = TagMuCharge[maxPtMuonIndex];
		if(MuCharge!=-1 && MuCharge!=1){
			MuCharge = 0;
			CounterCharge++;
			tree.Fill();  
			continue;
		}		
	
		 
		if(MuIP3D[maxPtMuonIndex] < 0  || TagMuPt[maxPtMuonIndex] < 0){ 
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


		if(BplusEta == -999 || BplusPhi == -999){
			CounterMudR++;
			MuCharge =0;
			tree.Fill();  
			continue;
		}

		DeltaR = deltaR(BplusEta,BplusPhi,MuEta[maxPtMuonIndex],MuPhi[maxPtMuonIndex]);
		if(DeltaR<0){
			CounterMudR++;
			MuCharge =0;
			tree.Fill();  
			continue;
		}
	
		PFMuon = isPFTagMuon[maxPtMuonIndex];

		if(PFMuon !=1){
			CounterMaxInd++;
			MuCharge =0;
			tree.Fill();  
			continue;
		}

		MuIP = MuIP3D[maxPtMuonIndex];		
		MuPt = TagMuPt[maxPtMuonIndex];	
		PFMuon = isPFTagMuon[maxPtMuonIndex];

//		cout << "before oppo muon " << endl;
		OppoBpMuon = isOppoBpMuon(GenNumberOfBdecays_, GenNumberOfDaughters_,BDauIdMC_,BmesonsId_);	
		
//		if(TMath::Abs(MCMuonMom) == 511 && MuCharge == -1){ RtagBdcounter++;}
//		cout << "jet tracks " << JetTrackCharge->size() << " " << i<<endl;
//		cout << "SimuMuon " << SimuMuon << endl;
//		cout << "SimuMuonmom " << MCMuonMom << endl;

		tree.Fill();
//		cout << "AFTER filling tree  " << endl;				
	}
	tree.Print();
	tree.Write(); 

	
	cout << "stopped by Muarraysize " << couterMuArraySize << endl;
	cout << "stopped by dR " << CounterMudR << endl;
	cout << "stopped by IP and Pt " << CounterMuIP << endl;
	cout << "stopped by noPF " << CounterMaxInd << endl;
	cout << "stopped by MuCharge " << CounterCharge << endl;
	cout << "entries " << nentries << endl;
	cout << "rtag Bd mom counter " << RtagBdcounter << endl;

        cout << "B Events with more than 5 daughter particles" << DauCounter << endl;
}
