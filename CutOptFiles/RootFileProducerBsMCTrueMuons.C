#include "TMath.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TROOT.h"
#include "TStyle.h"

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
	
	if(maxIndex > 24){cout <<"max index > 24!" << endl;}	

	if(PFMuExists == true){		
	return maxIndex;                // return index of max value  in array
	}
	
	else return -1;

}

Int_t isOppoBsMuon(Int_t GenNumberOfBdecays_, Int_t GenNumberOfDaughters_[],Int_t BDauIdMC_[10][15], Int_t BmesonsId_[]){


 bool JpsiFound = false;
 bool PhiFound = false;

 int signalflag = 0;
 int BsSigCounter =0;
 int BsOppoCounter=0;


//	cout << " inside oppo method 1 " << endl; 

	for(int x=0; x < GenNumberOfBdecays_; x++){
	    int B_id = abs(BmesonsId_[x]);
            int BpId = BmesonsId_[x]; 

//	    cout << " B id " <<BpId << endl; 	
//	    if(GenNumberOfBdecays_ >2) cout << "More than two b quarks!" << endl; 
			
	    JpsiFound = false;	
	    PhiFound =false;	
	   		
		for(int j=0; j < GenNumberOfDaughters_[x]; j++){
			int Dau_id = BDauIdMC_[x][j];
      			//int absDau_id = abs(Dau_id);
			
			if(Dau_id == 443) JpsiFound = true;
			if(Dau_id == 333) PhiFound = true;
		} 

	  if( B_id == 531 &&  JpsiFound  == true && PhiFound == true && GenNumberOfDaughters_[x] == 2 && signalflag ==0){ 
	BsSigCounter++; 
	signalflag =1;  //cout << "Bs channel!" << endl;

	  }
	
	 else if(B_id == 531){
	   BsOppoCounter++;	   
	  }
			
	}

//	cout << " BsOppoCounter " <<BsOppoCounter << endl; 

	if(BsOppoCounter!=0){return 1;}
	else return 0;
}

std::vector<float> JetCharges( vector<float> *PtJetTrkPt_, vector<float> *BJetTrkPt_, vector<int> *PtJetTrkCharge_,vector<int> *BJetTrkCharge_ ){

	float PtSum_PtJet=0;
 	float ChargePtSum_PtJet=0; 

 	float PtSum_BJet=0;
 	float ChargePtSum_BJet=0; 

	int BJetVecSize = BJetTrkPt_->size();
	int PtJetVecSize= PtJetTrkPt_->size();

	if(BJetVecSize > PtJetVecSize){

	  for(int k =0; k<BJetTrkPt_->size(); k++){
		PtSum_BJet = PtSum_BJet + BJetTrkPt_->at(k);
		ChargePtSum_BJet = ChargePtSum_BJet + BJetTrkPt_->at(k)*BJetTrkCharge_->at(k);
		
	
		if(k<PtJetVecSize){
			PtSum_PtJet = PtSum_PtJet + PtJetTrkPt_->at(k);
			ChargePtSum_PtJet = ChargePtSum_PtJet + PtJetTrkPt_->at(k)*PtJetTrkCharge_->at(k);
			
		}
	  }
	}

	else{ 

//		cout << "BJetVecSize " << BJetVecSize << " PtJetVecSize " << PtJetVecSize << endl;
//		myBsTree->GetEntry(810);
//		cout << "BsJetPx " << BsJetPx << " BsJetPy " << BsJetPy << endl; 
//		cout << "Jet Pt " << sqrt(BsJetPx*BsJetPx + BsJetPy*BsJetPy) << endl;

	  for(int k =0; k<PtJetTrkPt_->size(); k++){
		PtSum_PtJet = PtSum_PtJet + PtJetTrkPt_->at(k);
		ChargePtSum_PtJet = ChargePtSum_PtJet + PtJetTrkPt_->at(k)*PtJetTrkCharge_->at(k);
/*		cout << "track " << k << endl; 
		cout << "pt " << PtJetTrkPt_->at(k) << endl;
		cout << "charge " << PtJetTrkCharge_->at(k) << endl;
		cout << " ChargePtSum_PtJet " << ChargePtSum_PtJet <<endl;
		cout << " PtSum_PtJet " <<PtSum_PtJet  <<endl;
*/
		
		if(k<BJetVecSize){
			PtSum_BJet = PtSum_BJet + BJetTrkPt_->at(k);
			ChargePtSum_BJet = ChargePtSum_BJet + BJetTrkPt_->at(k)*BJetTrkCharge_->at(k);

/*		cout << "pt " << BJetTrkPt_->at(k) << endl;
		cout << "charge " << BJetTrkCharge_->at(k) << endl;
		cout << " ChargePtSum_BJet " << ChargePtSum_BJet<<endl;
		cout << " PtSum_BJet " <<PtSum_BJet  <<endl;
*/
		}
	  }
//	if(ChargePtSum_PtJet/PtSum_PtJet == 1.0){ cout << "ChargePtSum_PtJet/PtSum_PtJet " << ChargePtSum_PtJet/PtSum_PtJet << " " << endl; }		
      }

	vector<float> JetCharges;
	JetCharges.push_back(ChargePtSum_BJet/PtSum_BJet);
	JetCharges.push_back(ChargePtSum_PtJet/PtSum_PtJet);

	return JetCharges;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
	

	Double_t deltaEta = eta1-eta2;
	Double_t deltaPhi = phi1-phi2;

 	while (deltaPhi > TMath::Pi()) deltaPhi -= TMath::TwoPi();
	while (deltaPhi <= -TMath::Pi()) deltaPhi += TMath::TwoPi();

	Double_t DeltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi );

	return DeltaR;
}


void RootFileProducerBsMCTrueMuons(){

//TFile *MCPUFile = TFile::Open("/tmp/terhi/BsMCCleanedNew.root"); 
TFile *MCPUFile = TFile::Open("/tmp/terhi/BsMC/BsMCRun2Cleaned.root");
TTree* myBsTree = (TTree*)MCPUFile->Get("BsTree");
//myBsTree->Print(); 

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

 Int_t MuArraySize;
 Double_t BsEta,BsFitPt;
 Double_t BsPhi, BsFitM,BsPtJetdR,BsBJetdR;
 Int_t BsIniFlavour; 
 Int_t nentries = myBsTree->GetEntries();
 
 std::vector<float> *PtJetTrkPt= new std::vector<float>(); 
 std::vector<int> *PtJetTrkCharge= new std::vector<int>();
 std::vector<float> *BJetTrkPt= new std::vector<float>(); 
 std::vector<int> *BJetTrkCharge= new std::vector<int>();


 myBsTree->SetBranchAddress("TagMuRecoPt", TagMuPt);
 myBsTree->SetBranchAddress("TagMuRecoIP3D", MuIP3D); 
 myBsTree->SetBranchAddress("TagMuListSize", &MuArraySize);
 myBsTree->SetBranchAddress("TagMuRecoChg", TagMuCharge);
 myBsTree->SetBranchAddress("GlobalTagMuon", isGlobalTagMuon);
 myBsTree->SetBranchAddress("TrackerTagMuon", isTrackerTagMuon);
 myBsTree->SetBranchAddress("PFTagMuon", isPFTagMuon);
 myBsTree->SetBranchAddress("TagMuRecoEta",MuEta);
 myBsTree->SetBranchAddress("TagMuRecoPhi",MuPhi);
 myBsTree->SetBranchAddress("BsIniFlavour", &BsIniFlavour);
 myBsTree->SetBranchAddress("BsFitEta",&BsEta);
 myBsTree->SetBranchAddress("BsFitPhi",&BsPhi);
 myBsTree->SetBranchAddress("BsFitM",&BsFitM);
 myBsTree->SetBranchAddress("BsFitPt",&BsFitPt);
 myBsTree->SetBranchAddress("TagMuSimu",MuSimu);
 myBsTree->SetBranchAddress("MuRecoMCmother",MuSimuMom);
 myBsTree->SetBranchAddress("MuRecoMCMotherMother",MuSimuGMom);
/*
 myBsTree->SetBranchAddress("TagMuRecoPtRel",TagMuPtRel);

 myBsTree->SetBranchAddress("PtJetTrkPt", &PtJetTrkPt);
 myBsTree->SetBranchAddress("PtJetTrkCharge", &PtJetTrkCharge);

 myBsTree->SetBranchAddress("BJetTrkPt", &BJetTrkPt);
 myBsTree->SetBranchAddress("BJetTrkCharge", &BJetTrkCharge);

 myBsTree->SetBranchAddress("BsPtJetdR",&BsPtJetdR);
 myBsTree->SetBranchAddress("BsBJetdR",&BsBJetdR);

 Double_t BsJetPx, BsJetPy;
// myBsTree->SetBranchAddress("BsJetPx", &BsJetPx);
// myBsTree->SetBranchAddress("BsJetPy", &BsJetPy); 
*/

 Int_t BDauIdMC_[10][15];
 Int_t BmesonsId_[10];
 Int_t GenNumberOfBdecays_;
 Int_t GenNumberOfDaughters_[10];
 myBsTree->SetBranchAddress("BDauIdMC", BDauIdMC_ );
 myBsTree->SetBranchAddress("GenNumberOfBdecays", &GenNumberOfBdecays_ );
 myBsTree->SetBranchAddress("BmesonsId", BmesonsId_  );
 myBsTree->SetBranchAddress("GenNumberOfDaughters", GenNumberOfDaughters_  );

	Double_t DeltaR, MuPt, MuIP,SimuMuon,MuPtRel,BJetCharge,PtJetCharge,BsPtJetPt;
	Int_t MuCharge,BsFlavour,PFMuon,TrackerMuon,GlobalMuon,MCMuonGMom,MCMuonMom,OppoBsMuon,OppoBsMeson;

		//CutOptFileMCBsRun2.root
 	TFile f("","RECREATE");
        TTree tree("cutTree","cutTree");
        tree.Branch("DeltaR",&DeltaR,"DeltaR/D");  
	tree.Branch("MuPt",&MuPt,"MuPt/D");  
	tree.Branch("MuIP",&MuIP,"MuIP/D");  
	tree.Branch("MuCharge",&MuCharge,"MuCharge/I");
	tree.Branch("BsFlavour",&BsIniFlavour,"BsFlavour/I");
	tree.Branch("PFMuon", &PFMuon,"PFMuon/I");
//	tree.Branch("GlobalMuon", &GlobalMuon,"GlobalMuon/I");
//	tree.Branch("TrackerMuon", &TrackerMuon,"TrackerMuon/I");
	tree.Branch("BpMass", &BsFitM,"BpMass/D");
        tree.Branch("SimuMuon", &SimuMuon,"SimuMuon/D");
	tree.Branch("BsPt", &BsFitPt,"BsPt/D");
	tree.Branch("MCMuonMom", &MCMuonMom,"MCMuonMom/I");
	tree.Branch("MCMuonGMom", &MCMuonGMom,"MCMuonGMom/I");
//	tree.Branch("MuPtRel",&MuPtRel,"MuPtRel/D");
	
	tree.Branch("OppoBsMuon",&OppoBsMuon,"OppoBsMuon/I");
	tree.Branch("OppoBsMeson",&OppoBsMeson,"OppoBsMeson/I"); 
/*
	tree.Branch("BJetCharge",&BJetCharge,"BJetCharge/D"); 
	tree.Branch("PtJetCharge",&PtJetCharge,"PtJetCharge/D"); 
 

	tree.Branch("PtJetTrkCharge","vector<int>",PtJetTrkCharge);
	tree.Branch("PtJetTrkPt","vector<float>",PtJetTrkPt);
	tree.Branch("BJetTrkCharge","vector<int>",BJetTrkCharge);
	tree.Branch("BJetTrkPt","vector<float>",BJetTrkPt);

	tree.Branch("BsPtJetPt", &BsPtJetPt,"BsPtJetPt/D");
	tree.Branch("BsPtJetdR", &BsPtJetdR,"BsPtJetdR/D");
	tree.Branch("BsBJetdR", &BsBJetdR,"BsBJetdR/D");
*/
	Int_t couterMuArraySize=0, CounterCharge=0, CounterMudR =0, CounterMuIP=0, CounterMaxInd =0, CounterFlavour=0;

			//nentries
	for(Int_t i = 0; i<nentries; i++){
		
		vector<float> myJetCharge;
		myBsTree->GetEntry(i);

		DeltaR= -999;
		MuPt= -999;
		PFMuon = 0;
		MuCharge = 0;
		MuIP = -999;
		SimuMuon = -999;	
		MCMuonGMom = -999;
		MCMuonMom = -999;
		MuPtRel = -999;
		BsPtJetPt = -999;
		OppoBsMeson = 0;
		OppoBsMuon =0;

		OppoBsMeson = isOppoBsMuon(GenNumberOfBdecays_, GenNumberOfDaughters_,BDauIdMC_,BmesonsId_);

/*		
		BsPtJetPt = TMath::Sqrt(BsJetPx*BsJetPx + BsJetPy*BsJetPy);
		myJetCharge = JetCharges(PtJetTrkPt,BJetTrkPt,PtJetTrkCharge,BJetTrkCharge);
		BJetCharge = myJetCharge[0];
		PtJetCharge = myJetCharge[1];
*/
		if(MuArraySize <= 0){
			couterMuArraySize++;
			tree.Fill(); 			
			continue;
		} 
		
		if(MuArraySize > 25) cout <<"Mu array size exceeded! " <<i << endl; 			
		Int_t maxPtMuonIndex = maximumIndex(TagMuPt,isPFTagMuon,MuArraySize);
		if(maxPtMuonIndex == -1) {
			CounterMaxInd++;
			tree.Fill();  
			continue;
		}
 		
		SimuMuon =  MuSimu[maxPtMuonIndex];
		MCMuonMom = MuSimuMom[maxPtMuonIndex];
		MCMuonGMom = MuSimuGMom[maxPtMuonIndex];	
		MuPtRel = TagMuPtRel[maxPtMuonIndex];	

		if(BsIniFlavour!=-1 && BsIniFlavour!=1){
			MuCharge = 0;
			CounterFlavour++;
			tree.Fill();  
			continue;
		}

		MuCharge = TagMuCharge[maxPtMuonIndex];
		if(MuCharge!=-1 && MuCharge!=1){
			MuCharge = 0;
			CounterCharge++;
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

		if(BsEta == -999 || BsPhi == -999) {
			CounterMudR++;
			MuCharge =0;
			tree.Fill();  
			continue;
		}
//		cout << "simu muon" << SimuMuon << endl;
		DeltaR = deltaR(BsEta,BsPhi,MuEta[maxPtMuonIndex],MuPhi[maxPtMuonIndex]);
//		TrackerMuon = isTrackerTagMuon[maxPtMuonIndex];
//		GlobalMuon = isGlobalTagMuon[maxPtMuonIndex]; 
		PFMuon = isPFTagMuon[maxPtMuonIndex];
		MuIP = MuIP3D[maxPtMuonIndex];		
		MuPt = TagMuPt[maxPtMuonIndex];	
		PFMuon = isPFTagMuon[maxPtMuonIndex];
		OppoBsMuon = isOppoBsMuon(GenNumberOfBdecays_, GenNumberOfDaughters_,BDauIdMC_,BmesonsId_);	

		tree.Fill();
		
		
	}
	tree.Print();
	tree.Write(); 

	cout << "stopped by Muarraysize " << couterMuArraySize << endl;
	cout << "stopped by dR " << CounterMudR << endl;
	cout << "stopped by IP and Pt " << CounterMuIP << endl;
	cout << "stopped by noPF " << CounterMaxInd << endl;
	cout << "stopped by MuCharge " << CounterCharge << endl;
	cout << "stopped by MC flavour " << CounterFlavour << endl;
	cout << "entries " << nentries << endl;

}
