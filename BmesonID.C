#include "TEfficiency.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>

using namespace std;


void BmesonID(){ 

//  TFile *Bfile = TFile::Open("/tmp/terhi/BsMCCleanedNew.root");
//  TFile *Bfile = TFile::Open("/tmp/terhi/BsMC/BsMCRun2Cleaned.root");

//  TFile *Bfile = TFile::Open("/tmp/terhi/BplusMCAllNew.root");
  TFile *Bfile = TFile::Open("/tmp/terhi/BdMCPreCleaned.root");
  TTree *BsTree = (TTree*)Bfile->Get("BsTree");

  bool isBsMC = false;
  bool isBpMC = false;
  bool isBdMC = true;

  int nentries = BsTree->GetEntries();
  cout<<"Entries: " <<nentries << endl;
 
  int GenNumberOfBdecays_;
  int BmesonsId_[10];
  int BDauIdMC_[10][15];
  int BDauDauIdMC_[10][15][10];
  int GenNumberOfDaughters_[10];
  int GenNumberOfDaughtersDaughters_[10][15];

  BsTree->SetBranchAddress(  "GenNumberOfBdecays"               , &GenNumberOfBdecays_  );
  BsTree->SetBranchAddress(  "BmesonsId"                   , BmesonsId_  );
  BsTree->SetBranchAddress(  "BDauIdMC"                    , BDauIdMC_  );
  BsTree->SetBranchAddress(  "BDauDauIdMC"                         , BDauDauIdMC_  );
  BsTree->SetBranchAddress(  "GenNumberOfDaughters"        , GenNumberOfDaughters_  );
  BsTree->SetBranchAddress(  "GenNumberOfDaughtersDaughters" ,   GenNumberOfDaughtersDaughters_  );

 int NoBdecays =0;
 int TotBdecays =0; 

 int BplusOppoCounter=0; 
 int BsOppoCounter=0;
 int BdOppoCounter=0;

 int BpSigCounter=0; 
 int BsSigCounter=0;
 int BdSigCounter=0;
 int bbpaircounter =0, bbpairthreecounter =0, threebcounter =0;

 bool JpsiFound = false;
 bool KplusFound = false;
 bool KstarFound = false;
 bool PhiFound = false;
 int signalflag = 0,nosigChannel=0;
 int trueflagcount = 0;
 vector<int> BsoppoDau;  
//nentries

 int BminusCounter =0;
 int BplusCounter =0;

 for(int l =0; l <nentries; l++){	
	BsTree->GetEntry(l);

	if(GenNumberOfBdecays_ == 0) {NoBdecays++;}	

//	cout << "number of B decays in event: "<< l << ": " << GenNumberOfBdecays_ << endl;	        
//	if(GenNumberOfBdecays_ == 4){bbpaircounter++; }
//	if(GenNumberOfBdecays_ > 2){bbpaircounter++; continue; }
//	if(GenNumberOfBdecays_ > 2 ) {threebcounter++; continue;}
	if(GenNumberOfBdecays_ == 6) { bbpairthreecounter++;}

	TotBdecays = TotBdecays + GenNumberOfBdecays_;
	signalflag = 0;
	trueflagcount =0;

	
  	
	
	for(int x=0; x < GenNumberOfBdecays_; x++){
	    int B_id = abs(BmesonsId_[x]);
            int BpId = BmesonsId_[x]; 



		if(BpId == -521) {BminusCounter++; }
		if(BpId == 521) {BplusCounter++;}
	
	    //if(GenNumberOfBdecays_ > 4) {	
//	    cout << " B id " << BpId <<endl; 
	    //}	
	    JpsiFound = false;	
	    KplusFound =false;	
	    KstarFound = false;
	    PhiFound = false;

		BsoppoDau.clear();
		
		   for(int j=0; j < GenNumberOfDaughters_[x]; j++){
			int Dau_id = BDauIdMC_[x][j];
      			//int absDau_id = abs(Dau_id);
			
			if(abs(Dau_id) == 443) JpsiFound = true;
			if(Dau_id == 321) KplusFound = true;
			if(abs(Dau_id) == 313) KstarFound = true;
			if(abs(Dau_id) == 333 ) PhiFound = true;
//			cout <<"Daughter id " << Dau_id << endl;
/*			if(abs(B_id) == 511 && (KstarFound == false || JpsiFound == false) ){ 
			   BsoppoDau.push_back(Dau_id);			   	
			}*/
		} 

	   if( B_id ==531 && JpsiFound  == true && PhiFound == true && GenNumberOfDaughters_[x] == 
2 && isBsMC == true && signalflag ==0){ BsSigCounter++; signalflag = 1; trueflagcount++; //cout << "Bs channel!" << endl;
	}

	  else if( BpId == 521 &&  JpsiFound  == true && KplusFound == true && GenNumberOfDaughters_[x] == 2 && isBpMC == true && signalflag ==0){ BpSigCounter++; signalflag =1;  //cout << "B+ channel!" << endl;
	}

	   else if( B_id == 511 && JpsiFound  == true && KstarFound == true && GenNumberOfDaughters_[x] == 2 && isBdMC == true && signalflag ==0){ BdSigCounter++; signalflag = 1; //cout << "Bd channel!" << endl;
}
	else{ 
	//	cout << "oppo meson" << endl;
		if(B_id == 521) BplusOppoCounter++;
		if(B_id == 531) BsOppoCounter++;		
		if(B_id == 511) {BdOppoCounter++; 
/*			cout << "Event " << l <<endl;
			for(int k =0; k<BsoppoDau.size(); k++){
				cout << BsoppoDau[k] << endl; 
*/			
			} 

		}
	    }
	}

/*	if(signalflag == 0){
//cout << "no signal channel B decay in the event " << l << endl;
nosigChannel++;
}
	if( trueflagcount >= 2 ) {cout << "event more than one accepted signals " << l << endl; }
 }
*/
  cout <<  "total number of B decays " << TotBdecays <<endl;
  cout << "no B decays " << NoBdecays << endl;
  cout << "Decays with no B signal channel " << nosigChannel << endl; 
  cout << "decays with 2 bb pairs " << bbpaircounter << endl;
  cout << "decays with 3 bb pairs " << bbpairthreecounter << endl;
  cout << "three b quarks " << threebcounter <<endl;
  cout << "BsSigCounter " <<  BsSigCounter << endl;
  cout << "BpSigCounter " <<  BpSigCounter << endl;
  cout << "BdSigCounter " <<  BdSigCounter << endl;
  cout << "Bd opposite side mesons " << BdOppoCounter << endl;
  cout << "Bs opposite side mesons " << BsOppoCounter << endl;
  cout << "Bplus opposite side mesons " << BplusOppoCounter << endl;

  cout << "Bplus mesons in general " << BplusCounter << endl;
    cout << "Bminus mesons in general " << BminusCounter << endl;	
  cout <<"" <<endl;

   Double_t xAxis[7] = {-0.05, 0, 0.0055, 0.015, 0.035, 0.1, 0.3}; 
   
   TEfficiency *BmesonEff = new TEfficiency("BmesonEff","",6,xAxis);
     
  //efficiency calculation for signal mesons

   BmesonEff->SetTotalEvents(0,nentries); 
   BmesonEff->SetPassedEvents(0,BsSigCounter);   
   
   double dBs = (BmesonEff->GetEfficiencyErrorLow(0) + BmesonEff->GetEfficiencyErrorUp(0) )/2.0;

   BmesonEff->SetTotalEvents(1,nentries );  
   BmesonEff->SetPassedEvents(1,BpSigCounter);   
   
   double dBp = (BmesonEff->GetEfficiencyErrorLow(1) + BmesonEff->GetEfficiencyErrorUp(1) )/2.0;

   BmesonEff->SetTotalEvents(2,nentries );
   BmesonEff->SetPassedEvents(2,BdSigCounter);   
  
   double dBd = (BmesonEff->GetEfficiencyErrorLow(2) + BmesonEff->GetEfficiencyErrorUp(2) )/2.0;

   //efficiency calculation for opposite side mesons
  
	Int_t signalMesonEvents;
 	
	if( isBdMC == true){signalMesonEvents =BdSigCounter; }
	if( isBsMC == true){signalMesonEvents =BsSigCounter; } 
	if( isBpMC == true){signalMesonEvents =BpSigCounter; } 

   BmesonEff->SetTotalEvents(3,TotBdecays-signalMesonEvents); 
   BmesonEff->SetPassedEvents(3,BsOppoCounter);   
   
   double dOppoBs = (BmesonEff->GetEfficiencyErrorLow(3) + BmesonEff->GetEfficiencyErrorUp(3) )/2.0;

   BmesonEff->SetTotalEvents(4,TotBdecays-signalMesonEvents );  
   BmesonEff->SetPassedEvents(4,BplusOppoCounter);   
   
   double dOppoBp = (BmesonEff->GetEfficiencyErrorLow(4) + BmesonEff->GetEfficiencyErrorUp(4) )/2.0;

   BmesonEff->SetTotalEvents(5,TotBdecays-signalMesonEvents );
   BmesonEff->SetPassedEvents(5,BdOppoCounter);   
  
   double dOppoBd = (BmesonEff->GetEfficiencyErrorLow(5) + BmesonEff->GetEfficiencyErrorUp(5) )/2.0;
   
//  cout <<  "total number of B decays " << TotBdecays <<endl;
  
  cout << "Bs signal mesons from all events " << BmesonEff->GetEfficiency(0) << "  +- " << dBs  << endl;
  cout << "Bplus signal mesons from all events " << BmesonEff->GetEfficiency(1) << "  +- " << dBp << endl;
  cout << "Bd signal mesons from all events " << BmesonEff->GetEfficiency(2)  << "  +- " << dBd << endl;

  cout << "" << endl;

  cout << "Bs opposite side mesons from all OST particles " << BmesonEff->GetEfficiency(3) << "  +- " << dOppoBs  << endl;
  cout << "Bplus opposite side mesons from all OST particles " << BmesonEff->GetEfficiency(4) << "  +- " << dOppoBp << endl;
  cout << "Bd opposite side mesons from all OST particles " << BmesonEff->GetEfficiency(5)  << "  +- " << dOppoBd << endl;

}
