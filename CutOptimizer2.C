#include "TMath.h"
#include <vector>
#include <cmath>
int maximumIndex(vector<double> array, int size){
 
     int maxIndex = 0;       // start with max = first element

     for(int i = 0; i<size; i++)
     {
          if(array[i] > array[maxIndex]){ 
		maxIndex = i;
	 }	
     }	
     return maxIndex;                // return index of max value  in array
}

void CutOptimizer2(){
gROOT->SetStyle("Plain");
//gStyle->SetFillColor(0);  //white fill color
gStyle->SetFrameBorderMode(0);  //no frame border
gStyle->SetCanvasBorderMode(0);  //no canvas border
gStyle->SetOptTitle(0);  //no title
gStyle->SetOptStat(1111111);  //no statistics _or_
gStyle->SetCanvasColor(0);
gStyle->SetStatBorderSize(0);

//TFile *MCPUFile = TFile::Open("/tmp/terhi/BsMCCleanedNew.root");
//TFile *MCPUFile = TFile::Open("/tmp/terhi/BplusMCFiles/BplusMCCleaned2.root");
//TTree* myBsTree = (TTree*)MCPUFile->Get("BsTree");


 
 vector<double> WrongTagFrac; 
 vector<double> TagEfficiency; 
 vector<double> TagPower; 
 vector<TString> Cutvec; 

   //TFile* myCutFile = new TFile("CutOptFileMCBsNewTrig2.root"); (old one)
  TFile* myCutFile = TFile::Open("CutOptFiles/CutOptFileBsMCFinal2.root");

   TTree* mycutTree = (TTree*)myCutFile->Get("cutTree");		
	
    Int_t nentries = mycutTree->GetEntries("BpMass > 5.265 && BpMass < 5.475");
   TEventList *DR = new TEventList("DR"); 
 //  mycutTree->Draw(">>DR"); 
  // cout << nentries << endl;
/*     TCanvas *c = new TCanvas("c","c");
     c->Divide(2,2);
     c->cd(1); 
     TH1F* delta = new TH1F("delta","delta",4,-2,2);	
     mycutTree->Draw("BsFlavour>>delta(4,-2,2)"); 

   
*/
  vector<double> Wrong; 
  vector<double> WandR;

/*
 //  TString Wtagcut = "&& PFMuon == 1 && ((BsFlavour == 1 && MuCharge == 1) || (BsFlavour == -1 && MuCharge == -1) )";
//   TString Rtagcut = "&& PFMuon == 1 && ((BsFlavour == -1 && MuCharge == 1) || (BsFlavour == 1 && MuCharge == -1) )"; 	

 TString Wtagcut = "&& (( BsFlavour == 1 && MuCharge == 1 ) || ( BsFlavour == -1 && MuCharge == -1 )  )";
   TString Rtagcut = "DeltaR>0.300000 && MuIP<0.100000 && MuPt>1.900000 && PFMuon==1 && ((BsFlavour == -1 && MuCharge == 1) || (BsFlavour == 1 && MuCharge == -1) )"; 
	
	//cout << Wtagcut << endl;
	//cout << Rtagcut << endl;

// mycutTree->Draw(">>DR",Wtagcut);
 // int Wtags = DR->GetN();

 mycutTree->Draw(">>DR",Rtagcut);
  int Rtags = DR->GetN();

  mycutTree->Draw(">>DR","DeltaR>0.300000 && MuIP<0.100000 && MuPt>1.900000 && PFMuon==1 && MuCharge !=0");
  int Ttags = DR->GetN();

	double dR = 0.3; double IP = 0.1; double pt = 1.9;

	TString cut = Form("DeltaR>%f && MuIP<%f && MuPt>%f ",dR,IP,pt); 
	 TString cutWtag = cut + Wtagcut;	
	cout << " My wtag cut: "<< cutWtag << endl;
*/  

 TString Wtagcut = "&& TMath::Abs(SimuMuon) == 13 && BpMass > 5.265 && BpMass < 5.475 && PFMuon == 1 && ((BsFlavour == 1 && MuCharge == 1) || (BsFlavour == -1 && MuCharge == -1) )";

 TString Rtagcut = "&& BpMass > 5.265 && BpMass < 5.475 && PFMuon == 1 && ((BsFlavour == -1 && MuCharge == 1) || (BsFlavour == 1 && MuCharge == -1) )"; 	


   TH1F *WTagF = new TH1F("WTagF","WTagF", 100,0.2,0.7); 
   TH1F *Efficiency = new TH1F("Efficiency","Efficiency",100,-0.02,0.2);
   TH1F *TagginPower = new TH1F("TagginPower","TagginPower",100,-0.01,0.01);
      
   TH3F *Scatter = new TH3F("Scatter","scatter",5,0,0.5,50,0,5,100,0,0.01);

   double WTagFrac, effi, tagpow;
   double tagpowBest=-1; 
   double dR =0;
   int counter =0; 	
 //CHECK THE LIMITS!
 
  for(double IP=0; IP<0.5; IP = IP+0.1){	// IP<0.3
    for(double dR=0; dR<5; dR = dR+0.1){     // dR <0.1 
	for(double pt=0; pt<5; pt = pt+0.1){  // Pt < 2

   TString cut = Form("DeltaR>%f && MuIP<%f && MuPt>%f ",dR,IP,pt); 
  	
   TString cutWtag = cut + Wtagcut;
   TString cutRtag = cut + Rtagcut; 



    mycutTree->Draw(">>DR",cutWtag);
    int Wtags = DR->GetN();
	
 //  cout << cutWtag << endl;
 //  cout << cutRtag << endl;
	

    mycutTree->Draw(">>DR",cutRtag);
    int Rtags = DR->GetN();
	
    int Totaltags = Rtags + Wtags;
    if(Totaltags != 0){	
	WTagFrac = double(Wtags)/double(Totaltags);	
    }
    else WTagFrac = 0.5;
		    
    effi = double(Totaltags)/double(nentries);
    tagpow = effi*(1.0-2.0*WTagFrac)*(1.0-2.0*WTagFrac);
     
    Wrong.push_back(Wtags);
    WandR.push_back(Totaltags);

    if(tagpow > tagpowBest){
	tagpowBest = tagpow;
	TString BestCuts = cut;
    } 		

    WTagF->Fill(WTagFrac);    
    Efficiency->Fill(effi);    
    TagginPower->Fill(tagpow); 

//	Scatter->Fill(IP,pt,tagpow);

    WrongTagFrac.push_back(WTagFrac); 
    TagEfficiency.push_back(effi); 
    TagPower.push_back(tagpow);

	counter++;				
	}		
     }		
  }
	cout << "count " << counter << endl;
		
 /*
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,2);
  c1->cd(1); 
  WTagF->Draw();
  c1->cd(2);
  Efficiency->Draw();
  c1->cd(3); 
  TagginPower->Draw();
  c1->cd(4); 	 
  Scatter->Draw("surf3");	
*/
  int maxIndex = maximumIndex(TagPower,TagPower.size());
  cout << "max tagging power " << TagPower[maxIndex] << endl;
  cout << "efficiency " << TagEfficiency[maxIndex] << endl; 
  cout << "wrong tag frac " << WrongTagFrac[maxIndex] << endl; 
  cout << "Another best cut " << BestCuts << endl;		 
  cout << "number of wrongs "<< Wrong[maxIndex] << endl;
  cout << "number of totals "<< WandR[maxIndex] << endl;
  cout << "entries "<< nentries << endl;

   TEfficiency* eff_err = new TEfficiency("eff","eff",10,0,2);
   
    eff_err->SetTotalEvents(1,nentries);  
    eff_err->SetPassedEvents(1, WandR[maxIndex] ); // WandR[maxIndex]
    
   cout <<"eff error low " <<eff_err->GetEfficiencyErrorLow(1) << endl;
   cout <<"eff error up " <<eff_err->GetEfficiencyErrorUp(1) << endl;
   cout <<"eff  " <<eff_err->GetEfficiency(1) << endl;

   TEfficiency* wrongtag_eff= new TEfficiency("w_eff","w_eff",10,0,1);
   wrongtag_eff->SetTotalEvents(1, WandR[maxIndex]  ); //   
   wrongtag_eff->SetPassedEvents(1,Wrong[maxIndex]); //Wrong[maxIndex]
   
   cout <<"wfrac error low " <<wrongtag_eff->GetEfficiencyErrorLow(1) << endl;
   cout <<"wfrac error up " <<wrongtag_eff->GetEfficiencyErrorUp(1) << endl;
   cout <<"wfrac " <<wrongtag_eff->GetEfficiency(1) << endl;
   
   cout << "max tagging power " << eff_err->GetEfficiency(1)*(1-2*wrongtag_eff->GetEfficiency(1))*(1-2*wrongtag_eff->GetEfficiency(1)) << endl;  

   double tagpowerBs = eff_err->GetEfficiency(1)*(1-2*wrongtag_eff->GetEfficiency(1))*(1-2*wrongtag_eff->GetEfficiency(1));
   double deff = (eff_err->GetEfficiencyErrorLow(1) + eff_err->GetEfficiencyErrorUp(1))/2.0; 
   double dwtag = (wrongtag_eff->GetEfficiencyErrorLow(1) + wrongtag_eff->GetEfficiencyErrorUp(1) )/2.0;

   double term1 = pow((1.0-2.0*wrongtag_eff->GetEfficiency(1)),4.0)*deff*deff;
   double term2 = 16*pow(eff_err->GetEfficiency(1),2)*pow( (1.0-2.0*wrongtag_eff->GetEfficiency(1)),2)*dwtag*dwtag;

   double dtagpow = sqrt(term1+term2);		
   cout << "tag power maximum error " << dtagpow<< endl;

  cout << " " << endl; 
  cout << "Bp MC sample" << endl; 
  cout << " " << endl; 


 // TFile *BpcutFile = TFile::Open("CutOptFileMCBpNewTrig.root"); (old one) 
  TFile *BpcutFile = TFile::Open("CutOptFiles/CutOptFileBpMCFinal.root");
  TTree* BpcutTree = (TTree*)BpcutFile->Get("cutTree");
  int entriesBp = BpcutTree->GetEntries("BpMass > 5.103 && BpMass < 5.455");
/*
   TString WtagcutBp = "DeltaR>0.300000 && MuIP<0.100000 && MuPt>1.900000 && MuCharge == 1 && PFMuon==1"; 
   	
   TString OptimizedCut = "DeltaR>0.300000 && MuIP<0.100000 && MuPt>1.900000 && PFMuon==1 && MuCharge !=0"; 
*/

   TString WtagcutBp = BestCuts + "&& BpMass > 5.103 && BpMass < 5.455 && MuCharge == 1 && PFMuon==1"; 
   TString OptimizedCut = BestCuts + "&& BpMass > 5.103 && BpMass < 5.455 && PFMuon==1 && MuCharge !=0"; 
	
  cout << WtagcutBp << endl;
  cout << OptimizedCut << endl;  

  TEventList *lista = new TEventList("lista"); 
  BpcutTree->Draw(">>lista",WtagcutBp);   
  Int_t WtagsBp = lista->GetN();

  BpcutTree->Draw(">>lista",OptimizedCut);
  Int_t TotalTagsBp = lista->GetN();

  cout << "total Bp " << TotalTagsBp << " " << "WtagsBp " << WtagsBp << endl;

  TEfficiency* eff_errBp = new TEfficiency("effBp","effBp",10,0,2);
   
    eff_errBp->SetTotalEvents(1,entriesBp);  
    eff_errBp->SetPassedEvents(1, TotalTagsBp );

   cout << "eff  "<< eff_errBp->GetEfficiency(1) << endl;
   cout <<"eff error low " <<eff_errBp->GetEfficiencyErrorLow(1) << endl;
   cout <<"eff error up " <<eff_errBp->GetEfficiencyErrorUp(1) << endl;

   TEfficiency* wrongtag_effBp = new TEfficiency("w_effBp","w_effBp",10,0,1);
     
   wrongtag_effBp->SetTotalEvents(1,TotalTagsBp);   
   wrongtag_effBp->SetPassedEvents(1,WtagsBp ); 
 
   cout <<"wfrac " << wrongtag_effBp->GetEfficiency(1) << endl;
   cout <<"wfrac error low " <<wrongtag_effBp->GetEfficiencyErrorLow(1) << endl;
   cout <<"wfrac error up " <<wrongtag_effBp->GetEfficiencyErrorUp(1) << endl;

   cout << "tagging power " << eff_errBp->GetEfficiency(1)*(1-2*wrongtag_effBp->GetEfficiency(1))*(1-2*wrongtag_effBp->GetEfficiency(1)) << endl; 

 
   double tagpowerBp = eff_errBp->GetEfficiency(1)*(1-2*wrongtag_effBp->GetEfficiency(1))*(1-2*wrongtag_effBp->GetEfficiency(1));

 double  deff = (eff_errBp->GetEfficiencyErrorLow(1) + eff_errBp->GetEfficiencyErrorUp(1))/2.0;
 
 double  dwtag = (wrongtag_effBp->GetEfficiencyErrorLow(1) + wrongtag_effBp->GetEfficiencyErrorUp(1) )/2.0;

   cout << "deff " << deff << endl;
   cout << "dwtag " << dwtag << endl;
double   term1 = pow((1.0-2.0*wrongtag_effBp->GetEfficiency(1)),4.0)*deff*deff;
   cout << "termi1 " << term1 << endl;	
 double  term2 = 16*pow(eff_errBp->GetEfficiency(1),2)*pow( (1.0-2.0*wrongtag_effBp->GetEfficiency(1)),2)*dwtag*dwtag;

    cout << "termi2 " << term2 << endl;
//   term3 = 8*eff_errBp->GetEfficiency(1)*pow((1.0-2.0*wrongtag_effBp->GetEfficiency(1)),3.0)*deff*dwtag;

double   dtagpow = sqrt(term1+term2);		
   cout << "tag power maximum error " << dtagpow<< endl;


//---------------------------- Bd MC sample

cout << endl;
cout << "Bd MC sample " << endl;
cout << endl;

 //TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileMCBdNewTrig.root"); //old one
TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileBdMCFinal.root"); 
 TTree* BdTree = (TTree*)BdFile->Get("cutTree");

   TString WtagcutBd = BestCuts + "&& BpMass > 5.039 && BpMass < 5.52 && PFMuon==1 && ( ( BdFlavour == 1 && MuCharge == 1 ) || ( MuCharge == -1 && BdFlavour == -1  ) )"; 
   TString OptimizedCut = BestCuts + "&& BpMass > 5.039 && BpMass < 5.52 && PFMuon==1 && MuCharge!=0"; 
  

  cout << WtagcutBd << endl;
  cout << OptimizedCut << endl;  

  TEventList *list = new TEventList("list"); 
  BdTree->Draw(">>list",WtagcutBd);   
  Int_t WtagsBd = list->GetN();


  BdTree->Draw(">>list",OptimizedCut);
  Int_t TotalTagsBd = list->GetN();

  cout << "total Bd " << TotalTagsBd << " " << "WtagsBd " << WtagsBd << endl;
  int entriesBd = BdTree->GetEntries("BpMass > 5.039 && BpMass < 5.52");


 TEfficiency* eff_errBd = new TEfficiency("effBd","effBd",10,0,2);
   
    eff_errBd->SetTotalEvents(1,entriesBd);  
    eff_errBd->SetPassedEvents(1, TotalTagsBd );

   cout << "eff  "<< eff_errBd->GetEfficiency(1) << endl;
   cout <<"eff error low " <<eff_errBd->GetEfficiencyErrorLow(1) << endl;
   cout <<"eff error up " <<eff_errBd->GetEfficiencyErrorUp(1) << endl;

   TEfficiency* wrongtag_effBd = new TEfficiency("w_effBd","w_effBd",10,0,1);
     
   wrongtag_effBd->SetTotalEvents(1,TotalTagsBd);   
   wrongtag_effBd->SetPassedEvents(1,WtagsBd ); 
 
   cout <<"wfrac " << wrongtag_effBd->GetEfficiency(1) << endl;
   cout <<"wfrac error low " <<wrongtag_effBd->GetEfficiencyErrorLow(1) << endl;
   cout <<"wfrac error up " <<wrongtag_effBd->GetEfficiencyErrorUp(1) << endl;

   cout << "tagging power " << eff_errBd->GetEfficiency(1)*(1-2*wrongtag_effBd->GetEfficiency(1))*(1-2*wrongtag_effBd->GetEfficiency(1)) << endl; 


   double tagpowerBd = eff_errBd->GetEfficiency(1)*(1-2*wrongtag_effBd->GetEfficiency(1))*(1-2*wrongtag_effBd->GetEfficiency(1));

   deff = (eff_errBd->GetEfficiencyErrorLow(1) + eff_errBd->GetEfficiencyErrorUp(1))/2.0;
 
   dwtag = (wrongtag_effBd->GetEfficiencyErrorLow(1) + wrongtag_effBd->GetEfficiencyErrorUp(1) )/2.0;

   cout << "deff " << deff << endl;
   cout << "dwtag " << dwtag << endl;

   term1 = pow((1.0-2.0*wrongtag_effBd->GetEfficiency(1)),4.0)*deff*deff;
   cout << "termi1 " << term1 << endl;	

   term2 = 16*pow(eff_errBd->GetEfficiency(1),2)*pow( (1.0-2.0*wrongtag_effBd->GetEfficiency(1)),2)*dwtag*dwtag;

    cout << "termi2 " << term2 << endl;

  double dtagpowBd = sqrt(term1+term2);		

   cout << "tag power maximum error " << dtagpowBd<< endl;

}
