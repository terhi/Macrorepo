void taggingWeighted(){

   TFile *McBpFile = TFile::Open("CutOptFiles/MuMcOppo.root");
   TTree* myBpTree = (TTree*)McBpFile->Get("cutTree");	

   TFile *McBpFile = TFile::Open("CutOptFiles/MuMcOppoBs.root");
   TTree* myBsTree = (TTree*)McBpFile->Get("cutTree");	

    TFile *McBdFile = TFile::Open("CutOptFiles/MuMcOppoBd.root");
   TTree* myBdTree = (TTree*)McBdFile->Get("cutTree");			

   TString GeneralCut = "DeltaR > 0.3 && MuIP < 0.3 && MuPt > 2.9 ";

   TString WtagcutBp_OppoBp = "&& MuCharge == 1 && PFMuon==1 && OppoBpMuon == 1";
   TString WtagcutBp = "&& MuCharge == 1 && PFMuon==1";

   TString TotaltagcutBp_OppoBp = "&& MuCharge!=0 && OppoBpMuon == 1"; 
   TString TotaltagcutBp ="&& MuCharge!=0"; 

//   cout << GeneralCut + TotaltagcutBp << endl;
//   cout << GeneralCut + WtagcutBp << endl;

    TEventList *Bplista = new TEventList("Bplista"); 
    myBpTree->Draw(">>Bplista",GeneralCut + WtagcutBp);   
    Int_t WtagsBpTot = Bplista->GetN();

    myBpTree->Draw(">>Bplista",GeneralCut + WtagcutBp_OppoBp);   
    Int_t WtagsBpOppoBp = Bplista->GetN();	

    cout << "WtagsBpTot " << WtagsBpTot  << endl;
    cout << "WtagsBpOppoBp " << WtagsBpOppoBp  << endl;

    // Wtag muons assumed to come from opposite side B meson must be multiplied by factor 2. this means that WtagsBpTot will be increased by amount of WtagsBpOppoBp's

    
    WtagsBpTot = WtagsBpOppoBp + WtagsBpTot;
    cout << "WtagsBpTot scaled " << WtagsBpTot  << endl;	 	

    myBpTree->Draw(">>Bplista",GeneralCut + TotaltagcutBp);   
    Int_t TtagsBpTot = Bplista->GetN();

    myBpTree->Draw(">>Bplista",GeneralCut + TotaltagcutBp_OppoBp);   
    Int_t TtagsBpOppoBp = Bplista->GetN();

    cout << "TtagsBpTot " << TtagsBpTot  << endl;
    cout << "TtagsBpOppoBp " << TtagsBpOppoBp  << endl;	

    myBpTree->Draw(">>Bplista"," OppoBpMeson == 1" ); 
    Int_t nentiresOppoBp = Bplista->GetN();

    int nentriesBp = myBpTree->GetEntries() + nentiresOppoBp;
// int nentriesBp = myBpTree->GetEntries();
     cout << " nentriesBp "<< nentriesBp << endl; 

   TtagsBpTot = TtagsBpOppoBp + TtagsBpTot;
	 
    cout << "TtagsBpTot scaled " << TtagsBpTot  << endl;
	
    TEfficiency* wrongtag_effBp = new TEfficiency("w_effBp","w_effBp",10,0,1);

  //WTAG FRAC
   
    wrongtag_effBp->SetTotalEvents(1,TtagsBpTot);   
    wrongtag_effBp->SetPassedEvents(1,WtagsBpTot); 
  
    double dwtag = (wrongtag_effBp->GetEfficiencyErrorLow(1) + wrongtag_effBp->GetEfficiencyErrorUp(1) )/2.0;	

   cout <<"wfrac Bp " << wrongtag_effBp->GetEfficiency(1) << endl;
   cout <<"wfrac error Bp " <<dwtag << endl;

  //EFFICIENCY 

    wrongtag_effBp->SetTotalEvents(2,nentriesBp ); 
    wrongtag_effBp->SetPassedEvents(2,TtagsBpTot);   
    
    double deff = (wrongtag_effBp->GetEfficiencyErrorLow(2) + wrongtag_effBp->GetEfficiencyErrorUp(2) )/2.0;	

   cout <<"efficiency  Bp " << wrongtag_effBp->GetEfficiency(2) << endl;
   cout <<"efficiency error Bp " <<deff << endl;
	
   //TAGGING POWER
 double eff = wrongtag_effBp->GetEfficiency(2);
 double wtag = wrongtag_effBp->GetEfficiency(1);

 double tagpowerBp = eff*(1-2*wtag)*(1-2*wtag);

   cout << "tagging power Bp " << tagpowerBp << endl; 	

  double term1 = pow((1.0-2.0*wtag),4.0)*deff*deff;
//   cout << "termi1 " << term1 << endl;
	
   double term2 = 16*pow(eff,2)*pow( (1.0-2.0*wtag),2)*dwtag*dwtag;
//   cout << "termi2 " << term2 << endl;

   double dtagpow = sqrt(term1+term2);		
   cout << "tag power error Bp " << dtagpow << endl;	
  
  
  //////////////// ---------- Bs --------------- //////////////

	cout << "" <<endl;
	cout << "Bs " << endl;

   TString WtagcutBs_OppoBs = "&& OppoBsMuon == 1 && PFMuon == 1 && ((BsFlavour == 1 && MuCharge == 1) || (BsFlavour == -1 && MuCharge == -1) )";

   TString WtagcutBs = "&& PFMuon == 1 && ((BsFlavour == 1 && MuCharge == 1) || (BsFlavour == -1 && MuCharge == -1) )";
      
    TString TotaltagcutBs_OppoBs = "&& OppoBsMuon == 1 && PFMuon == 1 && (BsFlavour == -1 || BsFlavour == 1) ";
 
   TString TotaltagcutBs = "&& PFMuon == 1 && (BsFlavour == -1 || BsFlavour == 1) "; 

    TEventList *Bslista = new TEventList("Bslista"); 
    myBsTree->Draw(">>Bslista",GeneralCut + WtagcutBs);   
    Int_t WtagsBsTot = Bslista->GetN();

    myBsTree->Draw(">>Bslista",GeneralCut + WtagcutBs_OppoBs);   
    Int_t WtagsBsOppoBs = Bslista->GetN();	

    cout << "WtagsBsTot " << WtagsBsTot  << endl;
    cout << "WtagsBsOppoBs " << WtagsBsOppoBs  << endl;

   WtagsBsTot = WtagsBsOppoBs + WtagsBsTot;
    cout << "WtagsBsTot scaled " << WtagsBsTot  << endl;	 	

    myBsTree->Draw(">>Bslista",GeneralCut + TotaltagcutBs);   
    Int_t TtagsBsTot = Bslista->GetN();

    myBsTree->Draw(">>Bslista",GeneralCut + TotaltagcutBs_OppoBs );   
    Int_t TtagsBsOppoBs = Bslista->GetN();

    cout << "TtagsBsTot " << TtagsBsTot  << endl;
    cout << "TtagsBsOppoBs " << TtagsBsOppoBs  << endl;	


 //   cout << " nentriesBs "<< nentriesBs << endl;

    TtagsBsTot = TtagsBsOppoBs + TtagsBsTot;
	 
    cout << "TtagsBsTot scaled " << TtagsBsTot  << endl;
	
    TEfficiency* wrongtag_effBs = new TEfficiency("w_effBs","w_effBs",10,0,1);
     
  //WTAG FRAC

    wrongtag_effBs->SetTotalEvents(1,TtagsBsTot);   
    wrongtag_effBs->SetPassedEvents(1,WtagsBsTot ); 
  
    double dwtag = (wrongtag_effBs->GetEfficiencyErrorLow(1) + wrongtag_effBs->GetEfficiencyErrorUp(1) )/2.0;	

   cout <<"wfrac Bs " << wrongtag_effBs->GetEfficiency(1) << endl;
   cout <<"wfrac error Bs " <<dwtag << endl;

  //EFFICIENCY 
  myBsTree->Draw(">>Bslista"," OppoBsMeson == 1" ); 
  Int_t nentiresOppoBs = Bslista->GetN();

  int nentriesBs = myBsTree->GetEntries() + nentiresOppoBs;
//int nentriesBs = myBsTree->GetEntries() ;
     cout << " nentriesBs "<< nentriesBs << endl; 


    wrongtag_effBs->SetTotalEvents(2,nentriesBs ); 
    wrongtag_effBs->SetPassedEvents(2,TtagsBsTot);   
    
    double deff = (wrongtag_effBs->GetEfficiencyErrorLow(2) + wrongtag_effBs->GetEfficiencyErrorUp(2) )/2.0;	

   cout <<"efficiency  Bs " << wrongtag_effBs->GetEfficiency(2) << endl;
   cout <<"efficiency error Bs " <<deff << endl;


  //TAGGING POWER

   
 double eff = wrongtag_effBs->GetEfficiency(2);
 double wtag = wrongtag_effBs->GetEfficiency(1);

 double tagpowerBs = eff*(1-2*wtag)*(1-2*wtag);

   cout << "tagging power Bs " << tagpowerBs << endl; 	

  double term1 = pow((1.0-2.0*wtag),4.0)*deff*deff;
//   cout << "termi1 " << term1 << endl;
	
   double term2 = 16*pow(eff,2)*pow( (1.0-2.0*wtag),2)*dwtag*dwtag;
//   cout << "termi2 " << term2 << endl;

   double dtagpow = sqrt(term1+term2);		
   cout << "tag power error Bs " << dtagpow << endl;	


 //////////////// ---------- Bd --------------- //////////////

	cout << "" <<endl;
	cout << "Bd " << endl;

   TString WtagcutBd_OppoBd = "&& OppoBdMeson == 1 && PFMuon == 1 && ((BdFlavour == 1 && MuCharge == 1) || (BdFlavour == -1 && MuCharge == -1) )";

   TString WtagcutBd = "&& PFMuon == 1 && ((BdFlavour == 1 && MuCharge == 1) || (BdFlavour == -1 && MuCharge == -1) )";
      
    TString TotaltagcutBd_OppoBd = "&& OppoBdMeson==1 && PFMuon == 1 && (BdFlavour == -1 || BdFlavour == 1) ";
 
   TString TotaltagcutBd = "&& PFMuon == 1 && (BdFlavour == -1 || BdFlavour == 1) "; 

    TEventList *Bdlista = new TEventList("Bdlista"); 
    myBdTree->Draw(">>Bdlista",GeneralCut + WtagcutBd);   
    Int_t WtagsBdTot = Bdlista->GetN();

    myBdTree->Draw(">>Bdlista",GeneralCut + WtagcutBd_OppoBd);   
    Int_t WtagsBdOppoBd = Bdlista->GetN();	

    cout << "WtagsBdTot " << WtagsBdTot  << endl;
    cout << "WtagsBdOppoBd " << WtagsBdOppoBd  << endl;

    WtagsBdTot = WtagsBdOppoBd + WtagsBdTot;
    cout << "WtagsBdTot scaled " << WtagsBdTot  << endl;	 	

    myBdTree->Draw(">>Bdlista",GeneralCut + TotaltagcutBd);   
    Int_t TtagsBdTot = Bdlista->GetN();

    myBdTree->Draw(">>Bdlista",GeneralCut + TotaltagcutBd_OppoBd );   
    Int_t TtagsBdOppoBd = Bdlista->GetN();

    cout << "TtagsBdTot " << TtagsBdTot  << endl;
    cout << "TtagsBdOppoBd " << TtagsBdOppoBd  << endl;	

    TtagsBdTot = TtagsBdOppoBd + TtagsBdTot;
	 
    cout << "TtagsBdTot scaled " << TtagsBdTot  << endl;

    TEfficiency* wrongtag_effBd = new TEfficiency("w_effBd","w_effBd",10,0,1);
     
  //WTAG FRAC

    wrongtag_effBd->SetTotalEvents(1,TtagsBdTot);   
    wrongtag_effBd->SetPassedEvents(1,WtagsBdTot ); 
  
    double dwtag = (wrongtag_effBd->GetEfficiencyErrorLow(1) + wrongtag_effBd->GetEfficiencyErrorUp(1) )/2.0;	

   cout <<"wfrac Bd " << wrongtag_effBd->GetEfficiency(1) << endl;
   cout <<"wfrac error Bd " <<dwtag << endl;

  //EFFICIENCY 
  myBdTree->Draw(">>Bdlista"," OppoBdMeson == 1" ); 
  Int_t nentiresOppoBd = Bdlista->GetN();

	cout <<"total oppoBd events " << nentiresOppoBd << endl;
  int nentriesBd = myBdTree->GetEntries() + nentiresOppoBd;
     cout << " nentriesBd "<< nentriesBd << endl; 


    wrongtag_effBd->SetTotalEvents(2,nentriesBd ); 
    wrongtag_effBd->SetPassedEvents(2,TtagsBdTot);   
    
    double deff = (wrongtag_effBd->GetEfficiencyErrorLow(2) + wrongtag_effBd->GetEfficiencyErrorUp(2) )/2.0;	

   cout <<"efficiency  Bd " << wrongtag_effBd->GetEfficiency(2) << endl;
   cout <<"efficiency error Bd " <<deff << endl;


  //TAGGING POWER

   
 double eff = wrongtag_effBd->GetEfficiency(2);
 double wtag = wrongtag_effBd->GetEfficiency(1);

 double tagpowerBd = eff*(1-2*wtag)*(1-2*wtag);

   cout << "tagging power Bd " << tagpowerBd << endl; 	

  double term1 = pow((1.0-2.0*wtag),4.0)*deff*deff;
//   cout << "termi1 " << term1 << endl;
	
   double term2 = 16*pow(eff,2)*pow( (1.0-2.0*wtag),2)*dwtag*dwtag;
//   cout << "termi2 " << term2 << endl;

   double dtagpow = sqrt(term1+term2);		
   cout << "tag power error Bd " << dtagpow << endl;	




}
