#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TStyle.h"
#include "Roo1DTable.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "TRandom3.h"

#include "RooConstVar.h"
#include "RooGlobalFunc.h"
#include "RooGaussian.h"


using namespace RooFit;
using namespace std;


void Bdmixing(){
//TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileBdCuts.root");
//TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileBdDataTag.root");
TFile *BdFile = TFile::Open("CutOptFiles/CutOptFileBdDataTagNew.root");
TTree* myBdTree = (TTree*)BdFile->Get("cutTree");

  RooRealVar BdCt("BdCt2D","B^{0} 2D ct [cm/c]",0.04,0.3);
  RooRealVar BpMass("BpMass","B^{0} mass [GeV/c^{2}]",5.0,5.5);
  
  RooRealVar MuCharge("MuCharge","MuCharge",-1,1);
  RooRealVar BdecayFlavour("BdecayFlavour","BdecayFlavour",-1,1);

  RooRealVar MuPt("MuPt","MuPt",0,20);
  RooRealVar MuIP("MuIP","MuIP",0,0.5);
  RooRealVar DeltaR("DeltaR","DeltaR",0,5); 


  RooCategory mixState("mixState","B0/B0bar mixing state") ;
  mixState.defineType("mixed",-1) ;
  mixState.defineType("unmixed",1) ;

  RooCategory tagFlav("tagFlav","Flavour of the tagged B0") ;
  tagFlav.defineType("B0",1) ;
  tagFlav.defineType("B0bar",-1) ;

  TCut signalcut = "BpMass>5.17 && BpMass < 5.37"; 
  TCut datacut = "DeltaR > 0.3 && MuPt > 2.9 && MuIP < 0.3 && BdCt2D > 0.02"; 
  TCut bkgcut1 = "BpMass < 5.17 && BpMass > 5.0"; 
  TCut bkgcut2 = "BpMass > 5.37 && BpMass< 5.5";
 
  RooDataSet data("Alldata","Alldata",RooArgSet(tagFlav,mixState,BpMass,BdCt,MuPt,MuIP,DeltaR),Import(*myBdTree),Cut(datacut));
  data.Print("v");
   
   RooDataSet* signaldata = (RooDataSet*)data.reduce(RooArgSet(tagFlav,mixState,BdCt),signalcut);
  // signaldata->Print("v") ;


   RooDataSet* bkgdatalow = (RooDataSet*)data.reduce(RooArgSet(tagFlav,BdCt,BpMass),bkgcut1); //onko ok?
   RooDataSet* bkgdataup = (RooDataSet*)data.reduce(RooArgSet(tagFlav,BdCt,BpMass),bkgcut2);
	
   bkgdataup->append(*bkgdatalow);
 //  bkgdataup->Print("v") ;	

   cout << "total events bkg " << bkgdataup->numEntries() << endl;

   Roo1DTable* btable = data.table(tagFlav);
   btable->Print("v") ;  
   Roo1DTable* table = data.table(mixState);
   table->Print("v") ;  

   //1st component total bkg PDF  

   // two exponentials as a form of Roodecay to describe the ct bkg and gaussian resolution
   RooRealVar mean_ct("mean_ct","mean_ct",0);
   RooRealVar sigma_ct("sigma_ct","sigma_ct",0.0022) ; //0.0022
   RooGaussModel gm_ct("gm_ct","gauss model",BdCt,mean_ct,sigma_ct);

   RooRealVar TauBkg1("TauBkg1","Tau 1 from the SideBands",0.001,2);
   RooDecay Bck_Ctau1("Bck_Ctau1","ctau bck",BdCt,TauBkg1,gm_ct,RooDecay::SingleSided);

   RooRealVar TauBkg2("TauBkg2","Tau 2 from the SideBands",0.001,2);
   RooDecay Bck_Ctau2("Bck_Ctau2","ctau bck",BdCt,TauBkg2,gm_ct,RooDecay::SingleSided);
   RooRealVar fracBkgCt("fracBkg","fraction of bkg events",0.5,0.,1.);

   RooAddPdf doubleRooDecay("doubleRooDecay","RooDecay1+RooDecay2", RooArgList(Bck_Ctau2,Bck_Ctau1),fracBkgCt);

   doubleRooDecay.fitTo(*bkgdataup);

   TCanvas *bkg = new TCanvas("bkg","bkg");
   bkg->SetLogy();
   RooPlot* bkgframe = BdCt.frame(Title("Bkg events"));	
   bkgdataup->plotOn(bkgframe);
   doubleRooDecay.plotOn(bkgframe);  
   doubleRooDecay.plotOn(bkgframe,Components("Bck_Ctau1"),LineStyle(kDashed),LineColor(6));
   doubleRooDecay.plotOn(bkgframe,Components("Bck_Ctau2"),LineStyle(kDashed),LineColor(2));
   bkgframe->Draw(); 

   TauBkg2.setConstant(kTRUE);
   TauBkg1.setConstant(kTRUE);
   fracBkgCt.setConstant(kTRUE); 
   mean_ct.setConstant(kTRUE); 
   sigma_ct.setConstant(kTRUE);
    
	//PDF for Bd mass

   TCanvas *BdMassCanv = new TCanvas("BdMass","BdMass");
   RooPlot* massframe = BpMass.frame(Title("Bd mass events"));	

        // Bkg for the Bd mass
   RooRealVar coeff1("coeff1","coeff1",-10.,0.);
   RooExponential BdMassBkg("BdMassBkg","BdMassBkg",BpMass,coeff1);
 
   RooRealVar nbkg("nbkg","background events BpMass",0,data.numEntries());
   RooRealVar nsig("nsig","signal events BpMass",0,data.numEntries());
  
      //gaussians for the Bd mass
   RooRealVar mean("mean","mean of two gaussians",5.25,5.2,5.35);
   RooRealVar sigma1("sigma1","width of gaussian R",0.02,0.007,0.045);
   RooGaussian gauss1("gauss1","gaussian PDF R",BpMass,mean,sigma1);
   RooRealVar fGG("fGGR","fraction of signal events",0.5,0.,1.);

   RooRealVar sigma2("sigma2","width of 2nd gaussian R",0.03,0.007,0.035); 
   RooGaussian gauss2("gauss2","gaussian PDF 2",BpMass,mean,sigma2);

   RooAddPdf DoubleGauss("DoubleGauss","gauss1+gauss2", RooArgList(gauss1,gauss2),fGG);

   //Bd mass PDF tot 
   RooAddPdf PDFBdMass("PDFBdMass","Signal + Bck Pdf for the mass",RooArgList(DoubleGauss,BdMassBkg),RooArgList(nsig,nbkg));

  PDFBdMass.fitTo(data,Extended(kTRUE));
  data.plotOn(massframe) ;
  PDFBdMass.plotOn(massframe) ;
  PDFBdMass.plotOn(massframe,Components("gauss1"),LineStyle(kDashed),LineColor(6));
  PDFBdMass.plotOn(massframe,Components("gauss2"),LineStyle(kDashed),LineColor(7));
  PDFBdMass.plotOn(massframe,Components("BdMassBkg"),LineStyle(kDashed),LineColor(3)); 
  massframe->Draw(); 

  nsig.Print();
  nbkg.Print();

//    BpMass.setRange("SB1",5.0,5.17) ;
//    BpMass.setRange("SB2",5.37,5.5) ;
//    BpMass.setRange("Signal",5.17,5.37) ;

   mean.setConstant(kTRUE);
   sigma1.setConstant(kTRUE);
//   sigma2.setConstant(kTRUE);
   fGG.setConstant(kTRUE);

   //product PDF for total bkg
   RooProdPdf PDF_Bck( "PDF_Bck" ,"BdMassBkg*doubleRooDecay", RooArgSet(doubleRooDecay,BdMassBkg));

  //product PDF for total signal

  // Construct Bdecay with mixing
  RooRealVar dm("dm","delta m(B0)",10.,20.); // dm = 0.472ps^-1 -> dm/ct = 15.89 cm^-1 ! 
  RooRealVar tau("tau","tau (B0)",0.,0.1); //  tau = 1.547ps -> ctau = 0.0462 cm!
  RooRealVar w("w","flavour mistag rate",0.308);
  RooRealVar dw("dw","delta mistag rate for B0/B0bar",0.01);
  RooBMixDecay bmix("bmix","decay",BdCt,mixState,tagFlav,tau,dm,w,dw,gm_ct,RooBMixDecay::SingleSided);

   RooProdPdf PDF_Signal( "PDF_Signal" ,"DoubleGauss*bmix", RooArgSet(DoubleGauss,bmix));

  //total 2D PDF  
   RooAddPdf PDF_Total("PDF_Total","Signal + Bck Pdf",RooArgList(PDF_Signal,PDF_Bck),RooArgList(nsig,nbkg));
 
  RooFitResult *Result = PDF_Total.fitTo(data,Extended(kTRUE),Save(),NumCPU(8));

  //Projection of 2D PDF on mass and ct axis
  RooPlot* ProjframeCt = BdCt.frame(Title("Projection of total PDF on BdCt"));
  data.plotOn(ProjframeCt);
  PDF_Total.plotOn(ProjframeCt);
  PDF_Total.plotOn(ProjframeCt,Components("bmix"),ProjWData(BdCt,data),LineStyle(kDashed),LineColor(9));
  PDF_Total.plotOn(ProjframeCt,Components("Bck_Ctau1"),ProjWData(BdCt,data),LineStyle(kDashed),LineColor(6));
  PDF_Total.plotOn(ProjframeCt,Components("Bck_Ctau2"),ProjWData(BdCt,data),LineStyle(kDashed),LineColor(2));  

  RooPlot* ProjframeBdM = BpMass.frame(Title("Projection of total PDF on Bd Mass"));  
  data.plotOn(ProjframeBdM);
  PDF_Total.plotOn(ProjframeBdM);
  PDF_Total.plotOn(ProjframeBdM ,Components("gauss1"),ProjWData(BpMass,data),LineStyle(kDashed),LineColor(6));
  PDF_Total.plotOn(ProjframeBdM ,Components("gauss2"),ProjWData(BpMass,data),LineStyle(kDashed),LineColor(7));
  PDF_Total.plotOn(ProjframeBdM ,Components("BdMassBkg"),ProjWData(BpMass,data),LineStyle(kDashed),LineColor(3));
   
  dm.Print();
  tau.Print();

 /*
  TCanvas* c = new TCanvas("Bdmixing","Bdmixing",800,600) ;  
  c->Divide(1,2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; gPad->SetLogy();  ProjframeCt->Draw();
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; ProjframeBdM->Draw();
	 
  //frame with mixing B0 -> B0bar and  B0bar -> B0
  RooPlot* frame1 = BdCt.frame(Title("B0 ct distribution of mixed events")) ;

  //how does the projection work? should I projectom ct or not?

  TCanvas *mixed = new TCanvas("mixed","mixed");
  data.plotOn(frame1,Cut("tagFlav==tagFlav::B0 && mixState==mixState::mixed"));
//  data.plotOn(frame1,Cut("tagFlav==tagFlav::B0bar && mixState==mixState::mixed"),MarkerColor(3));
  //  data.plotOn(frame1,Cut("mixState==mixState::mixed"));
    PDF_Total.plotOn(frame1,Slice(mixState,"mixed"),Slice(tagFlav,"B0"));
    PDF_Total.plotOn(frame1,Slice(mixState,"mixed"),Slice(tagFlav,"B0"),Components("bmix"));
 //   PDF_Total.plotOn(frame1,Slice(mixState,"mixed"),Components("bmix"),LineStyle(kDashed),ProjWData(mixState,data),LineColor(9));
 //   PDF_Total.plotOn(ProjframeCt,Components("Bck_Ctau1"),ProjWData(BdCt,data),LineStyle(kDashed),LineColor(6));
//  PDF_Total.plotOn(ProjframeCt,Components("Bck_Ctau2"),ProjWData(BdCt,data),LineStyle(kDashed),LineColor(2));  

    frame1->Draw();
*/
/*
  RooPlot* frame2 = BdCt.frame(Title("B0 ct distribution of unmixed events")) ;  

  TCanvas *unmixed = new TCanvas("unmixed","unmixed");
//  data->plotOn(frame2,Cut("tagFlav==tagFlav::B0 && mixState==mixState::unmixed")); 
//  data->plotOn(frame2,Cut("tagFlav==tagFlav::B0bar && mixState==mixState::unmixed"), MarkerColor(3)); 
   data.plotOn(frame2,Cut("mixState==mixState::unmixed"));
PDF_Total.plotOn(frame2,Slice(mixState,"unmixed")); 
  
  frame2->Draw(); 
   
*/

  // Construct signal, total likelihood observables and likelihood ratio
  RooAbsPdf* sigProj = PDF_Signal.createProjection(BdCt) ; 
  RooAbsPdf* totalProj= PDF_Total.createProjection(BdCt) ; 
  RooFormulaVar llratio_func("llratio","log10(@0)-log10(@1)",RooArgList(*sigProj,*totalProj)) ;

 RooRealVar* llratio = (RooRealVar*)data.addColumn(llratio_func) ;
 data.Print("v");
  // Calculate likelihood ratio for each event, define subset of events with high signal likelihood

  RooDataSet* dataSel = (RooDataSet*) data.reduce(Cut("llratio>0.2")) ;
 
  // Make plot frame and plot data
//  RooPlot* frame = llratio->frame(-0.1,0.4,50) ;
//  data.plotOn(frame);
 
/////////// simulation

//RooRealVar ct("ct", "ct", 0.02,0.3);
  RooRealVar dmSim("dm","delta m(B0)",15.98) ;  // dm = 0.472ps^-1 -> dm/ct = 15.89 cm^-1 !
  RooRealVar tauSim("tau","tau (B0)",0.047) ; //  tau = 1.547ps -> ctau = 0.0462 cm!
//  RooRealVar wSim("w","flavour mistag rate",0.308) ;
 // RooRealVar dwSim("dw","delta mistag rate for B0/B0bar",0.01) ;
 RooBMixDecay bmixSimu("bmixSim","decaySim",BdCt,mixState,tagFlav,tauSim,dmSim,w,dw,gm_ct,RooBMixDecay::SingleSided);

 RooDataSet* dataSimu = bmixSimu.generate(RooArgSet(BdCt,mixState,tagFlav),4700) ;
 RooDataSet* dataSimuMoremistag = new RooDataSet("dataSim","dataSim",RooArgList(mixState,tagFlav,BdCt));
	
Roo1DTable* tableSim = dataSimu->table(mixState);
   tableSim->Print("v") ;

  TRandom3 *rand = new TRandom3();

  Int_t Nentries = dataSimu->numEntries();
//  cout << nentries << endl; 	
 
    for(int i =0; i < Nentries; i++){
     const RooArgSet* obs = dataSimu->get(i);
    
     RooCategory* tagFlavCopy = (RooCategory*)obs->find("tagFlav") ;
     RooRealVar* dtCopy = (RooRealVar*)obs->find("BdCt2D");
     RooCategory* tagmix = (RooCategory*)obs->find("mixState") ;
  //   cout << "index ini "<< tagmix->getIndex() <<endl;
     double random = rand->Uniform(0,1);
//	cout << random << endl;
     if(random < 0.1){ 		   	
     tagmix->setIndex(-1*tagmix->getIndex());     	
     }	
	dataSimuMoremistag->add(RooArgSet(*tagmix,*tagFlavCopy,*dtCopy));
 //    cout << "index swapped "<<tagmix->getIndex() <<endl;	

    }

 Roo1DTable* tablemis = dataSimuMoremistag->table(mixState);
   tablemis->Print("v") ;

///////////

 //Mixing with B0bar -> B0
  RooPlot* frameB0bar = BdCt.frame(Title("Mixed B^{0} events")) ;
  dataSel->plotOn(frameB0bar,Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0bar"));
  
 PDF_Total.plotOn(frameB0bar,Components("bmix"),Slice(mixState,"mixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),NumCPU(4),LineStyle(kDashed),LineColor(9)) ;

 PDF_Total.plotOn(frameB0bar,Slice(mixState,"mixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),NumCPU(4)) ;

  //dataSimu->plotOn(frame,Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0"),MarkerColor(3),LineColor(3));
  //dataSimuMoremistag->plotOn(frame,Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0bar"),MarkerColor(2),LineColor(2));
  //bmixSimu.plotOn(frame,Slice(tagFlav,"B0bar"),Slice(mixState,"mixed"),LineColor(3));

 PDF_Total.plotOn(frameB0bar,Components("Bck_Ctau1"),Slice(mixState,"mixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(6));
  PDF_Total.plotOn(frameB0bar,Components("Bck_Ctau2"),Slice(mixState,"mixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(2));  

  TCanvas* B0barMix = new TCanvas("B0barMix","B0barMix"); 
  B0barMix->SetLogy();
  frameB0bar->Draw();
  B0barMix->SaveAs("B0barmix.pdf");   

  //Mixing B0-> B0bar 

  RooPlot* frameB0 = BdCt.frame(Title("Mixed #bar{B}^{0} events")) ;
  dataSel->plotOn(frameB0,Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0"));
  
  PDF_Total.plotOn(frameB0,Components("bmix"),Slice(mixState,"mixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),NumCPU(4),LineStyle(kDashed),LineColor(9)) ;

 PDF_Total.plotOn(frameB0,Slice(mixState,"mixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),NumCPU(4)) ;


//  dataSimuMoremistag->plotOn(frame,Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0"),MarkerColor(2),LineColor(2));
//  bmixSimu.plotOn(frame,Slice(tagFlav,"B0"),Slice(mixState,"mixed"),LineColor(3));

 PDF_Total.plotOn(frameB0,Components("Bck_Ctau1"),Slice(mixState,"mixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(6));
  PDF_Total.plotOn(frameB0,Components("Bck_Ctau2"),Slice(mixState,"mixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(2)); 

  TCanvas* B0Mix = new TCanvas("B0Mix","B0Mix"); 
  B0Mix->SetLogy();
  frameB0->Draw();
  B0Mix->SaveAs("B0mix.pdf");   

 //Unmixed B0: B0-> B0 

  RooPlot* frameB0nomix = BdCt.frame(Title("Unmixed B^{0} events")) ;
  dataSel->plotOn(frameB0nomix,Cut("mixState==mixState::unmixed && tagFlav==tagFlav::B0"));
  
  PDF_Total.plotOn(frameB0nomix,Components("bmix"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),NumCPU(4),LineStyle(kDashed),LineColor(9)) ;

 PDF_Total.plotOn(frameB0nomix,Slice(mixState,"unmixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),NumCPU(4)) ;

  PDF_Total.plotOn(frameB0nomix,Components("Bck_Ctau1"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(6));

  PDF_Total.plotOn(frameB0nomix,Components("Bck_Ctau2"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(2)); 

  TCanvas* B0Unmix = new TCanvas("B0Unmix","B0Unmix"); 
  B0Unmix->SetLogy();
  frameB0nomix->Draw();
  B0Unmix->SaveAs("B0Unmix.pdf");  

   //Unmixed B0bar: B0bar-> B0bar 

  RooPlot* frameB0barnomix = BdCt.frame(Title("Unmixed #bar{B}^{0} events")) ;
  dataSel->plotOn(frameB0barnomix,Cut("mixState==mixState::unmixed && tagFlav==tagFlav::B0bar"));
  
  PDF_Total.plotOn(frameB0barnomix,Components("bmix"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),NumCPU(4),LineStyle(kDashed),LineColor(9)) ;

 PDF_Total.plotOn(frameB0barnomix,Slice(mixState,"unmixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),NumCPU(4)) ;

  PDF_Total.plotOn(frameB0barnomix,Components("Bck_Ctau1"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(6));

  PDF_Total.plotOn(frameB0barnomix,Components("Bck_Ctau2"),Slice(mixState,"unmixed"),Slice(tagFlav,"B0bar"),ProjWData(*dataSel),LineStyle(kDashed),LineColor(2)); 

  TCanvas* B0barUnmix = new TCanvas("B0barUnmix","B0barUnmix"); 
  B0barUnmix->SetLogy();
  frameB0barnomix->Draw();
  B0barUnmix->SaveAs("B0barUnmix.pdf");  

   TCanvas* Logratio = new TCanvas("LogRatio","LogRatio");
   Logratio->cd(); 
   RooPlot* LLframe = llratio->frame(0.,0.5,100) ;
   data.plotOn(LLframe);
   LLframe->Draw();	
/*
  // Perform parallel projection using MC integration of pdf using given input dataSet. 
  // In this mode the data-weighted average of the pdf is calculated by splitting the
  // input dataset in N equal pieces and calculating in parallel the weighted average
  // one each subset. The N results of those calculations are then weighted into the
  // final result
  
  // Use four processes
  PDF_Total.plotOn(frame,ProjWData(*dataSel),NumCPU(4)) ;
  frame->Draw();
*/

 Roo1DTable* bttable = dataSel->table(RooArgSet(tagFlav,mixState)) ;
 bttable->Print("v") ;
/*
 Int_t B0barmixed = dataSel->numEntries(Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0bar"));
  
  Int_t B0mixed = dataSel->numEntries(Cut("mixState==mixState::mixed && tagFlav==tagFlav::B0"));

    Int_t B0Unmixed = dataSel->numEntries(Cut("mixState==mixState::unmixed && tagFlav==tagFlav::B0"));
 //   Int_t B0barUnmixed = dataSel->numEntries(Cut("mixState==mixState::unmixed && tagFlav==tagFlav::B0bar"));

   cout << "mixed B0 " << B0mixed << endl;
   cout << "mixed B0bar " << B0barmixed << endl;
   cout << "unmixed B0 " << B0Unmixed << endl;
*/
}
