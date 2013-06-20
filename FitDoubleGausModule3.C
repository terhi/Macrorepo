#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCut.h"
#include "TAxis.h"
#include "TString.h"
#include "RooCategory.h"
#include "TEventList.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "Roo1DTable.h"
#include "TPaveText.h"
#include "TText.h"
#include <iomanip>  
using namespace RooFit;
using namespace std;

vector<int> FitDoubleGausModule3(TTree *tree, TString Precut,Double_t upperlim2R,Double_t lowerlimR,Int_t flag, Float_t binLow, Float_t binHigh,Int_t index){

//TFile *McBplusFile = TFile::Open("SlimmedBdata.root");
//TTree* tree = (TTree*)McBplusFile->Get("cutTree");
	
 //RooFormulaVar mycut("mycut",Precut); 
 RooRealVar MuCharge("MuCharge","MuCharge",-1,1);
 RooRealVar MuPt("MuPt","MuPt",0,20);
 RooRealVar MuIP("MuIP","MuIP",0,0.5);
 RooRealVar DeltaR("DeltaR","DeltaR",0,5); 
 RooRealVar PFMuon("PFMuon","PFMuon",0,1); 
//let's choose category to be the integer tree branch "MuCharge", this is name of the branch!
  RooCategory TagCat("MuCharge"," Tag muon Charge"); 

//let's cast the suitable values for MuCharge, wrong tag is the one in which MuCharge == 1 and right tag the one where MuCharge == -1
  TagCat.defineType("Wtag",1); 
  TagCat.defineType("Rtag",-1); 


  RooRealVar BpMass("BpMass","B^{+} mass [GeV/c^{2}]",5.18,5.5); //Sulkujen sisapuolella olevan nimen tulee olla sama kuin branchin nimi
  
 //import Bpmass events from the branch BpMass with condition, that MuCharge in tagging category has above defined values
 
  RooDataSet dataTotal("dataTotal","dataTotal",RooArgSet(TagCat,BpMass,MuIP,MuPt,DeltaR,PFMuon),Import(*tree),Cut(Precut)); 

	dataTotal.Print();
	//Double gaussian + Background PDF for RightTags 	

   Roo1DTable* btable = dataTotal.table(TagCat) ;
   btable->Print() ;
   btable->Print("v") ;

  // Retrieve number of events from table
  // Number can be non-integer if source dataset has weighed events
  Int_t Nwtag = btable->get("Wtag");
  Int_t Nrtag = btable->get("Rtag");
	
	RooRealVar meanR("meanR","mean of gaussian R",5.25,5.2,5.35);
	RooRealVar sigmaR("sigmaR","width of gaussian R",0.01,lowerlimR,0.052);
	RooGaussian gaussR("gaussR","gaussian PDF R",BpMass,meanR,sigmaR);
 	RooRealVar nsigR("nsigR","signal events R",0,Nrtag);
	RooRealVar fGGR("fGGR","fraction of signal events R",0.5,0.,1.);
 
	RooRealVar sigma2R("sigma2R","width of 2nd gaussian R",0.03,0.009,upperlim2R); 
//sigma2 changed many times lower limit 0.009-0.1 works for IP first bin2
//upper limit 0.009-0.07 works for second,third and 4th IP bin

	RooGaussian gauss2R("gauss2R","gaussian PDFR",BpMass,meanR,sigma2R);	
	RooAddPdf doubleGaussRight("DoubleGaussianR","gauss1+gauss2", RooArgList(gaussR,gauss2R),fGGR);

	 RooRealVar* coeff1R;
       
	if(binHigh < 3.16 && binHigh > 3.14 && binLow > 2.8 && binLow < 2.82 && flag == 3){
	cout <<"agfad adfafdafd asdfasfadsfad sdfafda"<< endl;	
	 coeff1R = new RooRealVar("coeff1R","coeff1",-3.,-2.);  //-2.. -1.65 ok arvo -1.81	
	}
	else coeff1R = new RooRealVar("coeff1R","coeff1",-4,0.); 
	

	RooExponential expoR("expoR","expoR",BpMass,*coeff1R);
	RooRealVar nbkgR("nbkgR","background events R",0,Nrtag);
/*	RooRealVar fEER("fEER","fEER",0.,1);
	RooRealVar div("div","div",0.035,0.01,0.12);
	RooGenericPdf erfPdf("erf","(TMath::Erf((-BpMass+5.15)/div)+1)",RooArgSet(div,BpMass)); 
	RooAddPdf bkgtot("BkgtotR","erf+expo R", RooArgList(erfPdf,expoR),fEER);
*/	
	RooAddPdf RightTagPDF("RtagPDF","DuobleGaussR + expo",RooArgList(expoR,doubleGaussRight),RooArgList(nbkgR,nsigR));

	//Double gaussian + Background PDF for WrongTags 
	//Wrong tag PDF uses the same BG function defined in Rtag section

	RooRealVar nsigW("nsigW","signal events W",0,Nwtag);
	RooRealVar nbkgW("nbkgW","background events W",0,Nwtag);
	RooAddPdf WrongTagPDF("WtagPDF","DoubleGausW + expo",RooArgList(expoR,doubleGaussRight),RooArgList(nbkgW,nsigW));

	//Total PDF for both right tags and wrong tags 
	RooSimultaneous simTagPDF("TaggingPDF","Tagging PDF",TagCat);
	simTagPDF.addPdf(RightTagPDF,"Rtag"); 
	simTagPDF.addPdf(WrongTagPDF,"Wtag"); 

	RooFitResult* Result = simTagPDF.fitTo(dataTotal,Extended(kTRUE),Save()); 
	TString CanvGood= "a";
	TString CanvBad= "b";	

	TString BintitleGood;
	TString BintitleBad;

	if(flag == 2 ){
	BintitleGood = Form("Correctly tagged subsample %f < p_{T} < %f",binLow, binHigh);
	CanvGood = Form("GoodtagsPtBdcuts_%d",index);	
	CanvBad = Form("BadtagsPtBdcuts_%d",index); 
	BintitleBad = Form("Wrongly tagged subsample %f < p_{T} < %f",binLow, binHigh); 	
	}

	if(flag == 1 ){
	BintitleGood = Form("Correctly tagged subsample %f < IP < %f",binLow, binHigh); 
	BintitleBad = Form("Wrongly tagged subsample %f < IP < %f",binLow, binHigh);
	CanvGood = Form("GoodtagsIPBdcuts_%d",index);	
	CanvBad = Form("BadtagsIPBdcuts_%d",index); 	
	}

	if(flag == 3 ){
	
	BintitleGood = Form("Correctly tagged subsample %f < #DeltaR < %f",binLow, binHigh); 
	BintitleBad = Form("Wrongly tagged subsample %f < #DeltaR < %f",binLow,binHigh);
	CanvGood = Form("GoodtagsdRBdcuts_%d",index);	
	CanvBad = Form("BadtagsdRBdcuts_%d",index);  	
	}

	TCanvas* c1 = new TCanvas(CanvGood,CanvGood);
//	c1->Divide(2,1);

	c1->cd();
	c1->SetLeftMargin(0.1457286);
	c1->SetBottomMargin(0.1541463);


	TString tags = Form("Ntags %d #pm %d",int(nsigR.getVal()),int(nsigR.getError()) );
	TString bkg = Form("Nbkg %d #pm %d",int(nbkgR.getVal()),int(nbkgR.getError()) );

	TLegend *leg1 = new TLegend(0.459799,0.6608392,0.7876884,0.8811189,NULL,"brNDC");
//TLegend(0.6168342,0.6503497,0.8668342,0.8706294,NULL,"brNDC");
	//TLegend *leg1 = new TLegend(0.6645729,0.6993007,0.8718593,0.8741259,NULL,"brNDC");
	leg1->SetTextSize(0.07);
	leg1->SetLineColor(0);
	leg1->SetFillStyle(0);
	leg1->AddEntry((TObject*)0, tags,"");
	leg1->AddEntry((TObject*)0, bkg,"");

	
	RooPlot* frame1 = BpMass.frame(Title(BintitleGood)); 
	frame1->SetTitleSize(0.07);
 	frame1->GetYaxis()->SetTitle("Events / (3.2 MeV/c^{2})");
 	dataTotal.plotOn(frame1,Cut("MuCharge==MuCharge::Rtag"));
	simTagPDF.plotOn(frame1,Slice(TagCat,"Rtag"),ProjWData(TagCat,dataTotal));
 
	//Plot components of the RtagPDF separately
	simTagPDF.plotOn(frame1,Slice(TagCat,"Rtag"),Components("gaussR"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(6));
	simTagPDF.plotOn(frame1,Slice(TagCat,"Rtag"),Components("gauss2R"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(7));
	simTagPDF.plotOn(frame1,Slice(TagCat,"Rtag"),Components("expoR"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(3));
//	simTagPDF.plotOn(frame1,Slice(TagCat,"Rtag"),Components("erf"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(2));

	frame1->Draw();
	leg1->Draw("same");

	TCanvas* c2 = new TCanvas(CanvBad,CanvBad);
	c2->cd();
	c2->SetLeftMargin(0.1457286);
	c2->SetBottomMargin(0.1541463); //0.1341463
	

	tags = Form("Ntags %d #pm %d",int(nsigW.getVal()),int(nsigW.getError()) );
	bkg = Form("Nbkg %d #pm %d",int(nbkgW.getVal()),int(nbkgW.getError()) );

	TLegend *leg2 = new TLegend(0.459799,0.6608392,0.7876884,0.8811189,NULL,"brNDC");
//new TLegend(0.6168342,0.6503497,0.8668342,0.8706294,NULL,"brNDC");
	
	leg2->SetLineColor(0);
	leg2->SetTextSize(0.07);
	leg2->SetFillStyle(0);
	leg2->AddEntry((TObject*)0, tags,"");
	leg2->AddEntry((TObject*)0, bkg,"");

	RooPlot* frame2 = BpMass.frame(Title(BintitleBad));
	frame2->GetYaxis()->SetTitle("Events / (3.2 MeV/c^{2})");
 	dataTotal.plotOn(frame2,Cut("MuCharge==MuCharge::Wtag"));
	simTagPDF.plotOn(frame2,Slice(TagCat,"Wtag"),ProjWData(TagCat,dataTotal)); 

	//Plot components of the WtagPDF separately
	simTagPDF.plotOn(frame2,Slice(TagCat,"Wtag"),Components("gaussR"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(6));
	simTagPDF.plotOn(frame2,Slice(TagCat,"Wtag"),Components("gauss2R"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(7));
	simTagPDF.plotOn(frame2,Slice(TagCat,"Wtag"),Components("expoR"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(3));
//	simTagPDF.plotOn(frame2,Slice(TagCat,"Wtag"),Components("erf"),ProjWData(TagCat,dataTotal),LineStyle(kDashed),LineColor(2));
	frame2->Draw();
	leg2->Draw("same");

	//Save fit results to see convergence of fits	
//	 TFile *output = new TFile("Convergence.root","recreate");
//	 Result->Write("Convergence");
//    	 output->Close();
	CanvGood = "Binmassplots/"+CanvGood +".pdf";
	CanvBad = "Binmassplots/"+CanvBad +".pdf";

	cout << CanvGood << endl;
	cout << CanvBad << endl;
//	c1->SaveAs(CanvGood);
//	c2->SaveAs(CanvBad);

// Print final value of parameters
meanR.Print();
sigmaR.Print();
//sigma2.Print();
nsigR.Print();
nbkgR.Print();

nsigW.Print();
nbkgW.Print();

sigma2R.Print();
coeff1R->Print();

cout << "bins " << binHigh << " " << binLow  << " histof "<< flag<< endl; 


vector<int> TagVector;
TagVector.push_back(nsigR.getVal()); //number right tags
TagVector.push_back(nsigW.getVal()); //number of wrong tags

//c1->Delete();
//c2->Delete();
return TagVector;
} 
