#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/BsToJpsiPhiRootTree.h"
#include<vector>
#include <iostream>
using namespace std;


BsToJpsiPhiRootTree::BsToJpsiPhiRootTree()
{
  resetEntries();
  bsTree_ = 0;
  bsFile_ = 0;

}


void BsToJpsiPhiRootTree::createTree(const std::string filename)
{

  // open root file
  bsFile_ = new TFile (filename.c_str(), "RECREATE" );
  int bufsize = 256000;
  // create tree structure
  bsTree_ = new TTree("BsTree","BsTree",bufsize);

PtJetTrkPt_= new std::vector<float>(); 
PtJetTrkCharge_= new std::vector<int>();

BJetTrkCharge_ = new std::vector<int>();
BJetTrkPt_ = new std::vector<float>();

PVTrkCharge_ = new std::vector<int>();
PVTrkPt_ = new std::vector<float>();
PVTrkEta_ = new std::vector<float>();
PVTrkPhi_ = new std::vector<float>();

/*
vector<int> &vr = *JetTrkCharge_; //Create a reference
cout << vr[0] << endl;
cout << vr[1] << endl;
cout << vr[2] << endl;
cout << vr[3] << endl;
bsTree_->Branch("JetTrkCharge","vector<int>",&vr);
*/

//bsTree_->Branch("JetTrkCharge","vector< vector<int> >",JetTrkCharge_);
//bsTree_->Branch("JetTrkPt","vector< vector<float> >",JetTrkPt_);

bsTree_->Branch("PtJetTrkCharge","vector<int>",PtJetTrkCharge_);
bsTree_->Branch("PtJetTrkPt","vector<float>",PtJetTrkPt_);

//bsTree_->Branch("BOldJetTrkCharge","vector<int>",BOldJetTrkCharge_);
//bsTree_->Branch("BOldJetTrkPt","vector<float>",BOldJetTrkPt_);


bsTree_->Branch("BJetTrkCharge","vector<int>",BJetTrkCharge_);
bsTree_->Branch("BJetTrkPt","vector<float>",BJetTrkPt_);

bsTree_->Branch("PVTrkCharge","vector<int>",PVTrkCharge_);
bsTree_->Branch("PVTrkPt","vector<float>",PVTrkPt_);
bsTree_->Branch("PVTrkEta","vector<float>",PVTrkEta_);
bsTree_->Branch("PVTrkPhi","vector<float>",PVTrkPhi_);


bsTree_->Branch("TagMuRecoIP3D"      , TagMuRecoIP3D_       ,"TagMuRecoIP3D_[25]/D");
bsTree_->Branch("TagMuRecoIP3DErr"   , TagMuRecoIP3DErr_    ,"TagMuRecoIP3DErr_[25]/D" );

bsTree_->Branch(  "TagMuRecoPt"             , TagMuRecoPt_ ,"TagMuRecoPt_[25]/D");   
bsTree_->Branch(  "TagMuRecoP"             , TagMuRecoP_ ,"TagMuRecoP_[25]/D"); 
bsTree_->Branch(  "TagMuRecoEta"             , TagMuRecoEta_ ,"TagMuRecoEta_[25]/D");
bsTree_->Branch(  "TagMuRecoPhi"             , TagMuRecoPhi_ ,"TagMuRecoPhi_[25]/D");
bsTree_->Branch(  "TagMuRecoChg"             , TagMuRecoChg_ ,"TagMuRecoChg_[25]/D");
bsTree_->Branch(  "TagMuSimu"             , TagMuSimu_ ,"TagMuSimu_[25]/D");
bsTree_->Branch(  "MuRecoMCmother"             , MuRecoMCmother_ ,"MuRecoMCmother_[25]/I");
bsTree_->Branch(  "MuRecoMCMotherMother", MuRecoMCMotherMother_ ,"MuRecoMCMotherMother_[25]/I");
bsTree_->Branch(  "TagMuRecoPtRel", TagMuRecoPtRel_ ,"TagMuRecoPtRel_[25]/D");
bsTree_->Branch(  "TagMuJetdR", TagMuJetdR_ ,"TagMuJetdR_[25]/D");
bsTree_->Branch(  "TagMuPtJetPtRatio", TagMuPtJetPtRatio_ ,"TagMuPtJetPtRatio_[25]/D");


bsTree_->Branch("GlobalTagMuon", GlobalTagMuon_, "GlobalTagMuon_[25]/I");
bsTree_->Branch("TrackerTagMuon", TrackerTagMuon_, "TrackerTagMuon_[25]/I");
bsTree_->Branch("PFTagMuon", PFTagMuon_, "PFTagMuon_[25]/I");


bsTree_->Branch(  "TagMuListSize"             , &TagMuListSize_ ,"TagMuListSize/I");
bsTree_->Branch("BsLxy3DMC", &BsLxy3DMC_, "BsLxy3DMC/D");
bsTree_->Branch("BsLxy2DMC", &BsLxy2DMC_, "BsLxy2DMC/D");
bsTree_->Branch("BsPMC", &BsPMC_, "BsPMC/D");

bsTree_->Branch("BsPVDist2d", &BsPVDist2d_, "BsPVDist2d/D");
bsTree_->Branch("BsPVDist3d", &BsPVDist3d_, "BsPVDist3d/D");

bsTree_->Branch(  "PVZpos"             , PVZpos_ ,"PVZpos_[30]/D");
bsTree_->Branch(  "NTracksInPV"             , NTracksInPV_ ,"NTracksInPV_[30]/I" ); 
bsTree_->Branch(  "PVAbsPt"             , PVAbsPt_ ,"PVAbsPt_[30]/D");     
bsTree_->Branch("SVZpos", SVZpos_, "SVZpos_[30]/D");
bsTree_->Branch("BsPtMC", &BsPtMC_, "BsPtMC/D");
bsTree_->Branch("FirstBsMCZpos", &FirstBsMCZpos_, "FirstBsMCZpos/D");
bsTree_->Branch("BsCosThetaReco", &BsCosThetaReco_, "BsCosThetaReco/D");


bsTree_->Branch("BsCt3DPVCosTheta", &BsCt3DPVCosTheta_, "BsCt3DPVCosTheta/D");
bsTree_->Branch("BsCt2DPVCosTheta", &BsCt2DPVCosTheta_, "BsCt2DPVCosTheta/D");

bsTree_->Branch("BsCt3DPVHighestPt", &BsCt3DPVHighestPt_, "BsCt3DPVHighestPt/D");
bsTree_->Branch("BsCt2DPVHighestPt", &BsCt2DPVHighestPt_, "BsCt2DPVHighestPt/D");

bsTree_->Branch("BsCosThetaReco", &BsCosThetaReco_, "BsCosThetaReco/D");
bsTree_->Branch("PVindex", &PVindex_, "PVindex/I");

// Jet variables
bsTree_->Branch("BsJetPx", &BsJetPx_, "BsJetPx/D");
bsTree_->Branch("BsJetPy", &BsJetPy_, "BsJetPy/D");
bsTree_->Branch("BsJetPz", &BsJetPz_, "BsJetPz/D");
bsTree_->Branch("BsJetdR", &BsJetdR_, "BsJetdR/D");

bsTree_->Branch("BdJetPx", &BsJetPx_, "BsJetPx/D");
bsTree_->Branch("BdJetPy", &BsJetPy_, "BsJetPy/D");
bsTree_->Branch("BdJetPz", &BsJetPz_, "BsJetPz/D");
bsTree_->Branch("BdJetdR", &BdJetdR_, "BdJetdR/D");

bsTree_->Branch("BplusJetPx", &BplusJetPx_,"BplusJetPx/D");
bsTree_->Branch("BplusJetPy", &BplusJetPy_,"BplusJetPy/D");
bsTree_->Branch("BplusJetPz", &BplusJetPz_,"BplusJetPz/D");
bsTree_->Branch("BplusPtJetdR", &BplusPtJetdR_,"BplusPtJetdR/D");
bsTree_->Branch("BplusBJetdR", &BplusBJetdR_,"BplusBJetdR/D");

bsTree_->Branch("BsPtJetdR", &BsPtJetdR_,"BsPtJetdR/D");
bsTree_->Branch("BsBJetdR", &BsBJetdR_,"BsBJetdR/D");
bsTree_->Branch("BsBOldJetdR", &BsBOldJetdR_,"BsBOldJetdR/D");

bsTree_->Branch("BJetParton", &BJetParton_, "BJetParton/I");
bsTree_->Branch("BJetEta", &BJetEta_, "BJetEta/D");
bsTree_->Branch("BJetPhi", &BJetPhi_,"BJetPhi/D");

bsTree_->Branch("PtJetPt", &PtJetPt_,"PtJetPt/D");
bsTree_->Branch("BJetPt", &BJetPt_,"BJetPt/D");
bsTree_->Branch("BOldJetPt", &BOldJetPt_,"BOldJetPt/D");

bsTree_->Branch("JetBTagProb", &JetBTagProb_,"JetBTagProb/D");
bsTree_->Branch("BsJetCharge", &BsJetCharge_ ,"BsJetCharge/D");
bsTree_->Branch("BdJetCharge", &BdJetCharge_ ,"BdJetCharge/D");
bsTree_->Branch("BplusJetCharge", &BplusJetCharge_ ,"BplusJetCharge/D");

//Bplus variables

bsTree_->Branch("BplusIsPV", &BplusIsPV_, "BplusIsPV/I");
bsTree_->Branch("BplusIsBS", &BplusIsBS_, "BplusIsBS/I");
bsTree_->Branch("BplusPVindex", &BplusPVindex_, "BplusPVindex/I");
bsTree_->Branch("BplusVtxProb"	 , 	 &BplusVtxProb_ ,  "BplusVtxProb/D");
bsTree_->Branch("BplusM_fit"	 ,	 &BplusM_fit_   ,  "BplusM_fit/D");
bsTree_->Branch("BplusChi2"      ,       &BplusChi2_    ,  "BplusChi2/D");
bsTree_->Branch("BplusPt"	 ,	 &BplusPt_      ,  "BplusPt/D");
bsTree_->Branch("KplusPt"        ,       &KplusPt_      ,  "KplusPt/D");
bsTree_->Branch("KplusPtot"      ,       &KplusPtot_    ,  "KplusPtot/D");
bsTree_->Branch("BplusMu1Pt"     ,       &BplusMu1Pt_    ,  "BplusMu1Pt/D");	
bsTree_->Branch("BplusMu2Pt"     ,       &BplusMu2Pt_    ,  "BplusMu2Pt/D");
bsTree_->Branch("BplusMu1Ptot"   ,       &BplusMu1Ptot_  ,  "BplusMu1Ptot/D");
bsTree_->Branch("BplusMu2Ptot"   ,       &BplusMu2Ptot_  ,  "BplusMu2Ptot/D");
bsTree_->Branch("BplusLxy"       ,       &BplusLxy_      ,  "BplusLxy/D");
bsTree_->Branch("BplusLxyz"      ,       &BplusLxyz_     ,  "BplusLxyz/D");
bsTree_->Branch("BplusLxyErr"    ,      &BplusLxyErr_    ,  "BplusLxyErr/D");
bsTree_->Branch("BplusLxyzErr"   ,      &BplusLxyzErr_   ,  "BplusLxyzErr/D");	
bsTree_->Branch("BplusCt3D"       ,      &BplusCt3D_     ,  "BplusCt3D/D");
bsTree_->Branch("BplusCtErr3D"    ,      &BplusCtErr3D_  ,  "BplusCtErr3D/D");
bsTree_->Branch("BplusPtot"   ,       &BplusPtot_ ,     "BplusPtot/D" );
bsTree_->Branch("BplusCt3D_2"   ,       &BplusCt3D_2_ , "BplusCt3D_2/D" );
bsTree_->Branch("BplusCt2D_2"   ,       &BplusCt2D_2_ , "BplusCt2D_2/D" );
bsTree_->Branch("BplusCt2D"   ,       &BplusCt2D_ ,     "BplusCt2D/D" );
bsTree_->Branch("BplusCtPerSigma3D"   ,  &BplusCtPerSigma3D_,    "BplusCtPerSigma3D/D" );
bsTree_->Branch("BplusKIP3D"   ,       &BplusKIP3D_ ,     "BplusKIP3D/D" );
bsTree_->Branch("BplusKIP3DErr"   ,       &BplusKIP3DErr_ ,     "BplusKIP3DErr/D" );
bsTree_->Branch("IP3DKandJpsiVtx"  ,       &IP3DKandJpsiVtx_ , "IP3DKandJpsiVtx/D"  );
bsTree_->Branch("IP3DKandJpsiVtxErr"  ,     &IP3DKandJpsiVtxErr_, "IP3DKandJpsiVtxErr/D" );
//HLT trigger matching
bsTree_->Branch("BpmatchDoubleMu01", &BpmatchDoubleMu01_, "BpmatchDoubleMu01/I");
bsTree_->Branch("BpmatchDoubleMu02", &BpmatchDoubleMu02_,"BpmatchDoubleMu02/I");


//R matching 
bsTree_->Branch("BplusKmcId"   ,       &BplusKmcId_  , "BplusKmcId/I"  ); 
bsTree_->Branch("BplusKmomId"   ,       &BplusKmomId_ , "BplusKmomId/I" );


bsTree_->Branch("BplusMu1mcId"   ,       &BplusMu1mcId_  , "BplusMu1mcId/I"  );
bsTree_->Branch("BplusMu1momId"   ,       &BplusMu1momId_ , "BplusMu1momId/I" );
bsTree_->Branch("BplusMu1gmomId" ,       &BplusMu1gmomId_, "BplusMu1gmomId/I");
bsTree_->Branch("BplusMu2mcId"   ,       &BplusMu2mcId_  , "BplusMu2mcId/I"  );
bsTree_->Branch("BplusMu2momId"   ,       &BplusMu2momId_ , "BplusMu2momId/I" );
bsTree_->Branch("BplusMu2gmomId" ,       &BplusMu2gmomId_, "BplusMu2gmomId/I");

bsTree_->Branch("isMatchedBplus" ,       &isMatchedBplus_, "isMatchedBplus/I");
bsTree_->Branch("BplusCtErr2D" ,       &BplusCtErr2D_, "BplusCtErr2D/D");
bsTree_->Branch("BplusCtPerSigma2D" ,       &BplusCtPerSigma2D_, "BplusCtPerSigma2D/D");
bsTree_->Branch("BplusDecayChannel" ,     &BplusDecayChannel_, "BplusDecayChannel/D");
bsTree_->Branch("JpsiMass_bplus" ,     &JpsiMass_bplus_   , "JpsiMass_bplus/D");
bsTree_->Branch("JpsiPt_bplus"   ,     &JpsiPt_bplus_      , "JpsiPt_bplus/D");
bsTree_->Branch("BplusCosTheta"  ,     &BplusCosTheta_     ,  "BplusCosTheta/D");
bsTree_->Branch("LxyPerSigma2D"  ,     &LxyPerSigma2D_     ,  "LxyPerSigma2D/D");
bsTree_->Branch("BplusEta" , &BplusEta_ , "BplusEta/D"); 
bsTree_->Branch("BplusPhi" , &BplusPhi_ , "BplusPhi/D"); 

 
//Bplus variables


// new Bd variables
bsTree_->Branch("BdCt3D" , &BdCt3D_ , "BdCt3D/D"); 
bsTree_->Branch("BdCtErr3D" , &BdCtErr3D_ , "BdCtErr3D/D"); 
bsTree_->Branch("BdCt2D" , &BdCt2D_ , "BdCt2D/D");
bsTree_->Branch("BdCtErr2D" , &BdCtErr2D_ , "BdCtErr2D/D");  
// new Bd variables

bsTree_->Branch("BsCt3DPVCosTheta", &BsCt3DPVCosTheta_, "BsCt3DPVCosTheta/D");
bsTree_->Branch("BsCt2DPVCosTheta", &BsCt2DPVCosTheta_, "BsCt2DPVCosTheta/D");

bsTree_->Branch("BsCt3DPVHighestPt", &BsCt3DPVHighestPt_, "BsCt3DPVHighestPt/D");
bsTree_->Branch("BsCt2DPVHighestPt", &BsCt2DPVHighestPt_, "BsCt2DPVHighestPt/D");

bsTree_->Branch("BsCosThetaReco", &BsCosThetaReco_, "BsCosThetaReco/D");

bsTree_->Branch(  "ihaveajpsi"             , &ihaveajpsi_,                "ihaveajpsi/I");                                         
bsTree_->Branch(  "BsCowboy"             , &BsCowboy_,                "BsCowboy/I");                                         
bsTree_->Branch(  "BdCowboy"             , &BdCowboy_,                "BdCowboy/I");                                         
bsTree_->Branch(  "BsPhiVtxProb"             , &BsPhiVtxProb_,                "BsPhiVtxProb/D");                                         
bsTree_->Branch(  "BsMu1QualityG"             , &BsMu1QualityG_,                "BsMu1QualityG/I");                                         
bsTree_->Branch(  "BsMu2QualityG"             , &BsMu2QualityG_,                "BsMu2QualityG/I");                                         
bsTree_->Branch(  "BsMu1QualityT"             , &BsMu1QualityT_,                "BsMu1QualityT/I");                                         
bsTree_->Branch(  "BsMu2QualityT"             , &BsMu2QualityT_,                "BsMu2QualityT/I");                                         
bsTree_->Branch(  "BdMu1QualityG"             , &BdMu1QualityG_,                "BdMu1QualityG/I");                                         
bsTree_->Branch(  "BdMu2QualityG"             , &BdMu2QualityG_,                "BdMu2QualityG/I");                                         
bsTree_->Branch(  "BdMu1QualityT"             , &BdMu1QualityT_,                "BdMu1QualityT/I");                                         
bsTree_->Branch(  "BdMu2QualityT"             , &BdMu2QualityT_,                "BdMu2QualityT/I");                                         
bsTree_->Branch(  "NVertices"             , &NVertices_,                "NVertices/I");                                         
bsTree_->Branch(  "triggerbit_HLTmu3Tk"             , &triggerbit_HLTmu3Tk_,                "triggerbit_HLTmu3Tk/I");                                         
bsTree_->Branch(  "triggerbit_HLTmu5"		  , &triggerbit_HLTmu5_,                "triggerbit_HLTmu5/I");                                      
bsTree_->Branch(  "triggerbit_HLTmu7"		  , &triggerbit_HLTmu7_,                "triggerbit_HLTmu7/I");
bsTree_->Branch(  "triggerbit_HLTdoubleIsoMu3"	  , &triggerbit_HLTdoubleIsoMu3_,       "triggerbit_HLTdoubleIsoMu3/I");                                
bsTree_->Branch(  "triggerbit_HLTdoubleMu3"	  , &triggerbit_HLTdoubleMu3_,          "triggerbit_HLTdoubleMu3/I");        
bsTree_->Branch(  "triggerbit_HLTdoubleMu0"	  , &triggerbit_HLTdoubleMu0_,          "triggerbit_HLTdoubleMu0/I");                             
bsTree_->Branch(  "triggerbit_HLTL1DoubleMuOpen"  , &triggerbit_HLTL1DoubleMuOpen_,     "triggerbit_HLTL1DoubleMuOpen/I");             
bsTree_->Branch(  "triggerbit_HLTMu0Track0Jpsi"  , &triggerbit_HLTMu0Track0Jpsi_,     "triggerbit_HLTMu0Track0Jpsi/I");             
bsTree_->Branch(  "triggerbit_HLTL1DoubleMuOpenTight"  , &triggerbit_HLTL1DoubleMuOpenTight_,     "triggerbit_HLTL1DoubleMuOpenTight/I");             
bsTree_->Branch(  "triggerbit_HLT_DoubleMu3_Jpsi_v2"  , &triggerbit_HLT_DoubleMu3_Jpsi_v2_,     "triggerbit_HLT_DoubleMu3_Jpsi_v2/I");             
bsTree_->Branch(  "triggerbit_HLT_DoubleMu3_Jpsi_v2MC"  , &triggerbit_HLT_DoubleMu3_Jpsi_v2MC_,     "triggerbit_HLT_DoubleMu3_Jpsi_v2MC/I");             
bsTree_->Branch(  "triggerbit_HLT_DoubleMu3_Quarkonium_v2"  , &triggerbit_HLT_DoubleMu3_Quarkonium_v2_,     "triggerbit_HLT_DoubleMu3_Quarkonium_v2/I");             
bsTree_->Branch(  "triggerbit_HLT_DoubleMu3_Quarkonium_v2MC"  , &triggerbit_HLT_DoubleMu3_Quarkonium_v2MC_,     "triggerbit_HLT_DoubleMu3_Quarkonium_v2MC/I");             
bsTree_->Branch(  "triggerbit_HLT_DoubleMu3_Quarkonium_v1"  , &triggerbit_HLT_DoubleMu3_Quarkonium_v1_,     "triggerbit_HLT_DoubleMu3_Quarkonium_v1/I");             
bsTree_->Branch(  "triggerbit_Jpsi_Displaced_v1"  , &triggerbit_Jpsi_Displaced_v1_,     "triggerbit_Jpsi_Displaced_v1/I");             
bsTree_->Branch(  "triggerbit_7Jpsi_Displaced_v1"  , &triggerbit_7Jpsi_Displaced_v1_,     "triggerbit_7Jpsi_Displaced_v1/I");             
bsTree_->Branch(  "triggerbit_Jpsi_Displaced_v1MC"  , &triggerbit_Jpsi_Displaced_v1MC_,     "triggerbit_Jpsi_Displaced_v1MC/I");             
bsTree_->Branch(  "triggerbit_7Jpsi_Displaced_v2"  , &triggerbit_7Jpsi_Displaced_v2_,     "triggerbit_7Jpsi_Displaced_v2/I");             
bsTree_->Branch(  "triggerbit_7Jpsi_Displaced_v3"  , &triggerbit_7Jpsi_Displaced_v3_,     "triggerbit_7Jpsi_Displaced_v3/I");             
bsTree_->Branch(  "triggerbit_3p5Jpsi_Displaced_v2"  , &triggerbit_3p5Jpsi_Displaced_v2_,     "triggerbit_3p5Jpsi_Displaced_v2/I");             
bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v1"  , &triggerbit_4Jpsi_Displaced_v1_,     "triggerbit_4Jpsi_Displaced_v1/I");             
bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v4"  , &triggerbit_4Jpsi_Displaced_v4_,     "triggerbit_4Jpsi_Displaced_v4/I");             

bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v5"  , &triggerbit_4Jpsi_Displaced_v5_,     "triggerbit_4Jpsi_Displaced_v5/I");

bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v9"  , &triggerbit_4Jpsi_Displaced_v9_,     "triggerbit_4Jpsi_Displaced_v9/I");
        
bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v10" , &triggerbit_4Jpsi_Displaced_v10_,     "triggerbit_4Jpsi_Displaced_v10/I");

bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v11" , &triggerbit_4Jpsi_Displaced_v11_,     "triggerbit_4Jpsi_Displaced_v11/I");

bsTree_->Branch(  "triggerbit_4Jpsi_Displaced_v12" , &triggerbit_4Jpsi_Displaced_v12_,     "triggerbit_4Jpsi_Displaced_v12/I");

     
bsTree_->Branch(  "triggerbit_5Jpsi_Displaced_v1"  , &triggerbit_5Jpsi_Displaced_v1_,     "triggerbit_5Jpsi_Displaced_v1/I");             
bsTree_->Branch(  "triggerbit_5Jpsi_Displaced_v2"  , &triggerbit_5Jpsi_Displaced_v2_,     "triggerbit_5Jpsi_Displaced_v2/I");             
bsTree_->Branch(  "triggerbit_5Jpsi_Displaced_v4"  , &triggerbit_5Jpsi_Displaced_v4_,     "triggerbit_5Jpsi_Displaced_v4/I");             
bsTree_->Branch(  "triggerbit_5Jpsi_Displaced_v5"  , &triggerbit_5Jpsi_Displaced_v5_,     "triggerbit_5Jpsi_Displaced_v5/I");


bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_Muon_v15"	  , &triggerbit_Dimuon0_Jpsi_Muon_v15_,     "triggerbit_Dimuon0_Jpsi_Muon_v15/I"  ); 
      
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_Muon_v16"	  , &triggerbit_Dimuon0_Jpsi_Muon_v16_,"triggerbit_Dimuon0_Jpsi_Muon_v16/I"); 

 
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_Muon_v17"	  , &triggerbit_Dimuon0_Jpsi_Muon_v17_ ,"triggerbit_Dimuon0_Jpsi_Muon_v17/I"); 

bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_Muon_v18"	  , &triggerbit_Dimuon0_Jpsi_Muon_v18_, "triggerbit_Dimuon0_Jpsi_Muon_v18/I"); 
       
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v1"	  , &triggerbit_Dimuon0_Jpsi_v1_,     "triggerbit_Dimuon0_Jpsi_v1/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v3"	  , &triggerbit_Dimuon0_Jpsi_v3_,     "triggerbit_Dimuon0_Jpsi_v3/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v5"	  , &triggerbit_Dimuon0_Jpsi_v5_,     "triggerbit_Dimuon0_Jpsi_v5/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v6"	  , &triggerbit_Dimuon0_Jpsi_v6_,     "triggerbit_Dimuon0_Jpsi_v6/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v9"	  , &triggerbit_Dimuon0_Jpsi_v9_,     "triggerbit_Dimuon0_Jpsi_v9/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon0_Jpsi_v10"	  , &triggerbit_Dimuon0_Jpsi_v10_,     "triggerbit_Dimuon0_Jpsi_v10/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon10_Barrel"	  , &triggerbit_Dimuon10_Barrel_,     "triggerbit_Dimuon10_Barrel/I"  );                        
bsTree_->Branch(  "triggerbit_Dimuon13_Barrel"	  , &triggerbit_Dimuon13_Barrel_,     "triggerbit_Dimuon13_Barrel/I"  );                        
bsTree_->Branch(  "BSx"				  , &BSx_,                              "BSx/D");                                                     
bsTree_->Branch(  "BSy"				  , &BSy_,                              "BSy/D");                                                    
bsTree_->Branch(  "BSz"				  , &BSz_,                              "BSz/D");                                                       
bsTree_->Branch(  "BSdx"                           , &BSdx_,                              "BSdx/D");
bsTree_->Branch(  "BSdy"                           , &BSdy_,                              "BSdy/D");
bsTree_->Branch(  "BSdz"                           , &BSdz_,                              "BSdz/D");
bsTree_->Branch(  "BSsigmaZ"                           , &BSsigmaZ_,                              "BSsigmaZ/D");
bsTree_->Branch(  "BSdsigmaZ"                           , &BSdsigmaZ_,                              "BSdsigmaZ/D");
bsTree_->Branch(  "PVx"				  , &PVx_,                              "PVx/D");                                                       
bsTree_->Branch(  "PVy"				  , &PVy_,                              "PVy/D");                                                       
bsTree_->Branch(  "PVz"				  , &PVz_,                              "PVz/D");                                                       
bsTree_->Branch(  "PVerrx"			  , &PVerrx_,                           "PVerrx/D");                                                    
bsTree_->Branch(  "PVerry"			  , &PVerry_,                           "PVerry/D");                                                    
bsTree_->Branch(  "PVerrz"			  , &PVerrz_,                           "PVerrz/D"); 

bsTree_->Branch(  "MuMuDCA"			  , &MuMuDCA_,                           "MuMuDCA/D"); 
bsTree_->Branch(  "MuMuDistance"			  , &MuMuDistance_,                           "MuMuDistance/D"); 
bsTree_->Branch(  "MuMuDistanceSigma"			  , &MuMuDistanceSigma_,                           "MuMuDistanceSigma/D"); 
bsTree_->Branch(  "MuDr1"			  , &MuDr1_,                           "MuDr1/D"); 
bsTree_->Branch(  "MuDr2"			  , &MuDr2_,                           "MuDr2/D"); 
bsTree_->Branch(  "BdMuMuDCA"			  , &BdMuMuDCA_,                           "BdMuMuDCA/D"); 
bsTree_->Branch(  "BdMuMuDistance"			  , &BdMuMuDistance_,                           "BdMuMuDistance/D"); 
bsTree_->Branch(  "BdMuMuDistanceSigma"			  , &BdMuMuDistanceSigma_,                           "BdMuMuDistanceSigma/D"); 
bsTree_->Branch(  "BdMuDr1"			  , &BdMuDr1_,                           "BdMuDr1/D"); 
bsTree_->Branch(  "BdMuDr2"			  , &BdMuDr2_,                           "BdMuDr2/D"); 

bsTree_->Branch( "PVx_refit"   ,         &PVx_refit_   ,         "PVx_refit/D"   ); 
bsTree_->Branch( "PVy_refit"   ,	 &PVy_refit_   ,	 "PVy_refit/D"   );
bsTree_->Branch( "PVz_refit"   ,	 &PVz_refit_   ,	 "PVz_refit/D"   );
bsTree_->Branch( "PVerrx_refit",	 &PVerrx_refit_,	 "PVerrx_refit/D");
bsTree_->Branch( "PVerry_refit",	 &PVerry_refit_,	 "PVerry_refit/D");
bsTree_->Branch( "PVerrz_refit",	 &PVerrz_refit_,	 "PVerrz_refit/D");

bsTree_->Branch( "isPV"   ,         &isPV_   ,         "isPV/I"   ); 
bsTree_->Branch( "isBS"   ,         &isBS_   ,         "isBS/I"   ); 

bsTree_->Branch( "runNumber"   ,         &runNumber_   ,         "runNumber/I"   ); 
bsTree_->Branch( "eventNumber"   ,       &eventNumber_   ,       "eventNumber/i"   );
bsTree_->Branch( "lumiSection"   ,       &lumiSection_   ,       "lumiSection/I"   ); 

bsTree_->Branch( "PUinteraction"   ,         &PUinteraction_   ,         "PUinteraction/I"   ); 

bsTree_->Branch( "PionDeDx"                   ,&PionDeDx_ ,                      "PionDeDx_[25]/D");                         
bsTree_->Branch( "KaonDeDx"                   ,&KaonDeDx_ ,                      "KaonDeDx_[25]/D");                         
bsTree_->Branch( "KaonPt"                   ,&KaonPt_ ,                      "KaonPt_[25]/D");                         
bsTree_->Branch( "PionPt"                   ,&PionPt_ ,                      "PionPt_[25]/D");                         
              
// bsTree_->Branch( "MuRecoPt2"                   ,&MuRecoPt2_ ,                      "MuRecoPt2[15]/D");       
// bsTree_->Branch( "MuRecoChg2"                   ,&MuRecoChg2_ ,                      "MuRecoChg2[15]/I");                            
// bsTree_->Branch( "MuRecoEta2"                   ,&MuRecoEta2_ ,                      "MuRecoEta2[15]/D");                            
// bsTree_->Branch( "MuRecoPhi2"                   ,&MuRecoPhi2_ ,                      "MuRecoPhi2[15]/D");                                    

bsTree_->Branch( "matchDoubleMu01DiMuon0_",&matchDoubleMu01DiMuon0_ ,                      "matchDoubleMu01DiMuon0_/I");            

bsTree_->Branch( "matchDoubleMu02DiMuon0_",&matchDoubleMu02DiMuon0_ ,                      "matchDoubleMu02DiMuon0_/I");            
         
bsTree_->Branch( "BpmatchDoubleMu01DiMuon0_",&BpmatchDoubleMu01DiMuon0_ ,                      "matchDoubleMu01DiMuon0_/I");            

bsTree_->Branch( "BpmatchDoubleMu02DiMuon0_",&BpmatchDoubleMu02DiMuon0_ ,                      "BpmatchDoubleMu02DiMuon0_/I");   

bsTree_->Branch( "BdmatchDoubleMu01DiMuon0_",&BdmatchDoubleMu01DiMuon0_ ,                      "matchDoubleMu01DiMuon0_/I");            

bsTree_->Branch( "BdmatchDoubleMu02DiMuon0_",&BdmatchDoubleMu02DiMuon0_ ,                      "BdmatchDoubleMu02DiMuon0_/I");   



bsTree_->Branch( "matchL11"                   ,&matchL11_,                      "matchL11/I");                                                   
bsTree_->Branch( "matchL12"                   ,&matchL12_,                      "matchL12/I");                                                   
bsTree_->Branch( "match2mu01"                   ,&match2mu01_,                      "match2mu01/I");                                                   
bsTree_->Branch( "match2mu02"                   ,&match2mu02_,                      "match2mu02/I");                                                   
bsTree_->Branch( "match2mu31"                   ,&match2mu31_,                      "match2mu31/I");                                                   
bsTree_->Branch( "match2mu32"                   ,&match2mu32_,                      "match2mu32/I");                                                   
bsTree_->Branch( "match1mu01"                   ,&match1mu01_,                      "match1mu01/I");                                                   
bsTree_->Branch( "match1mu02"                   ,&match1mu02_,                      "match1mu02/I");                                                   
bsTree_->Branch( "matchDoubleMu31J"                   ,&matchDoubleMu31J_,                      "matchDoubleMu31J/I");                                                   
bsTree_->Branch( "matchDoubleMu32J"                   ,&matchDoubleMu32J_,                      "matchDoubleMu32J/I");                                                   
bsTree_->Branch( "matchDoubleMu31Q"                   ,&matchDoubleMu31Q_,                      "matchDoubleMu31Q/I");                                                   
bsTree_->Branch( "matchDoubleMu32Q"                   ,&matchDoubleMu32Q_,                      "matchDoubleMu32Q/I");                                                   
bsTree_->Branch( "matchDoubleMu71"                   ,&matchDoubleMu71_,                      "matchDoubleMu71/I");                                                   
bsTree_->Branch( "matchDoubleMu72"                   ,&matchDoubleMu72_,                      "matchDoubleMu72/I");                                                   
bsTree_->Branch( "matchDoubleMu41"                   ,&matchDoubleMu41_,                      "matchDoubleMu41/I");                                                   
bsTree_->Branch( "matchDoubleMu42"                   ,&matchDoubleMu42_,                      "matchDoubleMu42/I");                                                   
bsTree_->Branch( "matchDoubleMu51"                   ,&matchDoubleMu51_,                      "matchDoubleMu51/I");                                                   
bsTree_->Branch( "matchDoubleMu52"                   ,&matchDoubleMu52_,                      "matchDoubleMu52/I");                                                   
bsTree_->Branch( "matchDoubleMu01"                   ,&matchDoubleMu01_,                      "matchDoubleMu01/I");                                                   
bsTree_->Branch( "matchDoubleMu02"                   ,&matchDoubleMu02_,                      "matchDoubleMu02/I");                                                   
bsTree_->Branch( "matchDoubleMu101"                   ,&matchDoubleMu101_,                      "matchDoubleMu101/I");                                                   
bsTree_->Branch( "matchDoubleMu102"                   ,&matchDoubleMu102_,                      "matchDoubleMu102/I");                                                   
bsTree_->Branch( "matchDoubleMu131"                   ,&matchDoubleMu131_,                      "matchDoubleMu131/I");                                                   
bsTree_->Branch( "matchDoubleMu132"                   ,&matchDoubleMu132_,                      "matchDoubleMu132/I");                                                   
bsTree_->Branch( "matchmu0tk01"                   ,&matchmu0tk01_,                      "matchmu0tk01/I");                                                   
bsTree_->Branch( "matchmu0tk02"                   ,&matchmu0tk02_,                      "matchmu0tk02/I");                                                   
bsTree_->Branch( "matchFilterJpsi1"                   ,&matchFilterJpsi1_,                      "matchFilterJpsi1/I");                                                   
bsTree_->Branch( "matchFilterJpsi2"                   ,&matchFilterJpsi2_,                      "matchFilterJpsi2/I");                                                   
bsTree_->Branch( "BdmatchL11"                   ,&BdmatchL11_,                      "BdmatchL11/I");                                                   
bsTree_->Branch( "BdmatchL12"                   ,&BdmatchL12_,                      "BdmatchL12/I");                                                   
bsTree_->Branch( "Bdmatch2mu01"                   ,&Bdmatch2mu01_,                      "Bdmatch2mu01/I");                                                   
bsTree_->Branch( "Bdmatch2mu02"                   ,&Bdmatch2mu02_,                      "Bdmatch2mu02/I");                                                   
bsTree_->Branch( "Bdmatch2mu31"                   ,&Bdmatch2mu31_,                      "Bdmatch2mu31/I");                                                   
bsTree_->Branch( "Bdmatch2mu32"                   ,&Bdmatch2mu32_,                      "Bdmatch2mu32/I");                                                   
bsTree_->Branch( "Bdmatch1mu01"                   ,&Bdmatch1mu01_,                      "Bdmatch1mu01/I");                                                   
bsTree_->Branch( "Bdmatch1mu02"                   ,&Bdmatch1mu02_,                      "Bdmatch1mu02/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu31Q"                   ,&BdmatchDoubleMu31Q_,                      "BdmatchDoubleMu31Q/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu32Q"                   ,&BdmatchDoubleMu32Q_,                      "BdmatchDoubleMu32Q/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu71"                   ,&BdmatchDoubleMu71_,                      "BdmatchDoubleMu71/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu72"                   ,&BdmatchDoubleMu72_,                      "BdmatchDoubleMu72/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu41"                   ,&BdmatchDoubleMu41_,                      "BdmatchDoubleMu41/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu42"                   ,&BdmatchDoubleMu42_,                      "BdmatchDoubleMu42/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu51"                   ,&BdmatchDoubleMu51_,                      "BdmatchDoubleMu51/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu52"                   ,&BdmatchDoubleMu52_,                      "BdmatchDoubleMu52/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu01"                   ,&BdmatchDoubleMu01_,                      "BdmatchDoubleMu01/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu02"                   ,&BdmatchDoubleMu02_,                      "BdmatchDoubleMu02/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu101"                   ,&BdmatchDoubleMu101_,                      "BdmatchDoubleMu101/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu102"                   ,&BdmatchDoubleMu102_,                      "BdmatchDoubleMu102/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu131"                   ,&BdmatchDoubleMu131_,                      "BdmatchDoubleMu131/I");                                                   
bsTree_->Branch( "BdmatchDoubleMu132"                   ,&BdmatchDoubleMu132_,                      "BdmatchDoubleMu132/I");                                                   

                                
bsTree_->Branch( "MuonType", &MuonType_, "MuonType/I");

bsTree_->Branch(  "JpsiVtxProb"			  , &JpsiVtxProb_,                      "JpsiVtxProb/D");                                               
bsTree_->Branch(  "BdJpsiVtxProb"			  , &BdJpsiVtxProb_,                      "BdJpsiVtxProb/D");                                               
bsTree_->Branch(  "JpsiM_alone"			  , &JpsiM_alone_,                      "JpsiM_alone/D");                                               
bsTree_->Branch(  "JpsiPhi_alone"		  , &JpsiPhi_alone_,                    "JpsiPhi_alone/D");                                             
bsTree_->Branch(  "JpsiEta_alone"		  , &JpsiEta_alone_,                    "JpsiEta_alone/D");                                             
bsTree_->Branch(  "JpsiPt_alone"		  , &JpsiPt_alone_,                     "JpsiPt_alone/D");                                              
bsTree_->Branch(  "JpsiMu1Pt_alone"		  , &JpsiMu1Pt_alone_,                  "JpsiMu1Pt_alone/D");                                           
bsTree_->Branch(  "JpsiMu2Pt_alone"		  , &JpsiMu2Pt_alone_,                  "JpsiMu2Pt_alone/D");                                           
bsTree_->Branch(  "JpsiMu1Phi_alone"		  , &JpsiMu1Phi_alone_,                 "JpsiMu1Phi_alone/D");                                          
bsTree_->Branch(  "JpsiMu2Phi_alone"		  , &JpsiMu2Phi_alone_,                 "JpsiMu2Phi_alone/D");                                          
bsTree_->Branch(  "JpsiMu1Eta_alone"		  , &JpsiMu1Eta_alone_,                 "JpsiMu1Eta_alone/D");                                          
bsTree_->Branch(  "JpsiMu2Eta_alone"		  , &JpsiMu2Eta_alone_,                 "JpsiMu2Eta_alone/D");                                          
bsTree_->Branch(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_,               "JpsiMuon1Cat_alone/I");                                     
bsTree_->Branch(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_,               "JpsiMuon2Cat_alone/I");                                     
bsTree_->Branch(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_,               "BdJpsiMuon1Cat_alone/I");                                     
bsTree_->Branch(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_,               "BdJpsiMuon2Cat_alone/I");                                     
bsTree_->Branch(  "JpsiMu1d0_alone"               , &JpsiMu1d0_alone_,                  "JpsiMu1d0_alone/D");
bsTree_->Branch(  "JpsiMu2d0_alone"               , &JpsiMu2d0_alone_,                  "JpsiMu2d0_alone/D");
bsTree_->Branch(  "JpsiMu1dz_alone"               , &JpsiMu1dz_alone_,                  "JpsiMu1dz_alone/D");
bsTree_->Branch(  "JpsiMu2dz_alone"               , &JpsiMu2dz_alone_,                  "JpsiMu2dz_alone/D");
bsTree_->Branch(  "JpsiMu1chi2_alone"             , &JpsiMu1chi2_alone_,                  "JpsiMu1chi2_alone/D");
bsTree_->Branch(  "JpsiMu2chi2_alone"             , &JpsiMu2chi2_alone_,                  "JpsiMu2chi2_alone/D");
bsTree_->Branch(  "JpsiMu1ndof_alone"             , &JpsiMu1ndof_alone_,                  "JpsiMu1ndof_alone/I");
bsTree_->Branch(  "JpsiMu2ndof_alone"             , &JpsiMu2ndof_alone_,                  "JpsiMu2ndof_alone/I");
bsTree_->Branch(  "JpsiMu1nHits_alone"            , &JpsiMu1nHits_alone_,                  "JpsiMu1nHits_alone/I");
bsTree_->Branch(  "JpsiMu2nHits_alone"            , &JpsiMu2nHits_alone_,                  "JpsiMu2nHits_alone/I");
bsTree_->Branch(  "JpsiMu1nPixHits_alone"         , &JpsiMu1nPixHits_alone_,               "JpsiMu1nPixHits_alone/I");
bsTree_->Branch(  "JpsiMu2nPixHits_alone"         , &JpsiMu2nPixHits_alone_,               "JpsiMu2nPixHits_alone/I");

bsTree_->Branch(  "K1Pt_beffit"			  , &K1Pt_beffit_,                       "K1Pt_beffit/D");                                                
bsTree_->Branch(  "K1Pz_beffit"			  , &K1Pz_beffit_,                       "K1Pz_beffit/D");                                                
bsTree_->Branch(  "K1Eta_beffit"			  , &K1Eta_beffit_,                      "K1Eta_beffit/D");                
bsTree_->Branch(  "K1Phi_beffit"			  , &K1Phi_beffit_,                      "K1Phi_beffit/D");     
                                         
bsTree_->Branch(  "K2Pt_beffit"			  , &K2Pt_beffit_,                       "K2Pt_beffit/D");                                                
bsTree_->Branch(  "K2Pz_beffit"			  , &K2Pz_beffit_,                       "K2Pz_beffit/D");                                                
bsTree_->Branch(  "K2Eta_beffit"			  , &K2Eta_beffit_,                      "K2Eta_beffit/D");                                  
bsTree_->Branch(  "K2Phi_beffit"			  , &K2Phi_beffit_,                      "K2Phi_beffit/D");                                               

bsTree_->Branch(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_,                       "Mu1Pt_beffit/D");                                                
bsTree_->Branch(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_,                       "Mu1Pz_beffit/D");                                                
bsTree_->Branch(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_,                      "Mu1Eta_beffit/D");                
bsTree_->Branch(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_,                      "Mu1Phi_beffit/D");     
bsTree_->Branch(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_,                       "Mu2Pt_beffit/D");                                                
bsTree_->Branch(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_,                       "Mu2Pz_beffit/D");                                                
bsTree_->Branch(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_,                      "Mu2Eta_beffit/D");                                  
bsTree_->Branch(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_,                      "Mu2Phi_beffit/D");                                               

bsTree_->Branch(  "BdMu1Pt_beffit"			  , &BdMu1Pt_beffit_,                       "BdMu1Pt_beffit/D");                                                
bsTree_->Branch(  "BdMu1Pz_beffit"			  , &BdMu1Pz_beffit_,                       "BdMu1Pz_beffit/D");                                                
bsTree_->Branch(  "BdMu1Eta_beffit"			  , &BdMu1Eta_beffit_,                      "BdMu1Eta_beffit/D");                
bsTree_->Branch(  "BdMu1Phi_beffit"			  , &BdMu1Phi_beffit_,                      "BdMu1Phi_beffit/D");     
bsTree_->Branch(  "BdMu2Pt_beffit"			  , &BdMu2Pt_beffit_,                       "BdMu2Pt_beffit/D");                                                
bsTree_->Branch(  "BdMu2Pz_beffit"			  , &BdMu2Pz_beffit_,                       "BdMu2Pz_beffit/D");                                                
bsTree_->Branch(  "BdMu2Eta_beffit"			  , &BdMu2Eta_beffit_,                      "BdMu2Eta_beffit/D");                                  
bsTree_->Branch(  "BdMu2Phi_beffit"			  , &BdMu2Phi_beffit_,                      "BdMu2Phi_beffit/D");                                               

bsTree_->Branch(  "BsFitChi2"			  , &BsFitChi2_,                        "BsFitChi2/D");                                                 
bsTree_->Branch(  "BsFitNdof"			  , &BsFitNdof_,                        "BsFitNdof/I");                                                
bsTree_->Branch(  "BsFitVtxProb"		  , &BsFitVtxProb_,                     "BsFitVtxProb/D");                                              
        
bsTree_->Branch(  "BsFitM"			  , &BsFitM_,                           "BsFitM/D");                                                 
bsTree_->Branch(  "K1Pt_fit"			  , &K1Pt_fit_,                           "K1Pt_fit/D");                                                 
bsTree_->Branch(  "K2Pt_fit"			  , &K2Pt_fit_,                           "K2Pt_fit/D");                                                 
bsTree_->Branch(  "PhiM_fit"			  , &PhiM_fit_,                           "PhiM_fit/D");                                                 
bsTree_->Branch(  "BsFitEta"			  , &BsFitEta_,                         "BsFitEta/D");                                                  
bsTree_->Branch(  "BsFitPt"			  , &BsFitPt_,                          "BsFitPt/D");                                                   
bsTree_->Branch(  "BsFitPz"			  , &BsFitPz_,                          "BsFitPz/D");                                                   
bsTree_->Branch(  "BsFitPhi"			  , &BsFitPhi_,                         "BsFitPhi/D");                                                  
bsTree_->Branch(  "BsFitVtx_x"			  , &BsFitVtx_x_,                       "BsFitVtx_x/D");                                                
bsTree_->Branch(  "BsFitVtx_y"			  , &BsFitVtx_y_,                       "BsFitVtx_y/D");                                                
bsTree_->Branch(  "BsFitVtx_z"			  , &BsFitVtx_z_,                       "BsFitVtx_z/D");                                                
bsTree_->Branch(  "BsM_nofit"			  , &BsM_nofit_,                        "BsM_nofit/D");                                                 
bsTree_->Branch(  "BsPhi_nofit"			  , &BsPhi_nofit_,                      "BsPhi_nofit/D");                                               
bsTree_->Branch(  "BsEta_nofit"			  , &BsEta_nofit_,                      "BsEta_nofit/D");                                               
bsTree_->Branch(  "BsPt_nofit"			  , &BsPt_nofit_,                       "BsPt_nofit/D");                                                
bsTree_->Branch(  "BsPz_nofit"			  , &BsPz_nofit_,                       "BsPz_nofit/D");                                                
bsTree_->Branch(  "JpsiM_nofit"			  , &JpsiM_nofit_,                      "JpsiM_nofit/D");                                               
bsTree_->Branch(  "JpsiPhi_nofit"		  , &JpsiPhi_nofit_,                    "JpsiPhi_nofit/D");                                             
bsTree_->Branch(  "JpsiEta_nofit"		  , &JpsiEta_nofit_,                    "JpsiEta_nofit/D");                                             
bsTree_->Branch(  "JpsiPt_nofit"		  , &JpsiPt_nofit_,                     "JpsiPt_nofit/D");                                              
bsTree_->Branch(  "JpsiPz_nofit"		  , &JpsiPz_nofit_,                     "JpsiPz_nofit/D");                                              
bsTree_->Branch(  "BdJpsiM_nofit"		  , &BdJpsiM_nofit_,                      "BdJpsiM_nofit/D");                                               
bsTree_->Branch(  "BdJpsiPhi_nofit"		  , &BdJpsiPhi_nofit_,                    "BdJpsiPhi_nofit/D");                                             
bsTree_->Branch(  "BdJpsiEta_nofit"		  , &BdJpsiEta_nofit_,                    "BdJpsiEta_nofit/D");                                             
bsTree_->Branch(  "BdJpsiPt_nofit"		  , &BdJpsiPt_nofit_,                     "BdJpsiPt_nofit/D");                                              
bsTree_->Branch(  "BdJpsiPz_nofit"		  , &BdJpsiPz_nofit_,                     "BdJpsiPz_nofit/D");                                              
bsTree_->Branch(  "PhiM_nofit"			  , &PhiM_nofit_,                       "PhiM_nofit/D");                                                
bsTree_->Branch(  "PhiPhi_nofit"		  , &PhiPhi_nofit_,                     "PhiPhi_nofit/D");                                              
bsTree_->Branch(  "PhiEta_nofit"		  , &PhiEta_nofit_,                     "PhiEta_nofit/D");                                              
bsTree_->Branch(  "PhiPt_nofit"			  , &PhiPt_nofit_,                      "PhiPt_nofit/D");                                               
bsTree_->Branch(  "PhiPz_nofit"			  , &PhiPz_nofit_,                      "PhiPz_nofit/D");                                               
bsTree_->Branch(  "K1Pt_nofit"			  , &K1Pt_nofit_,                       "K1Pt_nofit/D");                                                
bsTree_->Branch(  "K1Pz_nofit"			  , &K1Pz_nofit_,                       "K1Pz_nofit/D");                                                
bsTree_->Branch(  "K1Eta_nofit"			  , &K1Eta_nofit_,                      "K1Eta_nofit/D");                                               
bsTree_->Branch(  "K1Phi_nofit"			  , &K1Phi_nofit_,                      "K1Phi_nofit/D");                                               
bsTree_->Branch(  "K1Key_nofit"			  , &K1Key_nofit_,                       "K1Key_nofit/I");                                                
bsTree_->Branch(  "K2Eta_nofit"			  , &K2Eta_nofit_,                      "K2Eta_nofit/D");                                               
bsTree_->Branch(  "K2Pt_nofit"			  , &K2Pt_nofit_,                       "K2Pt_nofit/D");                                                
bsTree_->Branch(  "K2Pz_nofit"			  , &K2Pz_nofit_,                       "K2Pz_nofit/D");                                                
bsTree_->Branch(  "K2Phi_nofit"			  , &K2Phi_nofit_,                      "K2Phi_nofit/D"); 
bsTree_->Branch(  "K2Key_nofit"			  , &K2Key_nofit_,                       "K2Key_nofit/I");                                              
bsTree_->Branch(  "K1Chi2"			  , &K1Chi2_,                           "K1Chi2/D");                                                    
bsTree_->Branch(  "K1nHits"			  , &K1nHits_,                          "K1nHits/I");                                                
bsTree_->Branch(  "K2Chi2"			  , &K2Chi2_,                           "K2Chi2/D");                                                    
bsTree_->Branch(  "K2nHits"			  , &K2nHits_,                          "K2nHits/I");                                                   
bsTree_->Branch(  "K1pixH"			  , &K1pixH_,                           "K1pixH/I");                                                  
bsTree_->Branch(  "K1trkH"			  , &K1trkH_,                           "K1trkH/I");                                                 
bsTree_->Branch(  "K2pixH"			  , &K2pixH_,                           "K2pixH/I");                                                 
bsTree_->Branch(  "K2trkH"			  , &K2trkH_,                           "K2trkH/I");                                                 
bsTree_->Branch(  "Mu1Chi2"			  , &Mu1Chi2_,                          "Mu1Chi2/D");                                                   
bsTree_->Branch(  "Mu1nHits"			  , &Mu1nHits_,                         "Mu1nHits/I");                                                 
bsTree_->Branch(  "Mu2Chi2"			  , &Mu2Chi2_,                          "Mu2Chi2/D");                                                   
bsTree_->Branch(  "Mu2nHits"			  , &Mu2nHits_,                         "Mu2nHits/I");                                               
bsTree_->Branch(  "Mu1pixH"			  , &Mu1pixH_,                          "Mu1pixH/I");                                               
bsTree_->Branch(  "Mu1trkH"			  , &Mu1trkH_,                          "Mu1trkH/I");                                                
bsTree_->Branch(  "Mu2pixH"			  , &Mu2pixH_,                          "Mu2pixH/I");                                                
bsTree_->Branch(  "Mu2trkH"			  , &Mu2trkH_,                          "Mu2trkH/I");                                                

bsTree_->Branch(  "Mu1d0"               , &Mu1d0_,                  "Mu1d0/D");
bsTree_->Branch(  "Mu2d0"               , &Mu2d0_,                  "Mu2d0/D");
bsTree_->Branch(  "Mu1dz"               , &Mu1dz_,                  "Mu1dz/D");
bsTree_->Branch(  "Mu2dz"               , &Mu2dz_,                  "Mu2dz/D");

bsTree_->Branch(  "costheta"			  , &costheta_,                         "costheta/D");                                                  
bsTree_->Branch(  "phi"				  , &phi_,                              "phi/D");                                                       
bsTree_->Branch(  "cospsi"			  , &cospsi_,                           "cospsi/D");                                                    
bsTree_->Branch(  "Bdcostheta"			  , &Bdcostheta_,                         "Bdcostheta/D");                                                  
bsTree_->Branch(  "Bdphi"		          , &Bdphi_,                              "Bdphi/D");                                                       
bsTree_->Branch(  "Bdcospsi"			  , &Bdcospsi_,                           "Bdcospsi/D");                                                    
bsTree_->Branch(  "BdcosthetaMC"		  , &BdcosthetaMC_,                         "BdcosthetaMC/D");                                                  
bsTree_->Branch(  "BdphiMC"		          , &BdphiMC_,                              "BdphiMC/D");                                                       
bsTree_->Branch(  "BdcospsiMC"			  , &BdcospsiMC_,                           "BdcospsiMC/D");                                                    

bsTree_->Branch(  "AngleBsDecayLength"		  , &AngleBsDecayLength_,               "AngleBsDecayLength/D");                                        
bsTree_->Branch(  "CosDeltaAlpha"		  , &CosDeltaAlpha_,               "CosDeltaAlpha/D");                                        
bsTree_->Branch(  "BdCosDeltaAlpha"		  , &BdCosDeltaAlpha_,               "BdCosDeltaAlpha/D");                                        
bsTree_->Branch(  "isMatched"			  , &isMatched_,                        "isMatched/I");                                                
bsTree_->Branch(  "isMatchedBd"			  , &isMatchedBd_,                      "isMatchedBd/I");                                              
bsTree_->Branch(  "BsLxy"			  , &BsLxy_,                            "BsLxy/D");                                                     
bsTree_->Branch(  "BsCt"                          , &BsCt_,                             "BsCt/D");
bsTree_->Branch(  "BsCtErr"                       , &BsCtErr_,                             "BsCtErr/D");
bsTree_->Branch(  "BsCt3D"                        , &BsCt3D_,                             "BsCt3D/D");
bsTree_->Branch(  "BsCt2D"                        , &BsCt2D_,                             "BsCt2D/D");
bsTree_->Branch(  "BsCt2DBS"                        , &BsCt2DBS_,                             "BsCt2DBS/D");
bsTree_->Branch(  "BdCt2DBS"                        , &BdCt2DBS_,                             "BdCt2DBS/D");
bsTree_->Branch(  "BdCt2DMC"                        , &BdCt2DMC_,                             "BdCt2DMC/D");
bsTree_->Branch(  "BdCt3DMC"                        , &BdCt3DMC_,                             "BdCt3DMC/D");
bsTree_->Branch(  "BsCtMPV"                       , &BsCtMPV_,                             "BsCtMPV/D");
bsTree_->Branch(  "BsCtErr3D"                     , &BsCtErr3D_,                             "BsCtErr3D/D");
bsTree_->Branch(  "BsCtErr2D"                     , &BsCtErr2D_,                             "BsCtErr2D/D");
bsTree_->Branch(  "BsCtErr2DBS"                     , &BsCtErr2DBS_,                             "BsCtErr2DBS/D");
bsTree_->Branch(  "BdCtErr2DBS"                     , &BdCtErr2DBS_,                             "BdCtErr2DBS/D");
bsTree_->Branch(  "BsCtErrMPV"                    , &BsCtErrMPV_,                             "BsCtErrMPV/D");
bsTree_->Branch(  "BsCt3Drefit"                   , &BsCt3Drefit_,                             "BsCt3Drefit/D");
bsTree_->Branch(  "BsCt2Drefit"                   , &BsCt2Drefit_,                             "BsCt2Drefit/D");
bsTree_->Branch(  "BsCtMPVrefit"                  , &BsCtMPVrefit_,                             "BsCtMPVrefit/D");
bsTree_->Branch(  "BsCtErr2D2"                    , &BsCtErr2D2_,                             "BsCtErr2D2/D");
bsTree_->Branch(  "BsCtErrMPV"                    , &BsCtErrMPV_,                             "BsCtErrMPV/D");
bsTree_->Branch(  "BsCtErr3Drefit"                , &BsCtErr3Drefit_,                            "BsCtErr3Drefit/D");
bsTree_->Branch(  "BsCtErr2Drefit"                , &BsCtErr2Drefit_,                            "BsCtErr2Drefit/D");
bsTree_->Branch(  "BsCtErrMPVrefit"               , &BsCtErrMPVrefit_,                             "BsCtErrMPVrefit/D");

bsTree_->Branch(  "BsLxyErr"			  , &BsLxyErr_,                            "BsLxyErr/D");                                                   
bsTree_->Branch(  "JpsiNumberOfCandidates"          , &JpsiNumberOfCandidates_,            "JpsiNumberOfCandidates/I");
bsTree_->Branch(  "PhiNumberOfCandidatesBeforeFit"          , &PhiNumberOfCandidatesBeforeFit_,            "PhiNumberOfCandidatesBeforeFit/I");
bsTree_->Branch(  "BsNumberOfCandidatesBeforeFit"          , &BsNumberOfCandidatesBeforeFit_,            "BsNumberOfCandidatesBeforeFit/I");
bsTree_->Branch(  "BsNumberOfCandidatesAfterFit"          , &BsNumberOfCandidatesAfterFit_,            "BsNumberOfCandidatesAfterFit/I");
bsTree_->Branch(  "BsNumberOfCandidatesAfterBestFit"          , &BsNumberOfCandidatesAfterBestFit_,            "BsNumberOfCandidatesAfterBestFit/I");

bsTree_->Branch(  "BsErrX"			  , &BsErrX_,                           "BsErrX/D");                                                    
bsTree_->Branch(  "BsErrY"			  , &BsErrY_,                           "BsErrY/D");                                                    
bsTree_->Branch(  "BsErrXY"			  , &BsErrXY_,                          "BsErrXY/D");                                                   
                                             
bsTree_->Branch(  "K1trkLay"			  , &K1trkLay_,                         "K1trkLay/I");                                                  
bsTree_->Branch(  "K1muDTh"			  , &K1muDTh_,                          "K1muDTh/I");                                               
bsTree_->Branch(  "K1muCSCh"			  , &K1muCSCh_,                         "K1muCSCh/I");                                                
bsTree_->Branch(  "K1muRPCh"			  , &K1muRPCh_,                         "K1muRPCh/I");                                               
bsTree_->Branch(  "K2trkLay"			  , &K2trkLay_,                         "K2trkLay/I");                                               
bsTree_->Branch(  "K2muDTh"			  , &K2muDTh_,                          "K2muDTh/I");                                               
bsTree_->Branch(  "K2muCSCh"			  , &K2muCSCh_,                         "K2muCSCh/I");                                                
bsTree_->Branch(  "K2muRPCh"			  , &K2muRPCh_,                         "K2muRPCh/I");                                               
bsTree_->Branch(  "Mu1trkLay"			  , &Mu1trkLay_,                        "Mu1trkLay/I");                                               
bsTree_->Branch(  "Mu1muDTh"			  , &Mu1muDTh_,                         "Mu1muDTh/I");                                              
bsTree_->Branch(  "Mu1muCSCh"			  , &Mu1muCSCh_,                        "Mu1muCSCh/I");                                               
bsTree_->Branch(  "Mu1muRPCh"			  , &Mu1muRPCh_,                        "Mu1muRPCh/I");                                              
bsTree_->Branch(  "Mu2trkLay"			  , &Mu2trkLay_,                        "Mu2trkLay/I");                                              
bsTree_->Branch(  "Mu2muDTh"			  , &Mu2muDTh_,                         "Mu2muDTh/I");                                              
bsTree_->Branch(  "Mu2muCSCh"			  , &Mu2muCSCh_,                        "Mu2muCSCh/I");                                               
bsTree_->Branch(  "Mu2muRPCh"			  , &Mu2muRPCh_,                        "Mu2muRPCh/I");                                              
bsTree_->Branch(  "K1mcId"			  , &K1mcId_,                           "K1mcId/I");                                                  
bsTree_->Branch(  "K1momId"			  , &K1momId_,                          "K1momId/I");                                                 
bsTree_->Branch(  "K1gmomId"			  , &K1gmomId_,                         "K1gmomId/I");                                                
bsTree_->Branch(  "K2mcId"			  , &K2mcId_,                           "K2mcId/I");                                               
bsTree_->Branch(  "K2momId"			  , &K2momId_,                          "K2momId/I");                                                 
bsTree_->Branch(  "K2gmomId"			  , &K2gmomId_,                         "K2gmomId/I");                                                
bsTree_->Branch(  "Mu1mcId"			  , &Mu1mcId_,                          "Mu1mcId/I");                                               
bsTree_->Branch(  "Mu1momId"			  , &Mu1momId_,                         "Mu1momId/I");                                                
bsTree_->Branch(  "Mu1gmomId"			  , &Mu1gmomId_,                        "Mu1gmomId/I");                                               
bsTree_->Branch(  "Mu2mcId"			  , &Mu2mcId_,                          "Mu2mcId/I");                                              
bsTree_->Branch(  "Mu2momId"			  , &Mu2momId_,                         "Mu2momId/I");                                                
bsTree_->Branch(  "Mu2gmomId"			  , &Mu2gmomId_,                        "Mu2gmomId/I");                                               
bsTree_->Branch(  "Mu1GlobalMuonPromptTight"      , &Mu1GlobalMuonPromptTight_,         "Mu1GlobalMuonPromptTight/I");
bsTree_->Branch(  "Mu2GlobalMuonPromptTight"      , &Mu2GlobalMuonPromptTight_,         "Mu2GlobalMuonPromptTight/I");
bsTree_->Branch(  "Mu1TrackerMuonArbitrated"      , &Mu1TrackerMuonArbitrated_,         "Mu1TrackerMuonArbitrated/I");
bsTree_->Branch(  "Mu1TMLastStationTight"         , &Mu1TMLastStationTight_,            "Mu1TMLastStationTight/I");
bsTree_->Branch(  "Mu1TMOneStationTight"          , &Mu1TMOneStationTight_,             "Mu1TMOneStationTight/I");
bsTree_->Branch(  "Mu1TMLastStationOptimizedLowPtTight", &Mu1TMLastStationOptimizedLowPtTight_, "Mu1TMLastStationOptimizedLowPtTight/I");
bsTree_->Branch(  "Mu1TMLastStationAngTight"      , &Mu1TMLastStationAngTight_,         "Mu1TMLastStationAngTight/I");
bsTree_->Branch(  "Mu1TMOneStationAngTight"       , &Mu1TMOneStationAngTight_,          "Mu1TMOneStationAngTight/I");
bsTree_->Branch(  "Mu1TMLastStationOptimizedBarrelLowPtTight", &Mu1TMLastStationOptimizedBarrelLowPtTight_, "Mu1TMLastStationOptimizedBarrelLowPtTight/I");
bsTree_->Branch(  "Mu2TrackerMuonArbitrated"      , &Mu2TrackerMuonArbitrated_,         "Mu2TrackerMuonArbitrated/I");
bsTree_->Branch(  "Mu2TMLastStationTight"         , &Mu2TMLastStationTight_,            "Mu2TMLastStationTight/I");
bsTree_->Branch(  "Mu2TMOneStationTight"          , &Mu2TMOneStationTight_,             "Mu2TMOneStationTight/I");
bsTree_->Branch(  "Mu2TMLastStationOptimizedLowPtTight", &Mu2TMLastStationOptimizedLowPtTight_, "Mu2TMLastStationOptimizedLowPtTight/I");
bsTree_->Branch(  "Mu2TMLastStationAngTight"      , &Mu2TMLastStationAngTight_,         "Mu2TMLastStationAngTight/I");
bsTree_->Branch(  "Mu2TMOneStationAngTight"       , &Mu2TMOneStationAngTight_,          "Mu2TMOneStationAngTight/I");
bsTree_->Branch(  "Mu2TMLastStationOptimizedBarrelLowPtTight", &Mu2TMLastStationOptimizedBarrelLowPtTight_, "Mu2TMLastStationOptimizedBarrelLowPtTight/I");
              
              
bsTree_->Branch(  "BsDist3d"			  , &BsDist3d_,                         "BsDist3d/D");                                                  
bsTree_->Branch(  "BsDist3dErr"			  , &BsDist3dErr_,                      "BsDist3dErr/D");                                               
bsTree_->Branch(  "BsTime3d"			  , &BsTime3d_,                         "BsTime3d/D");                                                  
bsTree_->Branch(  "BsTime3dErr"			  , &BsTime3dErr_,                      "BsTime3dErr/D");                                               
bsTree_->Branch(  "BsDist2d"			  , &BsDist2d_,                         "BsDist2d/D");                                                  
bsTree_->Branch(  "BsDist2dErr"			  , &BsDist2dErr_,                      "BsDist2dErr/D");                                               
bsTree_->Branch(  "BsTime2d"			  , &BsTime2d_,                         "BsTime2d/D");                                                  
bsTree_->Branch(  "BsTime2dErr"			  , &BsTime2dErr_,                      "BsTime2dErr/D");                                               
bsTree_->Branch(  "dedxTrk"			  , &dedxTrk_,                          "dedxTrk/D");                                                   
bsTree_->Branch(  "errdedxTrk"			  , &errdedxTrk_,                       "errdedxTrk/D");                                                
bsTree_->Branch(  "numdedxTrk"			  , &numdedxTrk_,                       "numdedxTrk/I");                                              
bsTree_->Branch(  "iPassedCutIdent"		  , &iPassedCutIdent_,                  "iPassedCutIdent/I");                                         
bsTree_->Branch(  "iPassedCutIdentBd"		  , &iPassedCutIdentBd_,                "iPassedCutIdentBd/I");                                        
bsTree_->Branch(  "BdTrack1Charge"		  , &BdTrack1Charge_,                "BdTrack1Charge/I");                                        
bsTree_->Branch(  "K1Fit_par"			  , K1Fit_par_,                         "K1Fit_par[7]/D");                                          
bsTree_->Branch(  "K2Fit_par"			  , K2Fit_par_,                         "K2Fit_par[7]/D");                                      
bsTree_->Branch(  "K1Fit_sigX"			  , &K1Fit_sigX_,                       "K1Fit_sigX/D");                                                
bsTree_->Branch(  "K1Fit_sigY"			  , &K1Fit_sigY_,                       "K1Fit_sigY/D");                                                
bsTree_->Branch(  "K1Fit_sigZ"			  , &K1Fit_sigZ_,                       "K1Fit_sigZ/D");                                                
bsTree_->Branch(  "K2Fit_sigX"			  , &K2Fit_sigX_,                       "K2Fit_sigX/D");                                                
bsTree_->Branch(  "K2Fit_sigY"			  , &K2Fit_sigY_,                       "K2Fit_sigY/D");                                                
bsTree_->Branch(  "K2Fit_sigZ"			  , &K2Fit_sigZ_,                       "K2Fit_sigZ/D");                                                
bsTree_->Branch(  "K1Fit_sigPX"			  , &K1Fit_sigPX_,                      "K1Fit_sigPX/D");                                               
bsTree_->Branch(  "K1Fit_sigPY"			  , &K1Fit_sigPY_,                      "K1Fit_sigPY/D");                                               
bsTree_->Branch(  "K1Fit_sigPZ"			  , &K1Fit_sigPZ_,                      "K1Fit_sigPZ/D");                                               
bsTree_->Branch(  "K2Fit_sigPX"			  , &K2Fit_sigPX_,                      "K2Fit_sigPX/D");                                               
bsTree_->Branch(  "K2Fit_sigPY"			  , &K2Fit_sigPY_,                      "K2Fit_sigPY/D");                                               
bsTree_->Branch(  "K2Fit_sigPZ"			  , &K2Fit_sigPZ_,                      "K2Fit_sigPZ/D");                                               

bsTree_->Branch(  "BBprod", &BBprod_, "BBprod/I");   
          
bsTree_->Branch(  "GenNumberOfBdecays"		  , &GenNumberOfBdecays_,               "GenNumberOfBdecays/I");                                        
bsTree_->Branch(  "BmesonsId"			  , BmesonsId_,                         "BmesonsId[10]/I");                                
bsTree_->Branch(  "BDauIdMC"			  , BDauIdMC_,                          "BDauIdMC[10][15]/I");                                 
bsTree_->Branch(  "BDauDauIdMC"			  , BDauDauIdMC_,                     	"BDauDauIdMC[10][15][10]/I");                        
bsTree_->Branch(  "GenNumberOfDaughters"	  , GenNumberOfDaughters_,              "GenNumberOfDaughters[10]/I");                           
bsTree_->Branch(  "GenNumberOfDaughtersDaughters" , GenNumberOfDaughtersDaughters_,     "GenNumberOfDaughtersDaughters[10][15]/I");            
bsTree_->Branch(  "BDauMMC"			  , BDauMMC_,                           "BDauMMC[10][15]/D");         
bsTree_->Branch(  "BDauPtMC"			  , BDauPtMC_,                          "BDauPtMC[10][15]/D");                                 
bsTree_->Branch(  "BDauPzMC"			  , BDauPzMC_,                          "BDauPzMC[10][15]/D");                             
bsTree_->Branch(  "BDauEtaMC"			  , BDauEtaMC_,                         "BDauEtaMC[10][15]/D");                             
bsTree_->Branch(  "BDauPhiMC"			  , BDauPhiMC_,                         "BDauPhiMC[10][15]/D");                            
bsTree_->Branch(  "BDauDauMMC"			  , BDauDauMMC_,                      	"BDauDauMMC[10][15][10]/D");                         
bsTree_->Branch(  "BDauDauPtMC"			  , BDauDauPtMC_,                     	"BDauDauPtMC[10][15][10]/D");                     
bsTree_->Branch(  "BDauDauPzMC"			  , BDauDauPzMC_,                     	"BDauDauPzMC[10][15][10]/D");                    
bsTree_->Branch(  "BDauDauEtaMC"		  , BDauDauEtaMC_,                    	"BDauDauEtaMC[10][15][10]/D");                    
bsTree_->Branch(  "BDauDauPhiMC"		  , BDauDauPhiMC_,                    	"BDauDauPhiMC[10][15][10]/D");                   
bsTree_->Branch(  "BMMC"			  , BMMC_,                              "BMMC[10]/D");                         
bsTree_->Branch(  "BsCt3DMC"			  , &BsCt3DMC_,                              "BsCt3DMC/D");                         
bsTree_->Branch(  "BsCt2DMC"			  , &BsCt2DMC_,                              "BsCt2DMC/D");                         
bsTree_->Branch(  "BsIniFlavour"			  , &BsIniFlavour_,                              "BsIniFlavour/I");                         
bsTree_->Branch(  "BdIniFlavour"			  , &BdIniFlavour_,                              "BdIniFlavour/I");                         
bsTree_->Branch(  "ChannelID"			  , &ChannelID_,                              "ChannelID_/I");                         
bsTree_->Branch(  "BdChannelID"			  , &BdChannelID_,                              "BdChannelID_/I");                         
bsTree_->Branch(  "BPtMC"			  , BPtMC_,                             "BPtMC[10]/D");                                         
bsTree_->Branch(  "BPxMC"			  , BPxMC_,                             "BPxMC[10]/D");                                        
bsTree_->Branch(  "BPyMC"			  , BPyMC_,                             "BPyMC[10]/D");                                        
bsTree_->Branch(  "BPzMC"			  , BPzMC_,                             "BPzMC[10]/D");                                        
bsTree_->Branch(  "BEtaMC"			  , BEtaMC_,                            "BEtaMC[10]/D");                                        
bsTree_->Branch(  "BPhiMC"			  , BPhiMC_,                            "BPhiMC[10]/D");    

 bsTree_->Branch(  "costhetaMC"			  , costhetaMC_,                        "costhetaMC[10]/D");    
 bsTree_->Branch(  "phiMC"			  , phiMC_,                             "phiMC[10]/D");    
 bsTree_->Branch(  "cospsiMC"			  , cospsiMC_,                          "cospsiMC[10]/D");    
 bsTree_->Branch(  "BscosthetaMC"			  , &BscosthetaMC_,                        "BscosthetaMC/D");    
 bsTree_->Branch(  "BsphiMC"			  , &BsphiMC_,                             "BsphiMC/D");    
 bsTree_->Branch(  "BscospsiMC"			  , &BscospsiMC_,                          "BscospsiMC/D");    
 
bsTree_->Branch(  "BVtxMC_x" ,   BVtxMC_x_  , "BVtxMC_x[10]/D" );
bsTree_->Branch(  "BVtxMC_y" ,	 BVtxMC_y_  , "BVtxMC_y[10]/D" );
bsTree_->Branch(  "BVtxMC_z" ,	 BVtxMC_z_  , "BVtxMC_z[10]/D" );
bsTree_->Branch(  "BSVtxMC_x",	 BSVtxMC_x_ , "BSVtxMC_x[10]/D");
bsTree_->Branch(  "BSVtxMC_y",	 BSVtxMC_y_ , "BSVtxMC_y[10]/D");
bsTree_->Branch(  "BSVtxMC_z",	 BSVtxMC_z_ , "BSVtxMC_z[10]/D");
bsTree_->Branch(  "BLxy_MC"  ,	 BLxy_MC_   , "BLxy_MC[10]/D"  );
bsTree_->Branch(  "BCt_MC"   ,	 BCt_MC_    , "BCt_MC[10]/D"   );   
bsTree_->Branch(  "BCt_MC2D"   ,	 BCt_MC2D_    , "BCt_MC2D[10]/D"   );   
bsTree_->Branch(  "BCt_MC3D"   ,	 BCt_MC3D_    , "BCt_MC3D[10]/D"   );   
                                   
bsTree_->Branch(  "genBsVtx_z"			  , &genBsVtx_z_,                       "genBsVtx_z/D");                                                
bsTree_->Branch(  "genBsVtx_y"			  , &genBsVtx_y_,                       "genBsVtx_y/D");                                                
bsTree_->Branch(  "genBsVtx_x"			  , &genBsVtx_x_,                       "genBsVtx_x/D");                                                
bsTree_->Branch(  "genBsSVtx_z"			  , &genBsSVtx_z_,                      "genBsSVtx_z/D");                                               
bsTree_->Branch(  "genBsSVtx_y" 		  , &genBsSVtx_y_,                      "genBsSVtx_y/D");     
                                       
bsTree_->Branch(  "genBsSVtx_x"			  , &genBsSVtx_x_,                      "genBsSVtx_x/D");                                               
bsTree_->Branch(  "isGenJpsiEvent"		  , &isGenJpsiEvent_,                   "isGenJpsiEvent/I");                                           
bsTree_->Branch(  "BdFitChi2_Hyp1"		  , &BdFitChi2_Hyp1_,                   "BdFitChi2_Hyp1/D");                                            
bsTree_->Branch(  "BdFitNdof_Hyp1"		  , &BdFitNdof_Hyp1_,                   "BdFitNdof_Hyp1/I");                                           
bsTree_->Branch(  "BdFitVtxProb_Hyp1"		  , &BdFitVtxProb_Hyp1_,                "BdFitVtxProb_Hyp1/D");                                         
bsTree_->Branch(  "BdFitM_Hyp1"			  , &BdFitM_Hyp1_,                      "BdFitM_Hyp1/D");                                               
bsTree_->Branch(  "BdFitEta_Hyp1"		  , &BdFitEta_Hyp1_,                    "BdFitEta_Hyp1/D");                                             
bsTree_->Branch(  "BdFitPt_Hyp1"		  , &BdFitPt_Hyp1_,                     "BdFitPt_Hyp1/D");                                              
bsTree_->Branch(  "BdFitPz_Hyp1"		  , &BdFitPz_Hyp1_,                     "BdFitPz_Hyp1/D");                                              
bsTree_->Branch(  "BdFitPhi_Hyp1"		  , &BdFitPhi_Hyp1_,                    "BdFitPhi_Hyp1/D");                                             
bsTree_->Branch(  "BdFitVtx_x_Hyp1"		  , &BdFitVtx_x_Hyp1_,                  "BdFitVtx_x_Hyp1/D");                                           
bsTree_->Branch(  "BdFitVtx_y_Hyp1"		  , &BdFitVtx_y_Hyp1_,                  "BdFitVtx_y_Hyp1/D");                                           
bsTree_->Branch(  "BdFitVtx_z_Hyp1"		  , &BdFitVtx_z_Hyp1_,                  "BdFitVtx_z_Hyp1/D");                                           
bsTree_->Branch(  "BdM_nofit"		  , &BdM_nofit_,                   "BdM_nofit/D");                                            
bsTree_->Branch(  "BdPhi_nofit"		  , &BdPhi_nofit_,                 "BdPhi_nofit/D");                                          
bsTree_->Branch(  "BdEta_nofit"		  , &BdEta_nofit_,                 "BdEta_nofit/D");                                          
bsTree_->Branch(  "BdPt_nofit"		  , &BdPt_nofit_,                  "BdPt_nofit/D");                                           
bsTree_->Branch(  "BdPz_nofit"		  , &BdPz_nofit_,                  "BdPz_nofit/D");                                           
bsTree_->Branch(  "KstarMass_nofit_Hyp1"	  , &KstarMass_nofit_Hyp1_,             "KstarMass_nofit_Hyp1/D"); 
bsTree_->Branch(  "KstarMass_nofit_Hyp2"	  , &KstarMass_nofit_Hyp2_,             "KstarMass_nofit_Hyp2/D");                                      
bsTree_->Branch(  "BdK1_kpi_par_Hyp1"		  , BdK1_kpi_par_Hyp1_,                 "BdK1_kpi_par_Hyp1[7]/D");                                  
bsTree_->Branch(  "BdK2_kpi_par_Hyp1"		  , BdK2_kpi_par_Hyp1_,                 "BdK2_kpi_par_Hyp1[7]/D");                              
bsTree_->Branch(  "BdK1_kpi_sigX_Hyp1"		  , &BdK1_kpi_sigX_Hyp1_,               "BdK1_kpi_sigX_Hyp1/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigY_Hyp1"		  , &BdK1_kpi_sigY_Hyp1_,               "BdK1_kpi_sigY_Hyp1/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigZ_Hyp1"		  , &BdK1_kpi_sigZ_Hyp1_,               "BdK1_kpi_sigZ_Hyp1/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigX_Hyp1"		  , &BdK2_kpi_sigX_Hyp1_,               "BdK2_kpi_sigX_Hyp1/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigY_Hyp1"		  , &BdK2_kpi_sigY_Hyp1_,               "BdK2_kpi_sigY_Hyp1/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigZ_Hyp1"		  , &BdK2_kpi_sigZ_Hyp1_,               "BdK2_kpi_sigZ_Hyp1/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigPX_Hyp1"		  , &BdK1_kpi_sigPX_Hyp1_,              "BdK1_kpi_sigPX_Hyp1/D");                                       
bsTree_->Branch(  "BdK1_kpi_sigPY_Hyp1"		  , &BdK1_kpi_sigPY_Hyp1_,              "BdK1_kpi_sigPY_Hyp1/D");                                       
bsTree_->Branch(  "BdK1_kpi_sigPZ_Hyp1"		  , &BdK1_kpi_sigPZ_Hyp1_,              "BdK1_kpi_sigPZ_Hyp1/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPX_Hyp1"		  , &BdK2_kpi_sigPX_Hyp1_,              "BdK2_kpi_sigPX_Hyp1/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPY_Hyp1"		  , &BdK2_kpi_sigPY_Hyp1_,              "BdK2_kpi_sigPY_Hyp1/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPZ_Hyp1"		  , &BdK2_kpi_sigPZ_Hyp1_,              "BdK2_kpi_sigPZ_Hyp1/D");                                       
bsTree_->Branch(  "BdFitChi2_Hyp2"		  , &BdFitChi2_Hyp2_,                   "BdFitChi2_Hyp2/D");                                            
bsTree_->Branch(  "BdFitNdof_Hyp2"		  , &BdFitNdof_Hyp2_,                   "BdFitNdof_Hyp2/I");                                           
bsTree_->Branch(  "BdFitVtxProb_Hyp2"		  , &BdFitVtxProb_Hyp2_,                "BdFitVtxProb_Hyp2/D");                                         
bsTree_->Branch(  "BdFitM_Hyp2"			  , &BdFitM_Hyp2_,                      "BdFitM_Hyp2/D");                                               
bsTree_->Branch(  "BdFitEta_Hyp2"		  , &BdFitEta_Hyp2_,                    "BdFitEta_Hyp2/D");                                             
bsTree_->Branch(  "BdFitPt_Hyp2"		  , &BdFitPt_Hyp2_,                     "BdFitPt_Hyp2/D");                                              
bsTree_->Branch(  "BdFitPz_Hyp2"		  , &BdFitPz_Hyp2_,                     "BdFitPz_Hyp2/D");                                              
bsTree_->Branch(  "BdFitPhi_Hyp2"		  , &BdFitPhi_Hyp2_,                    "BdFitPhi_Hyp2/D");                                             
bsTree_->Branch(  "BdFitVtx_x_Hyp2"		  , &BdFitVtx_x_Hyp2_,                  "BdFitVtx_x_Hyp2/D");                                           
bsTree_->Branch(  "BdFitVtx_y_Hyp2"		  , &BdFitVtx_y_Hyp2_,                  "BdFitVtx_y_Hyp2/D");                                           
bsTree_->Branch(  "BdFitVtx_z_Hyp2"		  , &BdFitVtx_z_Hyp2_,                  "BdFitVtx_z_Hyp2/D");                                           
bsTree_->Branch(  "BdNumberOfCandidates"          , &BdNumberOfCandidates_,            "BdNumberOfCandidates/I");
bsTree_->Branch(  "BdK1_kpi_par_Hyp2"		  , BdK1_kpi_par_Hyp2_,                 "BdK1_kpi_par_Hyp2[7]/D");                                  
bsTree_->Branch(  "BdK2_kpi_par_Hyp2"		  , BdK2_kpi_par_Hyp2_,                 "BdK2_kpi_par_Hyp2[7]/D");                              
bsTree_->Branch(  "BdK1_kpi_sigX_Hyp2"		  , &BdK1_kpi_sigX_Hyp2_,               "BdK1_kpi_sigX_Hyp2/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigY_Hyp2"		  , &BdK1_kpi_sigY_Hyp2_,               "BdK1_kpi_sigY_Hyp2/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigZ_Hyp2"		  , &BdK1_kpi_sigZ_Hyp2_,               "BdK1_kpi_sigZ_Hyp2/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigX_Hyp2"		  , &BdK2_kpi_sigX_Hyp2_,               "BdK2_kpi_sigX_Hyp2/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigY_Hyp2"		  , &BdK2_kpi_sigY_Hyp2_,               "BdK2_kpi_sigY_Hyp2/D");                                        
bsTree_->Branch(  "BdK2_kpi_sigZ_Hyp2"		  , &BdK2_kpi_sigZ_Hyp2_,               "BdK2_kpi_sigZ_Hyp2/D");                                        
bsTree_->Branch(  "BdK1_kpi_sigPX_Hyp2"		  , &BdK1_kpi_sigPX_Hyp2_,              "BdK1_kpi_sigPX_Hyp2/D");                                       
bsTree_->Branch(  "BdK1_kpi_sigPY_Hyp2"		  , &BdK1_kpi_sigPY_Hyp2_,              "BdK1_kpi_sigPY_Hyp2/D");                                       
bsTree_->Branch(  "BdK1_kpi_sigPZ_Hyp2"		  , &BdK1_kpi_sigPZ_Hyp2_,              "BdK1_kpi_sigPZ_Hyp2/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPX_Hyp2"		  , &BdK2_kpi_sigPX_Hyp2_,              "BdK2_kpi_sigPX_Hyp2/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPY_Hyp2"		  , &BdK2_kpi_sigPY_Hyp2_,              "BdK2_kpi_sigPY_Hyp2/D");                                       
bsTree_->Branch(  "BdK2_kpi_sigPZ_Hyp2"		  , &BdK2_kpi_sigPZ_Hyp2_,              "BdK2_kpi_sigPZ_Hyp2/D");                                       
bsTree_->Branch(  "BdK1Pt_nofit" 		  , &BdK1Pt_nofit_,                     "BdK1Pt_nofit/D");                                              
bsTree_->Branch(  "BdK1Pz_nofit" 		  , &BdK1Pz_nofit_,                     "BdK1Pz_nofit/D");                                              
bsTree_->Branch(  "BdK1Eta_nofit" 		  , &BdK1Eta_nofit_,                    "BdK1Eta_nofit/D");                                             
bsTree_->Branch(  "BdK1Phi_nofit" 		  , &BdK1Phi_nofit_,                    "BdK1Phi_nofit/D");                                             
bsTree_->Branch(  "BdK1Key_nofit" 		  , &BdK1Key_nofit_,                     "BdK1Key_nofit/I");                                              
bsTree_->Branch(  "BdK2Pt_nofit" 		  , &BdK2Pt_nofit_,                     "BdK2Pt_nofit/D");                                              
bsTree_->Branch(  "BdK2Pz_nofit" 		  , &BdK2Pz_nofit_,                     "BdK2Pz_nofit/D");                                              
bsTree_->Branch(  "BdK2Eta_nofit" 		  , &BdK2Eta_nofit_,                    "BdK2Eta_nofit/D");                                             
bsTree_->Branch(  "BdK2Phi_nofit" 		  , &BdK2Phi_nofit_,                    "BdK2Phi_nofit/D");                                             
bsTree_->Branch(  "BdK2Key_nofit" 		  , &BdK2Key_nofit_,                     "BdK2Key_nofit/I"); 

bsTree_->Branch(  "BdPVx_refit" ,  &BdPVx_refit_ ,  "BdPVx_refit/D" );    
bsTree_->Branch(  "BdPVy_refit" ,&BdPVy_refit_ ,    "BdPVy_refit/D");   
bsTree_->Branch(  "BdPVz_refit" ,  &BdPVz_refit_ ,  "BdPVz_refit/D" );              
bsTree_->Branch(  "BdPVerrx_refit" , &BdPVerrx_refit_ , "BdPVerrx_refit/D" ); 
bsTree_->Branch(  "BdPVerry_refit" ,&BdPVerry_refit_ ,"BdPVerry_refit/D" ); 
bsTree_->Branch(  "BdPVerrz_refit" , &BdPVerrz_refit_ , "BdPVerrz_refit/D" ); 
                                             
bsTree_->Branch(  "BdLxy"			  , &BdLxy_,                            "BdLxy/D");                                                     
bsTree_->Branch(  "BdLxyErr"			  , &BdLxyErr_,                            "BdLxyErr/D");                                                     
bsTree_->Branch(  "BdErrX"			  , &BdErrX_,                           "BdErrX/D");                                                    
bsTree_->Branch(  "BdErrY"			  , &BdErrY_,                           "BdErrY/D");                                                    
bsTree_->Branch(  "BdErrXY"			  , &BdErrXY_,                          "BdErrXY/D");                                                   
bsTree_->Branch(  "BdCt"			  , &BdCt_,                             "BdCt/D");                                                      
bsTree_->Branch(  "BdCtErr"			  , &BdCtErr_,                           "BdCtErr/D");                                                      
bsTree_->Branch(  "BdDist3d"			  , &BdDist3d_,                         "BdDist3d/D");                                                  
bsTree_->Branch(  "BdDist3dErr"			  , &BdDist3dErr_,                      "BdDist3dErr/D");                                               
bsTree_->Branch(  "BdTime3d"			  , &BdTime3d_,                         "BdTime3d/D");                                                  
bsTree_->Branch(  "BdTime3dErr"			  , &BdTime3dErr_,                      "BdTime3dErr/D");                                               
bsTree_->Branch(  "BdDist2d"			  , &BdDist2d_,                         "BdDist2d/D");                                                  
bsTree_->Branch(  "BdDist2dErr"			  , &BdDist2dErr_,                      "BdDist2dErr/D");                                               
bsTree_->Branch(  "BdTime2d"			  , &BdTime2d_,                         "BdTime2d/D");                                                  
bsTree_->Branch(  "BdTime2dErr"                   , &BdTime2dErr_,                      "BdTime2dErr/D");                                               
 bsTree_->Branch(  "BdK1mcId"     ,       &BdK1mcId_     ,  "BdK1mcId/I"    );
 bsTree_->Branch(  "BdK1momId"	  ,	  &BdK1momId_    ,  "BdK1momId/I"   );
 bsTree_->Branch(  "BdK1gmomId"	  ,	  &BdK1gmomId_   ,  "BdK1gmomId/I"  );
 bsTree_->Branch(  "BdK2mcId"	  ,	  &BdK2mcId_     ,  "BdK2mcId/I"    );
 bsTree_->Branch(  "BdK2momId"	  ,	  &BdK2momId_    ,  "BdK2momId/I"   );
 bsTree_->Branch(  "BdK2gmomId"	  ,	  &BdK2gmomId_   ,  "BdK2gmomId/I"  );
 bsTree_->Branch(  "BdMu1mcId"	  ,	  &BdMu1mcId_    ,  "BdMu1mcId/I"   );
 bsTree_->Branch(  "BdMu1momId"	  ,	  &BdMu1momId_   ,  "BdMu1momId/I"  );
 bsTree_->Branch(  "BdMu1gmomId"  ,	  &BdMu1gmomId_  ,  "BdMu1gmomId/I" );
 bsTree_->Branch(  "BdMu2mcId"	  ,	  &BdMu2mcId_    ,  "BdMu2mcId/I"   );
 bsTree_->Branch(  "BdMu2momId"	  ,	  &BdMu2momId_   ,  "BdMu2momId/I"  );
 bsTree_->Branch(  "BdMu2gmomId"  ,	  &BdMu2gmomId_  ,  "BdMu2gmomId/I" );

}

BsToJpsiPhiRootTree::~BsToJpsiPhiRootTree()
{}

void BsToJpsiPhiRootTree::writeFile()
{
  bsFile_->Write();
  bsFile_->Close();
  
}

void BsToJpsiPhiRootTree::resetEntries()
{
  runNumber_ =  -9999;
  eventNumber_ =  -9999;
  lumiSection_ =  -9999;
 
  PUinteraction_ = -9999999;

  ihaveajpsi_=  -9999;
  BsCowboy_=  -9999;
  BdCowboy_=  -9999;
  BsPhiVtxProb_=  -9999;
  BsMu1QualityG_=  -9999;
  BsMu2QualityG_=  -9999;
  BsMu1QualityT_=  -9999;
  BsMu2QualityT_=  -9999;
  BdMu1QualityG_=  -9999;
  BdMu2QualityG_=  -9999;
  BdMu1QualityT_=  -9999;
  BdMu2QualityT_=  -9999;


  NVertices_ = -9999999;
  triggerbit_HLTmu3Tk_ = -9999999;
  triggerbit_HLTmu5_ = -9999999;
  triggerbit_HLTdoubleIsoMu3_ = -9999999;
  triggerbit_HLTdoubleMu3_ = -9999999;
  triggerbit_HLTdoubleMu0_ = -9999999;
  triggerbit_HLTL1DoubleMuOpen_ = -9999999;
  triggerbit_HLTmu7_ = -9999999;
  triggerbit_HLTMu0Track0Jpsi_ = -9999999;
  triggerbit_HLTL1DoubleMuOpenTight_ = -9999999;                                     
  triggerbit_HLT_DoubleMu3_Jpsi_v2_ = -9999999;                                     
  triggerbit_HLT_DoubleMu3_Jpsi_v2MC_ = -9999999;                                     
  triggerbit_HLT_DoubleMu3_Quarkonium_v2_ = -9999999;                                     
  triggerbit_HLT_DoubleMu3_Quarkonium_v2MC_ = -9999999;                                     
  triggerbit_HLT_DoubleMu3_Quarkonium_v1_ = -9999999;                                     
  triggerbit_Jpsi_Displaced_v1_ = -9999999;                                     
  triggerbit_7Jpsi_Displaced_v1_ = -9999999;                                     
  triggerbit_Jpsi_Displaced_v1MC_ = -9999999;                                     
  triggerbit_7Jpsi_Displaced_v2_ = -9999999;                                     
  triggerbit_7Jpsi_Displaced_v3_ = -9999999;                                     
  triggerbit_3p5Jpsi_Displaced_v2_ = -9999999;                                     
  triggerbit_4Jpsi_Displaced_v1_ = -9999999;                                     
  triggerbit_4Jpsi_Displaced_v4_ = -9999999;                                     

  triggerbit_4Jpsi_Displaced_v5_ = -9999999;
  triggerbit_4Jpsi_Displaced_v9_ = -9999999;
  triggerbit_4Jpsi_Displaced_v10_ = -9999999;
  triggerbit_4Jpsi_Displaced_v11_ = -9999999; 
  triggerbit_4Jpsi_Displaced_v12_ = -9999999;                                         
                                      

  triggerbit_5Jpsi_Displaced_v1_ = -9999999;                                     
  triggerbit_5Jpsi_Displaced_v2_ = -9999999;                                     
  triggerbit_5Jpsi_Displaced_v4_ = -9999999;                                     
  triggerbit_5Jpsi_Displaced_v5_ = -9999999;                                     
  triggerbit_Dimuon0_Jpsi_v1_= -9999999;   
  triggerbit_Dimuon0_Jpsi_v3_= -9999999;   
  triggerbit_Dimuon0_Jpsi_v5_= -9999999;   
  triggerbit_Dimuon0_Jpsi_v6_= -9999999;   
  triggerbit_Dimuon0_Jpsi_v9_= -9999999;   
  triggerbit_Dimuon0_Jpsi_v10_= -9999999;   
  triggerbit_Dimuon10_Barrel_= -9999999;
  triggerbit_Dimuon13_Barrel_= -9999999;



  BSx_ = -9999999;
  BSy_ = -9999999;
  BSz_ = -9999999;
  BSdx_ = -9999999;
  BSdy_ = -9999999;
  BSdz_ = -9999999;
  BSsigmaZ_ = -9999999;
  BSdsigmaZ_ = -9999999;

  PVx_ = -9999999;
  PVy_ = -9999999;
  PVz_ = -9999999;
  PVerrx_ = -9999999;
  PVerry_ = -9999999;
  PVerrz_ = -9999999;

  MuMuDCA_ = -9999999;
  MuMuDistance_ = -9999999;
  MuMuDistanceSigma_ = -9999999;
  MuDr1_ = -9999999;
  MuDr2_ = -9999999;
  BdMuMuDCA_ = -9999999;
  BdMuMuDistance_ = -9999999;
  BdMuMuDistanceSigma_ = -9999999;
  BdMuDr1_ = -9999999;
  BdMuDr2_ = -9999999;

  PVx_refit_= -9999999; 
  PVy_refit_= -9999999;
  PVz_refit_= -9999999;   
  PVerrx_refit_= -9999999;
  PVerry_refit_= -9999999;
  PVerrz_refit_= -9999999;

  BdMu1Pt_beffit_= -9999999;
  BdMu1Pz_beffit_= -9999999;
  BdMu1Eta_beffit_= -9999999;
  BdMu1Phi_beffit_= -9999999;
  BdMu2Pt_beffit_= -9999999;
  BdMu2Pz_beffit_= -9999999;
  BdMu2Eta_beffit_= -9999999;
  BdMu2Phi_beffit_= -9999999;
  BdJpsiM_nofit_= -9999999;
  BdJpsiEta_nofit_= -9999999;
  BdJpsiPhi_nofit_= -9999999;
  BdJpsiPt_nofit_= -9999999;
  BdJpsiPz_nofit_= -9999999;
  
  matchL11_ =  -9999999;
  matchL12_ =  -9999999;
  match2mu01_ =  -9999999;
  match2mu02_ =  -9999999;
  match2mu31_ =  -9999999;
  match2mu32_ =  -9999999;
  match1mu01_ =  -9999999;
  match1mu02_ =  -9999999;
  matchmu0tk01_ =  -9999999;
  matchmu0tk02_ =  -9999999;
  matchFilterJpsi1_ =  -9999999;
  matchFilterJpsi2_ =  -9999999;
  matchDoubleMu31J_ =  -9999999;
  matchDoubleMu32J_ =  -9999999;
  matchDoubleMu31Q_ =  -9999999;
  matchDoubleMu32Q_ =  -9999999;
  matchDoubleMu71_ =  -9999999;
  matchDoubleMu72_ =  -9999999;
  matchDoubleMu41_ =  -9999999;
  matchDoubleMu42_ =  -9999999;
  matchDoubleMu51_ =  -9999999;
  matchDoubleMu52_ =  -9999999;
  matchDoubleMu01_ =  -9999999;
  matchDoubleMu02_ =  -9999999;
  matchDoubleMu101_ =  -9999999;
  matchDoubleMu102_ =  -9999999;
  matchDoubleMu131_ =  -9999999;
  matchDoubleMu132_ =  -9999999;

matchDoubleMu01DiMuon0_ = -999;
matchDoubleMu02DiMuon0_ = -999;
BpmatchDoubleMu01DiMuon0_ = -999;
BpmatchDoubleMu02DiMuon0_ = -999;
BdmatchDoubleMu01DiMuon0_ = -999;
BdmatchDoubleMu02DiMuon0_ = -999;

  BdmatchL11_ =  -9999999;
  BdmatchL12_ =  -9999999;
  Bdmatch2mu01_ =  -9999999;
  Bdmatch2mu02_ =  -9999999;
  Bdmatch2mu31_ =  -9999999;
  Bdmatch2mu32_ =  -9999999;
  Bdmatch1mu01_ =  -9999999;
  Bdmatch1mu02_ =  -9999999;
  BdmatchDoubleMu31Q_ =  -9999999;
  BdmatchDoubleMu32Q_ =  -9999999;
  BdmatchDoubleMu71_ =  -9999999;
  BdmatchDoubleMu72_ =  -9999999;
  BdmatchDoubleMu41_ =  -9999999;
  BdmatchDoubleMu42_ =  -9999999;
  BdmatchDoubleMu51_ =  -9999999;
  BdmatchDoubleMu52_ =  -9999999;
  BdmatchDoubleMu01_ =  -9999999;
  BdmatchDoubleMu02_ =  -9999999;
  BdmatchDoubleMu101_ =  -9999999;
  BdmatchDoubleMu102_ =  -9999999;
  BdmatchDoubleMu131_ =  -9999999;
  BdmatchDoubleMu132_ =  -9999999;

  MuonType_ = - 9999999;
  


  JpsiVtxProb_ = -9999999;
  BdJpsiVtxProb_ = -9999999;
  JpsiM_alone_ = -9999999;
  JpsiPhi_alone_ = -9999999;
  JpsiEta_alone_ = -9999999;
  JpsiPt_alone_ = -9999999;
  JpsiMu1Pt_alone_ = -9999999;
  JpsiMu2Pt_alone_ = -9999999;
  JpsiMu1Phi_alone_ = -9999999;
  JpsiMu2Phi_alone_ = -9999999;
  JpsiMu1Eta_alone_ = -9999999;
  JpsiMu2Eta_alone_ = -9999999;
  JpsiMuon1Cat_alone_ = -9999999;
  JpsiMuon2Cat_alone_ = -9999999;
  BdJpsiMuon1Cat_alone_ = -9999999;
  BdJpsiMuon2Cat_alone_ = -9999999;

  JpsiMu1d0_alone_ = -9999999;
  JpsiMu2d0_alone_ = -9999999;
  JpsiMu1dz_alone_ = -9999999;
  JpsiMu2dz_alone_ = -9999999;
  JpsiMu1chi2_alone_ = -9999999;
  JpsiMu2chi2_alone_ = -9999999;
  JpsiMu1ndof_alone_ = -9999999;
  JpsiMu2ndof_alone_ = -9999999;
  JpsiMu1nHits_alone_ = -9999999;
  JpsiMu2nHits_alone_ = -9999999;
  JpsiMu1nPixHits_alone_ = -9999999;
  JpsiMu2nPixHits_alone_ = -9999999;

  K1Pt_beffit_ = -9999999;
  K1Pz_beffit_ = -9999999;
  K1Eta_beffit_ = -9999999;
  K1Phi_beffit_ = -9999999;
  K2Pt_beffit_ = -9999999;
  K2Pz_beffit_ = -9999999;
  K2Eta_beffit_ = -9999999;
  K2Phi_beffit_ = -9999999;

  Mu1Pt_beffit_ = -9999999;
  Mu1Pz_beffit_ = -9999999;
  Mu1Eta_beffit_ = -9999999;
  Mu1Phi_beffit_ = -9999999;
  Mu2Pt_beffit_ = -9999999;
  Mu2Pz_beffit_ = -9999999;
  Mu2Eta_beffit_ = -9999999;
  Mu2Phi_beffit_ = -9999999;

  BsFitChi2_ = -9999999;
  BsFitNdof_ = -9999999;
  BsFitVtxProb_ = -9999999;
  JpsiNumberOfCandidates_ = 0;
  PhiNumberOfCandidatesBeforeFit_ =  0;
  BsNumberOfCandidatesBeforeFit_ =  0;
  BsNumberOfCandidatesAfterFit_ =  0;
  BsNumberOfCandidatesAfterBestFit_ =  0;
  BsFitM_ = -999;
  K1Pt_fit_ = -9999999;
  K2Pt_fit_ = -9999999;
  PhiM_fit_ = -9999999;
  BsFitEta_ = -9999999;
  BsFitPt_ = -9999999;
  BsFitPz_ = -9999999;
  BsFitPhi_ = -9999999;
  BsFitVtx_x_ = -9999999;
  BsFitVtx_y_ = -9999999;
  BsFitVtx_z_ = -9999999;
  BsM_nofit_ = -9999999;
  BsPhi_nofit_ = -9999999;
  BsEta_nofit_ = -9999999;
  BsPt_nofit_ = -9999999;
  BsPz_nofit_ = -9999999;
  JpsiM_nofit_ = -9999999;
  JpsiPhi_nofit_ = -9999999;
  JpsiEta_nofit_ = -9999999;
  JpsiPt_nofit_ = -9999999;
  JpsiPz_nofit_ = -9999999;
  PhiM_nofit_ = -9999999;
  PhiPhi_nofit_ = -9999999;
  PhiEta_nofit_ = -9999999;
  PhiPt_nofit_ = -9999999;
  PhiPz_nofit_ = -9999999;
  K1Pt_nofit_ = -9999999;
  K1Pz_nofit_ = -9999999;
  K1Eta_nofit_ = -9999999;
  K1Phi_nofit_ = -9999999;
 K1Key_nofit_ = -9999999;
  K2Eta_nofit_ = -9999999;
  K2Pt_nofit_ = -9999999;
  K2Pz_nofit_ = -9999999;
  K2Phi_nofit_ = -9999999;
 K2Key_nofit_ = -9999999;
  K1Chi2_ = -9999999;
  K1nHits_ = -9999999;
  K2Chi2_ = -9999999;
  K2nHits_ = -9999999;
  K1pixH_ = -9999999;
  K1trkH_ = -9999999;
  K2pixH_ = -9999999;
  K2trkH_ = -9999999;
  Mu1Chi2_ = -9999999;
  Mu1nHits_ = -9999999;
  Mu2Chi2_ = -9999999;
  Mu2nHits_ = -9999999;
  Mu1pixH_ = -9999999;
  Mu1trkH_ = -9999999;
  Mu2pixH_ = -9999999;
  Mu2trkH_ = -9999999;
  costheta_ = -9999999;
  phi_ = -9999999;
  cospsi_ = -9999999;
  Bdcostheta_ = -9999999;
  Bdphi_ = -9999999;
  Bdcospsi_ = -9999999;
  BdcosthetaMC_ = -9999999;
  BdphiMC_ = -9999999;
  BdcospsiMC_ = -9999999;
  AngleBsDecayLength_ = -9999999;
  CosDeltaAlpha_ = -9999999;
  BdCosDeltaAlpha_ = -9999999;

  isPV_ = -9999999;
  isBS_ = -9999999;

  isMatched_ = -9999999;
  isMatchedBd_ = -9999999;
  BsLxy_ = -9999999;
  BsLxyErr_ = -9999999;
  BsErrX_ = -9999999;
  BsErrY_ = -9999999;
  BsErrXY_ = -9999999;
  BsCt_ = -9999999;
  BsCtErr_ = -9999999;
  BsCt3D_ = -9999999;
  BsCt2D_ = -9999999;
  BsCt2DBS_ = -9999999;
  BdCt2DBS_ = -9999999;
  BdCt2DMC_ = -9999999;
  BdCt3DMC_ = -9999999;
  BsCtMPV_ = -9999999;
  BsCt3Drefit_ = -9999999;
  BsCt2Drefit_ = -9999999;
  BsCtMPVrefit_ = -9999999;
  BsCtErr3D_ = -9999999;
  BsCtErr2D_ = -9999999;
  BsCtErr2DBS_ = -9999999;
  BdCtErr2DBS_ = -9999999;
  BsCtErr2D2_ = -9999999;
  BsCtErrMPV_ = -9999999;
  BsCtErr3Drefit_ = -9999999;
  BsCtErr2Drefit_ = -9999999;
  BsCtErrMPVrefit_ = -9999999;
  K1trkLay_ = -9999999;
  K1muDTh_ = -9999999;
  K1muCSCh_ = -9999999;
  K1muRPCh_ = -9999999;
  K2trkLay_ = -9999999;
  K2muDTh_ = -9999999;
  K2muCSCh_ = -9999999;
  K2muRPCh_ = -9999999;
  Mu1trkLay_ = -9999999;
  Mu1muDTh_ = -9999999;
  Mu1muCSCh_ = -9999999;
  Mu1muRPCh_ = -9999999;
  Mu2trkLay_ = -9999999;
  Mu2muDTh_ = -9999999;
  Mu2muCSCh_ = -9999999;
  Mu2muRPCh_ = -9999999;
  K1mcId_ = -9999999;
  K1momId_ = -9999999;
  K1gmomId_ = -9999999;
  K2mcId_ = -9999999;
  K2momId_ = -9999999;
  K2gmomId_ = -9999999;
  Mu1mcId_ = -9999999;
  Mu1momId_ = -9999999;
  Mu1gmomId_ = -9999999;
  Mu2mcId_ = -9999999;
  Mu2momId_ = -9999999;
  Mu2gmomId_ = -9999999;

  Mu1d0_ = -9999999;
  Mu2d0_ = -9999999;
  Mu1dz_ = -9999999;
  Mu2dz_ = -9999999;
  Mu1GlobalMuonPromptTight_=-9999999;
  Mu2GlobalMuonPromptTight_=-9999999;
  Mu1TrackerMuonArbitrated_=-9999999;
  Mu1TMLastStationTight_=-9999999;
  Mu1TMOneStationTight_=-9999999;
  Mu1TMLastStationOptimizedLowPtTight_=-9999999;
  Mu1TMLastStationAngTight_=-9999999;
  Mu1TMOneStationAngTight_=-9999999;
  Mu1TMLastStationOptimizedBarrelLowPtTight_=-9999999;
  Mu2TrackerMuonArbitrated_=-9999999;
  Mu2TMLastStationTight_=-9999999;
  Mu2TMOneStationTight_=-9999999;
  Mu2TMLastStationOptimizedLowPtTight_=-9999999;
  Mu2TMLastStationAngTight_=-9999999;
  Mu2TMOneStationAngTight_=-9999999;
  Mu2TMLastStationOptimizedBarrelLowPtTight_=-9999999;

  

  BsDist3d_ = -9999999;
  BsDist3dErr_ = -9999999;
  BsTime3d_ = -9999999;
  BsTime3dErr_ = -9999999;
  BsDist2d_ = -9999999;
  BsDist2dErr_ = -9999999;
  BsTime2d_ = -9999999;
  BsTime2dErr_ = -9999999;
  dedxTrk_ = -9999999;
  errdedxTrk_ = -9999999;
  numdedxTrk_ = -9999999;
  iPassedCutIdent_ = -9999999;
  iPassedCutIdentBd_ = -9999999;
  BdTrack1Charge_ = -9999999;
 
  K1Fit_sigX_ = -9999999;
  K1Fit_sigY_ = -9999999;
  K1Fit_sigZ_ = -9999999;
  K2Fit_sigX_ = -9999999;
  K2Fit_sigY_ = -9999999;
  K2Fit_sigZ_ = -9999999;
  K1Fit_sigPX_ = -9999999;
  K1Fit_sigPY_ = -9999999;
  K1Fit_sigPZ_ = -9999999;
  K2Fit_sigPX_ = -9999999;
  K2Fit_sigPY_ = -9999999;
  K2Fit_sigPZ_ = -9999999;
 
  GenNumberOfBdecays_ = -9999999;

  genBsVtx_z_ = -9999999;
  genBsVtx_y_ = -9999999;
  genBsVtx_x_ = -9999999;
  genBsSVtx_z_ = -9999999;
  genBsSVtx_y_ = -9999999; 
  genBsSVtx_x_ = -9999999;
  isGenJpsiEvent_ = -9999999;
  BdFitChi2_Hyp1_ = -9999999;
  BdFitNdof_Hyp1_ = -9999999;
  BdFitVtxProb_Hyp1_ = -9999999;
  BdFitM_Hyp1_ = -9999999;
  BdFitEta_Hyp1_ = -9999999;
  BdFitPt_Hyp1_ = -9999999;
  BdFitPz_Hyp1_ = -9999999;
  BdFitPhi_Hyp1_ = -9999999;
  BdFitVtx_x_Hyp1_ = -9999999;
  BdFitVtx_y_Hyp1_ = -9999999;
  BdFitVtx_z_Hyp1_ = -9999999;
  BdM_nofit_ = -9999999;
  BdPhi_nofit_ = -9999999;
  BdEta_nofit_ = -9999999;
  BdPt_nofit_ = -9999999;
  BdPz_nofit_ = -9999999;
  KstarMass_nofit_Hyp1_ = -9999999;
  KstarMass_nofit_Hyp2_ = -9999999;
  BdK1_kpi_sigX_Hyp1_ = -9999999;
  BdK1_kpi_sigY_Hyp1_ = -9999999;
  BdK1_kpi_sigZ_Hyp1_ = -9999999;
  BdK2_kpi_sigX_Hyp1_ = -9999999;
  BdK2_kpi_sigY_Hyp1_ = -9999999;
  BdK2_kpi_sigZ_Hyp1_ = -9999999;
  BdK1_kpi_sigPX_Hyp1_ = -9999999;
  BdK1_kpi_sigPY_Hyp1_ = -9999999;
  BdK1_kpi_sigPZ_Hyp1_ = -9999999;
  BdK2_kpi_sigPX_Hyp1_ = -9999999;
  BdK2_kpi_sigPY_Hyp1_ = -9999999;
  BdK2_kpi_sigPZ_Hyp1_ = -9999999;
  BdFitChi2_Hyp2_ = -9999999;
  BdFitNdof_Hyp2_ = -9999999;
  BdFitVtxProb_Hyp2_ = -9999999;
  BdFitM_Hyp2_ = -9999999;
  BdFitEta_Hyp2_ = -9999999;
  BdFitPt_Hyp2_ = -9999999;
  BdFitPz_Hyp2_ = -9999999;
  BdFitPhi_Hyp2_ = -9999999;
  BdFitVtx_x_Hyp2_ = -9999999;
  BdFitVtx_y_Hyp2_ = -9999999;
  BdFitVtx_z_Hyp2_ = -9999999;
  BdNumberOfCandidates_ =  0; 

BdPVx_refit_    = -9999999;
BdPVy_refit_    = -9999999;
BdPVz_refit_    = -9999999;
               
BdPVerrx_refit_ = -9999999;
BdPVerry_refit_ = -9999999;
BdPVerrz_refit_ = -9999999;

  BdK1_kpi_sigX_Hyp2_ = -9999999;
  BdK1_kpi_sigY_Hyp2_ = -9999999;
  BdK1_kpi_sigZ_Hyp2_ = -9999999;
  BdK2_kpi_sigX_Hyp2_ = -9999999;
  BdK2_kpi_sigY_Hyp2_ = -9999999;
  BdK2_kpi_sigZ_Hyp2_ = -9999999;
  BdK1_kpi_sigPX_Hyp2_ = -9999999;
  BdK1_kpi_sigPY_Hyp2_ = -9999999;
  BdK1_kpi_sigPZ_Hyp2_ = -9999999;
  BdK2_kpi_sigPX_Hyp2_ = -9999999;
  BdK2_kpi_sigPY_Hyp2_ = -9999999;
  BdK2_kpi_sigPZ_Hyp2_ = -9999999;
  BdK1Pt_nofit_ = -9999999; 
  BdK1Pz_nofit_ = -9999999; 
  BdK1Eta_nofit_ = -9999999; 
  BdK1Phi_nofit_ = -9999999; 
 BdK1Key_nofit_ = -9999999; 
  BdK2Pt_nofit_ = -9999999; 
  BdK2Pz_nofit_ = -9999999; 
  BdK2Eta_nofit_ = -9999999; 
  BdK2Phi_nofit_ = -9999999; 
  BdK2Key_nofit_ = -9999999; 
  BdLxy_ = -9999999;
  BdLxyErr_ = -9999999;
  BdErrX_ = -9999999;
  BdErrY_ = -9999999;
  BdErrXY_ = -9999999;
  BdCt_ = -999;
  BdCtErr_ = -999;
  BdDist3d_ = -9999999;
  BdDist3dErr_ = -9999999;
  BdTime3d_ = -9999999;
  BdTime3dErr_ = -9999999;
  BdDist2d_ = -9999999;
  BdDist2dErr_ = -9999999;
  BdTime2d_ = -9999999;
  BdTime2dErr_ = -9999999;                                    

BdK1mcId_    = -9999999;
BdK1momId_   = -9999999;
BdK1gmomId_  = -9999999;
BdK2mcId_    = -9999999;
BdK2momId_   = -9999999;
BdK2gmomId_  = -9999999;
BdMu1mcId_   = -9999999;
BdMu1momId_  = -9999999;
BdMu1gmomId_ = -9999999;
BdMu2mcId_   = -9999999;
BdMu2momId_  = -9999999;
BdMu2gmomId_ = -9999999;

//jet variables 
//JetTrkPt_= new std::vector< std::vector<float>* >(); //initialisation of std vector pointer
//JetTrkCharge_= new std::vector< std::vector<int>* >(); //initialisation of std vector pointer

BJetEta_ = -999;
BJetPhi_ = -999;
BJetParton_ = -999;

BsJetdR_ = -999;
BsJetPx_ = -999;
BsJetPy_ = -999;
BsJetPz_ = -999;

BdJetdR_ = -999;
BdJetPx_ = -999;
BdJetPy_ = -999;
BdJetPz_ = -999;

BsPtJetdR_ = -999;
BsBJetdR_ = -999;
BsBOldJetdR_ = -999;

PtJetPt_ = -999;
BJetPt_ = -999;
BOldJetPt_ = -999;

BplusPtJetdR_ = -999;
BplusBJetdR_ = -999;
BplusJetPx_ = -999;
BplusJetPy_ = -999;
BplusJetPz_ = -999;
JetBTagProb_ = -999;
JetBOldTagProb_ = -999;



BdJetCharge_ = -999;
BsJetCharge_ = -999;
BplusJetCharge_ = -999;

triggerbit_Dimuon0_Jpsi_Muon_v15_ = -999;
triggerbit_Dimuon0_Jpsi_Muon_v16_ = -999;
triggerbit_Dimuon0_Jpsi_Muon_v17_ = -999;
triggerbit_Dimuon0_Jpsi_Muon_v18_ = -999;
//jet variables

 PVindex_ = -999;

//Bplus variables
  BplusIsBS_ = -999;
  BplusIsPV_ = -999; 
  BplusPVindex_ = -999;
  BplusPhi_ = -99999999;    
  BplusEta_ = -99999999; 
triggerbit_7Jpsi_Displaced_v1_ =  -9999999;  
LxyPerSigma2D_ = -9999999;
  JpsiMass_bplus_ = -9999999;
  JpsiPt_bplus_ = -9999999;
  BplusCosTheta_ = -9999999;
  BplusCtPerSigma2D_ = -9999999;
  BplusCtErr2D_ = -9999999; 
  BplusCtPerSigma3D_ = -9999999;
  BplusCt2D_= -9999999;
  BplusCt2D_2_= -9999999; 
  BplusCt3D_2_ = -9999999; 
  BplusVtxProb_ = -9999999; 
  BplusM_fit_ = -9999999;
  BplusChi2_ = -9999999;
  KplusPt_ = -9999999;
  KplusPtot_ = -9999999;
  BplusPt_ = -9999999; 
  BplusMu1Pt_ = -9999999; 
  BplusMu2Pt_ = -9999999; 
  BplusMu1Ptot_ = -9999999;
  BplusMu2Ptot_ = -9999999; 
  BplusLxy_ = -9999999;
  BplusLxyErr_ = -9999999;
  BplusLxyz_ = -9999999;
  BplusLxyzErr_= -9999999;	
  BplusCt3D_ =	-9999999;
  BplusCtErr3D_ = -9999999;
  BplusPtot_ = -9999999;
  BplusDecayChannel_ =  -9999999;
  BplusKmcId_= -9999999;	
  BplusKmomId_ = -9999999;


  BplusMu1mcId_= -9999999;
  BplusMu1momId_= -9999999;
  BplusMu1gmomId_= -9999999;	
	
  BplusMu2mcId_= -9999999;
  BplusMu2momId_ = -9999999;
  BplusMu2gmomId_ = -9999999;

  isMatchedBplus_ = -9999999;

  BpmatchDoubleMu01_ = -999;
  BpmatchDoubleMu02_ = -999;

  IP3DKandJpsiVtx_ = -999;
  IP3DKandJpsiVtxErr_ = -999;  	 
  BplusKIP3D_ = -999;
  BplusKIP3DErr_ = -999; 	
//Bplus variables

//new Bd variables
BdCt3D_ = -999; 
BdCtErr3D_ =-999; 
BdCt2D_ =-999;
BdCtErr2D_ =-999; 
//new Bd variables


BsCt3DMC_ = -9999999;
BsCt2DMC_ = -9999999;
BsIniFlavour_ = -9999999;
BdIniFlavour_ = -9999999;
ChannelID_ = -9999999;
BdChannelID_ = -9999999;
BscosthetaMC_ = -9999999;
BsphiMC_ = -9999999;
BscospsiMC_ = -9999999;
BsCosThetaReco_ = -9999999;

BsCt3DPVCosTheta_ = -999;
BsCt2DPVCosTheta_ =-999;

BsCt3DPVHighestPt_ = -999;
BsCt2DPVHighestPt_ = -999;

BBprod_ = -999;

  for(int i=0; i<7; i++){
    K1Fit_par_[i] = -9999999;
    K2Fit_par_[i] = -9999999;

    BdK1_kpi_par_Hyp1_[i] = -9999999;
    BdK2_kpi_par_Hyp1_[i] = -9999999;
    BdK1_kpi_par_Hyp2_[i] = -9999999;
    BdK2_kpi_par_Hyp2_[i] = -9999999;
  }
    

	for(int i = 0; i<30; i++){

		SVZpos_[i] = -999.0;
		PVZpos_[i] = -999.0;		
		PVAbsPt_[i] = -999.0;
		NTracksInPV_[i] = -999;

		if(i<25){
			TagMuRecoPtRel_[i] = -999;
			TagMuRecoIP3D_[i] = -999.0;
			TagMuRecoIP3DErr_[i] = -999.0;
			TagMuRecoPt_[i] = -999.0;   
			TagMuRecoP_[i]= -999.0; 
			TagMuRecoEta_[i] = -999.0;
			TagMuRecoPhi_[i] = -999.0;
			TagMuRecoChg_[i] = -999.0;
			TagMuSimu_[i] = -999.0;
			MuRecoMCmother_[i] = -999;
			MuRecoMCMotherMother_[i] = -999;
			GlobalTagMuon_[i] = 0;
			TrackerTagMuon_[i] = 0;
			TagMuJetdR_[i]=-999.0;
			PFTagMuon_[i] = 0;
			TagMuPtJetPtRatio_[i] = -999;
		}
		
	}

	BsLxy3DMC_ = -999;
	BsLxy2DMC_ = -999;
	TagMuListSize_ = -999;
	BsPMC_ = -999;
	BsPtMC_ = -999;
	FirstBsMCZpos_ = -999;
	BsPVDist3d_= -999;
	BsPVDist2d_ = -99;

  for(int i=0; i<10; i++){ 
    BmesonsId_[i]  =  -9999999;

    costhetaMC_[i] = -9999999;
    phiMC_[i] = -9999999;
    cospsiMC_[i] = -9999999;
   
    BMMC_[i] =  -9999999;
    BPtMC_[i] =  -9999999;
    BPxMC_[i] =  -9999999;
    BPyMC_[i] =  -9999999;
    BPzMC_[i] =  -9999999;
    BEtaMC_[i] =  -9999999;
    BPhiMC_[i] =  -9999999;

    BVtxMC_x_[i] =  -9999999;
    BVtxMC_y_[i]  =  -9999999;
    BVtxMC_z_[i]  =  -9999999;
    BSVtxMC_x_[i] =  -9999999;
    BSVtxMC_y_[i] =  -9999999;
    BSVtxMC_z_[i] =  -9999999;
    BLxy_MC_[i]   =  -9999999;
    BCt_MC_[i]    =  -9999999;
    BCt_MC2D_[i]    =  -9999999;
    BCt_MC3D_[i]    =  -9999999;
    
    GenNumberOfDaughters_[i] =  -9999999;
	
    for(int j=0;j<15;j++){
      BDauIdMC_[i][j]= -9999999;       
	
      BDauMMC_[i][j]= -9999999;  
      BDauPtMC_[i][j]= -9999999; 
      BDauPzMC_[i][j]= -9999999; 
      BDauEtaMC_[i][j]= -9999999;
      BDauPhiMC_[i][j]= -9999999;
      
      GenNumberOfDaughtersDaughters_[i][j] =  -9999999;
      


      for(int k=0; k<10; k++){
	BDauDauIdMC_[i][j][k]= -9999999;
	
	BDauDauMMC_[i][j][k]= -9999999;  
	BDauDauPtMC_[i][j][k]= -9999999; 
	BDauDauPzMC_[i][j][k]= -9999999; 
	BDauDauEtaMC_[i][j][k]= -9999999;
	BDauDauPhiMC_[i][j][k]= -9999999;
      }
    }
  }

 
 
    
  
} 

void BsToJpsiPhiRootTree::getDeDx(const double f1, const double f2, const int f3)
{
  dedxTrk_ = f1;
  errdedxTrk_ = f2;
  numdedxTrk_ = f3;
}





void BsToJpsiPhiRootTree::getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff,
				 const double gg, const double hh, const double ii)
{
  BSx_ = aa;
  BSy_ = bb;
  BSz_ = cc;
  PVx_ = dd;
  PVy_ = ee;
  PVz_ = ff;
  PVerrx_ = gg;
  PVerry_ = hh;
  PVerrz_ = ii;
}





void BsToJpsiPhiRootTree::getAngles(const double aa, const double bb, const double cc, const double dd)
{
  costheta_ = aa;
  phi_ = bb;
  cospsi_ = cc;
  AngleBsDecayLength_ = dd;
}







void BsToJpsiPhiRootTree::fill()
{
  bsTree_->Fill();
}






void BsToJpsiPhiRootTree::readTree(const std::string filename)
{

  // open root file
  bsFile_ = new TFile (filename.c_str(), "READ" );
 
  // create tree structure
  bsTree_ =  (TTree*) bsFile_->Get("BsTree");
  
  setBranchAddresses();
}

void BsToJpsiPhiRootTree::readTree(std::vector<std::string> filenames){
  TChain * myChain = new TChain("BsTree");
//  for(int i=0;i<filenames.size();++i) {
//    std::string filenome = filenames[i] ;
//    myChain->Add(filenome,-1);
//  }

  for(std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++){
    myChain->Add( (*it).c_str());
  }

  bsTree_ = myChain;
  setBranchAddresses();
}

void BsToJpsiPhiRootTree::setBranchAddresses(){

bsTree_->SetBranchAddress("BJetParton", &BJetParton_);
bsTree_->SetBranchAddress("PtJetTrkPt", PtJetTrkPt_);
bsTree_->SetBranchAddress("PtJetTrkCharge", PtJetTrkCharge_);
bsTree_->SetBranchAddress("PVTrkPt", PVTrkPt_);

bsTree_->SetBranchAddress(  "TagMuRecoIP3D"             , &TagMuRecoIP3D_ );
bsTree_->SetBranchAddress(  "TagMuRecoIP3DErr"          , &TagMuRecoIP3DErr_ );

bsTree_->SetBranchAddress(  "TagMuRecoPtRel"             , &TagMuRecoPtRel_ ); 
bsTree_->SetBranchAddress(  "TagMuRecoPt"             , &TagMuRecoPt_ );   
bsTree_->SetBranchAddress(  "TagMuRecoP"              , &TagMuRecoP_ ); 
bsTree_->SetBranchAddress(  "TagMuRecoEta"            , &TagMuRecoEta_ );
bsTree_->SetBranchAddress(  "TagMuRecoPhi"            , &TagMuRecoPhi_ );
bsTree_->SetBranchAddress(  "TagMuRecoChg"            , &TagMuRecoChg_ );
bsTree_->SetBranchAddress(  "TagMuListSize"	      , &TagMuListSize_);
bsTree_->SetBranchAddress(  "TagMuSimu"	      	      , &TagMuSimu_);	
bsTree_->SetBranchAddress(  "MuRecoMCmother"          , &MuRecoMCmother_ );
bsTree_->SetBranchAddress(  "MuRecoMCMotherMother"    , &MuRecoMCMotherMother_ );

bsTree_->SetBranchAddress(  "TagMuJetdR", &TagMuJetdR_ );


//bsTree_->SetBranchAddress("GlobalTagMuon", &GlobalTagMuon_);
//bsTree_->SetBranchAddress("TrackerTagMuon", &TrackerTagMuon_);
bsTree_->SetBranchAddress("PFTagMuon", &PFTagMuon_);
bsTree_->SetBranchAddress("BsLxy3DMC", &BsLxy3DMC_);
bsTree_->SetBranchAddress("BsLxy2DMC", &BsLxy2DMC_);
bsTree_->SetBranchAddress("BsPMC", &BsPMC_);

bsTree_->SetBranchAddress("BsPVDist2d", &BsPVDist2d_);
bsTree_->SetBranchAddress("BsPVDist3d", &BsPVDist3d_);

bsTree_->SetBranchAddress(  "PVZpos"            , &PVZpos_ );
bsTree_->SetBranchAddress(  "NTracksInPV"       , &NTracksInPV_  ); 
bsTree_->SetBranchAddress( "PVAbsPt"           , &PVAbsPt_); 
bsTree_->SetBranchAddress(  "SVZpos"            , &SVZpos_ );
bsTree_->SetBranchAddress("BsPtMC", &BsPtMC_);
bsTree_->SetBranchAddress("FirstBsMCZpos", &FirstBsMCZpos_);
bsTree_->SetBranchAddress("BsCosThetaReco", &BsCosThetaReco_);


bsTree_->SetBranchAddress("BsCt3DPVCosTheta", &BsCt3DPVCosTheta_ );
bsTree_->SetBranchAddress("BsCt2DPVCosTheta", &BsCt2DPVCosTheta_ );

bsTree_->SetBranchAddress("BsCt3DPVHighestPt", &BsCt3DPVHighestPt_ );
bsTree_->SetBranchAddress("BsCt2DPVHighestPt", &BsCt2DPVHighestPt_ );

bsTree_->SetBranchAddress(  "runNumber"             , &runNumber_  );                    
bsTree_->SetBranchAddress(  "eventNumber"             , &eventNumber_  );                    
bsTree_->SetBranchAddress(  "lumiSection"             , &lumiSection_  );                    

bsTree_->SetBranchAddress(  "PUinteraction"             , &PUinteraction_  );                    

bsTree_->SetBranchAddress(  "ihaveajpsi"             , &ihaveajpsi_);
bsTree_->SetBranchAddress(  "BsCowboy"             , &BsCowboy_);      
bsTree_->SetBranchAddress(  "BdCowboy"             , &BdCowboy_);
bsTree_->SetBranchAddress(  "BsPhiVtxProb"             , &BsPhiVtxProb_);
bsTree_->SetBranchAddress(  "BsMu1QualityG"             , &BsMu1QualityG_);
bsTree_->SetBranchAddress(  "BsMu2QualityG"             , &BsMu2QualityG_);
bsTree_->SetBranchAddress(  "BsMu1QualityT"             , &BsMu1QualityT_);
bsTree_->SetBranchAddress(  "BsMu2QualityT"             , &BsMu2QualityT_);
bsTree_->SetBranchAddress(  "BdMu1QualityG"             , &BdMu1QualityG_);
bsTree_->SetBranchAddress(  "BdMu2QualityG"             , &BdMu2QualityG_);
bsTree_->SetBranchAddress(  "BdMu1QualityT"             , &BdMu1QualityT_);
bsTree_->SetBranchAddress(  "BdMu2QualityT"             , &BdMu2QualityT_); 


bsTree_->SetBranchAddress(  "NVertices"             , &NVertices_  ); 
bsTree_->SetBranchAddress("PVindex", &PVindex_);

bsTree_->SetBranchAddress(  "BBprod", &BBprod_);   
//Jet variables

bsTree_->SetBranchAddress("BsBOldJetdR", &BsBOldJetdR_);

bsTree_->SetBranchAddress("PtJetPt", &PtJetPt_);
bsTree_->SetBranchAddress("BJetPt", &BJetPt_);

bsTree_->SetBranchAddress("BJetEta", &BJetEta_);
bsTree_->SetBranchAddress("BJetPhi", &BJetPhi_);

bsTree_->SetBranchAddress("BsJetdR", &BsJetdR_);
bsTree_->SetBranchAddress("BsJetPx", &BsJetPx_);
bsTree_->SetBranchAddress("BsJetPy", &BsJetPy_);
bsTree_->SetBranchAddress("BsJetPz", &BsJetPz_);

bsTree_->SetBranchAddress("BdJetdR", &BdJetdR_);
bsTree_->SetBranchAddress("BdJetPx", &BdJetPx_);
bsTree_->SetBranchAddress("BdJetPy", &BdJetPy_);
bsTree_->SetBranchAddress("BdJetPz", &BdJetPz_);

bsTree_->SetBranchAddress("BplusPtJetdR", &BplusPtJetdR_);
bsTree_->SetBranchAddress("BplusBJetdR", &BplusBJetdR_);

bsTree_->SetBranchAddress("BsPtJetdR", &BsPtJetdR_);
bsTree_->SetBranchAddress("BsBJetdR", &BsBJetdR_);

bsTree_->SetBranchAddress("BplusJetPx", &BplusJetPx_);
bsTree_->SetBranchAddress("BplusJetPy", &BplusJetPy_);
bsTree_->SetBranchAddress("BplusJetPz", &BplusJetPz_);
bsTree_->SetBranchAddress("JetBTagProb", &JetBTagProb_);
bsTree_->SetBranchAddress("JetBOldTagProb", &JetBOldTagProb_);

bsTree_->SetBranchAddress("BsJetCharge", &BsJetCharge_ );
bsTree_->SetBranchAddress("BdJetCharge", &BdJetCharge_ );
bsTree_->SetBranchAddress("BplusJetCharge", &BplusJetCharge_ );

//Bplus variables


bsTree_->SetBranchAddress("BplusIsPV", &BplusIsPV_);
bsTree_->SetBranchAddress("BplusIsBS", &BplusIsBS_);
bsTree_->SetBranchAddress("BplusPVindex", &BplusPVindex_);
bsTree_->SetBranchAddress("BplusEta" , &BplusEta_ ); 
bsTree_->SetBranchAddress("JpsiMass_bplus" ,     &JpsiMass_bplus_ );
bsTree_->SetBranchAddress("JpsiPt_bplus"   ,     &JpsiPt_bplus_ );
bsTree_->SetBranchAddress("BplusCosTheta"  ,     &BplusCosTheta_ );
bsTree_->SetBranchAddress("BplusDecayChannel" ,       &BplusDecayChannel_);
bsTree_->SetBranchAddress("BplusCtPerSigma2D" ,       &BplusCtPerSigma2D_);
bsTree_->SetBranchAddress("BplusCtErr2D" ,       &BplusCtErr2D_);
bsTree_->SetBranchAddress("BplusCtPerSigma3D"	 , &BplusCtPerSigma3D_  );
bsTree_->SetBranchAddress("BplusVtxProb"	 , &BplusVtxProb_  );
bsTree_->SetBranchAddress("BplusM_fit"	 , &BplusM_fit_    );
bsTree_->SetBranchAddress("BplusChi2"		 , &BplusChi2_     );
bsTree_->SetBranchAddress("KplusPtot"      	 , &KplusPtot_     );
bsTree_->SetBranchAddress("KplusPt"      	 , &KplusPt_       );
bsTree_->SetBranchAddress("BplusLxyz"      	 , &BplusLxyz_     );
bsTree_->SetBranchAddress("BplusLxyzErr"       , &BplusLxyzErr_);
bsTree_->SetBranchAddress("BplusMu1Pt"     ,       &BplusMu1Pt_ );	
bsTree_->SetBranchAddress("BplusMu2Pt"     ,       &BplusMu2Pt_ );
bsTree_->SetBranchAddress("BplusLxy"       ,       &BplusLxy_   );
bsTree_->SetBranchAddress("BplusLxyErr"    ,       &BplusLxyErr_ );
bsTree_->SetBranchAddress("BplusPt"	   ,	   &BplusPt_    );
bsTree_->SetBranchAddress("BplusCt3D"      ,      &BplusCt3D_   );
bsTree_->SetBranchAddress("BplusCtErr3D"   ,      &BplusCtErr3D_);
bsTree_->SetBranchAddress("BplusMu1Ptot"   ,       &BplusMu1Ptot_  );
bsTree_->SetBranchAddress("BplusMu2Ptot"   ,       &BplusMu2Ptot_  );
bsTree_->SetBranchAddress("BplusPtot"      ,       &BplusPtot_  );
bsTree_->SetBranchAddress("BplusCt3D_2"   ,       &BplusCt3D_2_ );
bsTree_->SetBranchAddress("BplusCt2D_2"   ,       &BplusCt2D_2_ );
bsTree_->SetBranchAddress("BplusCt2D"     ,       &BplusCt2D_);
bsTree_->SetBranchAddress("BplusKmcId"     ,       &BplusKmcId_  );
bsTree_->SetBranchAddress("BplusKmomId"    ,       &BplusKmomId_ );
bsTree_->SetBranchAddress("BplusMu1mcId"   ,       &BplusMu1mcId_  );
bsTree_->SetBranchAddress("BplusMu1momId"  ,       &BplusMu1momId_);
bsTree_->SetBranchAddress("BplusMu1gmomId" ,       &BplusMu1gmomId_);
bsTree_->SetBranchAddress("BplusMu2mcId"   ,       &BplusMu2mcId_  );
bsTree_->SetBranchAddress("BplusMu2momId"  ,       &BplusMu2momId_);
bsTree_->SetBranchAddress("BplusMu2gmomId" ,       &BplusMu2gmomId_);
bsTree_->SetBranchAddress("isMatchedBplus" ,       &isMatchedBplus_);
bsTree_->SetBranchAddress("BplusKIP3D"     ,       &BplusKIP3D_    );
bsTree_->SetBranchAddress("BplusKIP3DErr"  ,       &BplusKIP3DErr_ );
bsTree_->SetBranchAddress("IP3DKandJpsiVtx"  ,       &IP3DKandJpsiVtx_ );
bsTree_->SetBranchAddress("IP3DKandJpsiVtxErr"  ,     &IP3DKandJpsiVtxErr_ );

bsTree_->SetBranchAddress("LxyPerSigma2D"  ,     &LxyPerSigma2D_ );
bsTree_->SetBranchAddress("BplusPhi" , &BplusPhi_ ); 

bsTree_->SetBranchAddress("BpmatchDoubleMu01", &BpmatchDoubleMu01_ );
bsTree_->SetBranchAddress("BpmatchDoubleMu02", &BpmatchDoubleMu02_ );
//Bplus variables
  

//new Bd variables 

bsTree_->SetBranchAddress("BdCt3D" , &BdCt3D_ ); 
bsTree_->SetBranchAddress("BdCtErr3D" , &BdCtErr3D_ ); 
bsTree_->SetBranchAddress("BdCt2D" , &BdCt2D_ );
bsTree_->SetBranchAddress("BdCtErr2D" , &BdCtErr2D_ );  
// new Bd variables                 

bsTree_->SetBranchAddress(  "NVertices"             , &NVertices_  ); 

bsTree_->SetBranchAddress(  "triggerbit_HLTmu3Tk"             , &triggerbit_HLTmu3Tk_  );                    
bsTree_->SetBranchAddress(  "triggerbit_HLTmu5"		  , &triggerbit_HLTmu5_  );                                   
bsTree_->SetBranchAddress(  "triggerbit_HLTmu7"		  , &triggerbit_HLTmu7_  );                                   
bsTree_->SetBranchAddress(  "triggerbit_HLTdoubleIsoMu3"	  , &triggerbit_HLTdoubleIsoMu3_  );                           
bsTree_->SetBranchAddress(  "triggerbit_HLTdoubleMu3"	  , &triggerbit_HLTdoubleMu3_  );                 
bsTree_->SetBranchAddress(  "triggerbit_HLTdoubleMu0"	  , &triggerbit_HLTdoubleMu0_  );                        
bsTree_->SetBranchAddress(  "triggerbit_HLTL1DoubleMuOpen"  , &triggerbit_HLTL1DoubleMuOpen_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLTMu0Track0Jpsi"	  , &triggerbit_HLTMu0Track0Jpsi_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLTL1DoubleMuOpenTight"	  , &triggerbit_HLTL1DoubleMuOpenTight_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLT_DoubleMu3_Jpsi_v2"	  , &triggerbit_HLT_DoubleMu3_Jpsi_v2_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLT_DoubleMu3_Jpsi_v2MC"	  , &triggerbit_HLT_DoubleMu3_Jpsi_v2MC_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLT_DoubleMu3_Quarkonium_v2"	  , &triggerbit_HLT_DoubleMu3_Quarkonium_v2_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLT_DoubleMu3_Quarkonium_v2MC"	  , &triggerbit_HLT_DoubleMu3_Quarkonium_v2MC_  );               
bsTree_->SetBranchAddress(  "triggerbit_HLT_DoubleMu3_Quarkonium_v1"	  , &triggerbit_HLT_DoubleMu3_Quarkonium_v1_  );               
bsTree_->SetBranchAddress(  "triggerbit_Jpsi_Displaced_v1"	  , &triggerbit_Jpsi_Displaced_v1_  );               
bsTree_->SetBranchAddress(  "triggerbit_7Jpsi_Displaced_v1"	  , &triggerbit_7Jpsi_Displaced_v1_  );               
bsTree_->SetBranchAddress(  "triggerbit_Jpsi_Displaced_v1MC"	  , &triggerbit_Jpsi_Displaced_v1MC_  );               
bsTree_->SetBranchAddress(  "triggerbit_7Jpsi_Displaced_v2"	  , &triggerbit_7Jpsi_Displaced_v2_  );               
bsTree_->SetBranchAddress(  "triggerbit_7Jpsi_Displaced_v3"	  , &triggerbit_7Jpsi_Displaced_v3_  );               
bsTree_->SetBranchAddress(  "triggerbit_3p5Jpsi_Displaced_v2"	  , &triggerbit_3p5Jpsi_Displaced_v2_  );               
bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v1"	  , &triggerbit_4Jpsi_Displaced_v1_  );               
bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v4"	  , &triggerbit_4Jpsi_Displaced_v4_  );               

bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v5"	  , &triggerbit_4Jpsi_Displaced_v5_  );
bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v9"	  , &triggerbit_4Jpsi_Displaced_v9_  );      

bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v10"	  , &triggerbit_4Jpsi_Displaced_v10_  );              

bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v11"	  , &triggerbit_4Jpsi_Displaced_v11_  );    
bsTree_->SetBranchAddress(  "triggerbit_4Jpsi_Displaced_v12"	  , &triggerbit_4Jpsi_Displaced_v12_  ); 
      

bsTree_->SetBranchAddress(  "triggerbit_5Jpsi_Displaced_v1"	  , &triggerbit_5Jpsi_Displaced_v1_  );               
bsTree_->SetBranchAddress(  "triggerbit_5Jpsi_Displaced_v2"	  , &triggerbit_5Jpsi_Displaced_v2_  );               
bsTree_->SetBranchAddress(  "triggerbit_5Jpsi_Displaced_v4"	  , &triggerbit_5Jpsi_Displaced_v4_  );               
bsTree_->SetBranchAddress(  "triggerbit_5Jpsi_Displaced_v5"	  , &triggerbit_5Jpsi_Displaced_v5_  );               

bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_Muon_v15"	  , &triggerbit_Dimuon0_Jpsi_Muon_v15_); 

bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_Muon_v16"	  , &triggerbit_Dimuon0_Jpsi_Muon_v16_); 

bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_Muon_v17"	  , &triggerbit_Dimuon0_Jpsi_Muon_v17_); 

bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_Muon_v18"	  , &triggerbit_Dimuon0_Jpsi_Muon_v18_); 



bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v1"	  , &triggerbit_Dimuon0_Jpsi_v1_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v3"	  , &triggerbit_Dimuon0_Jpsi_v3_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v5"	  , &triggerbit_Dimuon0_Jpsi_v5_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v6"	  , &triggerbit_Dimuon0_Jpsi_v6_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v9"	  , &triggerbit_Dimuon0_Jpsi_v9_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon0_Jpsi_v10"	  , &triggerbit_Dimuon0_Jpsi_v10_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon10_Barrel"	  , &triggerbit_Dimuon10_Barrel_  );                        
bsTree_->SetBranchAddress(  "triggerbit_Dimuon13_Barrel"	  , &triggerbit_Dimuon13_Barrel_  );                        

bsTree_->SetBranchAddress(  "BSx"				  , &BSx_  );                                                
bsTree_->SetBranchAddress(  "BSy"				  , &BSy_  );                                               
bsTree_->SetBranchAddress(  "BSz"				  , &BSz_  );                                                  
bsTree_->SetBranchAddress(  "BSdx"                                , &BSdx_  );
bsTree_->SetBranchAddress(  "BSdy"                                , &BSdy_  );
bsTree_->SetBranchAddress(  "BSdz"                                , &BSdz_  );
bsTree_->SetBranchAddress(  "BSsigmaZ"                            , &BSsigmaZ_  );
bsTree_->SetBranchAddress(  "BSdsigmaZ"                           , &BSdsigmaZ_  );
bsTree_->SetBranchAddress(  "PVx"				  , &PVx_  );                                                  
bsTree_->SetBranchAddress(  "PVy"				  , &PVy_  );                                                  
bsTree_->SetBranchAddress(  "PVz"				  , &PVz_  );                                                  
bsTree_->SetBranchAddress(  "PVerrx"			  , &PVerrx_  );                                               
bsTree_->SetBranchAddress(  "PVerry"			  , &PVerry_  );                                               
bsTree_->SetBranchAddress(  "PVerrz"			  , &PVerrz_  ); 

bsTree_->SetBranchAddress(  "MuMuDCA"			  , &MuMuDCA_  ); 
bsTree_->SetBranchAddress(  "MuMuDistance"			  , &MuMuDistance_  ); 
bsTree_->SetBranchAddress(  "MuMuDistanceSigma"			  , &MuMuDistanceSigma_  ); 
bsTree_->SetBranchAddress(  "MuDr1"			  , &MuDr1_  ); 
bsTree_->SetBranchAddress(  "MuDr2"			  , &MuDr2_  ); 
bsTree_->SetBranchAddress(  "BdMuMuDCA"			  , &BdMuMuDCA_  ); 
bsTree_->SetBranchAddress(  "BdMuMuDistance"			  , &BdMuMuDistance_  ); 
bsTree_->SetBranchAddress(  "BdMuMuDistanceSigma"			  , &BdMuMuDistanceSigma_  ); 
bsTree_->SetBranchAddress(  "BdMuDr1"			  , &BdMuDr1_  ); 
bsTree_->SetBranchAddress(  "BdMuDr2"			  , &BdMuDr2_  ); 

bsTree_->SetBranchAddress(  "isPV"				  , &isPV_  );                                                  
bsTree_->SetBranchAddress(  "isBS"				  , &isBS_  );                                                  

bsTree_->SetBranchAddress( "PVx_refit"   ,         &PVx_refit_      ); 
bsTree_->SetBranchAddress( "PVy_refit"   ,	 &PVy_refit_        );
bsTree_->SetBranchAddress( "PVz_refit"   ,	 &PVz_refit_        );
bsTree_->SetBranchAddress( "PVerrx_refit",	 &PVerrx_refit_     );
bsTree_->SetBranchAddress( "PVerry_refit",	 &PVerry_refit_     );
bsTree_->SetBranchAddress( "PVerrz_refit",	 &PVerrz_refit_     );

bsTree_->SetBranchAddress(  "PionDeDx"                         , &PionDeDx_  );
bsTree_->SetBranchAddress(  "KaonDeDx"                         , &KaonDeDx_  );
bsTree_->SetBranchAddress(  "KaonPt"                         , &KaonPt_  );
bsTree_->SetBranchAddress(  "PionPt"                         , &PionPt_  );
//bsTree_->SetBranchAddress(  "MuRecoChg2"                         , &MuRecoChg2_  );
bsTree_->SetBranchAddress(  "MuRecoPhi1"                         , &MuRecoPhi1_  );
//bsTree_->SetBranchAddress(  "MuRecoPhi2"                         , &MuRecoPhi2_  );

bsTree_->SetBranchAddress(  "matchL11"                         , &matchL11_  );
bsTree_->SetBranchAddress(  "matchL12"                         , &matchL12_  );
bsTree_->SetBranchAddress(  "match2mu01"                         , &match2mu01_  );
bsTree_->SetBranchAddress(  "match2mu02"                         , &match2mu02_  );
bsTree_->SetBranchAddress(  "match2mu31"                         , &match2mu31_  );
bsTree_->SetBranchAddress(  "match2mu32"                         , &match2mu32_  );
bsTree_->SetBranchAddress(  "match1mu01"                         , &match1mu01_  );
bsTree_->SetBranchAddress(  "match1mu02"                         , &match1mu02_  );
bsTree_->SetBranchAddress(  "matchDoubleMu31J"                         , &matchDoubleMu31J_  );
bsTree_->SetBranchAddress(  "matchDoubleMu32J"                         , &matchDoubleMu32J_  );
bsTree_->SetBranchAddress(  "matchDoubleMu31Q"                         , &matchDoubleMu31Q_  );
bsTree_->SetBranchAddress(  "matchDoubleMu32Q"                         , &matchDoubleMu32Q_  );
bsTree_->SetBranchAddress(  "matchDoubleMu71"                         , &matchDoubleMu71_  );
bsTree_->SetBranchAddress(  "matchDoubleMu72"                         , &matchDoubleMu72_  );
bsTree_->SetBranchAddress(  "matchDoubleMu41"                         , &matchDoubleMu41_  );
bsTree_->SetBranchAddress(  "matchDoubleMu42"                         , &matchDoubleMu42_  );
bsTree_->SetBranchAddress(  "matchDoubleMu51"                         , &matchDoubleMu51_  );
bsTree_->SetBranchAddress(  "matchDoubleMu52"                         , &matchDoubleMu52_  );
bsTree_->SetBranchAddress(  "matchDoubleMu01"                         , &matchDoubleMu01_  );
bsTree_->SetBranchAddress(  "matchDoubleMu02"                         , &matchDoubleMu02_  );
bsTree_->SetBranchAddress(  "matchDoubleMu101"                         , &matchDoubleMu101_  );
bsTree_->SetBranchAddress(  "matchDoubleMu102"                         , &matchDoubleMu102_  );
bsTree_->SetBranchAddress(  "matchDoubleMu131"                         , &matchDoubleMu131_  );
bsTree_->SetBranchAddress(  "matchDoubleMu132"                         , &matchDoubleMu132_  );
bsTree_->SetBranchAddress(  "matchmu0tk01"                         , &matchmu0tk01_  );
bsTree_->SetBranchAddress(  "matchmu0tk02"                         , &matchmu0tk02_  );
bsTree_->SetBranchAddress(  "matchFilterJpsi1"                         , &matchFilterJpsi1_  );
bsTree_->SetBranchAddress(  "matchDoubleMu01DiMuon0", &matchDoubleMu01DiMuon0_  );
bsTree_->SetBranchAddress(  "matchDoubleMu02DiMuon0", &matchDoubleMu02DiMuon0_  );

bsTree_->SetBranchAddress( "BpmatchDoubleMu01DiMuon0", &BpmatchDoubleMu01DiMuon0_ );            
bsTree_->SetBranchAddress( "BpmatchDoubleMu02DiMuon0",&BpmatchDoubleMu02DiMuon0_ ); 
bsTree_->SetBranchAddress("BdmatchDoubleMu01DiMuon0",&BdmatchDoubleMu01DiMuon0_ );            
bsTree_->SetBranchAddress( "BdmatchDoubleMu02DiMuon0",&BdmatchDoubleMu02DiMuon0_ );

bsTree_->SetBranchAddress(  "matchFilterJpsi2"                         , &matchFilterJpsi2_  );
bsTree_->SetBranchAddress(  "BdmatchL11"                         , &BdmatchL11_  );
bsTree_->SetBranchAddress(  "BdmatchL12"                         , &BdmatchL12_  );
bsTree_->SetBranchAddress(  "Bdmatch2mu01"                         , &Bdmatch2mu01_  );
bsTree_->SetBranchAddress(  "Bdmatch2mu02"                         , &Bdmatch2mu02_  );
bsTree_->SetBranchAddress(  "Bdmatch2mu31"                         , &Bdmatch2mu31_  );
bsTree_->SetBranchAddress(  "Bdmatch2mu32"                         , &Bdmatch2mu32_  );
bsTree_->SetBranchAddress(  "Bdmatch1mu01"                         , &Bdmatch1mu01_  );
bsTree_->SetBranchAddress(  "Bdmatch1mu02"                         , &Bdmatch1mu02_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu31Q"                         , &BdmatchDoubleMu31Q_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu32Q"                         , &BdmatchDoubleMu32Q_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu71"                         , &BdmatchDoubleMu71_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu72"                         , &BdmatchDoubleMu72_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu41"                         , &BdmatchDoubleMu41_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu42"                         , &BdmatchDoubleMu42_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu51"                         , &BdmatchDoubleMu51_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu52"                         , &BdmatchDoubleMu52_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu01"                         , &BdmatchDoubleMu01_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu02"                         , &BdmatchDoubleMu02_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu101"                         , &BdmatchDoubleMu101_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu102"                         , &BdmatchDoubleMu102_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu131"                         , &BdmatchDoubleMu131_  );
bsTree_->SetBranchAddress(  "BdmatchDoubleMu132"                         , &BdmatchDoubleMu132_  );




bsTree_->SetBranchAddress( "BdMu1Pt_beffit", &BdMu1Pt_beffit_);
bsTree_->SetBranchAddress( "BdMu1Pz_beffit", &BdMu1Pz_beffit_);
bsTree_->SetBranchAddress( "BdMu1Eta_beffit", &BdMu1Eta_beffit_);
bsTree_->SetBranchAddress( "BdMu1Phi_beffit", &BdMu1Phi_beffit_);
bsTree_->SetBranchAddress( "BdMu2Pt_beffit", &BdMu2Pt_beffit_);
bsTree_->SetBranchAddress( "BdMu2Pz_beffit", &BdMu2Pz_beffit_);
bsTree_->SetBranchAddress( "BdMu2Eta_beffit", &BdMu2Eta_beffit_);
bsTree_->SetBranchAddress( "BdMu2Phi_beffit", &BdMu2Phi_beffit_);
bsTree_->SetBranchAddress( "BdJpsiM_nofit", &BdJpsiM_nofit_);
bsTree_->SetBranchAddress( "BdJpsiEta_nofit", &BdJpsiEta_nofit_);
bsTree_->SetBranchAddress( "BdJpsiPhi_nofit", &BdJpsiPhi_nofit_);
bsTree_->SetBranchAddress( "BdJpsiPt_nofit", &BdJpsiPt_nofit_);
bsTree_->SetBranchAddress( "BdJpsiPz_nofit", &BdJpsiPz_nofit_);

 
bsTree_->SetBranchAddress( "MuonType", &MuonType_ );

bsTree_->SetBranchAddress(  "JpsiVtxProb"			  , &JpsiVtxProb_  );                                          
bsTree_->SetBranchAddress(  "BdJpsiVtxProb"			  , &BdJpsiVtxProb_  );                                          
bsTree_->SetBranchAddress(  "JpsiM_alone"			  , &JpsiM_alone_  );                                          
bsTree_->SetBranchAddress(  "JpsiPhi_alone"		  , &JpsiPhi_alone_  );                                        
bsTree_->SetBranchAddress(  "JpsiEta_alone"		  , &JpsiEta_alone_  );                                        
bsTree_->SetBranchAddress(  "JpsiPt_alone"		  , &JpsiPt_alone_  );                                         
bsTree_->SetBranchAddress(  "JpsiMu1Pt_alone"		  , &JpsiMu1Pt_alone_  );                                      
bsTree_->SetBranchAddress(  "JpsiMu2Pt_alone"		  , &JpsiMu2Pt_alone_  );                                      
bsTree_->SetBranchAddress(  "JpsiMu1Phi_alone"		  , &JpsiMu1Phi_alone_  );                                     
bsTree_->SetBranchAddress(  "JpsiMu2Phi_alone"		  , &JpsiMu2Phi_alone_  );                                     
bsTree_->SetBranchAddress(  "JpsiMu1Eta_alone"		  , &JpsiMu1Eta_alone_  );                                     
bsTree_->SetBranchAddress(  "JpsiMu2Eta_alone"		  , &JpsiMu2Eta_alone_  );                                     
bsTree_->SetBranchAddress(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_  );                                
bsTree_->SetBranchAddress(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_  );                                
bsTree_->SetBranchAddress(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_  );                                
bsTree_->SetBranchAddress(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_  );                                
bsTree_->SetBranchAddress(  "JpsiMu1d0_alone"             , &JpsiMu1d0_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2d0_alone"             , &JpsiMu2d0_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu1dz_alone"             , &JpsiMu1dz_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2dz_alone"             , &JpsiMu2dz_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu1chi2_alone"           , &JpsiMu1chi2_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2chi2_alone"           , &JpsiMu2chi2_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu1ndof_alone"           , &JpsiMu1ndof_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2ndof_alone"           , &JpsiMu2ndof_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu1nHits_alone"                  , &JpsiMu1nHits_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2nHits_alone"                  , &JpsiMu2nHits_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu1nPixHits_alone"               , &JpsiMu1nPixHits_alone_  );
bsTree_->SetBranchAddress(  "JpsiMu2nPixHits_alone"               , &JpsiMu2nPixHits_alone_  );

bsTree_->SetBranchAddress(  "BsFitChi2"			  , &BsFitChi2_  );                                            
bsTree_->SetBranchAddress(  "BsFitNdof"			  , &BsFitNdof_  );                                           
bsTree_->SetBranchAddress(  "BsFitVtxProb"		  , &BsFitVtxProb_  );                          
bsTree_->SetBranchAddress(  "JpsiNumberOfCandidates"        , &JpsiNumberOfCandidates_);               
bsTree_->SetBranchAddress(  "PhiNumberOfCandidatesBeforeFit"        , &PhiNumberOfCandidatesBeforeFit_);
bsTree_->SetBranchAddress(  "BsNumberOfCandidatesBeforeFit"        , &BsNumberOfCandidatesBeforeFit_);
bsTree_->SetBranchAddress(  "BsNumberOfCandidatesAfterFit"        , &BsNumberOfCandidatesAfterFit_);
bsTree_->SetBranchAddress(  "BsNumberOfCandidatesAfterBestFit"        , &BsNumberOfCandidatesAfterBestFit_);

bsTree_->SetBranchAddress(  "K1Pt_beffit"			  , &K1Pt_beffit_  );                                           
bsTree_->SetBranchAddress(  "K1Pz_beffit"			  , &K1Pz_beffit_  );                                           
bsTree_->SetBranchAddress(  "K1Eta_beffit"			  , &K1Eta_beffit_  );                                          
bsTree_->SetBranchAddress(  "K1Phi_beffit"			  , &K1Phi_beffit_  );                                          
bsTree_->SetBranchAddress(  "K2Pt_beffit"			  , &K2Pt_beffit_  );                                           
bsTree_->SetBranchAddress(  "K2Pz_beffit"			  , &K2Pz_beffit_  );                                           
bsTree_->SetBranchAddress(  "K2Eta_beffit"			  , &K2Eta_beffit_  );                                          
bsTree_->SetBranchAddress(  "K2Phi_beffit"			  , &K2Phi_beffit_  );                                          

bsTree_->SetBranchAddress(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_  );                                           
bsTree_->SetBranchAddress(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_  );                                           
bsTree_->SetBranchAddress(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_  );                                          
bsTree_->SetBranchAddress(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_  );                                          
bsTree_->SetBranchAddress(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_  );                                           
bsTree_->SetBranchAddress(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_  );                                           
bsTree_->SetBranchAddress(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_  );                                          
bsTree_->SetBranchAddress(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_  );                                          

bsTree_->SetBranchAddress(  "Mu1d0"             , &Mu1d0_  );
bsTree_->SetBranchAddress(  "Mu2d0"             , &Mu2d0_  );
bsTree_->SetBranchAddress(  "Mu1dz"             , &Mu1dz_  );
bsTree_->SetBranchAddress(  "Mu2dz"             , &Mu2dz_  );

bsTree_->SetBranchAddress(  "BsFitM"			  , &BsFitM_  );                                            
bsTree_->SetBranchAddress(  "K1Pt_fit"			  , &K1Pt_fit_  );                                            
bsTree_->SetBranchAddress(  "K2Pt_fit"			  , &K2Pt_fit_  );                                            
bsTree_->SetBranchAddress(  "PhiM_fit"			  , &PhiM_fit_  );                                            
bsTree_->SetBranchAddress(  "BsFitEta"			  , &BsFitEta_  );                                             
bsTree_->SetBranchAddress(  "BsFitPt"			  , &BsFitPt_  );                                              
bsTree_->SetBranchAddress(  "BsFitPz"			  , &BsFitPz_  );                                              
bsTree_->SetBranchAddress(  "BsFitPhi"			  , &BsFitPhi_  );                                             
bsTree_->SetBranchAddress(  "BsFitVtx_x"			  , &BsFitVtx_x_  );                                           
bsTree_->SetBranchAddress(  "BsFitVtx_y"			  , &BsFitVtx_y_  );                                           
bsTree_->SetBranchAddress(  "BsFitVtx_z"			  , &BsFitVtx_z_  );                                           
bsTree_->SetBranchAddress(  "BsM_nofit"			  , &BsM_nofit_  );                                            
bsTree_->SetBranchAddress(  "BsPhi_nofit"			  , &BsPhi_nofit_  );                                          
bsTree_->SetBranchAddress(  "BsEta_nofit"			  , &BsEta_nofit_  );                                          
bsTree_->SetBranchAddress(  "BsPt_nofit"			  , &BsPt_nofit_  );                                           
bsTree_->SetBranchAddress(  "BsPz_nofit"			  , &BsPz_nofit_  );                                           
bsTree_->SetBranchAddress(  "JpsiM_nofit"			  , &JpsiM_nofit_  );                                          
bsTree_->SetBranchAddress(  "JpsiPhi_nofit"		  , &JpsiPhi_nofit_  );                                        
bsTree_->SetBranchAddress(  "JpsiEta_nofit"		  , &JpsiEta_nofit_  );                                        
bsTree_->SetBranchAddress(  "JpsiPt_nofit"		  , &JpsiPt_nofit_  );                                         
bsTree_->SetBranchAddress(  "JpsiPz_nofit"		  , &JpsiPz_nofit_  );                                         
bsTree_->SetBranchAddress(  "PhiM_nofit"			  , &PhiM_nofit_  );                                           
bsTree_->SetBranchAddress(  "PhiPhi_nofit"		  , &PhiPhi_nofit_  );                                         
bsTree_->SetBranchAddress(  "PhiEta_nofit"		  , &PhiEta_nofit_  );                                         
bsTree_->SetBranchAddress(  "PhiPt_nofit"			  , &PhiPt_nofit_  );                                          
bsTree_->SetBranchAddress(  "PhiPz_nofit"			  , &PhiPz_nofit_  );                                          
bsTree_->SetBranchAddress(  "K1Pt_nofit"			  , &K1Pt_nofit_  );                                           
bsTree_->SetBranchAddress(  "K1Pz_nofit"			  , &K1Pz_nofit_  );                                           
bsTree_->SetBranchAddress(  "K1Eta_nofit"			  , &K1Eta_nofit_  );                                          
bsTree_->SetBranchAddress(  "K1Phi_nofit"			  , &K1Phi_nofit_  );                                          
bsTree_->SetBranchAddress(  "K1Key_nofit"			  , &K1Key_nofit_  );                                           
bsTree_->SetBranchAddress(  "K2Eta_nofit"			  , &K2Eta_nofit_  );                                          
bsTree_->SetBranchAddress(  "K2Pt_nofit"			  , &K2Pt_nofit_  );                                           
bsTree_->SetBranchAddress(  "K2Pz_nofit"			  , &K2Pz_nofit_  );                                           
bsTree_->SetBranchAddress(  "K2Phi_nofit"			  , &K2Phi_nofit_  );                                          
bsTree_->SetBranchAddress(  "K2Key_nofit"			  , &K2Key_nofit_  );                                           
bsTree_->SetBranchAddress(  "K1Chi2"			  , &K1Chi2_  );                                               
bsTree_->SetBranchAddress(  "K1nHits"			  , &K1nHits_  );                                           
bsTree_->SetBranchAddress(  "K2Chi2"			  , &K2Chi2_  );                                               
bsTree_->SetBranchAddress(  "K2nHits"			  , &K2nHits_  );                                              
bsTree_->SetBranchAddress(  "K1pixH"			  , &K1pixH_  );                                             
bsTree_->SetBranchAddress(  "K1trkH"			  , &K1trkH_  );                                            
bsTree_->SetBranchAddress(  "K2pixH"			  , &K2pixH_  );                                            
bsTree_->SetBranchAddress(  "K2trkH"			  , &K2trkH_  );                                            
bsTree_->SetBranchAddress(  "Mu1Chi2"			  , &Mu1Chi2_  );                                              
bsTree_->SetBranchAddress(  "Mu1nHits"			  , &Mu1nHits_  );                                            
bsTree_->SetBranchAddress(  "Mu2Chi2"			  , &Mu2Chi2_  );                                              
bsTree_->SetBranchAddress(  "Mu2nHits"			  , &Mu2nHits_  );                                          
bsTree_->SetBranchAddress(  "Mu1pixH"			  , &Mu1pixH_  );                                          
bsTree_->SetBranchAddress(  "Mu1trkH"			  , &Mu1trkH_  );                                           
bsTree_->SetBranchAddress(  "Mu2pixH"			  , &Mu2pixH_  );                                           
bsTree_->SetBranchAddress(  "Mu2trkH"			  , &Mu2trkH_  );                                           
bsTree_->SetBranchAddress(  "costheta"			  , &costheta_  ); 
bsTree_->SetBranchAddress(  "cospsi"			  , &cospsi_  );                                             
bsTree_->SetBranchAddress(  "phi"				  , &phi_  );                                                  
bsTree_->SetBranchAddress(  "Bdcospsi"			  , &Bdcospsi_  );                                               
bsTree_->SetBranchAddress(  "Bdcostheta"			  , &Bdcostheta_  );                                             
bsTree_->SetBranchAddress(  "Bdphi"				  , &Bdphi_  );                                                  
bsTree_->SetBranchAddress(  "Bdcospsi"			  , &Bdcospsi_  );                                               
bsTree_->SetBranchAddress(  "BdcospsiMC"			  , &BdcospsiMC_  );                                               
bsTree_->SetBranchAddress(  "BdcosthetaMC"			  , &BdcosthetaMC_  );                                             
bsTree_->SetBranchAddress(  "BdphiMC"				  , &BdphiMC_  );                                                  
bsTree_->SetBranchAddress(  "AngleBsDecayLength"		  , &AngleBsDecayLength_  );                                   
bsTree_->SetBranchAddress(  "CosDeltaAlpha"		  , &CosDeltaAlpha_  );                                   
bsTree_->SetBranchAddress(  "BdCosDeltaAlpha"		  , &BdCosDeltaAlpha_  );                                   
bsTree_->SetBranchAddress(  "isMatched"			  , &isMatched_  );                                           
bsTree_->SetBranchAddress(  "isMatchedBd"			  , &isMatchedBd_  );                                         
bsTree_->SetBranchAddress(  "BsLxy"			  , &BsLxy_  );                    
bsTree_->SetBranchAddress(  "BsLxyErr"			  , &BsLxyErr_  );                                                
bsTree_->SetBranchAddress(  "BsErrX"			  , &BsErrX_  );                                               
bsTree_->SetBranchAddress(  "BsErrY"			  , &BsErrY_  );                                               
bsTree_->SetBranchAddress(  "BsErrXY"			  , &BsErrXY_  );                                              
bsTree_->SetBranchAddress(  "BsCt"			  , &BsCt_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr"			  , &BsCtErr_  );                                                 
bsTree_->SetBranchAddress(  "BsCt3D"			  , &BsCt3D_  );                                                 
bsTree_->SetBranchAddress(  "BsCt2D"			  , &BsCt2D_  );                                                 
bsTree_->SetBranchAddress(  "BsCt2DBS"			  , &BsCt2DBS_  );                                                 
bsTree_->SetBranchAddress(  "BdCt2DBS"			  , &BdCt2DBS_  );                                                 
bsTree_->SetBranchAddress(  "BdCt2DMC"			  , &BdCt2DMC_  );                                                 
bsTree_->SetBranchAddress(  "BdCt3DMC"			  , &BdCt3DMC_  );                                                 
bsTree_->SetBranchAddress(  "BsCtMPV"			  , &BsCtMPV_  );                                                 
bsTree_->SetBranchAddress(  "BsCt3Drefit"		  , &BsCt3Drefit_  );                                                 
bsTree_->SetBranchAddress(  "BsCt2Drefit"		  , &BsCt2Drefit_  );                                                 
bsTree_->SetBranchAddress(  "BsCtMPVrefit"		  , &BsCtMPVrefit_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr3D"			  , &BsCtErr3D_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr2D"			  , &BsCtErr2D_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr2DBS"			  , &BsCtErr2DBS_  );                                                 
bsTree_->SetBranchAddress(  "BdCtErr2DBS"			  , &BdCtErr2DBS_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr2D2"		  , &BsCtErr2D2_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErrMPV"		  , &BsCtErrMPV_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr3Drefit"		  , &BsCtErr3Drefit_  );                                                 
bsTree_->SetBranchAddress(  "BsCtErr2Drefit"		  , &BsCtErr2Drefit_  );
bsTree_->SetBranchAddress(  "BsCtErrMPVrefit"		  , &BsCtErrMPVrefit_  );                                                 
bsTree_->SetBranchAddress(  "K1trkLay"			  , &K1trkLay_  );                                             
bsTree_->SetBranchAddress(  "K1muDTh"			  , &K1muDTh_  );                                          
bsTree_->SetBranchAddress(  "K1muCSCh"			  , &K1muCSCh_  );                                           
bsTree_->SetBranchAddress(  "K1muRPCh"			  , &K1muRPCh_  );                                          
bsTree_->SetBranchAddress(  "K2trkLay"			  , &K2trkLay_  );                                          
bsTree_->SetBranchAddress(  "K2muDTh"			  , &K2muDTh_  );                                          
bsTree_->SetBranchAddress(  "K2muCSCh"			  , &K2muCSCh_  );                                           
bsTree_->SetBranchAddress(  "K2muRPCh"			  , &K2muRPCh_  );                                          
bsTree_->SetBranchAddress(  "Mu1trkLay"			  , &Mu1trkLay_  );                                          
bsTree_->SetBranchAddress(  "Mu1muDTh"			  , &Mu1muDTh_  );                                         
bsTree_->SetBranchAddress(  "Mu1muCSCh"			  , &Mu1muCSCh_  );                                          
bsTree_->SetBranchAddress(  "Mu1muRPCh"			  , &Mu1muRPCh_  );                                         
bsTree_->SetBranchAddress(  "Mu2trkLay"			  , &Mu2trkLay_  );                                         
bsTree_->SetBranchAddress(  "Mu2muDTh"			  , &Mu2muDTh_  );                                         
bsTree_->SetBranchAddress(  "Mu2muCSCh"			  , &Mu2muCSCh_  );                                          
bsTree_->SetBranchAddress(  "Mu2muRPCh"			  , &Mu2muRPCh_  );                                         
bsTree_->SetBranchAddress(  "K1mcId"			  , &K1mcId_  );                                             
bsTree_->SetBranchAddress(  "K1momId"			  , &K1momId_  );                                            
bsTree_->SetBranchAddress(  "K1gmomId"			  , &K1gmomId_  );                                           
bsTree_->SetBranchAddress(  "K2mcId"			  , &K2mcId_  );                                          
bsTree_->SetBranchAddress(  "K2momId"			  , &K2momId_  );                                            
bsTree_->SetBranchAddress(  "K2gmomId"			  , &K2gmomId_  );                                           
bsTree_->SetBranchAddress(  "Mu1mcId"			  , &Mu1mcId_  );                                          
bsTree_->SetBranchAddress(  "Mu1momId"			  , &Mu1momId_  );                                           
bsTree_->SetBranchAddress(  "Mu1gmomId"			  , &Mu1gmomId_  );                                          
bsTree_->SetBranchAddress(  "Mu2mcId"			  , &Mu2mcId_  );                                         
bsTree_->SetBranchAddress(  "Mu2momId"			  , &Mu2momId_  );                                           
bsTree_->SetBranchAddress(  "Mu2gmomId"			  , &Mu2gmomId_  );                                          
bsTree_->SetBranchAddress(  "Mu1GlobalMuonPromptTight"    , &Mu1GlobalMuonPromptTight_   );
bsTree_->SetBranchAddress(  "Mu2GlobalMuonPromptTight"    , &Mu2GlobalMuonPromptTight_  );
bsTree_->SetBranchAddress(  "Mu1TrackerMuonArbitrated"    , &Mu1TrackerMuonArbitrated_  );
bsTree_->SetBranchAddress(  "Mu1TMLastStationTight"       , &Mu1TMLastStationTight_  );
bsTree_->SetBranchAddress(  "Mu1TMOneStationTight"         , &Mu1TMOneStationTight_  );
bsTree_->SetBranchAddress(  "Mu1TMLastStationOptimizedLowPtTight", &Mu1TMLastStationOptimizedLowPtTight_  );
bsTree_->SetBranchAddress(  "Mu1TMLastStationAngTight"    , &Mu1TMLastStationAngTight_  );
bsTree_->SetBranchAddress(  "Mu1TMOneStationAngTight"     , &Mu1TMOneStationAngTight_  );
bsTree_->SetBranchAddress(  "Mu1TMLastStationOptimizedBarrelLowPtTight", &Mu1TMLastStationOptimizedBarrelLowPtTight_  );
bsTree_->SetBranchAddress(  "Mu2TrackerMuonArbitrated"    , &Mu2TrackerMuonArbitrated_  );
bsTree_->SetBranchAddress(  "Mu2TMLastStationTight"       , &Mu2TMLastStationTight_  );
bsTree_->SetBranchAddress(  "Mu2TMOneStationTight"         , &Mu2TMOneStationTight_  );
bsTree_->SetBranchAddress(  "Mu2TMLastStationOptimizedLowPtTight", &Mu2TMLastStationOptimizedLowPtTight_  );
bsTree_->SetBranchAddress(  "Mu2TMLastStationAngTight"    , &Mu2TMLastStationAngTight_  );
bsTree_->SetBranchAddress(  "Mu2TMOneStationAngTight"     , &Mu2TMOneStationAngTight_  );
bsTree_->SetBranchAddress(  "Mu2TMLastStationOptimizedBarrelLowPtTight", &Mu2TMLastStationOptimizedBarrelLowPtTight_  );
                    
bsTree_->SetBranchAddress(  "BsDist3d"			  , &BsDist3d_  );                                             
bsTree_->SetBranchAddress(  "BsDist3dErr"			  , &BsDist3dErr_  );                                          
bsTree_->SetBranchAddress(  "BsTime3d"			  , &BsTime3d_  );                                             
bsTree_->SetBranchAddress(  "BsTime3dErr"			  , &BsTime3dErr_  );                                          
bsTree_->SetBranchAddress(  "BsDist2d"			  , &BsDist2d_  );                                             
bsTree_->SetBranchAddress(  "BsDist2dErr"			  , &BsDist2dErr_  );                                          
bsTree_->SetBranchAddress(  "BsTime2d"			  , &BsTime2d_  );                                             
bsTree_->SetBranchAddress(  "BsTime2dErr"			  , &BsTime2dErr_  );                                          
bsTree_->SetBranchAddress(  "dedxTrk"			  , &dedxTrk_  );                                              
bsTree_->SetBranchAddress(  "errdedxTrk"			  , &errdedxTrk_  );                                           
bsTree_->SetBranchAddress(  "numdedxTrk"			  , &numdedxTrk_  );                                         
bsTree_->SetBranchAddress(  "iPassedCutIdent"		  , &iPassedCutIdent_  );                                    
bsTree_->SetBranchAddress(  "iPassedCutIdentBd"		  , &iPassedCutIdentBd_  );                                   
bsTree_->SetBranchAddress(  "BdTrack1Charge"		  , &BdTrack1Charge_  );                                   
bsTree_->SetBranchAddress(  "K1Fit_par"			  , K1Fit_par_  );                                         
bsTree_->SetBranchAddress(  "K2Fit_par"			  , K2Fit_par_  );                                     
bsTree_->SetBranchAddress(  "K1Fit_sigX"			  , &K1Fit_sigX_  );                                           
bsTree_->SetBranchAddress(  "K1Fit_sigY"			  , &K1Fit_sigY_  );                                           
bsTree_->SetBranchAddress(  "K1Fit_sigZ"			  , &K1Fit_sigZ_  );                                           
bsTree_->SetBranchAddress(  "K2Fit_sigX"			  , &K2Fit_sigX_  );                                           
bsTree_->SetBranchAddress(  "K2Fit_sigY"			  , &K2Fit_sigY_  );                                           
bsTree_->SetBranchAddress(  "K2Fit_sigZ"			  , &K2Fit_sigZ_  );                                           
bsTree_->SetBranchAddress(  "K1Fit_sigPX"			  , &K1Fit_sigPX_  );                                          
bsTree_->SetBranchAddress(  "K1Fit_sigPY"			  , &K1Fit_sigPY_  );                                          
bsTree_->SetBranchAddress(  "K1Fit_sigPZ"			  , &K1Fit_sigPZ_  );                                          
bsTree_->SetBranchAddress(  "K2Fit_sigPX"			  , &K2Fit_sigPX_  );                                          
bsTree_->SetBranchAddress(  "K2Fit_sigPY"			  , &K2Fit_sigPY_  );                                          
bsTree_->SetBranchAddress(  "K2Fit_sigPZ"			  , &K2Fit_sigPZ_  );                                          
              
bsTree_->SetBranchAddress(  "GenNumberOfBdecays"		  , &GenNumberOfBdecays_  );                                   
bsTree_->SetBranchAddress(  "BmesonsId"			  , BmesonsId_  );                                
bsTree_->SetBranchAddress(  "BDauIdMC"			  , BDauIdMC_  );                                     
bsTree_->SetBranchAddress(  "BDauDauIdMC"			  , BDauDauIdMC_  );                                
bsTree_->SetBranchAddress(  "GenNumberOfDaughters"	  , GenNumberOfDaughters_  );                           
bsTree_->SetBranchAddress(  "GenNumberOfDaughtersDaughters" , GenNumberOfDaughtersDaughters_  );                
bsTree_->SetBranchAddress(  "BDauMMC"			  , BDauMMC_  );                          
bsTree_->SetBranchAddress(  "BDauPtMC"			  , BDauPtMC_  );                                     
bsTree_->SetBranchAddress(  "BDauPzMC"			  , BDauPzMC_  );                                 
bsTree_->SetBranchAddress(  "BDauEtaMC"			  , BDauEtaMC_  );                                 
bsTree_->SetBranchAddress(  "BDauPhiMC"			  , BDauPhiMC_  );                                
bsTree_->SetBranchAddress(  "BDauDauMMC"			  , BDauDauMMC_  );                                 
bsTree_->SetBranchAddress(  "BDauDauPtMC"			  , BDauDauPtMC_  );                             
bsTree_->SetBranchAddress(  "BDauDauPzMC"			  , BDauDauPzMC_  );                            
bsTree_->SetBranchAddress(  "BDauDauEtaMC"		  , BDauDauEtaMC_  );                            
bsTree_->SetBranchAddress(  "BDauDauPhiMC"		  , BDauDauPhiMC_  );                           
bsTree_->SetBranchAddress(  "BMMC"			  , BMMC_  );                             
bsTree_->SetBranchAddress(  "BsCt2DMC"			  , &BsCt2DMC_  );                             
bsTree_->SetBranchAddress(  "BsCt3DMC"			  , &BsCt3DMC_  );                             
bsTree_->SetBranchAddress(  "BsIniFlavour"			  , &BsIniFlavour_  );                             
bsTree_->SetBranchAddress(  "BdIniFlavour"			  , &BdIniFlavour_  );                             
bsTree_->SetBranchAddress(  "ChannelID"			  , &ChannelID_  );                             
bsTree_->SetBranchAddress(  "BdChannelID"			  , &BdChannelID_  );                             
bsTree_->SetBranchAddress(  "BPtMC"			  , BPtMC_  );                                         
bsTree_->SetBranchAddress(  "BPxMC"			  , BPxMC_  );                                        
bsTree_->SetBranchAddress(  "BPyMC"			  , BPyMC_  );
bsTree_->SetBranchAddress(  "BPzMC"			  , BPzMC_  );                                        
bsTree_->SetBranchAddress(  "BEtaMC"			  , BEtaMC_  );                                        
bsTree_->SetBranchAddress(  "BPhiMC"			  , BPhiMC_  );        
 bsTree_->SetBranchAddress(  "costhetaMC"		  , costhetaMC_);        
 bsTree_->SetBranchAddress(  "phiMC"			  , phiMC_);        
 bsTree_->SetBranchAddress(  "cospsiMC"			  , cospsiMC_);        
 bsTree_->SetBranchAddress(  "BscosthetaMC"		  , &BscosthetaMC_);        
 bsTree_->SetBranchAddress(  "BsphiMC"			  , &BsphiMC_);        
 bsTree_->SetBranchAddress(  "BscospsiMC"			  , &BscospsiMC_);        

bsTree_->SetBranchAddress(  "BVtxMC_x" , BVtxMC_x_); 
bsTree_->SetBranchAddress(  "BVtxMC_y" , BVtxMC_y_); 
bsTree_->SetBranchAddress(  "BVtxMC_z" , BVtxMC_z_); 
bsTree_->SetBranchAddress(  "BSVtxMC_x", BSVtxMC_x_);
bsTree_->SetBranchAddress(  "BSVtxMC_y", BSVtxMC_y_);
bsTree_->SetBranchAddress(  "BSVtxMC_z", BSVtxMC_z_);
bsTree_->SetBranchAddress(  "BLxy_MC"  , BLxy_MC_);  
bsTree_->SetBranchAddress(  "BCt_MC"   , BCt_MC_);
 bsTree_->SetBranchAddress(  "BCt_MC2D"   , BCt_MC2D_);   
 bsTree_->SetBranchAddress(  "BCt_MC3D"   , BCt_MC3D_);   

                               
bsTree_->SetBranchAddress(  "genBsVtx_z"			  , &genBsVtx_z_  );                                           
bsTree_->SetBranchAddress(  "genBsVtx_y"			  , &genBsVtx_y_  );                                           
bsTree_->SetBranchAddress(  "genBsVtx_x"			  , &genBsVtx_x_  );                                           
bsTree_->SetBranchAddress(  "genBsSVtx_z"			  , &genBsSVtx_z_  );                                          
bsTree_->SetBranchAddress(  "genBsSVtx_y" 		  , &genBsSVtx_y_  );                                          
bsTree_->SetBranchAddress(  "genBsSVtx_x"			  , &genBsSVtx_x_  );                                          
bsTree_->SetBranchAddress(  "isGenJpsiEvent"		  , &isGenJpsiEvent_  );                                      
bsTree_->SetBranchAddress(  "BdFitChi2_Hyp1"		  , &BdFitChi2_Hyp1_  );                                       
bsTree_->SetBranchAddress(  "BdFitNdof_Hyp1"		  , &BdFitNdof_Hyp1_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtxProb_Hyp1"		  , &BdFitVtxProb_Hyp1_  );                                    
bsTree_->SetBranchAddress(  "BdFitM_Hyp1"			  , &BdFitM_Hyp1_  );                                          
bsTree_->SetBranchAddress(  "BdFitEta_Hyp1"		  , &BdFitEta_Hyp1_  );                                        
bsTree_->SetBranchAddress(  "BdFitPt_Hyp1"		  , &BdFitPt_Hyp1_  );                                         
bsTree_->SetBranchAddress(  "BdFitPz_Hyp1"		  , &BdFitPz_Hyp1_  );                                         
bsTree_->SetBranchAddress(  "BdFitPhi_Hyp1"		  , &BdFitPhi_Hyp1_  );                                        
bsTree_->SetBranchAddress(  "BdFitVtx_x_Hyp1"		  , &BdFitVtx_x_Hyp1_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtx_y_Hyp1"		  , &BdFitVtx_y_Hyp1_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtx_z_Hyp1"		  , &BdFitVtx_z_Hyp1_  );                                      
bsTree_->SetBranchAddress(  "BdM_nofit"		  , &BdM_nofit_  );                                       
bsTree_->SetBranchAddress(  "BdPhi_nofit"		  , &BdPhi_nofit_  );                                     
bsTree_->SetBranchAddress(  "BdEta_nofit"		  , &BdEta_nofit_  );                                     
bsTree_->SetBranchAddress(  "BdPt_nofit"		  , &BdPt_nofit_  );                                      
bsTree_->SetBranchAddress(  "BdPz_nofit"		  , &BdPz_nofit_  );                                      
bsTree_->SetBranchAddress(  "KstarMass_nofit_Hyp1"	  , &KstarMass_nofit_Hyp1_  );  
bsTree_->SetBranchAddress(  "KstarMass_nofit_Hyp2"	  , &KstarMass_nofit_Hyp2_  );                                
bsTree_->SetBranchAddress(  "BdK1_kpi_par_Hyp1"		  , BdK1_kpi_par_Hyp1_  );                                 
bsTree_->SetBranchAddress(  "BdK2_kpi_par_Hyp1"		  , BdK2_kpi_par_Hyp1_  );                             
bsTree_->SetBranchAddress(  "BdK1_kpi_sigX_Hyp1"		  , &BdK1_kpi_sigX_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigY_Hyp1"		  , &BdK1_kpi_sigY_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigZ_Hyp1"		  , &BdK1_kpi_sigZ_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigX_Hyp1"		  , &BdK2_kpi_sigX_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigY_Hyp1"		  , &BdK2_kpi_sigY_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigZ_Hyp1"		  , &BdK2_kpi_sigZ_Hyp1_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPX_Hyp1"		  , &BdK1_kpi_sigPX_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPY_Hyp1"		  , &BdK1_kpi_sigPY_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPZ_Hyp1"		  , &BdK1_kpi_sigPZ_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPX_Hyp1"		  , &BdK2_kpi_sigPX_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPY_Hyp1"		  , &BdK2_kpi_sigPY_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPZ_Hyp1"		  , &BdK2_kpi_sigPZ_Hyp1_  );                                  
bsTree_->SetBranchAddress(  "BdFitChi2_Hyp2"		  , &BdFitChi2_Hyp2_  );                                       
bsTree_->SetBranchAddress(  "BdFitNdof_Hyp2"		  , &BdFitNdof_Hyp2_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtxProb_Hyp2"		  , &BdFitVtxProb_Hyp2_  );                                    
bsTree_->SetBranchAddress(  "BdFitM_Hyp2"			  , &BdFitM_Hyp2_  );                                          
bsTree_->SetBranchAddress(  "BdFitEta_Hyp2"		  , &BdFitEta_Hyp2_  );                                        
bsTree_->SetBranchAddress(  "BdFitPt_Hyp2"		  , &BdFitPt_Hyp2_  );                                         
bsTree_->SetBranchAddress(  "BdFitPz_Hyp2"		  , &BdFitPz_Hyp2_  );                                         
bsTree_->SetBranchAddress(  "BdFitPhi_Hyp2"		  , &BdFitPhi_Hyp2_  );                                        
bsTree_->SetBranchAddress(  "BdFitVtx_x_Hyp2"		  , &BdFitVtx_x_Hyp2_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtx_y_Hyp2"		  , &BdFitVtx_y_Hyp2_  );                                      
bsTree_->SetBranchAddress(  "BdFitVtx_z_Hyp2"		  , &BdFitVtx_z_Hyp2_  );                                      
                  
bsTree_->SetBranchAddress(  "BdPVx_refit",      &BdPVx_refit_   );     
bsTree_->SetBranchAddress(  "BdPVy_refit",      &BdPVy_refit_   );     
bsTree_->SetBranchAddress(  "BdPVz_refit",      &BdPVz_refit_   );     
bsTree_->SetBranchAddress(  "BdPVerrx_refit",   &BdPVerrx_refit_);     
bsTree_->SetBranchAddress(  "BdPVerry_refit",   &BdPVerry_refit_);     
bsTree_->SetBranchAddress(  "BdPVerrz_refit",   &BdPVerrz_refit_);     
bsTree_->SetBranchAddress(  "BdNumberOfCandidates"        , &BdNumberOfCandidates_);
bsTree_->SetBranchAddress(  "BdK1_kpi_par_Hyp2"		  , BdK1_kpi_par_Hyp2_  );                                 
bsTree_->SetBranchAddress(  "BdK2_kpi_par_Hyp2"		  , BdK2_kpi_par_Hyp2_  );                             
bsTree_->SetBranchAddress(  "BdK1_kpi_sigX_Hyp2"		  , &BdK1_kpi_sigX_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigY_Hyp2"		  , &BdK1_kpi_sigY_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigZ_Hyp2"		  , &BdK1_kpi_sigZ_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigX_Hyp2"		  , &BdK2_kpi_sigX_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigY_Hyp2"		  , &BdK2_kpi_sigY_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK2_kpi_sigZ_Hyp2"		  , &BdK2_kpi_sigZ_Hyp2_  );                                   
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPX_Hyp2"		  , &BdK1_kpi_sigPX_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPY_Hyp2"		  , &BdK1_kpi_sigPY_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK1_kpi_sigPZ_Hyp2"		  , &BdK1_kpi_sigPZ_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPX_Hyp2"		  , &BdK2_kpi_sigPX_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPY_Hyp2"		  , &BdK2_kpi_sigPY_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK2_kpi_sigPZ_Hyp2"		  , &BdK2_kpi_sigPZ_Hyp2_  );                                  
bsTree_->SetBranchAddress(  "BdK1Pt_nofit" 		  , &BdK1Pt_nofit_  );                                         
bsTree_->SetBranchAddress(  "BdK1Pz_nofit" 		  , &BdK1Pz_nofit_  );                                         
bsTree_->SetBranchAddress(  "BdK1Eta_nofit" 		  , &BdK1Eta_nofit_  );                                        
bsTree_->SetBranchAddress(  "BdK1Phi_nofit" 		  , &BdK1Phi_nofit_  );  
bsTree_->SetBranchAddress(  "BdK1Key_nofit" 		  , &BdK1Key_nofit_  );                                                  
bsTree_->SetBranchAddress(  "BdK2Pt_nofit" 		  , &BdK2Pt_nofit_  );                                         
bsTree_->SetBranchAddress(  "BdK2Pz_nofit" 		  , &BdK2Pz_nofit_  );                                         
bsTree_->SetBranchAddress(  "BdK2Eta_nofit" 		  , &BdK2Eta_nofit_  );                                        
bsTree_->SetBranchAddress(  "BdK2Phi_nofit" 		  , &BdK2Phi_nofit_  );                                        
bsTree_->SetBranchAddress(  "BdK2Key_nofit" 		  , &BdK2Key_nofit_  );                                                  
bsTree_->SetBranchAddress(  "BdLxy"			  , &BdLxy_  );                                                
bsTree_->SetBranchAddress(  "BdLxyErr"			  , &BdLxyErr_  );                                                
bsTree_->SetBranchAddress(  "BdErrX"			  , &BdErrX_  );                                               
bsTree_->SetBranchAddress(  "BdErrY"			  , &BdErrY_  );                                               
bsTree_->SetBranchAddress(  "BdErrXY"			  , &BdErrXY_  );                                              
bsTree_->SetBranchAddress(  "BdCt"			  , &BdCt_  );                                                 
bsTree_->SetBranchAddress(  "BdCtErr"			  , &BdCtErr_  );                                                 
bsTree_->SetBranchAddress(  "BdDist3d"			  , &BdDist3d_  );                                             
bsTree_->SetBranchAddress(  "BdDist3dErr"			  , &BdDist3dErr_  );                                          
bsTree_->SetBranchAddress(  "BdTime3d"			  , &BdTime3d_  );                                             
bsTree_->SetBranchAddress(  "BdTime3dErr"			  , &BdTime3dErr_  );                                          
bsTree_->SetBranchAddress(  "BdDist2d"			  , &BdDist2d_  );                                             
bsTree_->SetBranchAddress(  "BdDist2dErr"			  , &BdDist2dErr_  );                                          
bsTree_->SetBranchAddress(  "BdTime2d"			  , &BdTime2d_  );                                             
bsTree_->SetBranchAddress(  "BdTime2dErr"                   , &BdTime2dErr_  );                                             

 bsTree_->SetBranchAddress(  "BdK1mcId"           ,       &BdK1mcId_     );
 bsTree_->SetBranchAddress(  "BdK1momId"	  ,	  &BdK1momId_    );
 bsTree_->SetBranchAddress(  "BdK1gmomId"	  ,	  &BdK1gmomId_   );
 bsTree_->SetBranchAddress(  "BdK2mcId"	          ,	  &BdK2mcId_     );
 bsTree_->SetBranchAddress(  "BdK2momId"	  ,	  &BdK2momId_    );
 bsTree_->SetBranchAddress(  "BdK2gmomId"	  ,	  &BdK2gmomId_   );
 bsTree_->SetBranchAddress(  "BdMu1mcId"	  ,	  &BdMu1mcId_    );
 bsTree_->SetBranchAddress(  "BdMu1momId"	  ,	  &BdMu1momId_   );
 bsTree_->SetBranchAddress(  "BdMu1gmomId"        ,	  &BdMu1gmomId_  );
 bsTree_->SetBranchAddress(  "BdMu2mcId"	  ,	  &BdMu2mcId_    );
 bsTree_->SetBranchAddress(  "BdMu2momId"	  ,	  &BdMu2momId_   );
 bsTree_->SetBranchAddress(  "BdMu2gmomId"        ,	  &BdMu2gmomId_  );


}
