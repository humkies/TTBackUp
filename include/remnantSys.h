
// STL includes
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>

// Root includes
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"


class HistSummary
    {
    public:
      
      void bookHist(const char *, const int&, const char *);
      //In case histogram title need to be in certain sequence
      TString Title(TString spec, unsigned int i);

      TH1D *hRSysOldTag_dR, *hRSysNewTag_dR;
      TH1D *hRSysOldTag_dPt, *hRSysNewTag_dPt;
      TH1D *hRSysOldTag_dPhi, *hRSysNewTag_dPhi;      

      TH1D *hRSysOldTagLep_dR, *hRSysNewTagLep_dR;
      TH1D *hRSysOldTagLep_dPt, *hRSysNewTagLep_dPt;
      TH1D *hRSysOldTagLep_dPhi, *hRSysNewTagLep_dPhi;

      TH1D *hRSysOldTagHad_dR, *hRSysNewTagHad_dR;
      TH1D *hRSysOldTagHad_dPt, *hRSysNewTagHad_dPt;
      TH1D *hRSysOldTagHad_dPhi, *hRSysNewTagHad_dPhi;

      TH1D *hDrMinMonoJetNob, *hDrMinMonoJetb, *hDrMinDiJetNob, *hDrMinDiJetb;
      TH1D *hDPtMonoJetNob, *hDPtMonoJetb, *hDPtDiJetNob, *hDPtDiJetb;
      TH1D *hDPhiMonoJetNob, *hDPhiMonoJetb, *hDPhiDiJetNob, *hDPhiDiJetb;

      TH1D *hDrMinMonoJet_Old, *hDrMinDiJet_Old;
      TH1D *hDPhiMonoJet_Old, *hDPhiDiJet_Old;
      TH1D *hDPtMonoJet_Old, *hDPtDiJet_Old;
      //Matched and unmatcahed top quarks
      TH1D *h1MatchedDrMin_New, *h1MatchedDrMin_Old;
      TH1D *h1MatchedDPhi_New, *h1MatchedDPhi_Old;
      TH1D *h1MatchedDPt_New, *h1MatchedDPt_Old;

      TH1D *h1UnMatchedDrMin_New, *h1UnMatchedDrMin_Old;
      TH1D *h1UnMatchedDPhi_New, *h1UnMatchedDPhi_Old;
      TH1D *h1UnMatchedDPt_New, *h1UnMatchedDPt_Old;

      TH1D *h1MT2LepDecay_New, *h1MT2HadDecay_New, *h1MT2IncDecay_New;
      TH1D *h1MT2LepDecay_Old, *h1MT2HadDecay_Old, *h1MT2IncDecay_Old;
 
      TH1D *h1MT2Matched_New, *h1MT2UnMatched_New; 
      TH1D *h1MT2Matched_Old, *h1MT2UnMatched_Old;
      TH1D *h1MT2DiJetNoB, *h1MT2DiJetB;
      TH1D *h1MT2MonoJetNoB, *h1MT2MonoJetB;
      TH1D *h1MT2MonoJet_Old, *h1MT2DiJet_Old;

      TH1D *rSysConst_New_Pt, *rSysConst_Old_Pt;
      TFile *outFileName;



    };


void HistSummary::bookHist(const char * outfilename, const int& filerun, const char *spec)
{
  TString filename(outfilename);
  TString specT(spec);
  TString index(std::to_string(filerun));
  filename+= "_remnantSys_"+specT+"_"+index+".root";

  outFileName = new TFile(filename, "RECREATE");

  rSysConst_New_Pt = new TH1D("h1rSysConst_New_Pt", "rSysConst_New_Pt;P{T};#Events", 100, 0, 1000);
  rSysConst_New_Pt->Sumw2();
  rSysConst_Old_Pt = new TH1D("h1rSysConst_Old_Pt", "rSysConst_Old_Pt;P{T};#Events", 100, 0, 1000);
  rSysConst_Old_Pt->Sumw2();



  //Matched
  h1MatchedDrMin_New = new TH1D("h1MatchedDrMin_New", "dRmin(MatchedTop, RSys);#Delta R;Events", 45, 0, 4.5);
  h1MatchedDrMin_Old = new TH1D("h1MatchedDrMin_Old", "dRmin(MatchedTop, RSys);#Delta R;Events", 45, 0, 4.5);
  h1MatchedDrMin_New->Sumw2();
  h1MatchedDrMin_Old->Sumw2();
  h1MatchedDPt_New = new TH1D("h1MatchedDPt_New", "dPt(MatchedTop, RSys);#Delta R;Events", 100, -600, 600);
  h1MatchedDPt_Old = new TH1D("h1MatchedDPt_Old", "dPt(MatchedTop, RSys);#Delta R;Events", 100, -600, 600);
  h1MatchedDPt_New->Sumw2();
  h1MatchedDPt_Old->Sumw2();
  h1MatchedDPhi_New = new TH1D("h1MatchedDPhi_New", "dPhi(MatchedTop, RSys);#Delta R;Events",100, -3.2, 3.2);
  h1MatchedDPhi_Old = new TH1D("h1MatchedDPhi_Old", "dPhi(MatchedTop, RSys);#Delta R;Events", 100, -3.2, 3.2);
  h1MatchedDPhi_New->Sumw2();
  h1MatchedDPhi_Old->Sumw2();
  //Unmatched
  h1UnMatchedDrMin_New = new TH1D("h1UnMatchedDrMin_New", "dRmin(UnMatchedTop, RSys);#Delta R;Events", 45, 0, 4.5);
  h1UnMatchedDrMin_Old = new TH1D("h1UnMatchedDrMin_Old", "dRmin(UnMatchedTop, RSys);#Delta R;Events", 45, 0, 4.5);
  h1UnMatchedDrMin_New->Sumw2();
  h1UnMatchedDrMin_Old->Sumw2();
  h1UnMatchedDPt_New = new TH1D("h1UnMatchedDPt_New", "dPt(UnMatchedTop, RSys);#Delta R;Events", 100, -600, 600);
  h1UnMatchedDPt_Old = new TH1D("h1UnMatchedDPt_Old", "dPt(UnMatchedTop, RSys);#Delta R;Events", 100, -600, 600);
  h1UnMatchedDPt_New->Sumw2();
  h1UnMatchedDPt_Old->Sumw2();
  h1UnMatchedDPhi_New = new TH1D("h1UnMatchedDPhi_New", "dPhi(UnMatchedTop, RSys);#Delta R;Events",100, -3.2, 3.2);
  h1UnMatchedDPhi_Old = new TH1D("h1UnMatchedDPhi_Old", "dPhi(UnMatchedTop, RSys);#Delta R;Events", 100, -3.2, 3.2);
  h1UnMatchedDPhi_New->Sumw2();
  h1UnMatchedDPhi_Old->Sumw2();

  hRSysOldTag_dR =  new TH1D("h1RSysOldTag_dR", "dRmin(genTop, RSys);#Delta R;Events", 45, 0, 4.5);
  hRSysNewTag_dR =  new TH1D("h1RSysNewTag_dR", "dRMin(genTop, RSys);#Delta R;Events",45, 0, 4.5);
  hRSysOldTag_dR->Sumw2();
  hRSysNewTag_dR->Sumw2();
  hRSysOldTag_dPt =  new TH1D("h1RSysOldTag_dPt", "dPtMin(genTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysNewTag_dPt =  new TH1D("h1RSysNewTag_dPt", "dPtMin(genTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysOldTag_dPt->Sumw2();
  hRSysNewTag_dPt->Sumw2();
  hRSysOldTag_dPhi =  new TH1D("h1RSysOldTag_dPhi", "dPhiMin(genTop, RSys);#Delta#Phi;Events", 100, -3.2, 3.2);
  hRSysNewTag_dPhi =  new TH1D("h1RSysNewTag_dPhi", "dPhiMin(genTop, RSys);#Delta#Phi;Events",100, -3.2, 3.2);
  hRSysOldTag_dPhi->Sumw2();
  hRSysNewTag_dPhi->Sumw2();

  hRSysOldTagLep_dR=  new TH1D("h1RSysOldTagLep_dR", "dR(lepTop, RSys);#Delta R;Events", 45, 0, 4.5);
  hRSysNewTagLep_dR=  new TH1D("h1RSysNewTagLep_dR", "dR(lepTop, RSys);#Delta R;Events",45, 0, 4.5);
  hRSysOldTagLep_dR->Sumw2();
  hRSysNewTagLep_dR->Sumw2();
  hRSysOldTagLep_dPt =  new TH1D("h1RSysOldTagLep_dPt", "dPt(lepTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysNewTagLep_dPt =  new TH1D("h1RSysNewTagLep_dPt", "dPt(lepTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysOldTagLep_dPt->Sumw2();
  hRSysNewTagLep_dPt->Sumw2();
  hRSysOldTagLep_dPhi =  new TH1D("h1RSysOldTagLep_dPhi", "dPhi(lepTop, RSys);#Delta#Phi;Events", 100, -3.2, 3.2);
  hRSysNewTagLep_dPhi =  new TH1D("h1RSysNewTagLep_dPhi", "dPhi(lepTop, RSys);#Delta#Phi;Events",100, -3.2, 3.2);
  hRSysOldTagLep_dPhi->Sumw2();
  hRSysNewTagLep_dPhi->Sumw2();

  hRSysOldTagHad_dR=  new TH1D("h1RSysOldTagHad_dR", "dR(hadTop, RSys);#Delta R;Events", 45, 0, 4.5);
  hRSysNewTagHad_dR=  new TH1D("h1RSysNewTagHad_dR", "dR(hadTop, RSys);#Delta R;Events",45, 0, 4.5);
  hRSysOldTagHad_dR->Sumw2();
  hRSysNewTagHad_dR->Sumw2();
  hRSysOldTagHad_dPt =  new TH1D("h1RSysOldTagHad_dPt", "dPt(hadTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysNewTagHad_dPt =  new TH1D("h1RSysNewTagHad_dPt", "dPt(hadTop, RSys);#Delta P_{T};Events",100, -600, 600);
  hRSysOldTagHad_dPt->Sumw2();
  hRSysNewTagHad_dPt->Sumw2();
  hRSysOldTagHad_dPhi =  new TH1D("h1RSysOldTagHad_dPhi", "dPhi(hadTop, RSys);#Delta#Phi;Events", 100, -3.2, 3.2);
  hRSysNewTagHad_dPhi =  new TH1D("h1RSysNewTagHad_dPhi", "dPhi(hadTop, RSys);#Delta#Phi;Events",100, -3.2, 3.2);
  hRSysOldTagHad_dPhi->Sumw2();
  hRSysNewTagHad_dPhi->Sumw2();

  hDrMinMonoJetNob = new TH1D("h1dRMinMonoJetNob", "dRmin(genTop, MonoJetsNob); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinMonoJetNob->Sumw2();
  hDrMinDiJetNob = new TH1D("h1dRMinDiJetNob", "dRmin(genTop, DiJetsNob); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinDiJetNob->Sumw2();
  hDrMinMonoJetb = new TH1D("h1dRMinMonoJetb", "dRmin(genTop, MonoJetsb); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinMonoJetb->Sumw2();
  hDrMinDiJetb = new TH1D("h1dRMinDiJetb", "dRmin(genTop, DiJetsb); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinDiJetb->Sumw2();

  hDPtMonoJetNob = new TH1D("h1dPtMonoJetNob", "dPt(genTop, MonoJetsNob); #Delta P_{T}; Events", 100, -600, 600);
  hDPtMonoJetNob->Sumw2();
  hDPtDiJetNob = new TH1D("h1dPtDiJetNob", "dPt(genTop, DiJetsNob); #Delta P_{T}; Events", 100, -600, 600);
  hDPtDiJetNob->Sumw2();
  hDPtMonoJetb = new TH1D("h1dPtMonoJetb", "dPt(genTop, MonoJetsb); #Delta P_{T}; Events", 100, -600, 600);
  hDPtMonoJetb->Sumw2();
  hDPtDiJetb = new TH1D("h1dPtDiJetb", "dPt(genTop, DiJetsb); #Delta P_{T}; Events", 100, -600, 600);
  hDPtDiJetb->Sumw2(); 

  hDPhiMonoJetNob = new TH1D("h1dPhiMonoJetNob", "dPhi(genTop, MonoJetsNob); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiMonoJetNob->Sumw2();
  hDPhiDiJetNob = new TH1D("h1dPhiDiJetNob", "dPhi(genTop, DiJetsNob); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiDiJetNob->Sumw2();

  hDPhiMonoJetb = new TH1D("h1dPhiMonoJetb", "dPhi(genTop, MonoJetsb); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiMonoJetb->Sumw2();
  hDPhiDiJetb = new TH1D("h1dPhiDiJetb", "dPhi(genTop, DiJetsb); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiDiJetb->Sumw2();

  hDrMinMonoJet_Old = new TH1D("h1dRMinMonoJet_Old", "dRmin(genTop, MonoJets_Old); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinMonoJet_Old->Sumw2();
  hDrMinDiJet_Old = new TH1D("h1dRMinDiJet_Old", "dRmin(genTop, DiJets_Old); #Delta R_{min}; Events", 45, 0, 4.5);
  hDrMinDiJet_Old->Sumw2();

  hDPtMonoJet_Old = new TH1D("h1dPtMonoJet_Old", "dPt(genTop, MonoJets_Old); #Delta P_{T}; Events", 100, -600, 600);
  hDPtMonoJet_Old->Sumw2();
  hDPtDiJet_Old = new TH1D("h1dPtDiJet_Old", "dPt(genTop, DiJets_Old); #Delta P_{T}; Events", 100, -600, 600);
  hDPtDiJet_Old->Sumw2();

  hDPhiMonoJet_Old = new TH1D("h1dPhiMonoJet_Old", "dPhi(genTop, MonoJets_Old); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiMonoJet_Old->Sumw2();
  hDPhiDiJet_Old = new TH1D("h1dPhiDiJet_Old", "dPhi(genTop, DiJets_Old); #Delta #Phi; Events", 100, -3.2, 3.2);
  hDPhiDiJet_Old->Sumw2();

  h1MT2Matched_New = new TH1D("h1MT2Matched_New", "MT2_Matched_New;M_{T2};Events", 100, 0, 1000);
  h1MT2Matched_New->Sumw2();
  h1MT2UnMatched_New = new TH1D("h1MT2UnMatched_New", "MT2_UnMatched_New;M_{T2};Events", 100, 0, 1000);
  h1MT2UnMatched_New->Sumw2();
  h1MT2Matched_Old = new TH1D("h1MT2Matched_Old", "MT2_Matched_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2Matched_Old->Sumw2();
  h1MT2UnMatched_Old = new TH1D("h1MT2UnMatched_Old", "MT2_UnMatched_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2UnMatched_Old->Sumw2();
  h1MT2DiJetNoB = new TH1D("h1MT2DiJetNoB", "MT2_DiJetNoB;M_{T2};Events", 100, 0, 1000);
  h1MT2DiJetNoB->Sumw2();
  h1MT2DiJetB = new TH1D("h1MT2DiJetB", "MT2_DiJetB;M_{T2};Events", 100, 0, 1000);
  h1MT2DiJetB->Sumw2();
  h1MT2MonoJetNoB = new TH1D("h1MT2MonoJetNoB", "MT2_MonoJetNoB;M_{T2};Events", 100, 0, 1000);
  h1MT2MonoJetNoB->Sumw2();
  h1MT2MonoJetB = new TH1D("h1MT2MonoJetB", "MT2_MonoJetB;M_{T2};Events", 100, 0, 1000);
  h1MT2MonoJetB->Sumw2();

  h1MT2MonoJet_Old = new TH1D("h1MT2MonoJet_Old", "MT2_MonoJet_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2MonoJet_Old->Sumw2();
  h1MT2DiJet_Old= new TH1D("h1MT2DiJet_Old", "MT2_DiJet_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2DiJet_Old->Sumw2();


  h1MT2LepDecay_New = new TH1D("h1MT2Lep_New", "MT2_Lep_New;M_{T2};Events", 100, 0, 1000);
  h1MT2LepDecay_New->Sumw2();
  h1MT2HadDecay_New =new TH1D("h1MT2Had_New", "MT2_Had_New;M_{T2};Events", 100, 0, 1000);
  h1MT2HadDecay_New->Sumw2();
  h1MT2IncDecay_New = new TH1D("h1MT2Inc_New", "MT2_Inc_New;M_{T2};Events", 100, 0, 1000);
  h1MT2IncDecay_New->Sumw2();


  h1MT2LepDecay_Old = new TH1D("h1MT2Lep_Old", "MT2_Lep_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2LepDecay_Old->Sumw2();
  h1MT2HadDecay_Old =new TH1D("h1MT2Had_Old", "MT2_Had_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2HadDecay_Old->Sumw2();
  h1MT2IncDecay_Old = new TH1D("h1MT2Inc_Old", "MT2_Inc_Old;M_{T2};Events", 100, 0, 1000);
  h1MT2IncDecay_Old->Sumw2();
}


TString HistSummary::Title(TString spec, unsigned int idx)
{
  TString title = spec;
  title += "_";
  title += idx;
  return title;
}




//Gives hadronically decaying top quark
std::vector< TLorentzVector > GetHadTopLVec(std::vector<TLorentzVector> genDecayLVec,
                                         std::vector<int>genDecayPdgIdVec,
                                         std::vector<int>genDecayIdxVec,
                                         std::vector<int>genDecayMomIdxVec)
{

  std::vector<TLorentzVector> tLVec;
  for(unsigned it=0; it<genDecayLVec.size(); it++){
    int pdgId = genDecayPdgIdVec.at(it);
    if(abs(pdgId)==6){
      for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ){
	  int pdgId = genDecayPdgIdVec.at(ig);
	  if(abs(pdgId)==24){
	    int flag = 0;
	    for(unsigned iq=0; iq<genDecayLVec.size(); iq++){
	      if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) {
		int pdgid = genDecayPdgIdVec.at(iq);
		if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++;
	      }
	    }
	    if(!flag) tLVec.push_back(genDecayLVec.at(it));
	  }
	}
      }//dau. loop
    }//top cond
  }//genloop
  return tLVec;
}

//Gives leptoniconically decaying top quark
std::vector< TLorentzVector > GetLepTopLVec(std::vector<TLorentzVector> genDecayLVec,
					    std::vector<int>genDecayPdgIdVec,
					    std::vector<int>genDecayIdxVec,
					    std::vector<int>genDecayMomIdxVec,
					    std::vector<TLorentzVector>& neutrinoLVec)
{

  std::vector<TLorentzVector> tLVec;
  for(unsigned it=0; it<genDecayLVec.size(); it++)
    {
      int pdgId = genDecayPdgIdVec.at(it);
      if(abs(pdgId)==6)
	{
	  for(unsigned ig=0; ig<genDecayLVec.size(); ig++)
	    {
	      if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) )
		{
		  int pdgId = genDecayPdgIdVec.at(ig);
		  if(abs(pdgId)==24)
		    {
		      int flag = 0;
		      for(unsigned iq=0; iq<genDecayLVec.size(); iq++)
			{
			  if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) 
			    {
			      int pdgid = genDecayPdgIdVec.at(iq);
			      if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++;
			      if(abs(pdgid)== 12 || abs(pdgid)== 14 || abs(pdgid)== 16)
				{
				  neutrinoLVec.push_back(genDecayLVec.at(iq));
				}
			    }
			}
		      if(flag) tLVec.push_back(genDecayLVec.at(it));
		    }
		}
	    }//dau. loop                                                                                                                                         
	}//top cond                                                                                                                                            
    }//genloop                                                                                                                                               
  return tLVec;
}

// Performs matching and updates Matched top using  New 
// tagger  and returns true if matched 
bool getMatchedTop(std::vector<TopObject*> NTop, 
		   std::vector<TopObject*>&MachedNTop, 
		   std::vector<TLorentzVector> topLVec, 
		   std::vector<TLorentzVector>& MtopLVec)
  {
    bool match = false;
    if(topLVec.size() ==  0) return match;
    double deltaR = 0.4;
    //Loop over  NTops to check mathcing
    for(unsigned int nT = 0; nT < NTop.size(); nT++)
      {
	double drMin = 99999.;
	unsigned int topId = -1;
	// Now loop over topLVec to match
	for (unsigned int genT = 0; genT < topLVec.size(); genT++)
	  {
	    const double dr =  NTop[nT]->p().DeltaR(topLVec.at(genT));
	    if( dr < drMin){drMin = dr; topId = genT;}
	  }
	if(drMin < deltaR)
	  {
	    MachedNTop.push_back(NTop[nT]);
	    MtopLVec.push_back(topLVec[topId]);
	    match = true;
	  }
      }
    return match;

  }
bool getMatchedTopOld(std::vector<TLorentzVector> NTop_Old,
                   std::vector<TLorentzVector>& MachedNTop_Old,
                   std::vector<TLorentzVector> topLVec,
                   std::vector<TLorentzVector>& MtopLVec_Old)
  {
    bool match = false;
    if(topLVec.size() ==  0) return match;
    double deltaR = 0.4;
    //Loop over  NTops to check mathcing                                                                                      
    for(unsigned int nT = 0; nT < NTop_Old.size(); nT++)
      {
        double drMin = 99999.;
        unsigned int topId = -1;
        // Now loop over topLVec to match                                                                                     
        for (unsigned int genT = 0; genT < topLVec.size(); genT++)
          {
            const double dr =  NTop_Old.at(nT).DeltaR(topLVec.at(genT));
            if( dr < drMin){drMin = dr; topId = genT;}
          }
        if(drMin < deltaR)
          {
            MachedNTop_Old.push_back(NTop_Old[nT]);
            MtopLVec_Old.push_back(topLVec[topId]);
            match = true;
          }
      }
    return match;


  }

double GetdRMin(std::vector<TLorentzVector> genTopLVec, TLorentzVector rSysLVec){
  double dRMin = 9999;
  
  for(unsigned dr =0; dr<genTopLVec.size(); dr++){
    double dR = genTopLVec[dr].DeltaR(rSysLVec);
    if(dR < dRMin) dRMin = dR;
  }
  return dRMin;
}

double getDPhiForMindR(std::vector<TLorentzVector> genTopLVec, TLorentzVector rSysLVec)
  {
    double dRMin = 9999;
    double dPhiForMindR = 0.0;
    for(unsigned dr =0; dr<genTopLVec.size(); dr++){
     
      double dR = genTopLVec[dr].DeltaR(rSysLVec);
      double dPhi = genTopLVec[dr].DeltaPhi(rSysLVec);
      if(dR < dRMin)
	{
	  dRMin = dR;
	  dPhiForMindR = dPhi;
	}
    } 
    return dPhiForMindR;
  }

double getDPtForMindR(std::vector<TLorentzVector> genTopLVec, TLorentzVector rSysLVec)
{
  double dRMin = 9999;
  double dPtForMindR = 0.0;
  for(unsigned dr =0; dr<genTopLVec.size(); dr++){

    double dR = genTopLVec[dr].DeltaR(rSysLVec);
    double dPt = genTopLVec[dr].Pt() - rSysLVec.Pt();
    if(dR < dRMin)
      {
	dRMin = dR;
	dPtForMindR = dPt;
      }
  }
  return dPtForMindR;
}


void GetAllTopLVec(std::vector<TLorentzVector> genDecayLVec,
		   std::vector<int>genDecayPdgIdVec,
		   std::vector<int>genDecayIdxVec,
		   std::vector<int>genDecayMomIdxVec,
		   std::vector<std::vector<TLorentzVector> >& topDausLVec, 
		   std::vector<std::vector<int> >& topDausIdxVec,
		   std::vector<TLorentzVector>& genTopLVec,
		   std::vector<int>& topDecayTypeVec, // 0: hadronic  1: leptonic 
		   int& cntHadDecay, 
		   int& cntLeptDecay,
		   bool badGenInfo = false
		   
)
{

  std::vector<int> cachedgenWIdxVec;
  for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++)
    {
      const int pdgId = genDecayPdgIdVec.at(ig);
      if( std::abs(pdgId) == 24 ) cachedgenWIdxVec.push_back(ig);
    }

  /************************************************************************************/
  /************************************************************************************/
  
  // 0: hadronic  1: leptonic
  //std::vector<int> topDecayTypeVec;
  //bool badGenInfo = false;
  for(unsigned int iw=0; iw<cachedgenWIdxVec.size(); iw++) //1st cachedgenWIdxVec
    {
      int idxW = cachedgenWIdxVec[iw];
      int momIdx_odx = genDecayMomIdxVec[idxW];
      int W_odx = genDecayIdxVec[idxW];
      bool isFromTop = false;
      std::vector<int> perbIdxVec, perWIdxVec;
      std::vector<TLorentzVector> perbLVec, perWLVec;
      TLorentzVector pergenTopLVec;
      int momIdx = -1;
      int decayType = 0;
      for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++) // 2nd genDecayPdgIdVec
	{
	  const int pdgId = genDecayPdgIdVec.at(ig);
	  if( genDecayIdxVec[ig] == momIdx_odx && std::abs(pdgId) == 6 ){ isFromTop = true; momIdx = ig; }
	  if( isFromTop )
	    {
	      if( genDecayMomIdxVec[ig] == momIdx_odx && std::abs(pdgId) == 5 )
		{
		  perbIdxVec.push_back(ig); 
		  perbLVec.push_back(genDecayLVec.at(ig));
		}
	      if( genDecayMomIdxVec[ig] == W_odx )
		{
		  perWIdxVec.push_back(ig);
		  perWLVec.push_back(genDecayLVec.at(ig));
		  if( std::abs(pdgId) >= 11 && std::abs(pdgId) <=16 ) decayType =1;

		}
	    } 
	}// 2nd genDecayPdgIdVec ends
      if( perbIdxVec.size() != 1 || perWIdxVec.size() != 2 || momIdx == -1 )
	{
	  badGenInfo = true;
	}
      std::vector<int> perIdxVec = perbIdxVec;
      std::vector<TLorentzVector> perLVec = perbLVec;
      for(unsigned int ip=0; ip<perWIdxVec.size(); ip++){
	perIdxVec.push_back(perWIdxVec.at(ip));
	perLVec.push_back(perWLVec.at(ip));
      }
      topDausLVec.push_back(perLVec);
      topDausIdxVec.push_back(perIdxVec);
      topDecayTypeVec.push_back(decayType);
      if( momIdx != -1 ) pergenTopLVec = genDecayLVec.at(momIdx); 
      genTopLVec.push_back(pergenTopLVec);
    } //1st cachedgenWIdxVec ends
  /************************************************************************************/
  /************************************************************************************/
  
  //int cntHadDecay = 0, cntLeptDecay = 0;
  for(unsigned int id=0; id<topDecayTypeVec.size(); id++)
    {
      if( topDecayTypeVec[id] == 0 ) cntHadDecay ++;
      else if( topDecayTypeVec[id] == 1 ) cntLeptDecay ++;
    }
}

int find_idx(int genIdx, const std::vector<int> &genDecayIdxVec){
  for(int ig=0; ig<(int)genDecayIdxVec.size(); ig++){
    if( genDecayIdxVec[ig] == genIdx ) return ig;
  }
  return -1;
}


bool find_mother(int momIdx, int dauIdx, const std::vector<int> &genDecayIdxVec, const std::vector<int> &genDecayMomIdxVec){
  if( momIdx == -1 || dauIdx == -1 ) return false;

  if( dauIdx == momIdx ) return true;

  int thisIdx = dauIdx;
  while( thisIdx >=0 ){
    int momGenIdx = genDecayMomIdxVec[thisIdx];
    thisIdx = find_idx(momGenIdx, genDecayIdxVec);
    if( thisIdx == momIdx ) return true;
  }

  return false;
}
