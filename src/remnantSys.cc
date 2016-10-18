#include "../include/remnantSys.h"

#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>

#include "SusyAnaTools/Tools/SATException.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

 
BaselineVessel * ttTest =0;
void mypassBaselineFunc(NTupleReader &tr)
{
  (*ttTest)(tr);
}

const std::string spec = "topTag";
const char * spec1 = "match";

// === Main Function ===================================================
int main(int argc, char* argv[]) 
  {
  try
    {
    if (argc < 5)
      {
	std::cerr <<"Please give 4 arguments "<<"SubsampleName"<< " MaxEvent"<<" No. of Files to run" << "Startfile" <<std::endl;
	std::cerr <<" Valid configurations are " << std::endl;
	std::cerr <<" ./remnantSys TTbarInc 1000 1 0" << std::endl;
	return -1;
      }

    const char *subSampleName = argv[1];
    const char *maxevent = argv[2];
    
    int numFiles = -1;
    int startFile = 0;
    // Change string arg to int
    int  maxEvent =  std::atoi(maxevent);
    numFiles =  std::atoi(argv[3]);
    startFile =  std::atoi(argv[4]);

    // Prepare file list and finalize it
    TChain *fChain = 0;
    HistSummary histSummary;
    int stfile = int(startFile/numFiles);
    histSummary.bookHist(subSampleName, stfile, spec1);
    const string condor =  (argc == 6) ? argv[5]: "";
  
    AnaSamples::SampleSet ss = condor.empty()? AnaSamples::SampleSet():AnaSamples::SampleSet(argv[5]);
    AnaSamples::SampleCollection sc(ss);
                                   
    double ScaleMC = 1.;                                                                              
    if(ss[subSampleName] != ss.null())                                                                             
      {                                                                                                               
	fChain = new TChain(ss[subSampleName].treePath.c_str());                                                           
	ss[subSampleName].addFilesToChain(fChain, startFile, numFiles);

	ScaleMC = ss[subSampleName].getWeight();
      }          

    AnaFunctions::prepareForNtupleReader();
    AnaFunctions::prepareTopTagger();
    type3Ptr->setdebug(true);

    NTupleReader *tr =0;
    tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);

    ttTest = new BaselineVessel(spec);
    tr->registerFunction(&mypassBaselineFunc);

    // Jump in to new topTagger method.
    TopTagger tt;
    tt.setCfgFile("Example_TopTagger.cfg");

    TString sample(subSampleName);
    // --- Analyse events --------------------------------------------
    std::cout<<"First loop begin: "<<std::endl;
    int entries = tr->getNEntries();
    std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<ScaleMC<<std::endl; 
    cout<<"maxevent: "<<maxEvent<<endl;

    Mt2::ChengHanBisect_Mt2_332_Calculator mt2Calculator;

    // Loop over the events (tree entries)
    while(tr->getNextEvent())
      {
	if(maxEvent>=0 && tr->getEvtNum() > maxEvent ) break;
	// Add print out of the progress of looping
	if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) 
	  {
	    std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
	  }

	const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
	const vector<int> &genDecayIdxVec = tr->getVec<int>("genDecayIdxVec");
	const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
	const vector<int> &genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
	const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
	const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");


	const int nbJets = tr->getVar<int>("cntCSVS" + spec);
	const int nJets = tr->getVar<int>("cntNJetsPt30Eta24" + spec);

	const double MT2 = tr->getVar<double>("best_had_brJet_MT2" + spec);

	const double& met    = tr->getVar<double>("met");
	const double& metphi = tr->getVar<double>("metphi");
            
	TLorentzVector metLVec;
	metLVec.SetPtEtaPhiM(met, 0, metphi, 0);




	if(nbJets < 1 || nJets < 4 || met < 200.) continue;



	double EvtWt = tr->getVar<double>("evtWeight");


	// New lorentzvector for old tagger
	//*******************************************
	TLorentzVector rSysConst_Old;
	//*******************************************

	//Get constituents of old remnant system
	std::vector<TLorentzVector> rSysConst_Old_Vec = type3Ptr->rSysConstituentLVecs;
	//Add constituents to make single LVec of remnant system
	for(unsigned int i =0; i < rSysConst_Old_Vec.size(); i++)
	  {
	    rSysConst_Old += rSysConst_Old_Vec.at(i);
	  }
	

	
	double EventWeight = EvtWt;
	double scaleMC  = ScaleMC * EventWeight;



	// Create gentop lorentz vector
	//*******************************************
	std::vector<TLorentzVector>  neutrinoLVec;
	vector<TLorentzVector> hadtopLVec = GetHadTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
	vector<TLorentzVector> leptopLVec = GetLepTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec, neutrinoLVec);

	//std:: cout << "neutrinoLVec size:  " << neutrinoLVec.size() << std::endl;
	
	std::vector<std::vector<TLorentzVector> > topDausLVec;
	std::vector<std::vector<int> > topDausIdxVec;
	std::vector<TLorentzVector> genTopLVec;
	std::vector<TLorentzVector> genTopLVec_notUsed;
	std::vector<int> topDecayTypeVec; // 0: hadronic  1: leptonic
	int cntHadDecay = 0;
	int cntLeptDecay = 0;
	bool badGenInfo = false;
	GetAllTopLVec( genDecayLVec,genDecayPdgIdVec, genDecayIdxVec,genDecayMomIdxVec,
		       topDausLVec, topDausIdxVec, genTopLVec_notUsed, topDecayTypeVec,
		       cntHadDecay,  cntLeptDecay, badGenInfo);

	vector<TLorentzVector> visleptopLVec;
	/*
	TLorentzVector NuLVec;
	// Nu comes from W > enu munu and taunu and tau> tauNu and tau > w> nu 
	if( !badGenInfo && genTopLVec.size()>=2 && cntHadDecay ==1 && cntLeptDecay ==1 )
	  {
	    for(unsigned int id=0; id<topDecayTypeVec.size(); id++)
	      {
		std::vector<int> perDauIdxVec = topDausIdxVec.at(id);
		std::vector<TLorentzVector> perDauLVec = topDausLVec.at(id);
		TLorentzVector pergenTopLVec = genTopLVec_notUsed.at(id);
		
		if( topDecayTypeVec[id] == 1 )
		  { // If leptonic decay.
		    for(unsigned int idau=0; idau< perDauIdxVec.size(); idau++)
		      {
			int pdgId = genDecayPdgIdVec.at(perDauIdxVec.at(idau));
			
			if( std::abs(pdgId) == 12 || std::abs(pdgId) == 14 || std::abs(pdgId) == 16)
			  {
			    NuLVec = NuLVec +  perDauLVec.at(idau);
			  }
		      }
		  }
	      } // for(unsigned int id=0; id<topDecayTypeVec.size()
	  } // if( !badGenInfo )  
	*/
	if(cntHadDecay != 1 && cntLeptDecay != 1) continue;
	// Create a visible system on leptonic decay by removing neutrinos 
	// from it      .
	
	//Remove neutrion from leptonic decay
	if(neutrinoLVec.size() == leptopLVec.size()){
	  for(unsigned int t = 0; t < leptopLVec.size(); t++)
	    {
	      visleptopLVec.push_back(leptopLVec.at(t) - neutrinoLVec.at(t));
	    }
	}
    
	// Add leptonic to inclusive top L Vec
	for(unsigned int i = 0; i < visleptopLVec.size(); i++){

	  genTopLVec.push_back(visleptopLVec.at(i));
	}	 
	// Add hadronic to genTop L vec
	for(unsigned int i = 0;i < hadtopLVec.size(); i++){

          genTopLVec.push_back(hadtopLVec.at(i));
	}


	if((hadtopLVec.size() <= 0) && (visleptopLVec.size() <= 0)){ continue;}
	//*******************************************

	//construct vector of constituents 
	vector<Constituent> constituents = ttUtility::packageConstituents(jetsLVec, recoJetsBtag_0);
	//run tagger
	tt.runTagger(constituents);
	//get output of tagger
	const TopTaggerResults& ttr = tt.getResults();

	//Use result for top var
	vector<TopObject*> NTop = ttr.getTops();
	const int nTops_Old = tr->getVar<int>("nTopCandSortedCnt" + spec);


	const vector<TLorentzVector> NTop_Old = tr->getVec<TLorentzVector>("vTops"+spec);


	//vector<TopObject> NtopCand = ttr.getTopCandidates();    
	const TopObject rSysNew = ttr.getRsys();

	//*******************************************	
	TLorentzVector rSysConst_New = rSysNew.p();

	// Check the size r sys(No of constituents in two taggers                                                                                             
        const std::vector<Constituent const *> constituents_new  = rSysNew.getConstituents();
        const int cnstSize_New =  rSysNew.getNConstituents();


	//Calculate New MT2
	double newMT2 = -999.9;
	double deltaR = -999;
	if(ttr.getTops().size() >= 1)
	  {
	    const double massOfSystemA = ttr.getTops()[0]->P().M(); // GeV
	    const double pxOfSystemA   = ttr.getTops()[0]->P().Px(); // GeV
	    const double pyOfSystemA   = ttr.getTops()[0]->P().Py(); // GeV
            
	    double massOfSystemB =  rSysConst_New.M(); // GeV
	    double pxOfSystemB   =  rSysConst_New.Px(); // GeV
	    double pyOfSystemB   =  rSysConst_New.Py(); // GeV

	    deltaR = ttr.getTops()[0]->P().DeltaR(rSysConst_New);

	    if(ttr.getTops().size() >= 2)
	      {
		massOfSystemB = ttr.getTops()[1]->P().M(); // GeV 
		pxOfSystemB   = ttr.getTops()[1]->P().Px(); // GeV
		pyOfSystemB   = ttr.getTops()[1]->P().Py(); // GeV

		deltaR = ttr.getTops()[0]->P().DeltaR(ttr.getTops()[1]->P());
	      }
            
	    // The missing transverse momentum:
	    const double pxMiss        = metLVec.Px(); // GeV
	    const double pyMiss        = metLVec.Py(); // GeV
            
	    const double invis_mass    = metLVec.M(); // GeV
          
	    Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, pyOfSystemA), massOfSystemA);
	    Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB), massOfSystemB);
	    Mt2::TwoVector                pT_Miss(pxMiss, pyMiss);
          
	    newMT2 = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, invis_mass);
	  }


	if(NTop.size() >= 1)
	  
	  {
	    
	    double dRMinGenNew = 0; double dPhiMinGenNew = 0; double dPtMinGenNew = 0;
	    dRMinGenNew = GetdRMin(genTopLVec, rSysConst_New);
	    dPhiMinGenNew = getDPhiForMindR(genTopLVec, rSysConst_New);
	    dPtMinGenNew = getDPtForMindR(genTopLVec, rSysConst_New);
	    
	    
	    histSummary.hRSysNewTag_dR->Fill(dRMinGenNew,scaleMC);
	    histSummary.hRSysNewTag_dPt->Fill(dPtMinGenNew,scaleMC);
	    histSummary.hRSysNewTag_dPhi->Fill(dPhiMinGenNew,scaleMC);
	    histSummary.h1MT2IncDecay_New->Fill(newMT2, scaleMC);

	  }
	if (nTops_Old >=1)
	  {  
	    double dRMinGenOld = 0; double dPhiMinGenOld = 0; double dPtMinGenOld = 0;
	    dRMinGenOld = GetdRMin(genTopLVec, rSysConst_Old);
	    dPhiMinGenOld = getDPhiForMindR(genTopLVec, rSysConst_Old);
	    dPtMinGenOld = getDPtForMindR(genTopLVec, rSysConst_Old);
	
	    histSummary.hRSysOldTag_dR->Fill(dRMinGenOld,scaleMC);
	    histSummary.hRSysOldTag_dPt->Fill(dPtMinGenOld,scaleMC);
	    histSummary.hRSysOldTag_dPhi->Fill(dPhiMinGenOld,scaleMC);
	    histSummary.h1MT2IncDecay_Old->Fill(MT2, scaleMC);
	  }


	//New top tagged

	if(NTop.size() >= 1)
	  {

	    histSummary.rSysConst_New_Pt->Fill(rSysConst_New.Pt(), scaleMC);
	    
	    for(int i = 0; i < topDecayTypeVec.size(); i++)
	      {
		if(topDecayTypeVec[i] ==1) //Leptonic
		  {
		    double drLepNew = 0; double dPhiLepNew = 0; double dPtLepNew = 0;
		    
		    drLepNew = GetdRMin(visleptopLVec, rSysConst_New);                                                                
		    dPhiLepNew = getDPhiForMindR(visleptopLVec, rSysConst_New);
		    dPtLepNew = getDPtForMindR(visleptopLVec, rSysConst_New);
		    
		    histSummary.hRSysNewTagLep_dR->Fill(drLepNew, scaleMC );
		    histSummary.hRSysNewTagLep_dPt->Fill(dPtLepNew, scaleMC);
		    histSummary.hRSysNewTagLep_dPhi->Fill(dPhiLepNew,scaleMC);
		    histSummary.h1MT2LepDecay_New->Fill(newMT2, scaleMC);
		  }
		
		if(topDecayTypeVec[i] == 0) // Hadronic
		  {
		    double drHadNew = 0; double dPhiHadNew = 0; double dPtHadNew = 0;
		    drHadNew = GetdRMin(hadtopLVec, rSysConst_New);
		    dPhiHadNew = getDPhiForMindR(hadtopLVec, rSysConst_New);
		    dPtHadNew = getDPtForMindR(hadtopLVec, rSysConst_New);
		    
		    histSummary.hRSysNewTagHad_dR->Fill(drHadNew, scaleMC );
		    histSummary.hRSysNewTagHad_dPt->Fill(dPtHadNew, scaleMC);
		    histSummary.hRSysNewTagHad_dPhi->Fill(dPhiHadNew,scaleMC);
		    
		    histSummary.h1MT2HadDecay_New->Fill(newMT2, scaleMC);
		  }
	      }
	    //Get MAtched vectors
	    bool topmatch = false;
	    vector<TopObject*> MatchNtop;
	    vector<TLorentzVector> MatchGentop; 
	    
	    //Matched UnMatched Hadronically decaying top case.
	    if(getMatchedTop(NTop, MatchNtop, hadtopLVec, MatchGentop)) topmatch = true;//final top match	
	    
	    if(topmatch)
	      {
		double dRMinMatchNew = GetdRMin(visleptopLVec, rSysConst_New);
		double dPtMatchNew = getDPtForMindR(visleptopLVec, rSysConst_New);
		double dPhiMatchNew = getDPhiForMindR(visleptopLVec, rSysConst_New);
		histSummary.h1MatchedDrMin_New->Fill(dRMinMatchNew, scaleMC);
		histSummary.h1MatchedDPt_New->Fill(dPtMatchNew, scaleMC);
		histSummary.h1MatchedDPhi_New->Fill(dPhiMatchNew, scaleMC);
		histSummary.h1MT2Matched_New->Fill(newMT2, scaleMC);
	      }
	    
	    if(!topmatch)
	      {
		
		double dRMinUnMatchNew = GetdRMin(visleptopLVec, rSysConst_New);
		double dPtUnMatchNew = getDPtForMindR(visleptopLVec, rSysConst_New);
		double dPhiUnMatchNew = getDPhiForMindR(visleptopLVec, rSysConst_New);
		histSummary.h1UnMatchedDrMin_New->Fill(dRMinUnMatchNew, scaleMC);
		histSummary.h1UnMatchedDPt_New->Fill(dPtUnMatchNew, scaleMC);
		histSummary.h1UnMatchedDPhi_New->Fill(dPhiUnMatchNew, scaleMC);
		histSummary.h1MT2UnMatched_New->Fill(newMT2,scaleMC);
	      }
	    
	    
	    
	    // Catagorical 1b Nob dijet MonoJet case
	    bool monoJetWithb = false;
	    bool monoJetWithOutb = false;
	    bool diJetWithb = false;
	    bool diJetWithOutb =  false;
	    bool hasbJets_New = false;
	    
	    for(unsigned int ncst = 0; ncst < constituents_new.size(); ncst++)
	      {      // Check if there is b jet in remanant system 
		if (constituents_new.at(ncst)->getBTagDisc() > 0.8)  hasbJets_New = true;
	      }
		// If  remnant system has b tagged jets
	    if( hasbJets_New)
	      {
		if(constituents_new.size() ==1) monoJetWithb  = true;
		if(constituents_new.size() ==2) diJetWithb  = true;
	      }
	    // If remnant does't have a b jets
	    else
	      {
		if(constituents_new.size() ==1) monoJetWithOutb = true; // MonoJet
		if(constituents_new.size() ==2)   diJetWithOutb = true;//DiJet
	      }
	  
	    
	    // New Tagger (dijet, MonoJet) * (0 b, 1 b)
	    double dRMin_New = GetdRMin(visleptopLVec, rSysConst_New);
	    double dPhiForMindR_New = getDPhiForMindR(visleptopLVec, rSysConst_New);
	    double dPt_New = getDPtForMindR(visleptopLVec, rSysConst_New);	
	    
	    if(monoJetWithb){ histSummary.hDrMinMonoJetb->Fill(dRMin_New, scaleMC); histSummary.hDPhiMonoJetb->Fill(dPhiForMindR_New, scaleMC);
	      histSummary.hDPtMonoJetb->Fill(dPt_New, scaleMC); histSummary.h1MT2MonoJetB->Fill(newMT2, scaleMC);}
	    if(diJetWithb) {histSummary.hDrMinDiJetb->Fill(dRMin_New, scaleMC); histSummary.hDPhiDiJetb->Fill(dPhiForMindR_New, scaleMC);
	      histSummary.hDPtDiJetb->Fill(dPt_New, scaleMC);  histSummary.h1MT2DiJetB->Fill(newMT2, scaleMC);}
	    if(monoJetWithOutb){ histSummary.hDrMinMonoJetNob->Fill(dRMin_New, scaleMC); histSummary.hDPhiMonoJetNob->Fill(dPhiForMindR_New, scaleMC);
	      histSummary.hDPtMonoJetNob->Fill(dPt_New, scaleMC); histSummary.h1MT2MonoJetNoB->Fill(newMT2, scaleMC);}
	    if(diJetWithOutb) {histSummary.hDrMinDiJetNob->Fill(dRMin_New, scaleMC); histSummary.hDPhiDiJetNob->Fill(dPhiForMindR_New, scaleMC);
	      histSummary.hDPtDiJetNob->Fill(dPt_New, scaleMC); histSummary.h1MT2DiJetNoB->Fill(newMT2, scaleMC);}
	    
	  }
	
	//If Old Top Tagged
	if(nTops_Old >=1)
	  {
	
	    histSummary.rSysConst_Old_Pt->Fill(rSysConst_Old.Pt(), scaleMC);

	    for(int i = 0; i < topDecayTypeVec.size(); i++)
	      {
		if(topDecayTypeVec[i] == 1)
		  {
		    double drLepOld = 0; double dPhiLepOld = 0; double dPtLepOld = 0;
		    drLepOld = GetdRMin(visleptopLVec, rSysConst_Old);
		    dPhiLepOld = getDPhiForMindR(visleptopLVec, rSysConst_Old);
		    dPtLepOld = getDPtForMindR(visleptopLVec, rSysConst_Old);
		    
		    histSummary.hRSysOldTagLep_dR->Fill(drLepOld, scaleMC );
		    histSummary.hRSysOldTagLep_dPt->Fill(dPtLepOld, scaleMC);
		    histSummary.hRSysOldTagLep_dPhi->Fill(dPhiLepOld,scaleMC);
		    histSummary.h1MT2LepDecay_Old->Fill(MT2, scaleMC);
		  }

		if(topDecayTypeVec[i] == 0)
		  {
		    double drHadOld = 0; double dPhiHadOld = 0; double dPtHadOld = 0;
		    drHadOld = GetdRMin(hadtopLVec, rSysConst_Old);
		    dPhiHadOld = getDPhiForMindR(hadtopLVec, rSysConst_Old);
		    dPtHadOld = getDPtForMindR(hadtopLVec, rSysConst_Old);
		    
		    
		    histSummary.hRSysOldTagHad_dR->Fill(drHadOld, scaleMC ); 
		    histSummary.hRSysOldTagHad_dPt->Fill(dPtHadOld, scaleMC);
		    histSummary.hRSysOldTagHad_dPhi->Fill(dPhiHadOld,scaleMC);
		    histSummary.h1MT2HadDecay_Old->Fill(MT2, scaleMC);
		  }
	      }
	    
	    
	    //Get MAtched vectors                                                                                                                            
	    bool topmatch_Old = false;
	    vector<TLorentzVector> MatchNtop_Old;
	    vector<TLorentzVector> MatchGentop_Old;

	    //Matched UnMatched Hadronically decaying top case.                                                                                
 	    if(getMatchedTopOld(NTop_Old, MatchNtop_Old, hadtopLVec, MatchGentop_Old)) topmatch_Old = true;//final top match old                             

	
	    if(topmatch_Old)
	      {
		double dRMinMatchOld = GetdRMin(visleptopLVec, rSysConst_Old);
		double dPtMatchOld = getDPtForMindR(visleptopLVec, rSysConst_Old);
		double dPhiMatchOld = getDPhiForMindR(visleptopLVec, rSysConst_Old);
		histSummary.h1MatchedDrMin_Old->Fill(dRMinMatchOld, scaleMC);
		histSummary.h1MatchedDPt_Old->Fill(dPtMatchOld, scaleMC);
		histSummary.h1MatchedDPhi_Old->Fill(dPhiMatchOld, scaleMC);
		histSummary.h1MT2Matched_Old->Fill(MT2,scaleMC);
	      }

	    if(!topmatch_Old)
	      {
		double dRMinUnMatchOld = GetdRMin(visleptopLVec, rSysConst_Old);
		double dPtUnMatchOld = getDPtForMindR(visleptopLVec, rSysConst_Old);
		double dPhiUnMatchOld = getDPhiForMindR(visleptopLVec, rSysConst_Old);
		histSummary.h1UnMatchedDrMin_Old->Fill(dRMinUnMatchOld, scaleMC);
		histSummary.h1UnMatchedDPt_Old->Fill(dPtUnMatchOld, scaleMC);
		histSummary.h1UnMatchedDPhi_Old->Fill(dPhiUnMatchOld, scaleMC);
		histSummary.h1MT2UnMatched_Old->Fill(MT2,scaleMC);
	      }


	    // Catagorical 1b Nob dijet MonoJet case                                                                                                              
	   
	    bool monoJet_Old = false;
	    bool diJet_Old = false;
	   
	    if(rSysConst_Old_Vec.size() == 1) monoJet_Old = true;
	    if(rSysConst_Old_Vec.size() == 2) diJet_Old = true;
		

	      

	    // Old Tagger (dijet, MonoJet)                                                                                                                        
	    double dRMin_Old = GetdRMin(visleptopLVec, rSysConst_Old);
	    double dPhiForMindR_Old= getDPhiForMindR(visleptopLVec, rSysConst_Old);
	    double dPt_Old = getDPtForMindR(visleptopLVec, rSysConst_Old);

	    if(monoJet_Old){histSummary.hDrMinMonoJet_Old->Fill(dRMin_Old, scaleMC); histSummary.hDPhiMonoJet_Old->Fill(dPhiForMindR_Old, scaleMC);
	      histSummary.hDPtMonoJet_Old->Fill(dPt_Old, scaleMC); histSummary.h1MT2MonoJet_Old->Fill(MT2, scaleMC);}
	    if(diJet_Old) {histSummary.hDrMinDiJet_Old->Fill(dRMin_Old, scaleMC); histSummary.hDPhiDiJet_Old->Fill(dPhiForMindR_Old, scaleMC);
	      histSummary.hDPtDiJet_Old->Fill(dPt_Old, scaleMC); histSummary.h1MT2DiJet_Old->Fill(MT2, scaleMC);}


	  }
	//*******************************************
	/*
	cout << "**************************" << endl;
	cout << "drLepNew  " << drLepNew << "  dPhiLepNew  " << dPhiLepNew << "  dPtLepNew  " << dPtLepNew  << endl;  
	cout <<"drLepOld  " <<drLepOld << "  dPhiLepOld  " <<dPhiLepOld << "  dPtLepOld  " << dPtLepOld  << endl;
	cout << endl;
	cout << "drHadNew  " << drHadNew << "  dPhiHadNew  " << dPhiHadNew << "  dPtHadNew  " << dPtHadNew  << endl;
        cout <<"drHadOld  " <<drHadOld << "  dPhiHadOld  " <<dPhiHadOld << "  dPtHadOld  " << dPtHadOld  << endl;
	*/
       



      }



    histSummary.h1MatchedDrMin_New->Write();
    histSummary.h1MatchedDrMin_Old->Write();
    histSummary.h1MatchedDPt_New->Write();
    histSummary.h1MatchedDPt_Old->Write();
    histSummary.h1MatchedDPhi_New->Write();
    histSummary.h1MatchedDPhi_Old->Write();

    histSummary.h1UnMatchedDrMin_New->Write();
    histSummary.h1UnMatchedDrMin_Old->Write();
    histSummary.h1UnMatchedDPt_New->Write();
    histSummary.h1UnMatchedDPt_Old->Write();
    histSummary.h1UnMatchedDPhi_New->Write();
    histSummary.h1UnMatchedDPhi_Old->Write();


    histSummary.hRSysOldTag_dR->Write();
    histSummary.hRSysOldTag_dPt->Write();
    histSummary.hRSysOldTag_dPhi->Write();
    histSummary.hRSysNewTag_dR->Write();
    histSummary.hRSysNewTag_dPt->Write();
    histSummary.hRSysNewTag_dPhi->Write();
    
    histSummary.hRSysOldTagLep_dR->Write();
    histSummary.hRSysOldTagLep_dPt->Write();
    histSummary.hRSysOldTagLep_dPhi->Write();
    histSummary.hRSysNewTagLep_dR->Write();
    histSummary.hRSysNewTagLep_dPt->Write();
    histSummary.hRSysNewTagLep_dPhi->Write();

    histSummary.hRSysOldTagHad_dR->Write();
    histSummary.hRSysOldTagHad_dPt->Write();
    histSummary.hRSysOldTagHad_dPhi->Write();
    histSummary.hRSysNewTagHad_dR->Write();
    histSummary.hRSysNewTagHad_dPt->Write();
    histSummary.hRSysNewTagHad_dPhi->Write();
    
    histSummary.hDrMinMonoJetNob->Write();
    histSummary.hDrMinDiJetNob->Write();
    histSummary.hDrMinMonoJetb->Write();
    histSummary.hDrMinDiJetb->Write();

    histSummary.hDPtMonoJetNob->Write();
    histSummary.hDPtDiJetNob->Write();

    histSummary.hDPtMonoJetb->Write();
    histSummary.hDPtDiJetb->Write();

    histSummary.hDPhiMonoJetNob->Write();
    histSummary.hDPhiDiJetNob->Write();

    histSummary.hDPhiMonoJetb->Write();
    histSummary.hDPhiDiJetb->Write();


    histSummary.hDrMinMonoJet_Old->Write();
    histSummary.hDrMinDiJet_Old->Write();

    histSummary.hDPtMonoJet_Old->Write();
    histSummary.hDPtDiJet_Old->Write();

    histSummary.hDPhiMonoJet_Old->Write();
    histSummary.hDPhiDiJet_Old->Write();
   



    histSummary.h1MT2Matched_New->Write();
    histSummary.h1MT2UnMatched_New->Write();
    histSummary.h1MT2Matched_Old->Write();
    histSummary.h1MT2UnMatched_Old->Write();
    histSummary.h1MT2DiJetNoB->Write();
    histSummary.h1MT2DiJetB->Write();
    histSummary.h1MT2MonoJetNoB->Write();
    histSummary.h1MT2MonoJetB->Write();
    histSummary.h1MT2DiJet_Old->Write();
    histSummary.h1MT2MonoJet_Old->Write();

    histSummary.rSysConst_New_Pt->Write();
    histSummary.rSysConst_Old_Pt->Write();

    
    histSummary.h1MT2IncDecay_New->Write();
    histSummary.h1MT2IncDecay_Old->Write();

    histSummary.h1MT2LepDecay_New->Write();
    histSummary.h1MT2LepDecay_Old->Write();

    histSummary.h1MT2HadDecay_New->Write();
    histSummary.h1MT2HadDecay_Old->Write();


    histSummary.outFileName->Close();
    }
    catch(const SATException& e)
      {
	std::cout << e << std::endl;
      }
    catch(const TTException& e)
      {
	std::cout << e << std::endl;
      }
    return 0;
    
  }
