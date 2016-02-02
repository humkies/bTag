#include "BTagCorrector.h"
#include "NTupleReader.h"

#include "samples.h"
#include "customize.h"
#include "TH2Poly.h"
#include "searchBins.h"
#include "baselineDef.h"

#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

using namespace std;

BaselineVessel * SRblv =0;
void mypassBaselineFunc(NTupleReader &tr){
  (*SRblv)(tr);
}



const std::string spec = "MY";
int main(){
/*
int main(int argc, char* argv[])

{
  std::string condor = "";
  if(argc >= 5) condor = argv[4];
  AnaSamples::SampleSet        ss(argv[4], 3.0);
  AnaSamples::SampleCollection sc(ss);

  std::string ssName;
  int numFiles;
  int startFile;

  if(argc >= 4)
    {
      ssName = argv[1];
      numFiles = atoi(argv[2]);
      startFile = atoi(argv[3]);
    }

  const std::string samples = argv[1];

  TChain* ch = 0;
  if(ss[ssName] != ss.null())
    {
      ch = new TChain(ss[ssName].treePath.c_str());

      ss[ssName].addFilesToChain(ch, startFile, numFiles);
    }
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  type3Ptr->setdebug(true);
  print_searchBins();
  SRblv = new BaselineVessel(spec);
  NTupleReader tr(ch);
  tr.registerFunction(&mypassBaselineFunc);
*/
    const std::string samples = "SMS-T2tt_mStop-600-950_mLSP-1to450_baseline";
  const std::string nmStr = "_"+samples;
  const TString nmStrT = nmStr;
  //  size_t t0 = clock();  
  
  //char nBase[] = "root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Oct_2015_Ntp_v2X/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_Spring15DR74_Asympt25ns_Ntp_v2/150927_143844/0000/stopFlatNtuples_%d.root";
//  char nBase[] = "root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Oct_2015_Ntp_v2X/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_Spring15DR74_Asympt25ns_Ntp_v2/150927_143844/0000/stopFlatNtuples_%d.root";


  char nBase[]= "root://cmseos.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Dec_2015_Ntp_v4X//SMS-T2tt_mStop-600-950_mLSP-1to450_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Spring15_74X_Dec_2015_Ntp_v4p0_SMS-T2tt_FastSim_mStop-600-950_mLSP-1to450/160125_070532/0000/stopFlatNtuples_%d.root";



  TChain *f = new TChain("stopTreeMaker/AUX");                                          
                                                                                              
        size_t t0 = clock();                                                              
                                                                                              
      char fname[512];                                                                      
      for(int i = 1; i <= 10; ++i)                                                           
	{                                                                                  
	  sprintf(fname, nBase, i);                                                          
	  f->Add(fname);                                                                     
	}                                                                                  
                                                                                              
      print_searchBins();                                                                
                                                                                              
      NTupleReader tr(f); 

      AnaFunctions::prepareForNtupleReader();                                                                                                                                            
      AnaFunctions::prepareTopTagger();                                                                                                                                                                           
      type3Ptr->setdebug(true);   
      SRblv = new BaselineVessel(spec);
      tr.registerFunction(&mypassBaselineFunc);
     
      BTagCorrector btagcorr;
      tr.registerFunction(btagcorr);

      std::cout << "NEVTS: " << tr.getNEntries() << std::endl;
  /********************************************************************/
  // BTAG Corrector
  /********************************************************************/


  
  //  TFile* effFile = new TFile("bTagEfficiency"+nmStrT+".root", "READ");
  //  BTagCorrector btagcorr;
  //  btagcorr.SetEffs(effFile);
  //btagcorr.SetCalib("CSVv2_mod.csv");



  TFile* SFFile = new TFile("bTagSF"+nmStrT+".root", "RECREATE");   
  TH1D * h1_BTagSF_Cent = new TH1D("h1_BTagSF_Cent"+nmStrT, nmStrT+": BTags", 4, 0, 4);
  TH1D * h1_BTagSF_Up = new TH1D("h1_BTagSF_Up"+nmStrT, nmStrT+": BTags", 4, 0, 4);
  TH1D * h1_BTagSF_Down = new TH1D("h1_BTagSF_Down"+nmStrT, nmStrT+": BTags", 4, 0, 4);
  TH1D * h1_BTagSF_nBJets = new TH1D("h1_BTagSF_nBJets"+nmStrT, nmStrT+": BTags", 4, 0, 4);


  while(tr.getNextEvent())
    {
      if(tr.getEvtNum() == 1)
        {
	  tr.printTupleMembers();
	  FILE * fout = fopen("NTupleTypes.txt", "w");
	  tr.printTupleMembers(fout);
	  fclose(fout);
        }

      const double scaleFactor_up = tr.getVar<double>("bTagSF_EventWeightSimple_Up");
      const double scaleFactor_central = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
      const double scaleFactor_down = tr.getVar<double>("bTagSF_EventWeightSimple_Down");
      const vector<double> prob_up = tr.getVec<double>("bTagSF_EventWeightProb_Up");
      const vector<double> prob_cent = tr.getVec<double>("bTagSF_EventWeightProb_Central");
      const vector<double> prob_down = tr.getVec<double>("bTagSF_EventWeightProb_Down");
      //      const int nbJets = tr.getVar<int>("cntCSVS" +spec);
     const  int nBJets = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
     const bool passBaseline = tr.getVar<bool>("passBaseline" + spec);
      
      double Sumup = prob_up[0]+ prob_up[1]+ prob_up[2]+ prob_up[3];
      double Sumcent = prob_cent[0]+ prob_cent[1]+ prob_cent[2]+ prob_cent[3];
      double Sumdown = prob_down[0]+ prob_down[1]+ prob_down[2]+ prob_down[3];
      


      if(!passBaseline) continue;
      h1_BTagSF_Cent->Fill(0., prob_cent[0]);
      h1_BTagSF_Cent->Fill(1., prob_cent[1]);
      h1_BTagSF_Cent->Fill(2., prob_cent[2]);
      h1_BTagSF_Cent->Fill(3., prob_cent[3]);
      
      
      h1_BTagSF_Up->Fill(0., prob_up[0]);
      h1_BTagSF_Up->Fill(1., prob_up[1]);
      h1_BTagSF_Up->Fill(2., prob_up[2]);
      h1_BTagSF_Up->Fill(3., prob_up[3]);
      
      h1_BTagSF_Down->Fill(0., prob_down[0]);
      h1_BTagSF_Down->Fill(1., prob_down[1]);
      h1_BTagSF_Down->Fill(2., prob_down[2]);
      h1_BTagSF_Down->Fill(3., prob_down[3]);

      h1_BTagSF_nBJets->Fill(nBJets);


      









      /*             
       if(tr.getEvtNum()%1000 == 0)
	 {
	   // std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
	   cout << "ScaleFactor          :  " << scaleFactor << endl;
	   }*/
      if(tr.getEvtNum()%1000 == 0)
      {
      //  std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;     
	//	  printf("%6.5f\t\t %6.5f\t %6.5f\t %6.5f\n", prob[0], prob[1], prob[2], prob[3]);
	  cout << prob_up[0]<< ",    "<< prob_up[1]<<",    "<< prob_up[2]<<",     "<< prob_up[3]<< endl;
	  cout << prob_cent[0]<< ",    "<< prob_cent[1]<<",    "<< prob_cent[2]<<",     "<< prob_cent[3]<< endl;
	  cout << prob_down[0]<< ",    "<< prob_down[1]<<",    "<< prob_down[2]<<",     "<< prob_down[3]<< endl;
	  cout << "Sum:    " << Sumup << "    "<< Sumcent << "    " <<Sumdown  << endl;

	  cout <<"NBjets == " << nBJets << endl;
	 }/*
      if(tr.getEvtNum()%1000 == 0)
	{
	  //  std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;                                                                                                                                          
	  printf("%6.5f\t\t %6.5f\t %6.5f\t %6.5f\n", prob[0], prob[1], prob[2], prob[3]); 
	}
      */
       
    }




      h1_BTagSF_Cent->Write();
      h1_BTagSF_Up->Write();
      h1_BTagSF_Down->Write();
      h1_BTagSF_nBJets->Write();
      SFFile->Close();

      return 0;


}
