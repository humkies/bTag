//Tools Headers
#include "NTupleReader.h"
#include "samples.h"
#include "customize.h"
#include "searchBins.h"
#include "baselineDef.h"


//ROOT headers
#include "TFile.h"
#include "TStopwatch.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TChain.h"
//STL Headers
#include<iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;

const int nPtBins = 17;
const int nEtaBins = 3;

const double ptBins[]    =  {20,30,40,50,60,70,80,100,120,160,210,260,320,400,500,600,800,99999};
const double etaBins[]  =  {0.0,0.8,1.6,2.4};

BaselineVessel * SRblv =0;
void mypassBaselineFunc(NTupleReader &tr){
  (*SRblv)(tr);
}



const std::string spec = "MY";


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
      
      /*
  char nBase[] = "root://cmseos.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Dec_2015_Ntp_v4X//SMS-T2tt_mStop-600-950_mLSP-1to450_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Spring15_74X_Dec_2015_Ntp_v4p0_SMS-T2tt_FastSim_mStop-600-950_mLSP-1to450/160125_070532/0000/stopFlatNtuples_%d.root";
	TChain *f = new TChain("stopTreeMaker/AUX");

	//  size_t t0 = clock();

	char fname[512];
	for(int i = 1; i <= 10; ++i)
	   {
	   sprintf(fname, nBase, i);
	   f->Add(fname);
	   }

	   print_searchBins();
    
	   NTupleReader tr(f);
      
*/
      std::cout << "NEVTS: " << tr.getNEntries() << std::endl;
      
      //      const std::string samples = "SMS-T2tt_FastSim_mStop-600-950_mLSP-1to450";
      const std::string nmStr = "_"+samples;
      const TString nmStrT = nmStr;
      
      TFile *infile = new TFile("bTagEfficiency"+nmStrT+".root", "RECREATE");
      /*************************************************************/      
      //Declare efficiency histograms. n-> numerator d-> denominator
      // b-> bquark, c-> charm, usdg-> other light quarks.
      /*************************************************************/
      TH1::AddDirectory(kFALSE);
      TH2* n_eff_b =    new TH2D("h_eff_b", "bTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      TH2* n_eff_c =    new TH2D("h_eff_c", "cTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      TH2* n_eff_udsg = new TH2D("h_eff_udsg", "udsgTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      TH2* d_eff_b =    new TH2D("d_eff_b", "bTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      TH2* d_eff_c =    new TH2D("d_eff_c", "cTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      TH2* d_eff_udsg = new TH2D("d_eff_udsg", "udsgTag_Efficiency"+nmStrT, nPtBins, ptBins, nEtaBins, etaBins);
      
      n_eff_b->GetXaxis()->SetTitle("p_{T} [GeV]");
      n_eff_b->GetYaxis()->SetTitle("#eta");

      n_eff_c->GetXaxis()->SetTitle("p_{T} [GeV]");
      n_eff_c->GetYaxis()->SetTitle("#eta");
      
      n_eff_udsg->GetXaxis()->SetTitle("p_{T} [GeV]");
      n_eff_udsg->GetYaxis()->SetTitle("#eta");
      
      /*************************************************************/
      // Event loop begins                                                                                                            
      /*************************************************************/
      while(tr.getNextEvent())
	{

	  if(tr.getEvtNum() == 1)
	    {
	      tr.printTupleMembers();
	      FILE * fout = fopen("NTupleTypes.txt", "w");
	      tr.printTupleMembers(fout);
	      fclose(fout);
	    }

	  if( tr.getEvtNum()%100 == 0 ) cout<<"\n   Processing the "<<tr.getEvtNum()-1<<"th event ..."<<endl;
	  
	  
	  const  vector<TLorentzVector> inputJets = tr.getVec<TLorentzVector>("jetsLVec");
	  const vector<double> recoJetsBtag = tr.getVec<double>("recoJetsBtag_0");
	  const vector<int> recoJetsFlavor = tr.getVec<int>("recoJetsFlavor");
	  cout << "inputJets size : " << inputJets.size() << endl;
	  cout << "recoJetsBtag size : " << recoJetsBtag.size() << endl;
	  cout << "recoJetsFlavor size : " << recoJetsFlavor.size() << endl;
	  
	  
     
	  for(unsigned int ij=0; ij<inputJets.size(); ij++)
	    {
	      
	      double pt = inputJets[ij].Pt();
	      double eta = fabs(inputJets[ij].Eta());
	      if(pt< 30.0 || eta >2.4) continue; //HT Jets cut
	      double csv = recoJetsBtag.at(ij);
	      int flav =  abs(recoJetsFlavor.at(ij));
	      
	      if(flav==5) //b Jets
		{
		  d_eff_b->Fill(pt,eta);
		  if(csv > 0.890) n_eff_b->Fill(pt,eta);
		}
	      else if(flav==4) // c jets
		{
		  d_eff_c->Fill(pt,eta);
		  if(csv > 0.890) n_eff_c->Fill(pt,eta);
		}
	      else if(flav<4 || flav==21) // other flavours
		{
		  d_eff_udsg->Fill(pt,eta);
		  if(csv > 0.890) n_eff_udsg->Fill(pt,eta);
		}
	      
	    }
	}
      /*************************************************************/
      // End of Event loop
      /*************************************************************/
      
      n_eff_b->Divide(d_eff_b);
      n_eff_c->Divide(d_eff_c);
      n_eff_udsg->Divide(d_eff_udsg);
      n_eff_b->Write();
      n_eff_c->Write();
      n_eff_udsg->Write();  
      infile->Close();
     


      return 0;
 
}

/**************************************************************************/
// End of Main function
/**************************************************************************/





 





