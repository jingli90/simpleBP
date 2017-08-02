#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "simpleBP.h"
//#include "base.h"
//#include "base_plot.h"
//#include "base_function.h"

int main(int argc, char *argv[]){
	TString Filename_ttZ = "/afs/cern.ch/work/j/jing/private/FCPPL/learning_NN/regress/smallTree/small_ttZ.root";
	TString Filename_ttW = "/afs/cern.ch/work/j/jing/private/FCPPL/learning_NN/regress/smallTree/small_ttW.root";
	TString Filename_ttH = "/afs/cern.ch/work/j/jing/private/FCPPL/learning_NN/regress/smallTree/small_ttH.root";

	TString treeName = "Tree";
	TString outFilename = "output.root";

	simpleBP s(Filename_ttZ, treeName, outFilename, "ttZ");
	s.AddOtherFiles(Filename_ttW, "ttW");
	s.AddOtherFiles(Filename_ttH, "ttH");
	s.SetNLayer(1);
	//s.SetNNodes(1,10);
	//s.SetNNodes(1,24);
	s.SetNNodes(1,30);
	//s.SetNNodes(1,50);
	//s.SetNum_outputNN(5);
	s.SetNum_outputNN(1);
	s.SetEta(0.02);
	s.SetDecayRate(0.01);
	s.SetNEpochs(1000); //
	//s.SetNEpochs(250); //
	//s.SetNEpochs(10); //
	s.SetTestRate(5);
	s.SetNeuronType("tanh");
	s.SetCut("is_3l_TTH_SR==1 && mc_ttZhypAllowed==1 && catJets==0 && hasInputVars==1","Signal");
	s.SetCut("is_3l_TTH_SR==1 && mc_ttZhypAllowed==1 && catJets==0 && hasInputVars==1","Background");
	s.SetWeightExpression("weight");
	s.IsPrintEvolution(false);
	//s.SetTrainingMode("BP");
	s.SetTrainingMode("regress");

	s.AddVariable("multilepton_Bjet1_P4->E()");
	s.AddVariable("multilepton_Bjet1_P4->Theta()");
	s.AddVariable("multilepton_Bjet1_P4->Phi()");
	s.AddVariable("multilepton_Bjet2_P4->E()");
	s.AddVariable("multilepton_Bjet2_P4->Theta()");
	s.AddVariable("multilepton_Bjet2_P4->Phi()");
	s.AddVariable("multilepton_JetClosestMw1_P4->E()");
	s.AddVariable("multilepton_JetClosestMw1_P4->Theta()");
	s.AddVariable("multilepton_JetClosestMw1_P4->Phi()");
	s.AddVariable("multilepton_JetClosestMw2_P4->E()");
	s.AddVariable("multilepton_JetClosestMw2_P4->Theta()");
	s.AddVariable("multilepton_JetClosestMw2_P4->Phi()");
	s.AddVariable("multilepton_Lepton1_P4->E()");
	s.AddVariable("multilepton_Lepton1_P4->Theta()");
	s.AddVariable("multilepton_Lepton1_P4->Phi()");
	s.AddVariable("multilepton_Lepton2_P4->E()");
	s.AddVariable("multilepton_Lepton2_P4->Theta()");
	s.AddVariable("multilepton_Lepton2_P4->Phi()");
	s.AddVariable("multilepton_Lepton3_P4->E()");
	s.AddVariable("multilepton_Lepton3_P4->Theta()");
	s.AddVariable("multilepton_Lepton3_P4->Phi()");
	s.AddVariable("multilepton_mET->Pt()");
	s.AddVariable("multilepton_mET->Phi()");
	s.AddSpectator("is_3l_TTH_SR");
	s.AddSpectator("mc_ttZhypAllowed");
	s.AddSpectator("catJets");

	vector<double>* xL;
	vector<double>* xU;
	xU = new vector<double>;
	xL = new vector<double>;
	xU->push_back(1500);
	xU->push_back(1.563771);
	xU->push_back(1500);
	xU->push_back(6.283185);
	xU->push_back(1.563771);
	xL->push_back(4.7);
	xL->push_back(-1.545340);
	xL->push_back(4.7);
	xL->push_back(0);
	xL->push_back(-1.545340);
	s.SetXLXU(xL, xU);

	vector<double>* var_max_int;
	vector<double>* var_min_int;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;
	var_min_int->push_back(25); // multilepton_Bjet1_P4->E()
	var_max_int->push_back(2603);
	var_min_int->push_back(0); // multilepton_Bjet1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Bjet1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_Bjet2_P4->E()
	var_max_int->push_back(2203);
	var_min_int->push_back(0); // multilepton_Bjet2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Bjet2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_JetClosestMw1_P4->E()
	var_max_int->push_back(3170);
	var_min_int->push_back(0); // multilepton_JetClosestMw1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_JetClosestMw1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_JetClosestMw2_P4->E()
	var_max_int->push_back(2879);
	var_min_int->push_back(0); // multilepton_JetClosestMw2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_JetClosestMw2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(21); // multilepton_Lepton1_P4->E()
	var_max_int->push_back(2437);
	var_min_int->push_back(0); // multilepton_Lepton1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(10); // multilepton_Lepton2_P4->E()
	var_max_int->push_back(1343);
	var_min_int->push_back(0); // multilepton_Lepton2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_Lepton3_P4->E()
	var_max_int->push_back(755);
	var_min_int->push_back(0); // multilepton_Lepton3_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton3_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_mET->Pt()
	var_max_int->push_back(860);
	var_min_int->push_back(-4); // multilepton_mET->Phi() 
	var_max_int->push_back(4);
	s.SetMinMax(var_min_int, var_max_int);

	//s.SetBatch(-1); // -1 loop over all; 0 event by event; n per n events
	//s.SetBatch(1);
	s.SetBatch(500);

	s.myana();

}
