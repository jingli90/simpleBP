void MyRun(){
	TString sFilename = "/afs/cern.ch/work/j/jing/public/IPHCNtuple/fromXavier/trees_from_Xavier_20170302/ttHToNonbb.root";
	TString bFilename = "/afs/cern.ch/work/j/jing/public/IPHCNtuple/fromXavier/trees_from_Xavier_20170302/ttV.root";
	TString treeName = "Tree";
	TString outFilename = "output.root";
	gROOT->ProcessLine(".L simpleBP.C++");
	gROOT->ProcessLine(Form("simpleBP s(\"%s\",\"%s\",\"%s\",\"%s\")",sFilename.Data(),bFilename.Data(),treeName.Data(),outFilename.Data()));
	gROOT->ProcessLine(Form("s.AddVariable(\"nJet25_Recl\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"max_Lep_eta\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"MT_met_lep1\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"mindr_lep1_jet\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"mindr_lep2_jet\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"LepGood_conePt0\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"LepGood_conePt1\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddSpectator(\"is_3l_TTH_SR\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddSpectator(\"mc_ttZhypAllowed\",\"D\")")); //
	//gROOT->ProcessLine(Form("s.SetNLayer(1)")); //
	//gROOT->ProcessLine(Form("s.SetNNodes(1,8)")); //
	gROOT->ProcessLine(Form("s.SetNLayer(2)")); //
	gROOT->ProcessLine(Form("s.SetNNodes(1,7)")); //
	gROOT->ProcessLine(Form("s.SetNNodes(2,8)")); //
	gROOT->ProcessLine(Form("s.SetEta(0.02)"));
	gROOT->ProcessLine(Form("s.SetDecayRate(0.01)"));
	gROOT->ProcessLine(Form("s.SetNEpochs(350)")); //
	gROOT->ProcessLine(Form("s.SetTestRate(5)")); //
	//gROOT->ProcessLine(Form("s.SetCut(\"is_3l_TTH_SR==1\", \"Signal\")")); //
	//gROOT->ProcessLine(Form("s.SetCut(\"is_3l_TTH_SR==1\", \"Background\")")); //
	gROOT->ProcessLine(Form("s.SetCut(\"is_3l_TTH_SR==1 && mc_ttZhypAllowed==1\", \"Signal\")")); //
	gROOT->ProcessLine(Form("s.SetCut(\"is_3l_TTH_SR==1 && mc_ttZhypAllowed==1\", \"Background\")")); //
	gROOT->ProcessLine(Form("s.SetWeightExpression(\"weight\")"));
	gROOT->ProcessLine(Form("s.IsPrintEvolution(false)")); //
	gROOT->ProcessLine(Form("s.SetNeuronType(\"tanh\")"));
	//gROOT->ProcessLine(Form("s.SetNeuronType(\"sigmoid\")"));
	gROOT->ProcessLine(Form("s.SetThreshold(0.)"));
	gROOT->ProcessLine(Form("s.SetMinDeviate(0.1)"));
	gROOT->ProcessLine("s.myana()");
	//gROOT->ProcessLine(Form("s.~simpleBP()"));
}
