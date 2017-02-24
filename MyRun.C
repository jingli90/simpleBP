void MyRun(){
	TString sFilename = "signal_10000.root";
	TString bFilename = "background_10000.root";
	TString treeName = "t";
	TString outFilename = "output.root";
	gROOT->ProcessLine(".L base.C++");
	gROOT->ProcessLine(".L base_plot.h+");
	gROOT->ProcessLine(".L simpleBP.C+");
	gROOT->ProcessLine(Form("simpleBP s(\"%s\",\"%s\",\"%s\",\"%s\")",sFilename.Data(),bFilename.Data(),treeName.Data(),outFilename.Data()));
	gROOT->ProcessLine(Form("s.AddVariable(\"x\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"y\",\"D\")"));
	gROOT->ProcessLine(Form("s.AddVariable(\"z\",\"D\")"));
	gROOT->ProcessLine(Form("s.SetNLayer(1)"));
	//gROOT->ProcessLine(Form("s.SetNNodes(1,2)"));
	gROOT->ProcessLine(Form("s.SetNNodes(1,4)"));
	gROOT->ProcessLine(Form("s.SetEta(0.02)"));
	gROOT->ProcessLine(Form("s.SetDecayRate(0.01)"));
	gROOT->ProcessLine(Form("s.SetNEpochs(200)"));
	//gROOT->ProcessLine(Form("s.SetThreshold(0.1)"));
	gROOT->ProcessLine(Form("s.SetThreshold(0.)"));
	gROOT->ProcessLine(Form("s.SetMinDeviate(0.1)"));
	gROOT->ProcessLine("s.myana()");
	//gROOT->ProcessLine(Form("s.~simpleBP()"));
}
