#ifndef simpleBP_h
#define simpleBP_h

#include <iostream>
using namespace std;
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>

#include "base.h"
#include "base_plot.h"

class simpleBP{
	public:
		TString sFilename;
		TString bFilename;
		TString treeName;
		TString outFilename;
		TFile * fout;
		TChain * c_s;
		TChain * c_b;
		base * b_s;
		base * b_b;
		base * b_s_test;
		base * b_b_test;
		TTree * train_s;
		TTree * train_b;
		TTree * test_s;
		TTree * test_b;
		TTree * tout_train_s;
		TTree * tout_train_b;
		TTree * tout_test_s;
		TTree * tout_test_b;
		Long64_t nentries;
		Long64_t nentries_s;
		Long64_t nentries_b;
		Long64_t nentries_s_test;
		Long64_t nentries_b_test;
		Long64_t nentries_s_train;
		Long64_t nentries_b_train;
		vector<TString>* var_input;
		vector<TString>* type_input;
		vector<double>* var_max;
		vector<double>* var_min;
		vector<double>* var_max_int;
		vector<double>* var_min_int;

		int nEpochs;
		int nLayer;
		int nNodes[20];
		double eta; // learning rate
		double decay_rate;
		double threshold;
		double minDev;
		int TestRate;
		vector<double>* var;
		vector<double>* net[20];
		vector<double>* o[20];
		vector<vector<double>*>* weight[20];
		vector<vector<double>*>* weight_old[20];
		double netDis;
		double discriminant;
		double difference, difference_largest, difference_largest_old;
		double variance, variance_old;
		double deltaDis;
		vector<double>* delta[20];

		double sum_weight;
		double evt_weight;

		TCanvas * c;
		TGraph * g_difference;
		TGraph * g_variance;
		TH1D * g_variance_HistTrain;
		TH1D * g_variance_HistTest;
		vector<vector<TGraph*>*>* g_weight[20];

		simpleBP(TString SignalFile, TString BackgroundFile, TString treename, TString outfilename);
		virtual void myana();
		virtual void AddVariable(TString var, TString type);
		virtual void SetNLayer(int nlayer);
		virtual void SetNNodes(int i, int nnodes);
		virtual void SetEta(double Eta);
		void SetDecayRate(double rate){decay_rate=rate;}
		void SetNEpochs(int N){nEpochs=N;}
		void SetTestRate(int N){TestRate=N;}
		virtual void SetThreshold(double Th);
		virtual void SetMinDeviate(double Dev);
		virtual void InitParameters();
		virtual void InitTrees();
		virtual void doTraining();
		virtual void calculate(base * b);
		virtual void back_propogation(int isSig);
		virtual void Eval();
		virtual void plotHist();
		virtual void plotROC();
		virtual void Shuffle(Int_t* index, Int_t n); // from TMVA package
		void DecayLearningRate(Double_t rate){ eta *= (1-rate); } // from TMVA package
		virtual double CalculateEstimator(TString treeType, Int_t iEpoch);
};

#endif

#ifdef simpleBP_cxx
simpleBP::simpleBP(TString SignalFile, TString BackgroundFile, TString treename, TString outfilename){
	cout<<endl;
	cout<<"Welcome to the training ..."<<endl;
	cout<<"Signal file: "<<SignalFile<<endl;
	cout<<"Background file: "<<BackgroundFile<<endl;
	sFilename = SignalFile;
	bFilename = BackgroundFile;
	treeName = treename;
	outFilename = outfilename;
	c_s = new TChain(treeName.Data(),"c_s");
	c_b = new TChain(treeName.Data(),"c_b");
	c_s->Add(sFilename.Data());
	c_b->Add(bFilename.Data());
	b_s = new base();
	b_b = new base();
	b_s_test = new base();
	b_b_test = new base();

	TFile * fout = new TFile(outFilename.Data(),"RECREATE");

	var_input = new vector<TString>;
	type_input = new vector<TString>;
	var_max = new vector<double>;
	var_min = new vector<double>;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;

	for(int i=0;i<20;i++){
		nNodes[i]=0;
		net[i] = new vector<double>;
		o[i] = new vector<double>;
		delta[i] = new vector<double>;
		weight[i] = new vector<vector<double>*>;
		weight_old[i] = new vector<vector<double>*>;
		g_weight[i] = new vector<vector<TGraph*>*>;
	}
	var = new vector<double>;
	cout<<endl;

	g_difference = new TGraph();
	g_variance = new TGraph();

	c=new TCanvas("c","c",10,10,700,700);

	eta=0.02; 
	decay_rate=0.01;
	threshold=0.;
	minDev=0.1;
	TestRate=10;
}

void simpleBP::SetNLayer(int nlayer){
	cout<<endl;
	nLayer=nlayer;
	if(nlayer>18){
		cout<<"Too much inner layers. Set number of inner layers = 18"<<endl;
		nLayer=18;
	}
	nNodes[nLayer+1]=1;
	cout<<"Number of inner layers = "<<nLayer<<endl;
}

void simpleBP::SetNNodes(int i, int nnodes){
	nNodes[i]=nnodes;
	cout<<"Number of nodes in inner layer "<<i<<" is "<<nNodes[i]<<endl;
}

void simpleBP::AddVariable(TString var, TString type){
	cout<<"Init input variable "<<var_input->size()<<endl;
	var_input->push_back(var);
	type_input->push_back(type); 
	double var0Max=TMath::Max(c_s->GetMaximum(var.Data()),c_b->GetMaximum(var.Data()));
	double var0Min=TMath::Min(c_s->GetMinimum(var.Data()),c_b->GetMinimum(var.Data()));
	cout<<var<<":["<<var0Min<<","<<var0Max<<"]"<<endl;
	double var0Max_int= ceil(var0Max);
	double var0Min_int=floor(var0Min);
	cout<<var<<"_int:["<<var0Min_int<<","<<var0Max_int<<"]"<<endl;
	var_max->push_back(var0Max);
	var_min->push_back(var0Min);
	var_max_int->push_back(var0Max_int);
	var_min_int->push_back(var0Min_int);

	nNodes[0]=var_input->size();
}

void simpleBP::SetEta(double Eta){
	eta=Eta;
	cout<<"eta="<<eta<<endl;
}

void simpleBP::SetThreshold(double Th){
	threshold=Th;
	cout<<"threshold="<<threshold<<endl;
}

void simpleBP::SetMinDeviate(double Dev){
	minDev=Dev;
	cout<<"Min Deviate="<<minDev<<endl;
}

void simpleBP::InitParameters(){
	cout<<endl;
	cout<<"Init weights:"<<endl;
	gRandom = new TRandom3();
	gRandom->SetSeed(0);
	int i_color=1;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			vector<double>* tmp_old = new vector<double>;
			vector<TGraph*>* tmp_graph = new vector<TGraph*>;
			weight[i]->push_back(tmp);
			weight_old[i]->push_back(tmp_old);
			g_weight[i]->push_back(tmp_graph);
			delta[i]->push_back(0);
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0)
					weight[i]->at(j)->push_back(0);
				else
					weight[i]->at(j)->push_back(2*gRandom->Rndm()-1);
				cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";

				TGraph * tmp_graph_k=new TGraph();
				tmp_graph_k->SetTitle(Form("Layer %d to layer %d: w_{%d%d};Event;weight", i-1, i, j, k));
				//tmp_graph_k->SetMarkerColor(i_color);
				i_color++;
				g_weight[i]->at(j)->push_back(tmp_graph_k);
			}
			cout<<endl;
		}
	}
}

void simpleBP::InitTrees(){
	cout<<endl;
	cout<<"Init test and training trees ..."<<endl;
	nentries=0, nentries_s=0, nentries_b=0;
	nentries_s_train=0, nentries_b_train=0, nentries_s_test=0, nentries_b_test=0;
	nentries_s = c_s->GetEntries();
	nentries_b = c_b->GetEntries();

	if(nentries_s==nentries_b){
		nentries=nentries_s;
		cout<<"total: "<<nentries<<endl;
	}
	else{
		cout<<"Error! sig & bkg have different entries!"<<endl;
		nentries=TMath::Min(nentries_s,nentries_b);
		cout<<"signal: "<<nentries_s<<" background: "<<nentries_b;
		cout<<"total: "<<nentries<<endl;
	}

	train_s = new TTree();
	train_s = c_s->CloneTree(0);
	train_b = new TTree();	
	train_b = c_b->CloneTree(0);
	for (Long64_t jentry=0; jentry<nentries;jentry=jentry+2.){
		c_s->GetEntry(jentry);
		c_b->GetEntry(jentry);
		train_s->Fill();
		train_b->Fill();
		nentries_s_train++;
		nentries_b_train++;
	}

	test_s = new TTree();
	test_s = c_s->CloneTree(0);
	for (Long64_t jentry=1; jentry<nentries_s;jentry=jentry+2.){
		c_s->GetEntry(jentry);
		test_s->Fill();
		nentries_s_test++;
	}

	test_b = new TTree();
	test_b = c_b->CloneTree(0);
	for (Long64_t jentry=1; jentry<nentries_b;jentry=jentry+2.){
		c_b->GetEntry(jentry);
		test_b->Fill();
		nentries_b_test++;
	}

	b_s->Init(train_s);
	b_b->Init(train_b);
	b_s_test->Init(test_s);
	b_b_test->Init(test_b);

	cout<<"signal training: "<<nentries_s_train<<" background training: "<<nentries_b_train<<endl;
	cout<<"signal test: "<<nentries_s_test<<" background test: "<<nentries_b_test<<endl;

}

void simpleBP::doTraining(){
	cout<<endl;
	cout<<"start training ..."<<endl;

	//b_s->Init(train_s);
	//b_b->Init(train_b);
	//b_s_test->Init(test_s);
	//b_b_test->Init(test_b);

	//g_variance_HistTrain = new TH1D( "estimatorHistTrain", "training estimator", nEpochs, 0., nEpochs );
	//g_variance_HistTest = new TH1D( "estimatorHistTest", "test estimator", nEpochs, 0., nEpochs );
	Int_t nbinTest = Int_t(nEpochs/TestRate);
	g_variance_HistTrain = new TH1D( "estimatorHistTrain", "training estimator", nbinTest, Int_t(TestRate/2), nbinTest*TestRate+Int_t(TestRate/2));
	g_variance_HistTest = new TH1D( "estimatorHistTest", "test estimator", nbinTest, Int_t(TestRate/2), nbinTest*TestRate+Int_t(TestRate/2) );
	// estimators
	Double_t trainE = -1;
	Double_t testE  = -1;

	bool isStop=true;
	int N_loop=0;
	//variance_old=nentries_s_train*4.;
	variance_old=2.;
	difference_largest_old=2.;
	int iPoint=0;

	// from TMVA package
	// randomize the order events will be presented, important for sequential mode
	int nEvents = nentries_s_train+nentries_b_train;
	Int_t* index = new Int_t[nEvents];
	int isSig;

	Int_t lateEpoch = (Int_t)(nEpochs*0.95);

	do{
		isStop=false;
		N_loop++;
		cout<<"N_loop="<<N_loop<<endl;
		variance=0;
		difference_largest=0;
		sum_weight=0.;

		// from TMVA package
		// randomize the order events will be presented, important for sequential mode
		for (Int_t i = 0; i < nEvents; i++) index[i] = i;
		Shuffle(index, nEvents);

		for(Long64_t jentry=0; jentry<nEvents;jentry++)
			//for(Long64_t jentry=0; jentry<10;jentry++)
		{
			//if(jentry%1000==0)cout<<"entry: "<<jentry<<endl;
			//cout<<"entry: "<<jentry<<endl;
			int jentry_tmp=index[jentry];
			if(jentry_tmp<nentries_s_train){
				isSig=1;
				b_s->GetEntry(jentry_tmp);
				calculate(b_s);
				evt_weight=b_s->weight;
			}
			else{
				//isSig=-1;
				isSig=0; // test for TMVA
				b_b->GetEntry(jentry_tmp-nentries_s_train);
				calculate(b_b);
				evt_weight=b_b->weight;
			}

			sum_weight=sum_weight+evt_weight;
			//cout<<"discriminant="<<discriminant<<endl;

			difference=(discriminant-isSig)*(discriminant-isSig)/2.;
			//difference=(discriminant-isSig)*(discriminant-isSig); // test for TMVA
			difference=difference*evt_weight;
			//if(difference>2)cout<<"jentry_tmp="<<jentry_tmp<<" discriminant="<<discriminant<<" isSig="<<isSig<<endl;
			//cout<<"difference="<<difference<<endl;
			if(difference>difference_largest)difference_largest=difference;
			variance=variance+difference;
			if(difference>threshold){
				back_propogation(isSig);
			}

			for(int i=1;i<=nLayer+1;i++){
				for(int j=1;j<=nNodes[i];j++){
					for(int k=0;k<=nNodes[i-1];k++){
						g_weight[i]->at(j)->at(k)->SetPoint(iPoint,iPoint,weight[i]->at(j)->at(k));
					}
				}
			}
			g_difference->SetPoint(iPoint,iPoint,difference);
			iPoint++;
		}
		//cout<<"largest diff:"<<difference_largest<<endl;
		//cout<<"sum_weight:"<<sum_weight<<endl;
		variance=variance/sum_weight;
		//cout<<"variance="<<variance<<endl;
		g_variance->SetPoint(N_loop-1, N_loop, variance);
		if(N_loop%TestRate==0){
			trainE = CalculateEstimator( "Training", N_loop );
			testE  = CalculateEstimator( "Testing",  N_loop );
			g_variance_HistTrain->Fill( N_loop, trainE );
			g_variance_HistTest->Fill( N_loop, testE );
		}
		isStop=false;
		bool isStop1=TMath::Abs((variance-variance_old)/variance_old)<0.001;
		bool isStop2=TMath::Abs((difference_largest-difference_largest_old)/difference_largest_old)<0.001;
		//if(isStop1 && isStop2)isStop=true;
		if(N_loop>=nEpochs)isStop=true;
		variance_old=variance;
		difference_largest_old=difference_largest;
		if(N_loop >= lateEpoch)
			DecayLearningRate(TMath::Sqrt(decay_rate));
		else
			DecayLearningRate(decay_rate);
		//for(int i=1;i<=nLayer+1;i++){
		//	cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		//	for(int j=0;j<=nNodes[i];j++){
		//		for(int k=0;k<=nNodes[i-1];k++){
		//			cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";
		//		}
		//		cout<<endl;
		//	}
		//}
	}while(!isStop);
	cout<<"N_loop="<<N_loop<<endl;
	//cout<<"iPoint="<<iPoint<<endl;

	//TCanvas * c = new TCanvas("c","c",10,10,700,700);
	c->cd();
	TLatex * ltx = new TLatex();
	ltx->SetNDC(kTRUE);
	ltx->SetTextFont(22);
	ltx->SetTextSize(0.03);

	cout<<"Final weights:"<<endl;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";
			}
			cout<<endl;
		}
	}
	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=1;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				c->Clear();
				g_weight[i]->at(j)->at(k)->GetYaxis()->SetTitleOffset(1.5);
				g_weight[i]->at(j)->at(k)->Draw("AP");
				ltx->DrawLatex(0.15,0.85,Form("#splitline{nEntries=%d, nLoop=%d}{w_{%d%d}=%f}", nentries_s_train+nentries_b_train, N_loop, j, k, weight[i]->at(j)->at(k)));
				c->SaveAs(Form("layer%dtolayer%d_w%d%d.png", i-1, i, j, k));
			}
			cout<<endl;
		}
	}

	c->Clear();
	TPad * p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
	p1->SetFillColor(0);
	p1->Draw();
	p1->SetLeftMargin(0.14);
	p1->cd();
	g_difference->SetTitle(";Event;(#hat{y}-y)^{2}/2");
	g_difference->GetYaxis()->SetTitleOffset(2);
	g_difference->Draw("AP");
	ltx->DrawLatex(0.15,0.85,Form("nEntries=%d, nLoop=%d", nentries_s_train+nentries_b_train, N_loop));
	c->SaveAs("difference.png");

	c->Clear();
	TPad * p2 = new TPad("p2","p2", 0.00,0.00,1.00,0.97);
	p2->SetFillColor(0);
	p2->Draw();
	p2->SetLeftMargin(0.14);
	p2->cd();
	g_variance->SetTitle(";Epoch;#Sigma (w*(#hat{y}-y))^{2}/2/Sumw");
	g_variance->GetYaxis()->SetTitleOffset(2);
	g_variance->Draw("AP*");
	ltx->DrawLatex(0.15,0.85,Form("nEntries=%d, nLoop=%d", nentries_s_train+nentries_b_train, N_loop));
	c->SaveAs("variance.png");

	c->Clear();
	p2=new TPad("p2","p2", 0.00,0.00,1.00,0.97);
	p2->SetFillColor(0);
	p2->Draw();
	p2->SetLeftMargin(0.14);
	p2->SetTicks(1,1);
	p2->cd();
	gStyle->SetOptStat(0);
	g_variance_HistTrain->SetTitle("MLP Convergence Test; Epochs; Estimator");
	g_variance_HistTrain->SetMaximum((1.2*TMath::Max(g_variance_HistTrain->GetMaximum(), g_variance_HistTest->GetMaximum())));
	g_variance_HistTrain->GetYaxis()->SetTitleOffset(2);
	g_variance_HistTrain->SetLineColor(1);
	g_variance_HistTest->SetLineColor(2);
	g_variance_HistTrain->Draw();
	g_variance_HistTest->Draw("Same");
	TLegend* legend = new TLegend(0.7,0.75,0.9,0.9,"");
	legend->SetFillColor(kWhite);
	legend->AddEntry("estimatorHistTrain", "Training Sample", "l");
	legend->AddEntry("estimatorHistTest", "Test Sample", "l");
	legend->Draw();
	c->SaveAs("annconvergencetest.png");
}

void simpleBP::calculate(base * b){
	var->clear();
	net[0]->clear();
	o[0]->clear();

	net[0]->push_back(1);
	o[0]->push_back(1);

	for(int ivar=0;ivar<var_input->size();ivar++){
		double var_tmp=b->GetVal(var_input->at(ivar));
		var->push_back(2 * (var_tmp - var_min_int->at(ivar)) / (var_max_int->at(ivar) - var_min_int->at(ivar)) - 1);
		o[0]->push_back(var->at(ivar));
		//cout<<var_input->at(ivar)<<"="<<var_tmp<<endl;
		//cout<<var_input->at(ivar)<<"_nrm="<<var->at(ivar)<<endl;
	}

	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		net[i]->clear();
		o[i]->clear();
		net[i]->push_back(1);
		o[i]->push_back(1);
		for(int j=1;j<=nNodes[i];j++){
			double net_tmp=0;
			double o_tmp=0;
			for(int k=0;k<=nNodes[i-1];k++){
				net_tmp=net_tmp + weight[i]->at(j)->at(k) * o[i-1]->at(k);
				//cout<<"w"<<j<<k<<endl;
			}
			o_tmp=(TMath::Exp(net_tmp)-TMath::Exp(-net_tmp))/(TMath::Exp(net_tmp)+TMath::Exp(-net_tmp));
			net[i]->push_back(net_tmp);
			o[i]->push_back(o_tmp);
		}
	}

	netDis=net[nLayer+1]->at(1);
	discriminant=o[nLayer+1]->at(1);
	//cout<<"discriminant="<<discriminant<<endl;

}

void simpleBP::back_propogation(int isSig){
	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			weight_old[i]->at(j)->clear();
			for(int k=0;k<=nNodes[i-1];k++){
				double tmp = weight[i]->at(j)->at(k);
				weight_old[i]->at(j)->push_back(tmp);
				//cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" "<<weight_old[i]->at(j)->at(k)<<" ";
			}
			//cout<<endl;
		}
	}

	for(int i=1;i<=nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			weight[i]->at(j)->clear();
		}
	}

	deltaDis=-(isSig-discriminant)*(1+discriminant)*(1-discriminant);
	//cout<<"nLayer+1="<<nLayer+1<<endl;
	delta[nLayer+1]->clear();
	delta[nLayer+1]->push_back(0);
	delta[nLayer+1]->push_back(deltaDis);

	for(int i=nLayer; i>=1; i--){
		//cout<<"Inner layer "<<i+1<<" to "<<i<<endl;
		delta[i]->clear();
		for(int j=0;j<=nNodes[i];j++){
			double delta_tmp=0;
			for(int k=1;k<=nNodes[i+1];k++){
				double tmp_devNetj = ( 1 + o[i]->at(j) ) * ( 1 - o[i]->at(j) );
				double tmp_wkj = weight_old[i+1]->at(k)->at(j);
				delta_tmp = delta_tmp + tmp_devNetj  * tmp_wkj * delta[i+1]->at(k);
			}
			delta[i]->push_back(delta_tmp);
			//cout<<"delta"<<j<<" = "<<delta[i]->at(j)<<endl;
		}
		//cout<<endl;
	}

	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				double weight_tmp_change =  eta * delta[i]->at(j) * o[i-1]->at(k) * evt_weight;
				double weight_tmp = weight_old[i]->at(j)->at(k) - weight_tmp_change;
				weight[i]->at(j)->push_back(weight_tmp);
				//cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" "<<weight_old[i]->at(j)->at(k)<<" ";
			}
			//cout<<endl;
		}
	}
}

Double_t simpleBP::CalculateEstimator( TString treeType, Int_t iEpoch ){
	// calculate the estimator that training is attempting to minimize

	double estimator = 0.;

	// sanity check
	if (treeType!="Training" && treeType!="Testing") {
		cout<<"<CalculateEstimator> fatal error: wrong tree type: "<<treeType << endl;
	}
	base * b_s_tmp;
	base * b_b_tmp;
	Int_t nEvents;
	Int_t nentries_s_tmp, nentries_b_tmp;
	if(treeType=="Training"){
		b_s_tmp=b_s;
		b_b_tmp=b_b;
		nentries_s_tmp=nentries_s_train;
		nentries_b_tmp=nentries_b_train;
		nEvents=nentries_s_train+nentries_b_train;
	}
	if(treeType=="Testing"){
		b_s_tmp=b_s_test;
		b_b_tmp=b_b_test;
		nentries_s_tmp=nentries_s_test;
		nentries_b_tmp=nentries_b_test;
		nEvents=nentries_s_test+nentries_b_test;
	}

	double difference_tmp;
	double sum_weight_tmp = 0.;
	double evt_weight_tmp;
	int isSig;

	for(Long64_t jentry=0; jentry<nEvents;jentry++){
		if(jentry<nentries_s_tmp){
			isSig=1;
			b_s_tmp->GetEntry(jentry);
			calculate(b_s_tmp);
			evt_weight_tmp=b_s_tmp->weight;
		}
		else{
			//isSig=-1;
			isSig=0; // test for TMVA
			b_b_tmp->GetEntry(jentry-nentries_s_tmp);
			calculate(b_b_tmp);
			evt_weight_tmp=b_b_tmp->weight;
		}
		sum_weight_tmp=sum_weight_tmp+evt_weight_tmp;
		//difference_tmp=(discriminant-isSig)*(discriminant-isSig)/2.;
		difference_tmp=(discriminant-isSig)*(discriminant-isSig); // test for TMVA
		difference_tmp=difference_tmp*evt_weight;
		//if(jentry<10 && (iEpoch%100==0))
		//	cout<<"jentry="<<jentry<<" isSig="<<isSig<<" discriminant="<<discriminant<<endl;
		//if(jentry>=nentries_s_tmp && jentry<nentries_s_tmp+10 && (iEpoch%100==0))
		//	cout<<"jentry="<<jentry<<" isSig="<<isSig<<" discriminant="<<discriminant<<endl;
		estimator=estimator+difference_tmp;
	}
	estimator=estimator/sum_weight_tmp;

	return estimator;
}

void simpleBP::Eval(){
	//b_s->Init(train_s);
	//b_b->Init(train_b);
	tout_train_s = train_s->CloneTree(0);
	tout_train_s->Branch("discriminant",&discriminant,"discriminant/D");
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_train_s->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	tout_train_b = new TTree();	
	tout_train_b = train_b->CloneTree(0);
	tout_train_b->Branch("discriminant",&discriminant,"discriminant/D");	
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_train_b->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for (Long64_t jentry=0; jentry<nentries_s_train;jentry++){
		train_s->GetEntry(jentry);
		calculate(b_s);
		tout_train_s->Fill();
		train_b->GetEntry(jentry);
		calculate(b_b);
		tout_train_b->Fill();
	}

	//b_s->Init(test_s);
	tout_test_s = new TTree();
	tout_test_s = test_s->CloneTree(0);
	tout_test_s->Branch("discriminant",&discriminant,"discriminant/D"); 
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_test_s->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for (Long64_t jentry=0; jentry<nentries_s_test;jentry++){
		test_s->GetEntry(jentry);
		calculate(b_s_test);
		tout_test_s->Fill();
	}

	//b_b->Init(test_b);
	tout_test_b = new TTree();
	tout_test_b = test_b->CloneTree(0);
	tout_test_b->Branch("discriminant",&discriminant,"discriminant/D"); 
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_test_b->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for (Long64_t jentry=0; jentry<nentries_b_test;jentry++){
		test_b->GetEntry(jentry);
		calculate(b_b_test);
		tout_test_b->Fill();
	}

	tout_train_s->Write("train_s");
	tout_train_b->Write("train_b");
	tout_test_s->Write("test_s");
	tout_test_b->Write("test_b");
	//fout->Close();

}

void simpleBP::plotHist(){


	c->cd();
	gStyle->SetOptStat(0);

	c->Clear();
	TPad * p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
	p1->SetFillColor(0);
	p1->Draw();
	p1->SetLeftMargin(0.14);
	p1->cd();
	TH1D * h_s_train = GetHistoWeight(tout_train_s, "discriminant", 50, -1., 1., "1==1", "weight", "h_s_train");
	TH1D * h_b_train = GetHistoWeight(tout_train_b, "discriminant", 50, -1., 1., "1==1", "weight", "h_b_train");
	TH1D * h_s_test = GetHistoWeight(tout_test_s, "discriminant", 50, -1., 1., "1==1", "weight", "h_s_test");
	TH1D * h_b_test = GetHistoWeight(tout_test_b, "discriminant", 50, -1., 1., "1==1", "weight", "h_b_test");

	h_s_test->SetLineColor(2);
	h_b_test->SetLineColor(1);
	h_s_test->SetMaximum((1.5*TMath::Max(h_s_test->GetMaximum(), h_b_test->GetMaximum())));
	h_s_test->SetTitle(";Discriminant; Event / bin");
	h_s_test->GetYaxis()->SetTitleOffset(1.8);
	h_s_test->Draw();
	h_b_test->Draw("same");
	h_s_train->SetMarkerColor(2);
	h_s_train->SetMarkerStyle(7);
	h_b_train->SetMarkerColor(1);
	h_b_train->SetMarkerStyle(7);
	h_s_train->Draw("Psame");
	h_b_train->Draw("Psame");

	TLegend* legend_test = new TLegend(0.15,0.75,0.5,0.9,"");
	legend_test->SetFillColor(kWhite);
	legend_test->AddEntry("h_s_test", "Signal (test sample)", "l");
	legend_test->AddEntry("h_b_test", "Background (test sample)", "l");
	legend_test->Draw();

	TLegend* legend_train = new TLegend(0.5,0.75,0.85,0.9,"");
	legend_train->SetFillColor(kWhite);
	legend_train->AddEntry("h_s_train", "Signal (training sample)", "p");
	legend_train->AddEntry("h_b_train", "Background (training sample)", "p");
	legend_train->Draw();

	TLatex * ltx = new TLatex();
	ltx->SetNDC(kTRUE);
	ltx->SetTextFont(22);
	ltx->SetTextSize(0.03);

	double x_ltx=0.25, y_ltx=0.7;

	for(int i=1;i<=nLayer+1;i++){
		y_ltx=0.7;
		ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
		y_ltx=y_ltx-0.04;
		for(int j=1;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
				y_ltx=y_ltx-0.04;
			}
			cout<<endl;
		}
		x_ltx=x_ltx+0.2;
	}

	c->Print("discriminant.png");

	for(int i=0;i<nLayer+1;i++){
		for(int j=1;j<=nNodes[i];j++){
			c->Clear();
			p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
			p1->SetFillColor(0);
			p1->Draw();
			p1->SetLeftMargin(0.14);
			p1->cd();
			h_s_train = GetHistoWeight(tout_train_s, Form("l%dn%d",i,j), 50, -1., 1., "1==1", "weight", "h_s_train");
			h_b_train = GetHistoWeight(tout_train_b, Form("l%dn%d",i,j), 50, -1., 1., "1==1", "weight", "h_b_train");
			h_s_test = GetHistoWeight(tout_test_s,   Form("l%dn%d",i,j), 50, -1., 1., "1==1", "weight", "h_s_test");
			h_b_test = GetHistoWeight(tout_test_b,   Form("l%dn%d",i,j), 50, -1., 1., "1==1", "weight", "h_b_test");

			h_s_test->SetLineColor(2);
			h_b_test->SetLineColor(1);
			h_s_test->SetMaximum((1.5*TMath::Max(h_s_test->GetMaximum(),h_b_test->GetMaximum())));
			if(i==0)
				h_s_test->SetTitle(Form(";%s (norm);Event / bin",var_input->at(j-1).Data()));
			else
				h_s_test->SetTitle(Form(";Output of layer %d node %d;Event / bin",i,j));
			h_s_test->GetYaxis()->SetTitleOffset(1.8);
			h_s_test->Draw();
			h_b_test->Draw("same");
			h_s_train->SetMarkerColor(2);
			h_s_train->SetMarkerStyle(7);
			h_b_train->SetMarkerColor(1);
			h_b_train->SetMarkerStyle(7);
			h_s_train->Draw("Psame");
			h_b_train->Draw("Psame");

			legend_test->Clear();
			legend_test->SetFillColor(kWhite);
			legend_test->AddEntry("h_s_test", "Signal (test sample)", "l");
			legend_test->AddEntry("h_b_test", "Background (test sample)", "l");
			legend_test->Draw();

			legend_train->Clear();
			legend_train->SetFillColor(kWhite);
			legend_train->AddEntry("h_s_train", "Signal (training sample)", "p");
			legend_train->AddEntry("h_b_train", "Background (training sample)", "p");
			legend_train->Draw();

			c->Print(Form("O_l%dn%d.png",i,j));
		}
	}
}

void simpleBP::plotROC(){
	c->cd();

	c->Clear();
	TPad * p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
	p1->SetFillColor(0);
	p1->Draw();
	p1->SetTicks(1,1);
	p1->SetGridx();
	p1->SetGridy();
	p1->cd();
	TGraph * g_ROC_train=GetEffSvsEffB(tout_train_s, tout_train_b, "1==1", "1==1", "discriminant", -1, 1, "weight", 50, "g_ROC_train");
	TH2D* hGrid = new TH2D("Grid","Grid",1000,0,1,1000,0,1);
	hGrid->Draw();
	hGrid->SetTitle("training sample");
	hGrid->GetYaxis()->SetTitle("Signal Efficiency");
	hGrid->GetXaxis()->SetTitle("Background Efficiency");
	hGrid->GetYaxis()->SetTitleOffset(1.4);
	g_ROC_train->Draw("*same");
	TLatex * ltx = new TLatex();
	ltx->SetNDC(kTRUE);
	ltx->SetTextFont(22);
	ltx->SetTextSize(0.03);
	double x_ltx=0.5, y_ltx=0.5;
	for(int i=1;i<=nLayer+1;i++){
		y_ltx=0.5;
		ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
		y_ltx=y_ltx-0.04;
		for(int j=1;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
				y_ltx=y_ltx-0.04;
			}
			cout<<endl;
		}
		x_ltx=x_ltx+0.2;
	}
	c->SaveAs("ROC_train.png");

	c->Clear();
	TPad * p2 = new TPad("p2","p2", 0.00,0.00,1.00,0.97);
	p2->SetFillColor(0);
	p2->Draw();
	p2->SetTicks(1,1);
	p2->SetGridx();
	p2->SetGridy();
	p2->cd();
	TGraph * g_ROC_test=GetEffSvsEffB(tout_test_s, tout_test_b, "1==1", "1==1", "discriminant", -1, 1, "weight", 50, "g_ROC_test");
	hGrid->SetTitle("test sample");
	hGrid->Draw();
	g_ROC_test->Draw("*same");
	x_ltx=0.5, y_ltx=0.5;
	for(int i=1;i<=nLayer+1;i++){
		y_ltx=0.5;
		ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
		y_ltx=y_ltx-0.04;
		for(int j=1;j<=nNodes[i];j++){
			for(int k=0;k<=nNodes[i-1];k++){
				ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
				y_ltx=y_ltx-0.04;
			}
			cout<<endl;
		}
		x_ltx=x_ltx+0.2;
	}
	c->SaveAs("ROC_test.png");
}

void simpleBP::Shuffle(Int_t* index, Int_t n){
	// from TMVA package
	// Input:
	//   index: the array to shuffle
	//   n: the size of the array
	// Output:
	//   index: the shuffled indexes
	// This method is used for sequential training

	Int_t j, k;
	Int_t a = n - 1;
	gRandom = new TRandom3();
	gRandom->SetSeed(0);
	for (Int_t i = 0; i < n; i++) {
		j = (Int_t) (gRandom->Rndm() * a);
		k = index[j];
		index[j] = index[i];
		index[i] = k;
	}

}

#endif
