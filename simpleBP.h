#ifndef simpleBP_h
#define simpleBP_h

#include <iostream>
using namespace std;
#include <vector>
#include <time.h>

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
#include "base_function.h"

class simpleBP{
	public:
		TString Filename;
		TString treeName;
		TString label;
		TString Filename_other[10];
		TString Filelabel_other[10];
		TString outFilename;
		int nOther;

		TFile * fout;
		TChain * c_s;
		TChain * c_b;
		TTree * t_s_cut;
		TTree * t_b_cut;
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

		int nEpochs;
		int nLayer;
		double eta; // learning rate
		double decay_rate;
		double threshold;
		double minDev;
		int TestRate;
		TString NeuronType;
		int num_outputNN;
		TString cut_s;
		TString cut_b;
		TString weightExpression;
		bool isPrintEvolution;
		TString TrainingMode;

		vector<double>* var;
		vector<TString>* var_input;
		//vector<TString>* type_input;
		vector<TString>* spectator_var_input;
		//vector<TString>* spectator_type_input;
		//vector<double>* var_max;
		//vector<double>* var_min;
		vector<double>* var_max_int;
		vector<double>* var_min_int;
		vector<double>* xL;
		vector<double>* xU;
		vector<double>* var_input_NN_beforeNorm;
		vector<double>* var_kin_ttz;

		int nNodes[20];
		vector<double>* net[20];
		vector<double>* o[20];
		vector<vector<double>*>* weight[20];
		vector<vector<double>*>* weight_old[20];
		vector<vector<double>*>* weight_init[20];
		vector<vector<double>*>* weight_memory[20];
		double netDis;
		double discriminant;
		double difference, difference_largest, difference_largest_old;
		double variance, variance_old;
		double deltaDis;
		vector<double>* delta[20];
		vector<vector<TGraph*>*>* g_weight[20];

		int isBatch; // -1 loop over all; 0 event by event; n per n events
		int sizeBatch;
		int numBatch_tmp;

		double duration;

		double sum_weight;
		double evt_weight;
		double evt_weight_max;

		TCanvas * c;
		//TGraph * g_difference;
		//TGraph * g_variance;
		TH1D * g_variance_HistTrain;
		TH1D * g_variance_HistTest;

		simpleBP(TString File, TString treename, TString outfilename, TString label_);
		void AddOtherFiles(TString other, TString label){Filename_other[nOther]=other; Filelabel_other[nOther]=label; nOther++; cout<<"will eval also on "<<other<<endl;}
		virtual void myana();
		virtual void AddVariable(TString var){cout<<"Add input variable: "<<var<<endl; var_input->push_back(var);}
		virtual void AddSpectator(TString var){cout<<"Add spectator: "<<var<<endl; spectator_var_input->push_back(var);}
		virtual void SetNLayer(int nlayer);
		virtual void SetNNodes(int i, int nnodes);
		void SetNum_outputNN(int N){num_outputNN=N; nNodes[nLayer+1]=N; cout<<"number of input to MEM: "<< N <<endl;}
		void SetEta(double Eta){eta=Eta; cout<<"eta="<<eta<<endl;}
		void SetDecayRate(double rate){decay_rate=rate;cout<<"Decay rate = "<<rate<<endl;}
		void SetNEpochs(int N){nEpochs=N;cout<<"nEpochs = "<< N <<endl;}
		void SetTestRate(int N){TestRate=N;cout<<"TestRate = "<< N <<endl;}
		void SetThreshold(double Th){threshold=Th; cout<<"threshold="<<threshold<<endl;}
		void SetMinDeviate(double Dev){minDev=Dev;  cout<<"Min Deviate="<<minDev<<endl;}
		void SetNeuronType(TString myNeuronType){NeuronType=myNeuronType; cout<<"NeuronType="<<NeuronType<<endl;}
		virtual void SetCut(const TString& cut, const TString& className);
		void SetWeightExpression(TString expression){weightExpression=expression; cout<<"Weight expression is: \""<< weightExpression<<"\""<<endl;}
		void IsPrintEvolution(bool b){isPrintEvolution=b; cout<<"isPrintEvolution="<<isPrintEvolution<<endl;}
		void SetTrainingMode(TString mode){TrainingMode=mode; cout<<"TrainingMode="<<TrainingMode<<endl;}
		void SetXLXU(vector<double>* xL_, vector<double>* xU_){xL=xL_; xU=xU_;}
		void SetMinMax(vector<double>* min, vector<double>* max){var_max_int=max; var_min_int=min;}
		void SetBatch(int N){isBatch=N; cout<<"isBatch="<<isBatch<<endl;}
		virtual void InitParameters();
		virtual void InitTrees();
		virtual void doTraining();
		virtual void calculate(base * b);
		virtual void back_propogation(int isSig);
		virtual void Eval();
		virtual void Eval_regress();
		virtual void plotHist();
		virtual void plotHist_regress();
		virtual void plotROC();
		virtual void Shuffle(Int_t* index, Int_t n); // from TMVA package
		void DecayLearningRate(Double_t rate){ eta *= (1-rate); } // from TMVA package
		virtual double CalculateEstimator(TString treeType, Int_t iEpoch);
		//virtual void plotWeights();
};



simpleBP::simpleBP(TString File, TString treename, TString outfilename, TString label_){
	cout<<endl;
	cout<<"Welcome to the training ..."<<endl;
	cout<<"Filename: "<<File<<endl;
	Filename=File;
	treeName = treename;
	outFilename = outfilename;
	label=label_;
	nOther=0;

	for(int i=0;i<20;i++){
		nNodes[i]=0;
		net[i] = new vector<double>;
		o[i] = new vector<double>;
		delta[i] = new vector<double>;
		weight[i] = new vector<vector<double>*>;
		weight_old[i] = new vector<vector<double>*>;
		weight_init[i] = new vector<vector<double>*>;
		weight_memory[i] = new vector<vector<double>*>;
		g_weight[i] = new vector<vector<TGraph*>*>;
	}

	eta=0.02; 
	decay_rate=0.01;
	threshold=0.;
	minDev=0.1;
	TestRate=10;
	NeuronType="tanh";

	//c_s = new TChain(treeName.Data(),"c_s");
	//c_b = new TChain(treeName.Data(),"c_b");
	//c_s->Add(Filename.Data());
	//c_b->Add(bFilename.Data());
	//b_s = new base();
	//b_b = new base();
	//b_s_test = new base();
	//b_b_test = new base();

	fout = new TFile(outFilename.Data(),"RECREATE");

	var_input = new vector<TString>;
	//type_input = new vector<TString>;
	spectator_var_input = new vector<TString>;
	//spectator_type_input = new vector<TString>;
	//var_max = new vector<double>;
	//var_min = new vector<double>;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;
	xU = new vector<double>;
	xL = new vector<double>;
	var_input_NN_beforeNorm = new vector<double>;
	var_kin_ttz = new vector<double>;

	var = new vector<double>;
	cout<<endl;

	//g_difference = new TGraph();
	//g_variance = new TGraph();

	c=new TCanvas("c","c",10,10,700,700);

	isBatch=-1;
	sizeBatch=1;
	numBatch_tmp=0;
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
	if(i>nLayer){
		cout<<"Inner layer "<<i<<" does not exit"<<endl;
	}
	nNodes[i]=nnodes;
	cout<<"Number of nodes in inner layer "<<i<<" is "<<nNodes[i]<<endl;
}

void simpleBP::SetCut(const TString& cut, const TString& className){
	if(className=="Signal"){
		cut_s = cut;
		cout<<"Cut on signal sample: "<<cut_s<<endl;
	}
	else if(className=="Background"){
		cut_b = cut;
		cout<<"Cut on background sample: "<<cut_b<<endl;
	}
	else{
		cout<<"Error! Cut init failed!"<<endl;
	}
}

void simpleBP::InitParameters(){
	nNodes[0]=var_input->size();
	cout<<endl;

	if(TrainingMode=="BP")
		evt_weight_max=TMath::Max(t_s_cut->GetMaximum(weightExpression.Data()), t_b_cut->GetMaximum(weightExpression.Data()));
	if(TrainingMode=="regress")
		evt_weight_max=t_s_cut->GetMaximum(weightExpression.Data());
	if(isBatch==-1)
		sizeBatch=nentries_s_train;
	else if(isBatch==0)
		sizeBatch=1;
	else if(isBatch>0)
		sizeBatch=isBatch;
	else{
		cout<<"Error! isBatch="<<isBatch;
		sizeBatch=nentries_s_train;
	}

	cout<<"Init weights:"<<endl;
	gRandom = new TRandom3();
	gRandom->SetSeed(0);
	int i_color=1;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			vector<double>* tmp_old = new vector<double>;
			vector<double>* tmp_init = new vector<double>;
			vector<double>* tmp_memory = new vector<double>;
			vector<TGraph*>* tmp_graph = new vector<TGraph*>;
			weight[i]->push_back(tmp);
			weight_old[i]->push_back(tmp_old);
			weight_init[i]->push_back(tmp_init);
			weight_memory[i]->push_back(tmp_memory);
			g_weight[i]->push_back(tmp_graph);
			delta[i]->push_back(0);
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0)
					weight[i]->at(j)->push_back(0);
				else{
					weight[i]->at(j)->push_back(2*gRandom->Rndm()-1);
					cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";
				}
				weight_old[i]->at(j)->push_back(weight[i]->at(j)->at(k));
				weight_init[i]->at(j)->push_back(weight[i]->at(j)->at(k));
				weight_memory[i]->at(j)->push_back(weight[i]->at(j)->at(k));

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

	c_s = new TChain(treeName.Data(),"c_s");
	c_s->Add(Filename.Data());
	b_s = new base();
	b_s_test = new base();
	t_s_cut = c_s->CopyTree(cut_s.Data());
	nentries_s = t_s_cut->GetEntries();
	cout<<"signal: "<<nentries_s<<endl;

	//if(nentries_s==nentries_b){
	//	nentries=nentries_s;
	//	cout<<"total: "<<nentries<<endl;
	//}
	//else{
	//	cout<<"Attention! sig & bkg have different entries!"<<endl;
	//	nentries=TMath::Min(nentries_s,nentries_b);
	//	cout<<"signal: "<<nentries_s<<" background: "<<nentries_b<<endl;
	//	//cout<<"total: "<<nentries<<endl;
	//}

	train_s = new TTree();
	train_s = t_s_cut->CloneTree(0);
	for (Long64_t jentry=0; jentry<nentries_s;jentry=jentry+2.){
		t_s_cut->GetEntry(jentry);
		train_s->Fill();
		nentries_s_train++;
	}
	test_s = new TTree();
	test_s = t_s_cut->CloneTree(0);
	for (Long64_t jentry=1; jentry<nentries_s;jentry=jentry+2.){
		t_s_cut->GetEntry(jentry);
		test_s->Fill();
		nentries_s_test++;
	}

	b_s->Init(train_s);
	b_s_test->Init(test_s);

	cout<<"signal training: "<<nentries_s_train<<" signal test: "<<nentries_s_test<<endl;

	if(TrainingMode=="BP"){
		c_b = new TChain(treeName.Data(),"c_b");
		for(int i=0; i<nOther; i++){
			c_b->Add(Filename_other[i].Data());
		}
		b_b = new base();
		b_b_test = new base();
		t_b_cut = c_b->CopyTree(cut_b.Data());
		nentries_b = t_b_cut->GetEntries();
		cout<<"background: "<<nentries_b<<endl;

		train_b = new TTree();	
		train_b = t_b_cut->CloneTree(0);
		for (Long64_t jentry=0; jentry<nentries_b;jentry=jentry+2.){
			t_b_cut->GetEntry(jentry);
			train_b->Fill();
			nentries_b_train++;
		}
		test_b = new TTree();
		test_b = c_b->CloneTree(0);
		for (Long64_t jentry=1; jentry<nentries_b;jentry=jentry+2.){
			t_b_cut->GetEntry(jentry);
			test_b->Fill();
			nentries_b_test++;
		}

		b_b->Init(train_b);
		b_b_test->Init(test_b);
		cout<<"background training: "<<nentries_b_train<<" background test: "<<nentries_b_test<<endl;
	}

	//cout<<"signal training: "<<nentries_s_train<<" background training: "<<nentries_b_train<<endl;
	//cout<<"signal test: "<<nentries_s_test<<" background test: "<<nentries_b_test<<endl;
}

void simpleBP::doTraining(){
	cout<<endl;
	cout<<"start training ..."<<endl;

	clock_t start, finish;
	start = clock();

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
	if(TrainingMode=="regress")
		nEvents = nentries_s_train;
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
			//cout<<"jentry_tmp="<<jentry_tmp<<endl;
			if(jentry_tmp<nentries_s_train){
				isSig=1;
				b_s->GetEntry(jentry_tmp);
				calculate(b_s);  
				evt_weight=b_s->weight;
			}
			else if(jentry_tmp>=nentries_s_train && jentry_tmp<nEvents){
				//isSig=-1;
				isSig=0; // test for TMVA
				b_b->GetEntry(jentry_tmp-nentries_s_train);
				calculate(b_b);  
				evt_weight=b_b->weight;
			}

			sum_weight=sum_weight+evt_weight;
			//cout<<"discriminant="<<discriminant<<endl;

			difference=0;
			if(TrainingMode=="BP")
				difference=(discriminant-isSig)*(discriminant-isSig)/2.;
			if(TrainingMode=="regress"){
				for(int i_num_output=0; i_num_output<nNodes[nLayer+1]; i_num_output++){
					difference=difference+(o[nLayer+1]->at(i_num_output+1)-var_kin_ttz->at(i_num_output))*(o[nLayer+1]->at(i_num_output+1)-var_kin_ttz->at(i_num_output))/2/nNodes[nLayer+1];
				}
			}
			//difference=(discriminant-isSig)*(discriminant-isSig); // test for TMVA
			difference=difference*evt_weight;
			//if(difference>2)cout<<"jentry_tmp="<<jentry_tmp<<" discriminant="<<discriminant<<" isSig="<<isSig<<endl;
			//cout<<"difference="<<difference<<endl;
			if(difference>difference_largest)difference_largest=difference;
			variance=variance+difference;
			if(difference>threshold){
				back_propogation(isSig); 
			}

			if(isPrintEvolution){
				for(int i=1;i<=nLayer+1;i++){
					for(int j=1;j<=nNodes[i];j++){
						for(int k=0;k<=nNodes[i-1];k++){
							g_weight[i]->at(j)->at(k)->SetPoint(iPoint,iPoint,weight[i]->at(j)->at(k));
						}
					}
				}
				//g_difference->SetPoint(iPoint,iPoint,difference);
				iPoint++;
			}
		}
		//cout<<"largest diff:"<<difference_largest<<endl;
		//cout<<"sum_weight:"<<sum_weight<<endl;
		variance=variance/sum_weight;
		//cout<<"variance="<<variance<<endl;
		//g_variance->SetPoint(N_loop-1, N_loop, variance);
		if(N_loop%TestRate==0){
			trainE = CalculateEstimator( "Training", N_loop );
			testE  = CalculateEstimator( "Testing",  N_loop );
			g_variance_HistTrain->Fill( N_loop, trainE );
			g_variance_HistTest->Fill( N_loop, testE );
		}
		isStop=false;
		//bool isStop1=TMath::Abs((variance-variance_old)/variance_old)<0.001;
		//bool isStop2=TMath::Abs((difference_largest-difference_largest_old)/difference_largest_old)<0.001;
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
	if(isPrintEvolution){
		for(int i=1;i<=nLayer+1;i++){
			//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
			for(int j=1;j<=nNodes[i];j++){
				for(int k=0;k<=nNodes[i-1];k++){
					c->Clear();
					g_weight[i]->at(j)->at(k)->GetYaxis()->SetTitleOffset(1.5);
					g_weight[i]->at(j)->at(k)->Draw("AP");
					//ltx->DrawLatex(0.15,0.85,Form("#splitline{nEntries=%d, nLoop=%d}{w_{%d%d}=%f}", nentries_s_train+nentries_b_train, N_loop, j, k, weight[i]->at(j)->at(k)));
					ltx->DrawLatex(0.15,0.85,Form("#splitline{nEntries=%d, nLoop=%d}{w_{%d%d}=%f}", nEvents, N_loop, j, k, weight[i]->at(j)->at(k)));
					c->SaveAs(Form("layer%dtolayer%d_w%d%d.png", i-1, i, j, k));
				}
				//cout<<endl;
			}
		}
	}

	//if(isPrintEvolution){
	//	c->Clear();
	//	TPad * p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
	//	p1->SetFillColor(0);
	//	p1->Draw();
	//	p1->SetLeftMargin(0.14);
	//	p1->cd();
	//	g_difference->SetTitle(";Event;(#hat{y}-y)^{2}/2");
	//	g_difference->GetYaxis()->SetTitleOffset(2);
	//	g_difference->Draw("AP");
	//	ltx->DrawLatex(0.15,0.85,Form("nEntries=%d, nLoop=%d", nentries_s_train+nentries_b_train, N_loop));
	//	c->SaveAs("difference.png");
	//}

	//c->Clear();
	TPad * p2 = new TPad("p2","p2", 0.00,0.00,1.00,0.97);
	//p2->SetFillColor(0);
	//p2->Draw();
	//p2->SetLeftMargin(0.14);
	//p2->cd();
	//g_variance->SetTitle(";Epoch;#Sigma (w*(#hat{y}-y))^{2}/2/Sumw");
	//g_variance->GetYaxis()->SetTitleOffset(2);
	//g_variance->Draw("AP*");
	//ltx->DrawLatex(0.15,0.85,Form("nEntries=%d, nLoop=%d", nentries_s_train+nentries_b_train, N_loop));
	//c->SaveAs("variance.png");

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
	g_variance_HistTrain->Draw("HIST");
	g_variance_HistTest->Draw("HISTSame");
	TLegend* legend = new TLegend(0.7,0.75,0.9,0.9,"");
	legend->SetFillColor(kWhite);
	legend->AddEntry("estimatorHistTrain", "Training Sample", "l");
	legend->AddEntry("estimatorHistTest", "Test Sample", "l");
	legend->Draw();
	c->SaveAs("annconvergencetest.png");
	g_variance_HistTrain->Write();
	g_variance_HistTest->Write();

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;
	//cout<<"CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
	cout<<"duration="<<duration<<"sec"<<endl;
}

void simpleBP::calculate(base * b){
	var->clear();
	net[0]->clear();
	o[0]->clear();

	net[0]->push_back(1);
	o[0]->push_back(1);

	var_input_NN_beforeNorm->clear();
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet1_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet1_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet1_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet2_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet2_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_Bjet2_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw1_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw1_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw1_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw2_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw2_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_JetClosestMw2_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton1_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton1_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton1_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton2_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton2_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton2_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton3_P4->E());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton3_P4->Theta());
	var_input_NN_beforeNorm->push_back(b->multilepton_Lepton3_P4->Phi());
	var_input_NN_beforeNorm->push_back(b->multilepton_mET->Pt());
	var_input_NN_beforeNorm->push_back(b->multilepton_mET->Phi());
	//cout<<"var_input_NN_beforeNorm->size()="<<var_input_NN_beforeNorm->size()<<endl;

	for(unsigned int ivar=0;ivar<var_input->size();ivar++){
		double var_tmp=var_input_NN_beforeNorm->at(ivar);
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
			//o_tmp=(TMath::Exp(net_tmp)-TMath::Exp(-net_tmp))/(TMath::Exp(net_tmp)+TMath::Exp(-net_tmp));
			o_tmp=activation_function(net_tmp, NeuronType);
			net[i]->push_back(net_tmp);
			o[i]->push_back(o_tmp);
		}
	}

	netDis=net[nLayer+1]->at(1);
	discriminant=o[nLayer+1]->at(1);
	//cout<<"discriminant="<<discriminant<<endl;

	var_kin_ttz->clear();
	for(int i_num_output=0; i_num_output<nNodes[nLayer+1]; i_num_output++){
		double var_kin_ttz_tmp = 2 * (b->mc_kin_ttz_inputvars->at(i_num_output)-xL->at(i_num_output)) / (xU->at(i_num_output)-xL->at(i_num_output)) -1;
		var_kin_ttz->push_back(var_kin_ttz_tmp);
		//cout<<"var_kin_ttz->at(i_num_output)="<<var_kin_ttz->at(i_num_output)<<endl;
	}
}

void simpleBP::back_propogation(int isSig){
	numBatch_tmp++;

	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			weight_old[i]->at(j)->clear();
			for(int k=0;k<=nNodes[i-1];k++){
				//double tmp = weight[i]->at(j)->at(k);
				double tmp = weight_memory[i]->at(j)->at(k); //			--Jing Li @ 2017-08-01 
				weight_old[i]->at(j)->push_back(tmp);
				//cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" "<<weight_old[i]->at(j)->at(k)<<" ";
			}
			//cout<<endl;
		}
	}

	for(int i=1;i<=nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			//weight[i]->at(j)->clear();
			weight_memory[i]->at(j)->clear(); //		--Jing Li @ 2017-08-01 
		}
	}

	if(TrainingMode=="BP"){
		//deltaDis=-(isSig-discriminant)*(1+discriminant)*(1-discriminant);
		deltaDis=-(isSig-discriminant)*activation_function_dev(netDis, discriminant,NeuronType);
		//cout<<"nLayer+1="<<nLayer+1<<endl;
		delta[nLayer+1]->clear();
		delta[nLayer+1]->push_back(0);
		delta[nLayer+1]->push_back(deltaDis);
	}
	if(TrainingMode=="regress"){
		delta[nLayer+1]->clear();
		delta[nLayer+1]->push_back(0);
		for(int i_num_output=0; i_num_output<nNodes[nLayer+1]; i_num_output++){
			double deltaDis_tmp = -(var_kin_ttz->at(i_num_output)-o[nLayer+1]->at(i_num_output+1))*activation_function_dev(net[nLayer+1]->at(i_num_output+1), o[nLayer+1]->at(i_num_output+1),NeuronType);
			delta[nLayer+1]->push_back(deltaDis_tmp);
			//cout<<"delta[nLayer+1]->at(i_num_output+1)"<<delta[nLayer+1]->at(i_num_output+1)<<endl;
		}
	}

	for(int i=nLayer; i>=1; i--){
		//cout<<"Inner layer "<<i+1<<" to "<<i<<endl;
		delta[i]->clear();
		for(int j=0;j<=nNodes[i];j++){
			double delta_tmp=0;
			for(int k=1;k<=nNodes[i+1];k++){
				//double tmp_devNetj = ( 1 + o[i]->at(j) ) * ( 1 - o[i]->at(j) );
				double tmp_devNetj = activation_function_dev(net[i]->at(j), o[i]->at(j), NeuronType);
				//double tmp_wkj = weight_old[i+1]->at(k)->at(j);
				double tmp_wkj = weight[i+1]->at(k)->at(j); //			--Jing Li @ 2017-08-01 
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
				//double weight_tmp_change =  eta * delta[i]->at(j) * o[i-1]->at(k) * evt_weight;
				double weight_tmp_change =  eta * delta[i]->at(j) * o[i-1]->at(k) * evt_weight / evt_weight_max; // test			--Jing Li @ 2017-03-02 
				weight_tmp_change = weight_tmp_change / ((double)sizeBatch); //		--Jing Li @ 2017-08-01 
				double weight_tmp = weight_old[i]->at(j)->at(k) - weight_tmp_change;
				//weight[i]->at(j)->push_back(weight_tmp);
				weight_memory[i]->at(j)->push_back(weight_tmp); //		--Jing Li @ 2017-08-01 
				//cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" "<<weight_old[i]->at(j)->at(k)<<" ";
			}
			//cout<<endl;
		}
	}

	//			--Jing Li @ 2017-08-01 
	if(numBatch_tmp%sizeBatch == 0 ){
		for(int i=1;i<=nLayer+1;i++){
			for(int j=0;j<=nNodes[i];j++){
				weight[i]->at(j)->clear();
				for(int k=0;k<=nNodes[i-1];k++){
					double tmp = weight_memory[i]->at(j)->at(k); //       --Jing Li @ 2017-08-01 
					weight[i]->at(j)->push_back(tmp);
					//cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" "<<weight_memory[i]->at(j)->at(k)<<" ";
				}
			}
		}
	}
}

Double_t simpleBP::CalculateEstimator( TString treeType, Int_t iEpoch ){
	// calculate the estimator that training is attempting to minimize
	//cout<<"calculate the estimator that training is attempting to minimize ..."<<endl;

	double estimator = 0.;

	// sanity check
	if (treeType!="Training" && treeType!="Testing") {
		cout<<"<CalculateEstimator> fatal error: wrong tree type: "<<treeType << endl;
	}
	//base * b_s_tmp;
	//base * b_b_tmp;
	Int_t nEvents;
	Int_t nentries_s_tmp, nentries_b_tmp;
	if(treeType=="Training"){
		//b_s_tmp=b_s;
		//b_b_tmp=b_b;
		nentries_s_tmp=nentries_s_train;
		nentries_b_tmp=nentries_b_train;
		nEvents=nentries_s_train+nentries_b_train;
	}
	if(treeType=="Testing"){
		//b_s_tmp=b_s_test;
		//b_b_tmp=b_b_test;
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
			if(treeType=="Training"){
				b_s->GetEntry(jentry);
				calculate(b_s);
				evt_weight_tmp=b_s->weight;
			}
			if(treeType=="Testing"){
				b_s_test->GetEntry(jentry);
				calculate(b_s_test);
				evt_weight_tmp=b_s_test->weight;
			}
		}
		else{
			//isSig=-1;
			isSig=0; // test for TMVA
			if(treeType=="Training"){
				b_b->GetEntry(jentry);
				calculate(b_b);
				evt_weight_tmp=b_b->weight;
			}
			if(treeType=="Testing"){
				b_b_test->GetEntry(jentry);
				calculate(b_b_test);
				evt_weight_tmp=b_b_test->weight;
			}
		}
		sum_weight_tmp=sum_weight_tmp+evt_weight_tmp;
		//difference_tmp=(discriminant-isSig)*(discriminant-isSig)/2.;
		//difference_tmp=(discriminant-isSig)*(discriminant-isSig); // test for TMVA
		difference_tmp=0;
		if(TrainingMode=="BP")
			difference_tmp=(discriminant-isSig)*(discriminant-isSig)/2.;
		if(TrainingMode=="regress"){
			for(int i_num_output=0; i_num_output<nNodes[nLayer+1]; i_num_output++){
				difference_tmp=difference_tmp+(o[nLayer+1]->at(i_num_output+1)-var_kin_ttz->at(i_num_output))*(o[nLayer+1]->at(i_num_output+1)-var_kin_ttz->at(i_num_output))/2/nNodes[nLayer+1];
			}
		}
		difference_tmp=difference_tmp*evt_weight_tmp;
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
	fout->cd();
	tout_train_s = train_s->CloneTree(0);
	tout_train_s->Branch("discriminant",&discriminant,"discriminant/D");
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_train_s->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for (Long64_t jentry=0; jentry<nentries_s_train;jentry++){
		train_s->GetEntry(jentry);
		calculate(b_s);
		tout_train_s->Fill();
	}

	tout_train_b = new TTree();	
	tout_train_b = train_b->CloneTree(0);
	tout_train_b->Branch("discriminant",&discriminant,"discriminant/D");	
	for(int i=0;i<nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_train_b->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for (Long64_t jentry=0; jentry<nentries_b_train;jentry++){
		train_b->GetEntry(jentry);
		calculate(b_b);
		tout_train_b->Fill();
	}

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

void simpleBP::Eval_regress(){
	fout->cd();
	tout_train_s = train_s->CloneTree(0);
	tout_train_s->Branch("discriminant",&discriminant,"discriminant/D");
	for(int i=0;i<=nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_train_s->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for(int i=0; i<nNodes[nLayer+1]; i++){
		tout_train_s->Branch(Form("var_kin_ttz_%d",i), &(var_kin_ttz->at(i)), Form("var_kin_ttz_%d/D",i) );
	}
	for (Long64_t jentry=0; jentry<nentries_s_train;jentry++){
		train_s->GetEntry(jentry);
		calculate(b_s);
		tout_train_s->Fill();
	}

	tout_test_s = new TTree();
	tout_test_s = test_s->CloneTree(0);
	tout_test_s->Branch("discriminant",&discriminant,"discriminant/D"); 
	for(int i=0;i<=nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout_test_s->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}
	for(int i=0; i<nNodes[nLayer+1]; i++){
		tout_test_s->Branch(Form("var_kin_ttz_%d",i), &(var_kin_ttz->at(i)), Form("var_kin_ttz_%d/D",i) );
	}
	for (Long64_t jentry=0; jentry<nentries_s_test;jentry++){
		test_s->GetEntry(jentry);
		calculate(b_s_test);
		tout_test_s->Fill();
	}

	tout_train_s->Write("train_s");
	tout_test_s->Write("test_s");
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
	TH1D * h_s_train = GetHistoWeight(tout_train_s, "discriminant", 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_train");
	TH1D * h_b_train = GetHistoWeight(tout_train_b, "discriminant", 50, -1., 1., cut_b.Data(), weightExpression.Data(), "h_b_train");
	TH1D * h_s_test  = GetHistoWeight(tout_test_s,  "discriminant", 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_test");
	TH1D * h_b_test  = GetHistoWeight(tout_test_b,  "discriminant", 50, -1., 1., cut_b.Data(), weightExpression.Data(), "h_b_test");

	h_s_test->SetLineColor(2);
	h_b_test->SetLineColor(1);
	h_s_test->SetMaximum((1.5*TMath::Max(h_s_test->GetMaximum(), h_b_test->GetMaximum())));
	//h_s_test->SetTitle(";Discriminant; Event / bin");
	h_s_test->SetTitle(";Discriminant; Unit Area");
	h_s_test->GetYaxis()->SetTitleOffset(1.8);
	h_s_test->Sumw2();
	h_b_test->Sumw2();
	h_s_test->DrawNormalized("HIST");
	h_b_test->DrawNormalized("HISTsame");
	h_s_train->SetMarkerColor(2);
	h_s_train->SetMarkerStyle(20);
	h_b_train->SetMarkerColor(1);
	h_b_train->SetMarkerStyle(20);
	h_s_train->Sumw2();
	h_b_train->Sumw2();
	h_s_train->DrawNormalized("PEsame");
	h_b_train->DrawNormalized("PEsame");

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

	//for(int i=1;i<=nLayer+1;i++){
	//	y_ltx=0.7;
	//	ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
	//	y_ltx=y_ltx-0.04;
	//	for(int j=1;j<=nNodes[i];j++){
	//		for(int k=0;k<=nNodes[i-1];k++){
	//			ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
	//			y_ltx=y_ltx-0.04;
	//		}
	//		cout<<endl;
	//	}
	//	x_ltx=x_ltx+0.2;
	//}

	c->Print("discriminant.png");
	h_s_test->Write("discriminant_s_test");
	h_b_test->Write("discriminant_b_test");
	h_s_train->Write("discriminant_s_train");
	h_b_train->Write("discriminant_b_train");

	for(int i=0;i<nLayer+1;i++){
		for(int j=1;j<=nNodes[i];j++){
			c->Clear();
			p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
			p1->SetFillColor(0);
			p1->Draw();
			p1->SetLeftMargin(0.14);
			p1->cd();
			h_s_train = GetHistoWeight(tout_train_s, Form("l%dn%d",i,j), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_train");
			h_b_train = GetHistoWeight(tout_train_b, Form("l%dn%d",i,j), 50, -1., 1., cut_b.Data(), weightExpression.Data(), "h_b_train");
			h_s_test  = GetHistoWeight(tout_test_s,  Form("l%dn%d",i,j), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_test");
			h_b_test  = GetHistoWeight(tout_test_b,  Form("l%dn%d",i,j), 50, -1., 1., cut_b.Data(), weightExpression.Data(), "h_b_test");

			h_s_test->SetLineColor(2);
			h_b_test->SetLineColor(1);
			h_s_test->SetMaximum((1.5*TMath::Max(h_s_test->GetMaximum(),h_b_test->GetMaximum())));
			if(i==0)
				//h_s_test->SetTitle(Form(";%s (norm);Event / bin",var_input->at(j-1).Data()));
				h_s_test->SetTitle(Form(";%s (norm);Unit Area",var_input->at(j-1).Data()));
			else
				//h_s_test->SetTitle(Form(";Output of layer %d node %d;Event / bin",i,j));
				h_s_test->SetTitle(Form(";Output of layer %d node %d;Unit Area",i,j));
			h_s_test->GetYaxis()->SetTitleOffset(1.8);
			h_s_test->Sumw2();
			h_b_test->Sumw2();
			h_s_test->DrawNormalized("HIST");
			h_b_test->DrawNormalized("HISTsame");
			h_s_train->SetMarkerColor(2);
			h_s_train->SetMarkerStyle(20);
			h_b_train->SetMarkerColor(1);
			h_b_train->SetMarkerStyle(20);
			h_s_train->Sumw2();
			h_b_train->Sumw2();
			h_s_train->DrawNormalized("PEsame");
			h_b_train->DrawNormalized("PEsame");

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

void simpleBP::plotHist_regress(){

	c->cd();
	gStyle->SetOptStat(0);

	TH1D * h_s_train = new TH1D();
	TH1D * h_b_train = new TH1D();
	TH1D * h_s_test  = new TH1D();
	TH1D * h_b_test  = new TH1D();

	TLegend* legend_test = new TLegend(0.15,0.75,0.5,0.9,"");
	TLegend* legend_train = new TLegend(0.5,0.75,0.85,0.9,"");

	TPad * p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);

	int i=nLayer+1;
	for(int j=1;j<=nNodes[i];j++){
		c->Clear();
		p1 = new TPad("p1","p1", 0.00,0.00,1.00,0.97);
		p1->SetFillColor(0);
		p1->Draw();
		p1->SetLeftMargin(0.14);
		p1->cd();
		h_s_train = GetHistoWeight(tout_train_s, Form("l%dn%d",i,j), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_train");
		h_b_train = GetHistoWeight(tout_train_s, Form("var_kin_ttz_%d",j-1), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_b_train");
		h_s_test  = GetHistoWeight(tout_test_s,  Form("l%dn%d",i,j), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_s_test");
		h_b_test  = GetHistoWeight(tout_test_s,  Form("var_kin_ttz_%d",j-1), 50, -1., 1., cut_s.Data(), weightExpression.Data(), "h_b_test");

		h_s_test->SetLineColor(2);
		h_b_test->SetLineColor(1);
		h_s_test->SetMaximum((1.5*TMath::Max(h_s_test->GetMaximum(),h_b_test->GetMaximum())));
		if(i==0)
			//h_s_test->SetTitle(Form(";%s (norm);Event / bin",var_input->at(j-1).Data()));
			h_s_test->SetTitle(Form(";%s (norm);Unit Area",var_input->at(j-1).Data()));
		else
			//h_s_test->SetTitle(Form(";Output of layer %d node %d;Event / bin",i,j));
			h_s_test->SetTitle(Form(";Output of layer %d node %d;Unit Area",i,j));
		h_s_test->GetYaxis()->SetTitleOffset(1.8);
		h_s_test->Sumw2();
		h_b_test->Sumw2();
		h_s_test->DrawNormalized("HIST");
		h_b_test->DrawNormalized("HISTsame");
		h_s_train->SetMarkerColor(2);
		h_s_train->SetMarkerStyle(20);
		h_b_train->SetMarkerColor(1);
		h_b_train->SetMarkerStyle(20);
		h_s_train->Sumw2();
		h_b_train->Sumw2();
		h_s_train->DrawNormalized("PEsame");
		h_b_train->DrawNormalized("PEsame");

		legend_test->Clear();
		legend_test->SetFillColor(kWhite);
		legend_test->AddEntry("h_s_test", "NN output (test sample)", "l");
		legend_test->AddEntry("h_b_test", "mc_kin_ttz_inputvars_norm (test sample)", "l");
		legend_test->Draw();

		legend_train->Clear();
		legend_train->SetFillColor(kWhite);
		legend_train->AddEntry("h_s_train", "NN output (training sample)", "p");
		legend_train->AddEntry("h_b_train", "mc_kin_ttz_inputvars_norm (training sample)", "p");
		legend_train->Draw();

		c->Print(Form("O_l%dn%d.png",i,j));
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
	TGraph * g_ROC_train=GetEffSvsEffB(tout_train_s, tout_train_b, cut_s.Data(), cut_b.Data(), "discriminant", -1, 1, weightExpression.Data(), 50, "g_ROC_train");
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
	//for(int i=1;i<=nLayer+1;i++){
	//	y_ltx=0.5;
	//	ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
	//	y_ltx=y_ltx-0.04;
	//	for(int j=1;j<=nNodes[i];j++){
	//		for(int k=0;k<=nNodes[i-1];k++){
	//			ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
	//			y_ltx=y_ltx-0.04;
	//		}
	//		cout<<endl;
	//	}
	//	x_ltx=x_ltx+0.2;
	//}
	c->SaveAs("ROC_train.png");
	g_ROC_train->Write();

	c->Clear();
	TPad * p2 = new TPad("p2","p2", 0.00,0.00,1.00,0.97);
	p2->SetFillColor(0);
	p2->Draw();
	p2->SetTicks(1,1);
	p2->SetGridx();
	p2->SetGridy();
	p2->cd();
	TGraph * g_ROC_test=GetEffSvsEffB(tout_test_s, tout_test_b, cut_s.Data(), cut_b.Data(), "discriminant", -1, 1, weightExpression.Data(), 50, "g_ROC_test");
	hGrid->SetTitle("test sample");
	hGrid->Draw();
	g_ROC_test->Draw("*same");
	x_ltx=0.5, y_ltx=0.5;
	//for(int i=1;i<=nLayer+1;i++){
	//	y_ltx=0.5;
	//	ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
	//	y_ltx=y_ltx-0.04;
	//	for(int j=1;j<=nNodes[i];j++){
	//		for(int k=0;k<=nNodes[i-1];k++){
	//			ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
	//			y_ltx=y_ltx-0.04;
	//		}
	//		cout<<endl;
	//	}
	//	x_ltx=x_ltx+0.2;
	//}
	c->SaveAs("ROC_test.png");
	g_ROC_test->Write();
}
/*
	void simpleBP::plotWeights(){

	TCanvas * c_weights=new TCanvas("c_weights","c_weights",10,10,700,700*(nLayer+1));
	c_weights->cd();

	double x_ltx=0.1, y_ltx=0.9;
	double dx=0.2;
	double dy=0.04/double(nLayer+1);

	TLatex * ltx = new TLatex();
	ltx->SetNDC(kTRUE);
	ltx->SetTextFont(22);
	ltx->SetTextSize(0.03);

	ltx->DrawLatex(x_ltx,y_ltx,Form("#splitline{nEntries_sig=%d, nEntries_bkg=%d, nEntries=%d}{nEpochs=%d}",nentries_s_train, nentries_b_train, nentries_s_train+nentries_b_train, nEpochs));
	y_ltx = y_ltx-2*dy;

	for(int i=1;i<=nLayer+1;i++){
	ltx->DrawLatex(x_ltx,y_ltx,Form("Layer %d to %d", i-1, i));
	y_ltx=y_ltx-dy;
	for(int j=1;j<=nNodes[i];j++){
	for(int k=0;k<=nNodes[i-1];k++){
	ltx->DrawLatex(x_ltx,y_ltx,Form("w_{%d%d}=%f", j, k, weight[i]->at(j)->at(k)));
	y_ltx=y_ltx-dy;
	}
	if(j!=nNodes[i] && (j%4)!=0){
	x_ltx=x_ltx+dx;
	y_ltx=y_ltx+dy*(nNodes[i-1]+1);
	}
	else{
	x_ltx=0.1;
	y_ltx=y_ltx-dy;
	}
	}
	}

	ltx->DrawLatex(x_ltx,y_ltx,Form("duration = %.2f sec", duration));

	c_weights->SaveAs("myWeights.png");
	c_weights->SaveAs("myWeights.pdf");
	}
	*/
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

void simpleBP::myana(){
	InitTrees();
	InitParameters();
	doTraining();
	if(TrainingMode=="BP"){
		Eval();
		plotHist();
		plotROC();
	}
	if(TrainingMode=="regress"){
		Eval_regress();
		plotHist_regress();
	}
	// plotWeights();
}
#endif
