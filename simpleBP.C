#define simpleBP_cxx
#include "simpleBP.h"

void simpleBP::myana(){
	InitTrees();
	InitParameters();
	doTraining();
	Eval();
	plotHist();
	plotROC();
	plotWeights();
}
