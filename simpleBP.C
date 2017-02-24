#define simpleBP_cxx
#include "simpleBP.h"

void simpleBP::myana(){
	InitParameters();
	InitTrees();
	doTraining();
	Eval();
	plotHist();
	plotROC();
}
