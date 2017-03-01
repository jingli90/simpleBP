#ifdef simpleBP_cxx

#include "TMath.h"

double myTANH(double x){
	double y = (TMath::Exp(x)-TMath::Exp(-x))/(TMath::Exp(x)+TMath::Exp(-x));
	return y;
}

double myTANH_DEV(double y){
	double dev=(1+y)*(1-y);
	return dev;
}

double mySIGMOID(double x){
	double y = 1/(1+TMath::Exp(-x));
	return y;
}

double mySIGMOID_DEV(double y){
	double dev=y*(1-y);
	return dev;
}

double activation_function(double x, TString type){
	double y=-999.;
	if(type=="tanh") y=myTANH(x);
	else if (type=="sigmoid") y=mySIGMOID(x);
	else cout<<"Error! Invalidated neuron type!"<<endl;
	return y;
}

double activation_function_dev(double y, TString type){
	double dev=-999.;
	if(type=="tanh") dev=myTANH_DEV(y);
	else if (type=="sigmoid") dev=mySIGMOID_DEV(y);
	else cout<<"Error! Invalidated neuron type!"<<endl;
	return dev;
}

#endif
