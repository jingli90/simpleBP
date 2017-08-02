//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 26 11:44:47 2017 by ROOT version 6.06/01
// from TTree Tree/Tree
// found on file: ../smallTree/small_ttZ.root
//////////////////////////////////////////////////////////

#ifndef base_h
#define base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "vector"

class base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         weight;
   Int_t           catJets;
   Char_t          is_2lss_TTH_SR;
   Char_t          is_3l_TTH_SR;
   Int_t           mc_ttZhypAllowed;
   TLorentzVector  *multilepton_Bjet1_P4;
   TLorentzVector  *multilepton_Bjet2_P4;
   TLorentzVector  *multilepton_Lepton1_P4;
   TLorentzVector  *multilepton_Lepton2_P4;
   TLorentzVector  *multilepton_Lepton3_P4;
   TLorentzVector  *multilepton_JetClosestMw1_P4;
   TLorentzVector  *multilepton_JetClosestMw2_P4;
   TLorentzVector  *multilepton_mET;
   vector<double>  *mc_kin_tthfl_inputvars;
   vector<double>  *mc_kin_tthsl_inputvars;
   vector<double>  *mc_kin_ttz_inputvars;
   vector<double>  *mc_kin_ttw_inputvars;
   vector<double>  *mc_kin_ttwjj_inputvars;
   vector<double>  *mc_kin_ttbarfl_inputvars;
   vector<double>  *mc_kin_ttbarsl_inputvars;
   vector<double>  *mc_kin_tllj_inputvars;
   vector<double>  *mc_kin_wzjj_inputvars;
   vector<double>  *mc_kin_thj_inputvars;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_catJets;   //!
   TBranch        *b_is_2lss_TTH_SR;   //!
   TBranch        *b_is_3l_TTH_SR;   //!
   TBranch        *b_mc_ttZhypAllowed;   //!
   TBranch        *b_multilepton_Bjet1_P4;   //!
   TBranch        *b_multilepton_Bjet2_P4;   //!
   TBranch        *b_multilepton_Lepton1_P4;   //!
   TBranch        *b_multilepton_Lepton2_P4;   //!
   TBranch        *b_multilepton_Lepton3_P4;   //!
   TBranch        *b_multilepton_JetClosestMw1_P4;   //!
   TBranch        *b_multilepton_JetClosestMw2_P4;   //!
   TBranch        *b_multilepton_mET;   //!
   TBranch        *b_mc_kin_tthfl_inputvars;   //!
   TBranch        *b_mc_kin_tthsl_inputvars;   //!
   TBranch        *b_mc_kin_ttz_inputvars;   //!
   TBranch        *b_mc_kin_ttw_inputvars;   //!
   TBranch        *b_mc_kin_ttwjj_inputvars;   //!
   TBranch        *b_mc_kin_ttbarfl_inputvars;   //!
   TBranch        *b_mc_kin_ttbarsl_inputvars;   //!
   TBranch        *b_mc_kin_tllj_inputvars;   //!
   TBranch        *b_mc_kin_wzjj_inputvars;   //!
   TBranch        *b_mc_kin_thj_inputvars;   //!

   base(TTree *tree=0);
   virtual ~base();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};


base::base(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../smallTree/small_ttZ.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../smallTree/small_ttZ.root");
      }
      f->GetObject("Tree",tree);

   }
   Init(tree);
}

base::~base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void base::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   multilepton_Bjet1_P4 = 0;
   multilepton_Bjet2_P4 = 0;
   multilepton_Lepton1_P4 = 0;
   multilepton_Lepton2_P4 = 0;
   multilepton_Lepton3_P4 = 0;
   multilepton_JetClosestMw1_P4 = 0;
   multilepton_JetClosestMw2_P4 = 0;
   multilepton_mET = 0;
   mc_kin_tthfl_inputvars = 0;
   mc_kin_tthsl_inputvars = 0;
   mc_kin_ttz_inputvars = 0;
   mc_kin_ttw_inputvars = 0;
   mc_kin_ttwjj_inputvars = 0;
   mc_kin_ttbarfl_inputvars = 0;
   mc_kin_ttbarsl_inputvars = 0;
   mc_kin_tllj_inputvars = 0;
   mc_kin_wzjj_inputvars = 0;
   mc_kin_thj_inputvars = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("catJets", &catJets, &b_catJets);
   fChain->SetBranchAddress("is_2lss_TTH_SR", &is_2lss_TTH_SR, &b_is_2lss_TTH_SR);
   fChain->SetBranchAddress("is_3l_TTH_SR", &is_3l_TTH_SR, &b_is_3l_TTH_SR);
   fChain->SetBranchAddress("mc_ttZhypAllowed", &mc_ttZhypAllowed, &b_mc_ttZhypAllowed);
   fChain->SetBranchAddress("multilepton_Bjet1_P4", &multilepton_Bjet1_P4, &b_multilepton_Bjet1_P4);
   fChain->SetBranchAddress("multilepton_Bjet2_P4", &multilepton_Bjet2_P4, &b_multilepton_Bjet2_P4);
   fChain->SetBranchAddress("multilepton_Lepton1_P4", &multilepton_Lepton1_P4, &b_multilepton_Lepton1_P4);
   fChain->SetBranchAddress("multilepton_Lepton2_P4", &multilepton_Lepton2_P4, &b_multilepton_Lepton2_P4);
   fChain->SetBranchAddress("multilepton_Lepton3_P4", &multilepton_Lepton3_P4, &b_multilepton_Lepton3_P4);
   fChain->SetBranchAddress("multilepton_JetClosestMw1_P4", &multilepton_JetClosestMw1_P4, &b_multilepton_JetClosestMw1_P4);
   fChain->SetBranchAddress("multilepton_JetClosestMw2_P4", &multilepton_JetClosestMw2_P4, &b_multilepton_JetClosestMw2_P4);
   fChain->SetBranchAddress("multilepton_mET", &multilepton_mET, &b_multilepton_mET);
   fChain->SetBranchAddress("mc_kin_tthfl_inputvars", &mc_kin_tthfl_inputvars, &b_mc_kin_tthfl_inputvars);
   fChain->SetBranchAddress("mc_kin_tthsl_inputvars", &mc_kin_tthsl_inputvars, &b_mc_kin_tthsl_inputvars);
   fChain->SetBranchAddress("mc_kin_ttz_inputvars", &mc_kin_ttz_inputvars, &b_mc_kin_ttz_inputvars);
   fChain->SetBranchAddress("mc_kin_ttw_inputvars", &mc_kin_ttw_inputvars, &b_mc_kin_ttw_inputvars);
   fChain->SetBranchAddress("mc_kin_ttwjj_inputvars", &mc_kin_ttwjj_inputvars, &b_mc_kin_ttwjj_inputvars);
   fChain->SetBranchAddress("mc_kin_ttbarfl_inputvars", &mc_kin_ttbarfl_inputvars, &b_mc_kin_ttbarfl_inputvars);
   fChain->SetBranchAddress("mc_kin_ttbarsl_inputvars", &mc_kin_ttbarsl_inputvars, &b_mc_kin_ttbarsl_inputvars);
   fChain->SetBranchAddress("mc_kin_tllj_inputvars", &mc_kin_tllj_inputvars, &b_mc_kin_tllj_inputvars);
   fChain->SetBranchAddress("mc_kin_wzjj_inputvars", &mc_kin_wzjj_inputvars, &b_mc_kin_wzjj_inputvars);
   fChain->SetBranchAddress("mc_kin_thj_inputvars", &mc_kin_thj_inputvars, &b_mc_kin_thj_inputvars);
   Notify();
}

Bool_t base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_cxx
