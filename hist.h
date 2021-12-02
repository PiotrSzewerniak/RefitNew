//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep 11 19:34:58 2021 by ROOT version 6.24/02
// from TTree T/Particle data
// found on file: data.root
//////////////////////////////////////////////////////////

#ifndef hist_h
#define hist_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class hist {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Energy_1;
   Float_t         Energy_2;
   Float_t         Energy_3;
   Float_t         Energy_1_inp;
   Float_t         Energy_2_inp;
   Float_t         Energy_3_inp;
   Float_t         Energy_1_sm;
   Float_t         Energy_2_sm;
   Float_t         Energy_3_sm;
   Float_t         Energy_1_error;
   Float_t         Energy_2_error;
   Float_t         Energy_3_error;
   Float_t         Theta_1;
   Float_t         Theta_2;
   Float_t         Theta_3;
   Float_t         Theta_1_inp;
   Float_t         Theta_2_inp;
   Float_t         Theta_3_inp;
   Float_t         Theta_1_sm;
   Float_t         Theta_2_sm;
   Float_t         Theta_3_sm;
   Float_t         Theta_1_error;
   Float_t         Theta_2_error;
   Float_t         Theta_3_error;
   Float_t         Phi_1;
   Float_t         Phi_2;
   Float_t         Phi_3;
   Float_t         Phi_1_inp;
   Float_t         Phi_2_inp;
   Float_t         Phi_3_inp;
   Float_t         Phi_1_sm;
   Float_t         Phi_2_sm;
   Float_t         Phi_3_sm;
   Float_t         Phi_1_error;
   Float_t         Phi_2_error;
   Float_t         Phi_3_error;
   Float_t         Chi2;
   Float_t         Probability;
   Float_t         MissMassProt2Neut;
   Float_t         MissMassProt2NeutPre;

   // List of branches
   TBranch        *b_Energy_1;   //!
   TBranch        *b_Energy_2;   //!
   TBranch        *b_Energy_3;   //!
   TBranch        *b_Energy_1_inp;   //!
   TBranch        *b_Energy_2_inp;   //!
   TBranch        *b_Energy_3_inp;   //!
   TBranch        *b_Energy_1_sm;   //!
   TBranch        *b_Energy_2_sm;   //!
   TBranch        *b_Energy_3_sm;   //!
   TBranch        *b_Energy_1_error;   //!
   TBranch        *b_Energy_2_error;   //!
   TBranch        *b_Energy_3_error;   //!
   TBranch        *b_Theta_1;   //!
   TBranch        *b_Theta_2;   //!
   TBranch        *b_Theta_3;   //!
   TBranch        *b_Theta_1_inp;   //!
   TBranch        *b_Theta_2_inp;   //!
   TBranch        *b_Theta_3_inp;   //!
   TBranch        *b_Theta_1_sm;   //!
   TBranch        *b_Theta_2_sm;   //!
   TBranch        *b_Theta_3_sm;   //!
   TBranch        *b_Theta_1_error;   //!
   TBranch        *b_Theta_2_error;   //!
   TBranch        *b_Theta_3_error;   //!
   TBranch        *b_Phi_1;   //!
   TBranch        *b_Phi_2;   //!
   TBranch        *b_Phi_3;   //!
   TBranch        *b_Phi_1_inp;   //!
   TBranch        *b_Phi_2_inp;   //!
   TBranch        *b_Phi_3_inp;   //!
   TBranch        *b_Phi_1_sm;   //!
   TBranch        *b_Phi_2_sm;   //!
   TBranch        *b_Phi_3_sm;   //!
   TBranch        *b_Phi_1_error;   //!
   TBranch        *b_Phi_2_error;   //!
   TBranch        *b_Phi_3_error;   //!
   TBranch        *b_Chi2;   //!
   TBranch        *b_Probability;   //!
   TBranch        *b_MissMassProt2Neut;   //!
   TBranch        *b_MissMassProt2NeutPre;   //!

   hist(TTree *tree=0);
   virtual ~hist();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef hist_cxx
hist::hist(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Refit.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Refit.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

hist::~hist()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hist::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hist::LoadTree(Long64_t entry)
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

void hist::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // fChain->SetBranchAddress("Energy_1", &Energy_1, &b_Energy_1);
   // fChain->SetBranchAddress("Energy_2", &Energy_2, &b_Energy_2);
   // fChain->SetBranchAddress("Energy_3", &Energy_3, &b_Energy_3);
   // fChain->SetBranchAddress("Energy_1_inp", &Energy_1_inp, &b_Energy_1_inp);
   // fChain->SetBranchAddress("Energy_2_inp", &Energy_2_inp, &b_Energy_2_inp);
   // fChain->SetBranchAddress("Energy_3_inp", &Energy_3_inp, &b_Energy_3_inp);
   // fChain->SetBranchAddress("Energy_1_sm", &Energy_1_sm, &b_Energy_1_sm);
   // fChain->SetBranchAddress("Energy_2_sm", &Energy_2_sm, &b_Energy_2_sm);
   // fChain->SetBranchAddress("Energy_3_sm", &Energy_3_sm, &b_Energy_3_sm);
   // fChain->SetBranchAddress("Energy_1_error", &Energy_1_error, &b_Energy_1_error);
   // fChain->SetBranchAddress("Energy_2_error", &Energy_2_error, &b_Energy_2_error);
   // fChain->SetBranchAddress("Energy_3_error", &Energy_3_error, &b_Energy_3_error);
   // fChain->SetBranchAddress("Theta_1", &Theta_1, &b_Theta_1);
   // fChain->SetBranchAddress("Theta_2", &Theta_2, &b_Theta_2);
   // fChain->SetBranchAddress("Theta_3", &Theta_3, &b_Theta_3);
   // fChain->SetBranchAddress("Theta_1_inp", &Theta_1_inp, &b_Theta_1_inp);
   // fChain->SetBranchAddress("Theta_2_inp", &Theta_2_inp, &b_Theta_2_inp);
   // fChain->SetBranchAddress("Theta_3_inp", &Theta_3_inp, &b_Theta_3_inp);
   // fChain->SetBranchAddress("Theta_1_sm", &Theta_1_sm, &b_Theta_1_sm);
   // fChain->SetBranchAddress("Theta_2_sm", &Theta_2_sm, &b_Theta_2_sm);
   // fChain->SetBranchAddress("Theta_3_sm", &Theta_3_sm, &b_Theta_3_sm);
   // fChain->SetBranchAddress("Theta_1_error", &Theta_1_error, &b_Theta_1_error);
   // fChain->SetBranchAddress("Theta_2_error", &Theta_2_error, &b_Theta_2_error);
   // fChain->SetBranchAddress("Theta_3_error", &Theta_3_error, &b_Theta_3_error);
   // fChain->SetBranchAddress("Phi_1", &Phi_1, &b_Phi_1);
   // fChain->SetBranchAddress("Phi_2", &Phi_2, &b_Phi_2);
   // fChain->SetBranchAddress("Phi_3", &Phi_3, &b_Phi_3);
   // fChain->SetBranchAddress("Phi_1_inp", &Phi_1_inp, &b_Phi_1_inp);
   // fChain->SetBranchAddress("Phi_2_inp", &Phi_2_inp, &b_Phi_2_inp);
   // fChain->SetBranchAddress("Phi_3_inp", &Phi_3_inp, &b_Phi_3_inp);
   // fChain->SetBranchAddress("Phi_1_sm", &Phi_1_sm, &b_Phi_1_sm);
   // fChain->SetBranchAddress("Phi_2_sm", &Phi_2_sm, &b_Phi_2_sm);
   // fChain->SetBranchAddress("Phi_3_sm", &Phi_3_sm, &b_Phi_3_sm);
   // fChain->SetBranchAddress("Phi_1_error", &Phi_1_error, &b_Phi_1_error);
   // fChain->SetBranchAddress("Phi_2_error", &Phi_2_error, &b_Phi_2_error);
   // fChain->SetBranchAddress("Phi_3_error", &Phi_3_error, &b_Phi_3_error);
   // fChain->SetBranchAddress("Chi2", &Chi2, &b_Chi2);
   // fChain->SetBranchAddress("Probability", &Probability, &b_Probability);
   // fChain->SetBranchAddress("MissMassProt2Neut", &MissMassProt2Neut, &b_MissMassProt2Neut);
   // fChain->SetBranchAddress("MissMassProt2NeutPre", &MissMassProt2NeutPre, &b_MissMassProt2NeutPre);
   
   fChain->SetBranchAddress("Th1_f",&Theta_1);
   fChain->SetBranchAddress("Phi1_f",&Phi_1);
   fChain->SetBranchAddress("En1_f",&Energy_1);

   fChain->SetBranchAddress("Th2_f",&Theta_2);
   fChain->SetBranchAddress("Phi2_f",&Phi_2);
   fChain->SetBranchAddress("En2_f",&Energy_2);

   fChain->SetBranchAddress("Th1",&Theta_1_inp);
   fChain->SetBranchAddress("Phi1",&Phi_1_inp);
   fChain->SetBranchAddress("En1",&Energy_1_inp);

   fChain->SetBranchAddress("Th2",&Theta_2_inp);
   fChain->SetBranchAddress("Phi2",&Phi_2_inp);
   fChain->SetBranchAddress("En2",&Energy_2_inp);

   fChain->SetBranchAddress("Th1_sm",&Theta_1_sm);
   fChain->SetBranchAddress("Phi1_sm",&Phi_1_sm);
   fChain->SetBranchAddress("En1_sm",&Energy_1_sm);

   fChain->SetBranchAddress("Th2_sm",&Theta_2_sm);
   fChain->SetBranchAddress("Phi2_sm",&Phi_2_sm);
   fChain->SetBranchAddress("En2_sm",&Energy_2_sm);

   Notify();
}

Bool_t hist::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hist::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hist::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hist_cxx
