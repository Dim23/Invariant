//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 25 00:38:31 2019 by ROOT version 6.19/01
// from TTree mtree/Hadron EMC + TOF tree
// found on file: /home/dim2/NIR/310698.root
//////////////////////////////////////////////////////////

#ifndef hTANA_h
#define hTANA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define MAXNH 300
#define MAXRP 48


// Header file for the classes stored in the TTree if any.

class hTANA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Float_t         bbcz;
   Float_t         cent;
   Int_t           rh;
   Float_t         phir[MAXRP];   //[rh]
   Float_t         time[MAXRP];   //[rh]
   Float_t         qr0[MAXRP];   //[rh]
   Float_t         etar[MAXRP];   //[rh]
   Short_t         armr[MAXRP];   //[rh]
   Short_t         ring[MAXRP];   //[rh]
   Int_t           chid[MAXRP];   //[rh]
   Int_t           mh;
   Float_t         alpha[MAXNH];   //[mh]
   Short_t         dcarm[MAXNH];   //[mh]
   Float_t         p[MAXNH];   //[mh]
   Short_t         charge[MAXNH];   //[mh]
   Float_t         phi0[MAXNH];   //[mh]
   Float_t         the0[MAXNH];   //[mh]
   Float_t         phi[MAXNH];   //[mh]
   Float_t         ecore[MAXNH];   //[mh]
   Float_t         plemc[MAXNH];   //[mh]
   Float_t         ecent[MAXNH];   //[mh]
   Float_t         temc[MAXNH];   //[mh]
   Float_t         temcpi[MAXNH];   //[mh]
   Float_t         temcp[MAXNH];   //[mh]
   Float_t         temck[MAXNH];   //[mh]
   Short_t         sect[MAXNH];   //[mh]
   Float_t         isPiemc[MAXNH];   //[mh]
   Float_t         isPemc[MAXNH];   //[mh]
   Float_t         isKemc[MAXNH];   //[mh]
   Int_t           idtwr[MAXNH];   //[mh]
   Float_t         sigtof[MAXNH];   //[mh]
   Float_t         sigpc3[MAXNH];   //[mh]
   Float_t         sigemc[MAXNH];   //[mh]
   Float_t         res[MAXNH];   //[mh]
   Float_t         ttof[MAXNH];   //[mh]
   Int_t           slat[MAXNH];   //[mh]
   Float_t         pltof[MAXNH];   //[mh]
   Float_t         etof[MAXNH];   //[mh]
   Float_t         isPi[MAXNH];   //[mh]
   Float_t         isP[MAXNH];   //[mh]
   Float_t         isK[MAXNH];   //[mh]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_bbcz;   //!
   TBranch        *b_cent;   //!
   TBranch        *b_rh;   //!
   TBranch        *b_phir;   //!
   TBranch        *b_time;   //!
   TBranch        *b_qr0;   //!
   TBranch        *b_etar;   //!
   TBranch        *b_armr;   //!
   TBranch        *b_ring;   //!
   TBranch        *b_chid;   //!
   TBranch        *b_mh;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_dcarm;   //!
   TBranch        *b_p;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_phi0;   //!
   TBranch        *b_the0;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_ecore;   //!
   TBranch        *b_plemc;   //!
   TBranch        *b_ecent;   //!
   TBranch        *b_temc;   //!
   TBranch        *b_temcpi;   //!
   TBranch        *b_temcp;   //!
   TBranch        *b_temck;   //!
   TBranch        *b_sect;   //!
   TBranch        *b_isPiemc;   //!
   TBranch        *b_isPemc;   //!
   TBranch        *b_isKemc;   //!
   TBranch        *b_idtwr;   //!
   TBranch        *b_sigtof;   //!
   TBranch        *b_sigpc3;   //!
   TBranch        *b_sigemc;   //!
   TBranch        *b_res;   //!
   TBranch        *b_ttof;   //!
   TBranch        *b_slat;   //!
   TBranch        *b_pltof;   //!
   TBranch        *b_etof;   //!
   TBranch        *b_isPi;   //!
   TBranch        *b_isP;   //!
   TBranch        *b_isK;   //!

   hTANA(const char *file_adress);
   virtual ~hTANA();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


// additional member functions

   void add_file(const char *file);
   void ana_end(const char *file);
   void Loop_imp(void);
   void Book_Hist(void);
   void Fit_imp(void);
void INV(void);
void inv(void);


};

#endif

#ifdef hTANA_cxx
hTANA::hTANA(const char *file_adress) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tre


TFile *f = TFile::Open(file_adress);
  TTree *tree = (TTree*)f->Get("mtree");
      /*TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/dim2/NIR/310698.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/dim2/NIR/310698.root");
      }
      f->GetObject("mtree",tree);*/
 cout <<"tree found"<<endl;

  
   Init(tree);
}

hTANA::~hTANA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hTANA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hTANA::LoadTree(Long64_t entry)
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

void hTANA::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("bbcz", &bbcz, &b_bbcz);
   fChain->SetBranchAddress("cent", &cent, &b_cent);
   fChain->SetBranchAddress("rh", &rh, &b_rh);
   fChain->SetBranchAddress("phir", phir, &b_phir);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("qr0", qr0, &b_qr0);
   fChain->SetBranchAddress("etar", etar, &b_etar);
   fChain->SetBranchAddress("armr", armr, &b_armr);
   fChain->SetBranchAddress("ring", ring, &b_ring);
   fChain->SetBranchAddress("chid", chid, &b_chid);
   fChain->SetBranchAddress("mh", &mh, &b_mh);
   fChain->SetBranchAddress("alpha", alpha, &b_alpha);
   fChain->SetBranchAddress("dcarm", dcarm, &b_dcarm);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("phi0", phi0, &b_phi0);
   fChain->SetBranchAddress("the0", the0, &b_the0);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("ecore", ecore, &b_ecore);
   fChain->SetBranchAddress("plemc", plemc, &b_plemc);
   fChain->SetBranchAddress("ecent", ecent, &b_ecent);
   fChain->SetBranchAddress("temc", temc, &b_temc);
   fChain->SetBranchAddress("temcpi", temcpi, &b_temcpi);
   fChain->SetBranchAddress("temcp", temcp, &b_temcp);
   fChain->SetBranchAddress("temck", temck, &b_temck);
   fChain->SetBranchAddress("sect", sect, &b_sect);
   fChain->SetBranchAddress("isPiemc", isPiemc, &b_isPiemc);
   fChain->SetBranchAddress("isPemc", isPemc, &b_isPemc);
   fChain->SetBranchAddress("isKemc", isKemc, &b_isKemc);
   fChain->SetBranchAddress("idtwr", idtwr, &b_idtwr);
   fChain->SetBranchAddress("sigtof", sigtof, &b_sigtof);
   fChain->SetBranchAddress("sigpc3", sigpc3, &b_sigpc3);
   fChain->SetBranchAddress("sigemc", sigemc, &b_sigemc);
   fChain->SetBranchAddress("res", res, &b_res);
   fChain->SetBranchAddress("ttof", ttof, &b_ttof);
   fChain->SetBranchAddress("slat", slat, &b_slat);
   fChain->SetBranchAddress("pltof", pltof, &b_pltof);
   fChain->SetBranchAddress("etof", etof, &b_etof);
   fChain->SetBranchAddress("isPi", isPi, &b_isPi);
   fChain->SetBranchAddress("isP", isP, &b_isP);
   fChain->SetBranchAddress("isK", isK, &b_isK);
   Notify();
}

Bool_t hTANA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hTANA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hTANA::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hTANA_cxx
