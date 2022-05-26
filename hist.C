#define hist_cxx
#include "hist.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define Deg2Rad 0.0174533
#define Rad2Deg 57.2958
char name[200];
int iter=0;

void hist::Loop()
{
   // gStyle->SetOptFit();
   gStyle->SetHistLineWidth(2);
   gStyle->SetLabelSize(0.05, "xyz");
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(1101);
   gStyle->SetStatFontSize(0.07);
   TGaxis::SetMaxDigits(3);

   TFile *fout = new TFile("hist.root","RECREATE");

   TH1F* hConverged = new TH1F("hConverged","hConverged",200,-10,10);
   TH1F* hChiTheory = new TH1F("hChiTheory","hChiTheory",1000,0,2);

   TH1F* Res_spec_Ph1 = new TH1F("Res_spec_Ph1","Res_spec_Ph1",200,-0.25,0.25);
   TH1F* Res_spec_Ph2 = new TH1F("Res_spec_Ph2","Res_spec_Ph2",200,-0.25,0.25);
   // TH1F* Res_spec_Ph3 = new TH1F("Res_spec_Ph3","Res_spec_Ph3",200,-15,15);

   TH1F* QA_Res_Ph1 = new TH1F("QA_Res_Ph1","QA_Res_Ph1",200,-5,5);
   TH1F* QA_Res_Ph2 = new TH1F("QA_Res_Ph2","QA_Res_Ph2",200,-5,5);
   // TH1F* QA_Res_Ph3 = new TH1F("QA_Res_Ph3","QA_Res_Ph3",200,-40,40);

   TH1F* Res_spec_Th1 = new TH1F("Res_spec_Th1","Res_spec_Th1",200,-0.5,0.5);
   TH1F* Res_spec_Th2 = new TH1F("Res_spec_Th2","Res_spec_Th2",200,-0.5,0.5);
   // TH1F* Res_spec_Th3 = new TH1F("Res_spec_Th3","Res_spec_Th3",200,-15,15);

   TH1F* QA_Res_Th1 = new TH1F("QA_Res_Th1","QA_Res_Th1",200,-5,5);
   TH1F* QA_Res_Th2 = new TH1F("QA_Res_Th2","QA_Res_Th2",200,-5,5);
   // TH1F* QA_Res_Th3 = new TH1F("QA_Res_Th3","QA_Res_Th3",200,-50,50);

   TH1F* Res_spec_En1 = new TH1F("Res_spec_En1","Res_spec_En1",200,-16,16);
   TH1F* Res_spec_En2 = new TH1F("Res_spec_En2","Res_spec_En2",200,-16,16);
   // TH1F* Res_spec_En3 = new TH1F("Res_spec_En3","Res_spec_En3",200,-45,45);

   TH1F* QA_Res_En1 = new TH1F("QA_Res_En1","QA_Res_En1",200,-18,18);
   TH1F* QA_Res_En2 = new TH1F("QA_Res_En2","QA_Res_En2",200,-18,18);
   // TH1F* QA_Res_En3 = new TH1F("QA_Res_En3","QA_Res_En3",200,-120,120);

   TH1F* Pull_spec_Ph1 = new TH1F("Pull_spec_Ph1","Pull_spec_Ph1",200,-5,5);
   TH1F* Pull_spec_Ph2 = new TH1F("Pull_spec_Ph2","Pull_spec_Ph2",200,-5,5);
   // TH1F* Pull_spec_Ph3 = new TH1F("Pull_spec_Ph3","Pull_spec_Ph3",200,-30,30);

   TH1F* Pull_spec_Th1 = new TH1F("Pull_spec_Th1","Pull_spec_Th1",200,-5,5);
   TH1F* Pull_spec_Th2 = new TH1F("Pull_spec_Th2","Pull_spec_Th2",200,-5,5);
   // TH1F* Pull_spec_Th3 = new TH1F("Pull_spec_Th3","Pull_spec_Th3",200,-30,30);
   
   TH1F* Pull_spec_En1 = new TH1F("Pull_spec_En1","Pull_spec_En1",200,-8,8);
   TH1F* Pull_spec_En2 = new TH1F("Pull_spec_En2","Pull_spec_En2",200,-8,8);
   // TH1F* Pull_spec_En3 = new TH1F("Pull_spec_En3","Pull_spec_En3",200,-30,30);

   TH1F* MMProt2Neut = new TH1F("MMProt2Neut","MMProt2Neut",200,900,942);
   TH1F* MMProt2NeutPre = new TH1F("MMProt2NeutPre","MMProt2NeutPre",200,900,1000);

   TH2F* En1_smvsEn1_f = new TH2F("En1_smvsEn1_f","En1_smvsEn1_f",200,0,200,200,0,200);

   TH1F *Phi1_inp = new TH1F("Phi1_inp","Phi1_inp",200,-200,200);
   TH1F *Th1_inp = new TH1F("Th1_inp","Th1_inp",200,0,60);
   TH1F *En1_inp = new TH1F("En1_inp","En1_inp",200,-20,150);

   TH1F *Phi2_inp = new TH1F("Phi2_inp","Phi2_inp",200,-200,200);
   TH1F *Th2_inp = new TH1F("Th2_inp","Th2_inp",200,0,60);
   TH1F *En2_inp = new TH1F("En2_inp","En2_inp",200,-20,150);

   TH1F *Phi1_sm=new TH1F("Phi1_sm","Phi1_sm",200,-200,200);
   TH1F *Th1_sm=new TH1F("Th1_sm","Th1_sm",200,0,60);
   TH1F *En1_sm=new TH1F("En1_sm","En1_sm",200,-20,150);

   TH1F *Phi2_sm=new TH1F("Phi2_sm","Phi2_sm",200,-200,200);
   TH1F *Th2_sm=new TH1F("Th2_sm","Th2_sm",200,0,60);
   TH1F *En2_sm=new TH1F("En2_sm","En2_sm",200,0,150);

   // TH1F *Phi3_sm=new TH1F("Phi3_sm","Phi3_sm",200,-200,200);
   // TH1F *Th3_sm=new TH1F("Th3_sm","Th3_sm",200,0,60);
   // TH1F *En3_sm=new TH1F("En3_sm","En3_sm",200,0,150);

   TH1F *Phi1_fit=new TH1F("Phi1_fit","Phi1_fit",200,-200,200);
   TH1F *Th1_fit=new TH1F("Th1_fit","Th1_fit",200,0,60);
   TH1F *En1_fit=new TH1F("En1_fit","En1_fit",200,0,150);

   TH1F *Phi2_fit=new TH1F("Phi2_fit","Phi2_fit",200,-200,200);
   TH1F *Th2_fit=new TH1F("Th2_fit","Th2_fit",200,0,60);
   TH1F *En2_fit=new TH1F("En2_fit","En2_fit",200,0,150);

   // TH1F *Th3_fit=new TH1F("Th3_fit","Th3_fit",200,0,100);
   // TH1F *Phi3_fit=new TH1F("Phi3_fit","Phi3_fit",200,-200,200);
   // TH1F *En3_fit=new TH1F("En3_fit","En3_fit",200,0,150);

   TH1F *Phi1_f_Phi1_inp = new TH1F("Phi1_f_Phi1_inp","Phi1_f_Phi1_inp",200,-10,10);
   TH1F *Th1_f_Th1_inp = new TH1F("Th1_f_Th1_inp","Th1_f_Th1_inp",200,-10,10);
   TH1F *En1_f_En1_inp = new TH1F("En1_f_En1_inp","En1_f_En1_inp",200,-18,18);

   TH1F *Phi2_f_Phi2_inp = new TH1F("Phi2_f_Phi2_inp","Phi2_f_Phi2_inp",200,-10,10);
   TH1F *Th2_f_Th2_inp = new TH1F("Th2_f_Th2_inp","Th2_f_Th2_inp",200,-10,10);
   TH1F *En2_f_En2_inp = new TH1F("En2_f_En2_inp","En2_f_En2_inp",200,-18,18);

   // TH1F *Phi3_f_Phi3_inp = new TH1F("Phi3_f_Phi3_inp","Phi3_f_Phi3_inp",200,-10,10);
   // TH1F *Th3_f_Th3_inp = new TH1F("Th3_f_Th3_inp","Th3_f_Th3_inp",200,-10,10);
   // TH1F *En3_f_En3_inp = new TH1F("En3_f_En3_inp","En3_f_En3_inp",200,-10,10);

   TH1F *Phi1_inp_Phi1_smDivPhi1_sm = new TH1F("Phi1_inp_Phi1_smDivPhi1_sm", "Phi1_inp_Phi1_smDivPhi1_sm",200,-0.01,0.01);
   TH1F *Th1_inp_Th1_smDivTh1_sm = new TH1F("Th1_inp_Th1_smDivTh1_sm","Th1_inp_Th1_smDivTh1_sm",200,-0.5,0.5);
   TH1F *En1_inp_En1_smDivEn1_sm = new TH1F("En1_inp_En1_smDivEn1_sm","En1_inp_En1_smDivEn1_sm",200,-0.5,0.5);

   TH1F *Phi2_inp_Phi2_smDivPhi2_sm = new TH1F("Phi2_inp_Phi2_smDivPhi2_sm", "Phi2_inp_Phi2_smDivPhi2_sm",200,-0.01,0.01);
   TH1F *Th2_inp_Th2_smDivTh2_sm = new TH1F("Th2_inp_Th2_smDivTh2_sm","Th2_inp_Th2_smDivTh2_sm",200,-0.5,0.5);
   TH1F *En2_inp_En2_smDivEn2_sm = new TH1F("En2_inp_En2_smDivEn2_sm","En2_inp_En2_smDivEn2_sm",200,-0.5,0.5);

   TH1F *Phi1_inp_Phi1_fDivPhi1_f = new TH1F("Phi1_inp_Phi1_fDivPhi1_f","Phi1_inp_Phi1_fDivPhi1_f",200,-0.1,0.1);
   TH1F *Th1_inp_Th1_fDivTh1_f = new TH1F("Th1_inp_Th1_fDivTh1_f","Th1_inp_Th1_fDivTh1_f",200,-0.5,0.5);
   TH1F *En1_inp_En1_fDivEn1_f = new TH1F("En1_inp_En1_fDivEn1_f","En1_inp_En1_fDivEn1_f",200,-0.5,0.5);

   TH1F *Phi2_inp_Phi2_fDivPhi2_f = new TH1F("Phi2_inp_Phi2_fDivPhi2_f","Phi2_inp_Phi2_fDivPhi2_f",200,-0.1,0.1);
   TH1F *Th2_inp_Th2_fDivTh2_f = new TH1F("Th2_inp_Th2_fDivTh2_f","Th2_inp_Th2_fDivTh2_f",200,-0.5,0.5);
   TH1F *En2_inp_En2_fDivEn2_f = new TH1F("En2_inp_En2_fDivEn2_f","En2_inp_En2_fDivEn2_f",200,-0.5,0.5);

   TH1F *Phi1_sm_Phi1_fDivPhi1_f = new TH1F("Phi1_sm_Phi1_fDivPhi1_f","Phi1_sm_Phi1_fDivPhi1_f",200,-0.01,0.01);
   TH1F *Th1_sm_Th1_fDivTh1_f = new TH1F("Th1_sm_Th1_fDivTh1_f","Th1_sm_Th1_fDivTh1_f",200,-0.1,0.1);
   TH1F *En1_sm_En1_fDivEn1_f = new TH1F("En1_sm_En1_fDivEn1_f","En1_sm_En1_fDivEn1_f",200,-0.3,0.3);

   TH1F *Phi2_sm_Phi2_fDivPhi2_f = new TH1F("Phi2_sm_Phi2_fDivPhi2_f","Phi2_sm_Phi2_fDivPhi2_f",200,-0.01,0.01);
   TH1F *Th2_sm_Th2_fDivTh2_f = new TH1F("Th2_sm_Th2_fDivTh2_f","Th2_sm_Th2_fDivTh2_f",200,-0.1,0.1);
   TH1F *En2_sm_En2_fDivEn2_f = new TH1F("En2_sm_En2_fDivEn2_f","En2_sm_En2_fDivEn2_f",200,-0.3,0.3);
   
   TH1F *Phi1_inp_Phi1_sm = new TH1F("Phi1_inp_Phi1_sm","Phi1_inp_Phi1_sm",200,-6,6);
   TH1F *Th1_inp_Th1_sm = new TH1F("Th1_inp_Th1_sm","Th1_inp_Th1_sm",200,-6,6);
   TH1F *En1_inp_En1_sm = new TH1F("En1_inp_En1_sm","En1_inp_En1_sm",200,-18,18);

   TH1F *Phi2_inp_Phi2_sm = new TH1F("Phi2_inp_Phi2_sm","Phi2_inp_Phi2_sm",200,-6,6);
   TH1F *Th2_inp_Th2_sm = new TH1F("Th2_inp_Th2_sm","Th2_inp_Th2_sm",200,-6,6);
   TH1F *En2_inp_En2_sm = new TH1F("En2_inp_En2_sm","En2_inp_En2_sm",200,-18,18);

   TH1F *Phi1_inp_Phi1_f = new TH1F("Phi1_inp_Phi1_f","Phi1_inp_Phi1_f",200,-6,6);
   TH1F *Th1_inp_Th1_f = new TH1F("Th1_inp_Th1_f","Th1_inp_Th1_f",200,-6,6);
   TH1F *En1_inp_En1_f = new TH1F("En1_inp_En1_f","En1_inp_En1_f",200,-18,18);

   TH1F *Phi2_inp_Phi2_f = new TH1F("Phi2_inp_Phi2_f","Phi2_inp_Phi2_f",200,-6,6);
   TH1F *Th2_inp_Th2_f = new TH1F("Th2_inp_Th2_f","Th2_inp_Th2_f",200,-6,6);
   TH1F *En2_inp_En2_f = new TH1F("En2_inp_En2_f","En2_inp_En2_f",200,-18,18);

   TH1F *Phi1_sm_Phi1_f = new TH1F("Phi1_sm_Phi1_f","Phi1_sm_Phi1_f",200,-0.5,0.5);
   TH1F *Th1_sm_Th1_f = new TH1F("Th1_sm_Th1_f","Th1_sm_Th1_f",200,-0.5,0.5);
   TH1F *En1_sm_En1_f = new TH1F("En1_sm_En1_f","En1_sm_En1_f",200,-10,10);

   TH1F *Phi2_sm_Phi2_f = new TH1F("Phi2_sm_Phi2_f","Phi2_sm_Phi2_f",200,-0.5,0.5);
   TH1F *Th2_sm_Th2_f = new TH1F("Th2_sm_Th2_f","Th2_sm_Th2_f",200,-0.5,0.5);
   TH1F *En2_sm_En2_f = new TH1F("En2_sm_En2_f","En2_sm_En2_f",200,-10,10);

   Float_t sigma_En1, sigma_En2;
   Float_t sigma_Th1, sigma_Th2;
   Float_t sigma_Ph1, sigma_Ph2;

   double chitheory;
   for(double i=0;i<2;i=i+0.000001)
   {
      chitheory=ROOT::Math::chisquared_pdf(i,1,0);
      hChiTheory->Fill(chitheory);
   }

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(Converged==1 && Probability > 0.02) 
      {
         Res_spec_Ph1->Fill(Phi_1_sm-Phi_1);
         Res_spec_Ph2->Fill(Phi_2_sm-Phi_2);
         // Res_spec_Ph3->Fill(Phi_3_sm-Phi_3);

         Res_spec_Th1->Fill(Theta_1_sm-Theta_1);
         Res_spec_Th2->Fill(Theta_2_sm-Theta_2);
         // Res_spec_Th3->Fill(Theta_3_sm-Theta_3);

         Res_spec_En1->Fill(Energy_1_sm-Energy_1);
         Res_spec_En2->Fill(Energy_2_sm-Energy_2);
         // Res_spec_En3->Fill(Energy_3_sm-Energy_3);
         
         Phi1_f_Phi1_inp->Fill(Phi_1-Phi_1_inp);
         Th1_f_Th1_inp->Fill(Theta_1-Theta_1_inp);
         En1_f_En1_inp->Fill(Energy_1-Energy_1_inp);

         Phi2_f_Phi2_inp->Fill(Phi_2-Phi_2_inp);
         Th2_f_Th2_inp->Fill(Theta_2-Theta_2_inp);
         En2_f_En2_inp->Fill(Energy_2-Energy_2_inp);

         // Phi3_f_Phi3_inp->Fill(Phi_3-Phi_3_inp);
         // Th3_f_Th3_inp->Fill(Theta_3-Theta_3_inp);
         // // En3_f_En3_inp->Fill(En_3-En_3_inp);

         En1_smvsEn1_f->Fill(Energy_1_sm,Energy_1);

         Phi1_fit->Fill(Phi_1);
         Th1_fit->Fill(Theta_1);
         En1_fit->Fill(Energy_1);

         Phi2_fit->Fill(Phi_2);
         Th2_fit->Fill(Theta_2);
         En2_fit->Fill(Energy_2);

         // Phi3_fit->Fill(Phi_3);
         // Th3_fit->Fill(Theta_3);
         // En3_fit->Fill(Energy_3);

         Phi1_inp_Phi1_smDivPhi1_sm->Fill((Phi_1 - Phi_1_sm) / Phi_1_sm);
         Th1_inp_Th1_smDivTh1_sm->Fill((Theta_1_inp - Theta_1_sm) / Theta_1_sm); 
         En1_inp_En1_smDivEn1_sm->Fill((Energy_1_inp - Energy_1_sm) / Energy_1_sm);

         Phi2_inp_Phi2_smDivPhi2_sm->Fill((Phi_2 - Phi_2_sm) / Phi_2_sm);
         Th2_inp_Th2_smDivTh2_sm->Fill((Theta_2_inp - Theta_2_sm) / Theta_2_sm); 
         En2_inp_En2_smDivEn2_sm->Fill((Energy_2_inp - Energy_2_sm) / Energy_2_sm);

         Phi1_inp_Phi1_fDivPhi1_f->Fill((Phi_1_inp - Phi_1) / Phi_1);
         Th1_inp_Th1_fDivTh1_f->Fill((Theta_1_inp - Theta_1) / Theta_1); 
         En1_inp_En1_fDivEn1_f->Fill((Energy_1_inp - Energy_1) / Energy_1); 

         Phi2_inp_Phi2_fDivPhi2_f->Fill((Phi_2_inp - Phi_2) / Phi_2);
         Th2_inp_Th2_fDivTh2_f->Fill((Theta_2_inp - Theta_2) / Theta_2); 
         En2_inp_En2_fDivEn2_f->Fill((Energy_2_inp - Energy_2) / Energy_2); 

         Phi1_sm_Phi1_fDivPhi1_f->Fill((Phi_1_sm - Phi_1) / Phi_1); 
         Th1_sm_Th1_fDivTh1_f->Fill((Theta_1_sm - Theta_1) / Theta_1); 
         En1_sm_En1_fDivEn1_f->Fill((Energy_1_sm - Energy_1) / Energy_1); 

         Phi2_sm_Phi2_fDivPhi2_f->Fill((Phi_2_sm - Phi_2) / Phi_2); 
         Th2_sm_Th2_fDivTh2_f->Fill((Theta_2_sm - Theta_2) / Theta_2); 
         En2_sm_En2_fDivEn2_f->Fill((Energy_2_sm - Energy_2) / Energy_2); 

         Phi1_inp_Phi1_sm->Fill(Phi_1_inp - Phi_1_sm);
         Th1_inp_Th1_sm->Fill(Theta_1_inp - Theta_1_sm);
         En1_inp_En1_sm->Fill(Energy_1_inp - Energy_1_sm);

         Phi2_inp_Phi2_sm->Fill(Phi_2_inp - Phi_2_sm);
         Th2_inp_Th2_sm->Fill(Theta_2_inp - Theta_2_sm);
         En2_inp_En2_sm->Fill(Energy_2_inp - Energy_2_sm);

         Phi1_inp_Phi1_f->Fill(Phi_1_inp - Phi_1);
         Th1_inp_Th1_f->Fill(Theta_1_inp - Theta_1);
         En1_inp_En1_f->Fill(Energy_1_inp - Energy_1);

         Phi2_inp_Phi2_f->Fill(Phi_2_inp - Phi_2);
         Th2_inp_Th2_f->Fill(Theta_2_inp - Theta_2);
         En2_inp_En2_f->Fill(Energy_2_inp - Energy_2);

         Phi1_sm_Phi1_f->Fill(Phi_1_sm - Phi_1);
         Th1_sm_Th1_f->Fill(Theta_1_sm - Theta_1);
         En1_sm_En1_f->Fill(Energy_1_sm - Energy_1); 

         Phi2_sm_Phi2_f->Fill(Phi_2_sm - Phi_2);
         Th2_sm_Th2_f->Fill(Theta_2_sm - Theta_2);
         En2_sm_En2_f->Fill(Energy_2_sm - Energy_2); 
      }
      
      QA_Res_Ph1->Fill(Phi_1_sm-Phi_1_inp);
      QA_Res_Ph2->Fill(Phi_2_sm-Phi_2_inp);
      // QA_Res_Ph3->Fill(Phi_3_sm-Phi_3_inp);

      QA_Res_Th1->Fill(Theta_1_sm-Theta_1_inp);
      QA_Res_Th2->Fill(Theta_2_sm-Theta_2_inp);
      // QA_Res_Th3->Fill(Theta_3_sm-Theta_3_inp);

      QA_Res_En1->Fill(Energy_1_inp-Energy_1_sm);
      QA_Res_En2->Fill(Energy_2_inp-Energy_2_sm);
      // QA_Res_En3->Fill(Energy_3_inp-Energy_3_sm);

      hConverged->Fill(Converged);

      Phi1_inp->Fill(Phi_1_inp);
      Th1_inp->Fill(Theta_1_inp);
      En1_inp->Fill(Energy_1_inp);

      Phi2_inp->Fill(Phi_2_inp);
      Th2_inp->Fill(Theta_2_inp);
      En2_inp->Fill(Energy_2_inp);

      Phi1_sm->Fill(Phi_1_sm);
      Th1_sm->Fill(Theta_1_sm);
      En1_sm->Fill(Energy_1_sm);

      Phi2_sm->Fill(Phi_2_sm);
      Th2_sm->Fill(Theta_2_sm);
      En2_sm->Fill(Energy_2_sm);

      // Phi3_sm->Fill(Phi_3_sm);
      // Th3_sm->Fill(Theta_3_sm);
      // En3_sm->Fill(Energy_3_sm);
   }

   sigma_Ph1=Res_spec_Ph1->GetStdDev(1);
   sigma_Ph2=Res_spec_Ph2->GetStdDev(1);
   // sigma_Ph3=Res_spec_Ph3->GetStdDev(1);

   sigma_Th1=Res_spec_Th1->GetStdDev(1);
   sigma_Th2=Res_spec_Th2->GetStdDev(1);
   // sigma_Th3=Res_spec_Th3->GetStdDev(1);

   sigma_En1=Res_spec_En1->GetStdDev(1);
   sigma_En2=Res_spec_En2->GetStdDev(1);
   // sigma_En3=Res_spec_En3->GetStdDev(1);

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(Converged==1 && Probability > 0.02)
      {
         Pull_spec_Ph1->Fill((Phi_1_sm-Phi_1)/sigma_Ph1);
         Pull_spec_Ph2->Fill((Phi_2_sm-Phi_2)/sigma_Ph2);
         // Pull_spec_Ph3->Fill((Phi_3_sm*Rad2Deg-Phi_3*Rad2Deg)/sigma_Ph3);

         Pull_spec_Th1->Fill((Theta_1_sm-Theta_1)/sigma_Th1);
         Pull_spec_Th2->Fill((Theta_2_sm-Theta_2)/sigma_Th2);
         // Pull_spec_Th3->Fill((Theta_3_sm*Rad2Deg-Theta_3*Rad2Deg)/sigma_Th3);

         Pull_spec_En1->Fill((Energy_1_sm-Energy_1)/sigma_En1);
         Pull_spec_En2->Fill((Energy_2_sm-Energy_2)/sigma_En2);
         // Pull_spec_En3->Fill((Energy_3_sm*1000-Energy_3*1000)/sigma_En3);
      }
      
   }

   fout->cd();

   Res_spec_Ph1->Write();
   Res_spec_Ph2->Write();
   // Res_spec_Ph3->Write();

   QA_Res_Ph1->Write();
   QA_Res_Ph2->Write();
   // QA_Res_Ph3->Write();

   Res_spec_Th1->Write();
   Res_spec_Th2->Write();
   // Res_spec_Th3->Write();

   QA_Res_Th1->Write();
   QA_Res_Th2->Write();
   // QA_Res_Th3->Write();

   Res_spec_En1->Write();
   Res_spec_En2->Write();
   // Res_spec_En3->Write();

   QA_Res_En1->Write();
   QA_Res_En2->Write();
   // QA_Res_En3->Write();

   Pull_spec_Ph1->Write();
   Pull_spec_Ph2->Write();
   // Pull_spec_Ph3->Write();

   Pull_spec_Th1->Write();
   Pull_spec_Th2->Write();
   // Pull_spec_Th3->Write();

   Pull_spec_En1->Write();
   Pull_spec_En2->Write();
   // Pull_spec_En3->Write();

   Phi1_f_Phi1_inp->Write();
   Th1_f_Th1_inp->Write();
   En1_f_En1_inp->Write();

   Phi2_f_Phi2_inp->Write();
   Th2_f_Th2_inp->Write();
   En2_f_En2_inp->Write();

   hConverged->Write();
   hChiTheory->Write();

   En1_smvsEn1_f->Write();

   Phi1_inp->Write();
   Th1_inp->Write();
   En1_inp->Write();

   Phi2_inp->Write();
   Th2_inp->Write();
   En2_inp->Write();

   Phi1_sm->Write();
   Th1_sm->Write();
   En1_sm->Write();

   Phi1_fit->Write();
   Th1_fit->Write();
   En1_fit->Write();

   Phi1_inp_Phi1_smDivPhi1_sm->Write();
   Th1_inp_Th1_smDivTh1_sm->Write();
   En1_inp_En1_smDivEn1_sm->Write();

   Phi2_inp_Phi2_smDivPhi2_sm->Write();
   Th2_inp_Th2_smDivTh2_sm->Write();
   En2_inp_En2_smDivEn2_sm->Write();

   Phi1_inp_Phi1_fDivPhi1_f->Write();
   Th1_inp_Th1_fDivTh1_f->Write();
   En1_inp_En1_fDivEn1_f->Write();

   Phi2_inp_Phi2_fDivPhi2_f->Write();
   Th2_inp_Th2_fDivTh2_f->Write();
   En2_inp_En2_fDivEn2_f->Write();

   Phi1_sm_Phi1_fDivPhi1_f->Write();
   Th1_sm_Th1_fDivTh1_f->Write();
   En1_sm_En1_fDivEn1_f->Write();

   Phi2_sm_Phi2_fDivPhi2_f->Write();
   Th2_sm_Th2_fDivTh2_f->Write();
   En2_sm_En2_fDivEn2_f->Write();

   Phi1_inp_Phi1_sm->Write();
   Th1_inp_Th1_sm->Write();
   En1_inp_En1_sm->Write();

   Phi2_inp_Phi2_sm->Write();
   Th2_inp_Th2_sm->Write();
   En2_inp_En2_sm->Write();

   Phi1_inp_Phi1_f->Write();
   Th1_inp_Th1_f->Write();
   En1_inp_En1_f->Write();

   Phi2_inp_Phi2_f->Write();
   Th2_inp_Th2_f->Write();
   En2_inp_En2_f->Write();

   Phi1_sm_Phi1_f->Write();
   Th1_sm_Th1_f->Write();
   En1_sm_En1_f->Write();

   Phi2_sm_Phi2_f->Write();
   Th2_sm_Th2_f->Write();
   En2_sm_En2_f-> Write();


   TCanvas *Residuals = new TCanvas("Residuals","Residuals",3000,2000);
   Residuals->Divide(3,2);
   TCanvas *QA_res = new TCanvas("QA_res","QA_res",3000,2000);
   QA_res->Divide(3,2);
   TCanvas *Pull = new TCanvas("Pull","Pull",3000,2000);
   Pull->Divide(3,2);
   TCanvas *Compare = new TCanvas("Compare","Compare",3000,1000);
   Compare->Divide(3,1);
   TCanvas *Compare2 = new TCanvas("Compare2","Compare2",3000,1000);
   Compare2->Divide(3,1);
   TCanvas *QAandRes1 = new TCanvas("QAandRes1","QAandRes1",3000,2000);
   QAandRes1->Divide(3,2);
   TCanvas *QAandRes2 = new TCanvas("QAandRes2","QAandRes2",3000,2000);
   QAandRes2->Divide(3,2);
   TCanvas *FitInp = new TCanvas("FitInp","FitInp",3000,2000);
   FitInp->Divide(3,2);
   TCanvas *PhiSummary1 = new TCanvas("PhiSummary1","PhiSummary1",3000,2000);
   PhiSummary1->Divide(3,2);
   TCanvas *ThetaSummary1 = new TCanvas("ThetaSummary1","ThetaSummary1",3000,2000);
   ThetaSummary1->Divide(3,2);
   TCanvas *EnergySummary1 = new TCanvas("EnergySummary1","EnergySummary1",3000,2000);
   EnergySummary1->Divide(3,2);
   TCanvas *PhiSummary2 = new TCanvas("PhiSummary2","PhiSummary2",3000,2000);
   PhiSummary2->Divide(3,2);
   TCanvas *ThetaSummary2 = new TCanvas("ThetaSummary2","ThetaSummary2",3000,2000);
   ThetaSummary2->Divide(3,2);
   TCanvas *EnergySummary2 = new TCanvas("EnergySummary2","EnergySummary2",3000,2000);
   EnergySummary2->Divide(3,2);
   TCanvas *PhiSummary3 = new TCanvas("PhiSummary3","PhiSummary3",3000,2000);
   PhiSummary3->Divide(3,2);
   TCanvas *ThetaSummary3 = new TCanvas("ThetaSummary3","ThetaSummary3",3000,2000);
   ThetaSummary3->Divide(3,2);
   TCanvas *EnergySummary3 = new TCanvas("EnergySummary3","EnergySummary3",3000,2000);
   EnergySummary3->Divide(3,2);
   TCanvas *SmAndInp1 = new TCanvas("SmAndInp1", "SmAndInp1", 1300,1000);
   SmAndInp1->Divide(2,2);
   TCanvas *SmAndInp2 = new TCanvas("SmAndInp2", "SmAndInp2", 1300,1000);
   SmAndInp2->Divide(2,2);

   Residuals->cd(1);
   Res_spec_Ph1->Draw();
   Residuals->cd(4);
   Res_spec_Ph2->Draw();

   Residuals->cd(2);
   Res_spec_Th1->Draw();
   Residuals->cd(5);
   Res_spec_Th2->Draw();

   Residuals->cd(3);
   Res_spec_En1->Draw();
   Residuals->cd(6);
   Res_spec_En2->Draw();

   Residuals->Write();
   Residuals->SaveAs("Residuals.png");

   QA_res->cd(1);
   QA_Res_Ph1->Draw();
   QA_res->cd(4);
   QA_Res_Ph2->Draw();

   QA_res->cd(2);
   QA_Res_Th1->Draw();
   QA_res->cd(5);
   QA_Res_Th2->Draw();

   QA_res->cd(3);
   QA_Res_En1->Draw();
   QA_res->cd(6);
   QA_Res_En2->Draw();

   QA_res->Write();
   QA_res->SaveAs("QA_res.png");

   Pull->cd(1);
   Pull_spec_Ph1->Draw();
   Pull->cd(4);
   Pull_spec_Ph2->Draw();

   Pull->cd(2);
   Pull_spec_Th1->Draw();
   Pull->cd(5);
   Pull_spec_Th2->Draw();

   Pull->cd(3);
   Pull_spec_En1->Draw();
   Pull->cd(6);
   Pull_spec_En2->Draw();

   Pull->Write();
   Pull->SaveAs("Pull.png");

   FitInp->cd(1);
   Th1_f_Th1_inp->Draw();
   FitInp->cd(2);
   Phi1_f_Phi1_inp->Draw();
   FitInp->cd(3);
   En1_f_En1_inp->Draw();

   FitInp->cd(4);
   Th2_f_Th2_inp->Draw();
   FitInp->cd(5);
   Phi2_f_Phi2_inp->Draw();
   FitInp->cd(6);
   En2_f_En2_inp->Draw();

   FitInp->Write();
   FitInp->SaveAs("FitInp.png");

   gStyle->SetPalette(1);

   Compare->cd(1);
   Phi1_sm->Scale(1./Phi1_sm->Integral());
   Phi1_sm->Draw("PLC");

   Phi1_fit->Scale(1./Phi1_fit->Integral());
   Phi1_fit->Draw("SAME PLC");

   Compare->cd(2);
   Th1_sm->Scale(1./Th1_sm->Integral());
   Th1_sm->Draw("PLC");

   Th1_fit->Scale(1./Th1_fit->Integral());
   Th1_fit->Draw("SAME PLC");//first is violet

   Compare->cd(3);
   En1_sm->Scale(1./En1_sm->Integral());
   En1_sm->Draw("PLC");

   En1_fit->Scale(1./En1_fit->Integral());
   En1_fit->Draw("SAME PLC");

   Compare->Write();
   Compare->SaveAs("Compare.png");

   gStyle->SetPalette(1);

   Compare2->cd(1);
   Phi2_sm->Scale(1./Phi2_sm->Integral());
   Phi2_sm->Draw("PLC");

   Phi2_fit->Scale(1./Phi2_fit->Integral());
   Phi2_fit->Draw("SAME PLC");

   Compare2->cd(2);
   Th2_sm->Scale(1./Th2_sm->Integral());
   Th2_sm->Draw("PLC");

   Th2_fit->Scale(1./Th2_fit->Integral());
   Th2_fit->Draw("SAME PLC");//first is violet

   Compare2->cd(3);
   En2_sm->Scale(1./En2_sm->Integral());
   En2_sm->Draw("PLC");

   En2_fit->Scale(1./En2_fit->Integral());
   En2_fit->Draw("SAME PLC");

   Compare2->Write();
   Compare2->SaveAs("Compare2.png");

   QAandRes1->cd(1);
   QA_Res_Ph1->Draw();
   QAandRes1->cd(2);
   QA_Res_Th1->Draw();
   QAandRes1->cd(3);
   QA_Res_En1->Draw();
   QAandRes1->cd(4);
   Res_spec_Ph1->Draw();
   QAandRes1->cd(5);
   Res_spec_Th1->Draw();
   QAandRes1->cd(6);
   Res_spec_En1->Draw();

   QAandRes1->Write();
   QAandRes1->SaveAs("QAandRes1.png");

   QAandRes2->cd(1);
   QA_Res_Ph2->Draw();
   QAandRes2->cd(2);
   QA_Res_Th2->Draw();
   QAandRes2->cd(3);
   QA_Res_En2->Draw();
   QAandRes2->cd(4);
   Res_spec_Ph2->Draw();
   QAandRes2->cd(5);
   Res_spec_Th2->Draw();
   QAandRes2->cd(6);
   Res_spec_En2->Draw();

   QAandRes2->Write();
   QAandRes2->SaveAs("QAandRes2.png");

   PhiSummary1->cd(1);
   Res_spec_Ph1->Draw();
   PhiSummary1->cd(2);
   Pull_spec_Ph1->Draw();
   PhiSummary1->cd(3);
   Phi1_f_Phi1_inp->Draw();
   PhiSummary1->cd(4);
   Res_spec_Ph2->Draw();
   PhiSummary1->cd(5);
   Pull_spec_Ph2->Draw();
   PhiSummary1->cd(6);
   Phi2_f_Phi2_inp->Draw();

   PhiSummary1->Write();
   PhiSummary1->SaveAs("PhiSummary1.png");

   PhiSummary2->cd(1);

   TH1F *histogram;
   histogram = (TH1F*)fout->Get("Phi1_inp_Phi1_smDivPhi1_sm");
   double max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.00024,max+0.00018);
   Phi1_inp_Phi1_smDivPhi1_sm->Draw();
   PhiSummary2->cd(2);
   histogram = (TH1F*)fout->Get("Phi1_inp_Phi1_fDivPhi1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.017,max+0.017);
   Phi1_inp_Phi1_fDivPhi1_f->Draw();
   PhiSummary2->cd(3);
   histogram = (TH1F*)fout->Get("Phi1_sm_Phi1_fDivPhi1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.0002,max+0.0003);
   Phi1_sm_Phi1_fDivPhi1_f->Draw();
   PhiSummary2->cd(4);
   histogram = (TH1F*)fout->Get("Phi2_inp_Phi2_smDivPhi2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.00024,max+0.00018);
   Phi2_inp_Phi2_smDivPhi2_sm->Draw();
   PhiSummary2->cd(5);
   histogram = (TH1F*)fout->Get("Phi2_inp_Phi2_fDivPhi2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.017,max+0.017);
   Phi2_inp_Phi2_fDivPhi2_f->Draw();
   PhiSummary2->cd(6);
   histogram = (TH1F*)fout->Get("Phi2_sm_Phi2_fDivPhi2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.0002,max+0.0003);
   Phi2_sm_Phi2_fDivPhi2_f->Draw();

   PhiSummary2->Write();
   PhiSummary2->SaveAs("PhiSummary2.png");

   ThetaSummary1->cd(1);
   Res_spec_Th1->Draw();
   ThetaSummary1->cd(2);
   Pull_spec_Th1->Draw();
   ThetaSummary1->cd(3);
   Th1_f_Th1_inp->Draw();
   ThetaSummary1->cd(4);
   Res_spec_Th2->Draw();
   ThetaSummary1->cd(5);
   Pull_spec_Th2->Draw();
   ThetaSummary1->cd(6);
   Th2_f_Th2_inp->Draw();

   ThetaSummary1->Write();
   ThetaSummary1->SaveAs("ThetaSummary1.png");

   ThetaSummary2->cd(1);
   histogram = (TH1F*)fout->Get("Th1_inp_Th1_smDivTh1_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.08,max+0.08);
   Th1_inp_Th1_smDivTh1_sm->Draw();
   ThetaSummary2->cd(2);
   histogram = (TH1F*)fout->Get("Th1_inp_Th1_fDivTh1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.08,max+0.08);
   Th1_inp_Th1_fDivTh1_f->Draw();
   ThetaSummary2->cd(3);
   histogram = (TH1F*)fout->Get("Th1_sm_Th1_fDivTh1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.002,max+0.002);
   Th1_sm_Th1_fDivTh1_f->Draw();
   ThetaSummary2->cd(4);
   histogram = (TH1F*)fout->Get("Th2_inp_Th2_smDivTh2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.08,max+0.08);
   Th2_inp_Th2_smDivTh2_sm->Draw();
   ThetaSummary2->cd(5);
   histogram = (TH1F*)fout->Get("Th2_inp_Th2_fDivTh2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.08,max+0.08);
   Th2_inp_Th2_fDivTh2_f->Draw();
   ThetaSummary2->cd(6);
   histogram = (TH1F*)fout->Get("Th2_sm_Th2_fDivTh2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.0025,max+0.0025);
   Th2_sm_Th2_fDivTh2_f->Draw();

   ThetaSummary2->Write();
   ThetaSummary2->SaveAs("ThetaSummary2.png");

   EnergySummary1->cd(1);
   Res_spec_En1->Draw();
   EnergySummary1->cd(2);
   Pull_spec_En1->Draw();
   EnergySummary1->cd(3);
   En1_f_En1_inp->Draw();
   
   EnergySummary1->cd(4);
   Res_spec_En2->Draw();
   EnergySummary1->cd(5);
   Pull_spec_En2->Draw();
   EnergySummary1->cd(6);
   En2_f_En2_inp->Draw();

   EnergySummary1->Write();
   EnergySummary1->SaveAs("EnergySummary1.png");

   EnergySummary2->cd(1);
   histogram = (TH1F*)fout->Get("En1_inp_En1_smDivEn1_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.14,max+0.15);
   En1_inp_En1_smDivEn1_sm->Draw();
   EnergySummary2->cd(2);
   histogram = (TH1F*)fout->Get("En1_inp_En1_fDivEn1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.10,max+0.11);
   En1_inp_En1_fDivEn1_f->Draw();
   EnergySummary2->cd(3);
   histogram = (TH1F*)fout->Get("En1_sm_En1_fDivEn1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.010,max+0.010);
   En1_sm_En1_fDivEn1_f->Draw();
   EnergySummary2->cd(4);
   histogram = (TH1F*)fout->Get("En2_inp_En2_smDivEn2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.13,max+0.15);
   En2_inp_En2_smDivEn2_sm->Draw();
   EnergySummary2->cd(5);
   histogram = (TH1F*)fout->Get("En2_inp_En2_fDivEn2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.10,max+0.10);
   En2_inp_En2_fDivEn2_f->Draw();
   EnergySummary2->cd(6);
   histogram = (TH1F*)fout->Get("En2_sm_En2_fDivEn2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.010,max+0.010);
   En2_sm_En2_fDivEn2_f->Draw();

   EnergySummary2->Write();
   EnergySummary2->SaveAs("EnergySummary2.png");

   PhiSummary3->cd(1);
   histogram = (TH1F*)fout->Get("Phi1_inp_Phi1_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Phi1_inp_Phi1_sm->Draw();
   PhiSummary3->cd(2);
   histogram = (TH1F*)fout->Get("Phi1_inp_Phi1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Phi1_inp_Phi1_f->Draw();
   PhiSummary3->cd(3);
   histogram = (TH1F*)fout->Get("Phi1_sm_Phi1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.017,max+0.014);
   Phi1_sm_Phi1_f->Draw();
   PhiSummary3->cd(4);
   histogram = (TH1F*)fout->Get("Phi2_inp_Phi2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Phi2_inp_Phi2_sm->Draw();
   PhiSummary3->cd(5);
   histogram = (TH1F*)fout->Get("Phi2_inp_Phi2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Phi2_inp_Phi2_f->Draw();
   PhiSummary3->cd(6);
   histogram = (TH1F*)fout->Get("Phi2_sm_Phi2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.016,max+0.014);
   Phi2_sm_Phi2_f->Draw();

   PhiSummary3->Write();
   PhiSummary3->SaveAs("PhiSummary3.png");

   ThetaSummary3->cd(1);
   histogram = (TH1F*)fout->Get("Th1_inp_Th1_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Th1_inp_Th1_sm->Draw();
   ThetaSummary3->cd(2);
   histogram = (TH1F*)fout->Get("Th1_inp_Th1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Th1_inp_Th1_f->Draw();
   ThetaSummary3->cd(3);
   histogram = (TH1F*)fout->Get("Th1_sm_Th1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.015,max+0.018);
   Th1_sm_Th1_f->Draw();
   ThetaSummary3->cd(4);
   histogram = (TH1F*)fout->Get("Th2_inp_Th2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Th2_inp_Th2_sm->Draw();
   ThetaSummary3->cd(5);
   histogram = (TH1F*)fout->Get("Th2_inp_Th2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-2.3,max+2.3);
   Th2_inp_Th2_f->Draw();
   ThetaSummary3->cd(6);
   histogram = (TH1F*)fout->Get("Th2_sm_Th2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.017,max+0.014);
   Th2_sm_Th2_f->Draw();

   ThetaSummary3->Write();
   ThetaSummary3->SaveAs("ThetaSummary3.png");

   EnergySummary3->cd(1);
   histogram = (TH1F*)fout->Get("En1_inp_En1_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-10,max+10);
   En1_inp_En1_sm->Draw();
   EnergySummary3->cd(2);
   histogram = (TH1F*)fout->Get("En1_inp_En1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-8,max+8);
   En1_inp_En1_f->Draw();
   EnergySummary3->cd(3);
   histogram = (TH1F*)fout->Get("En1_sm_En1_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.3,max+0.3);
   En1_sm_En1_f->Draw();
   EnergySummary3->cd(4);
   histogram = (TH1F*)fout->Get("En2_inp_En2_sm");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-10,max+10);
   En2_inp_En2_sm->Draw();
   EnergySummary3->cd(5);
   histogram = (TH1F*)fout->Get("En2_inp_En2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-8,max+8);
   En2_inp_En2_f->Draw();
   EnergySummary3->cd(6);
   histogram = (TH1F*)fout->Get("En2_sm_En2_f");
   max = histogram->GetBinCenter(histogram->GetMaximumBin());
   histogram->Fit("gaus","WL","",max-0.3,max+0.3);
   En2_sm_En2_f->Draw();

   EnergySummary3->Write();
   EnergySummary3->SaveAs("EnergySummary3.png");

   SmAndInp1->cd(1);
   Phi1_sm->GetXaxis()->SetTitle("Azimuthal angle #phi [deg]");
   Phi1_sm->GetXaxis()->SetTitleSize(0.051);
   Phi1_sm->GetXaxis()->SetLabelOffset(0.01);
   Phi1_sm->GetYaxis()->SetTitle("Counts (normalized)");
   Phi1_sm->GetYaxis()->SetTitleSize(0.051);
   Phi1_sm->GetYaxis()->SetLabelOffset(0.01);
   Phi1_sm->GetXaxis()->SetNdivisions(9);
   Phi1_sm->Scale(1./Phi1_sm->Integral());
   Phi1_sm->Draw("PLC");
   Phi1_inp->Scale(1./Phi1_inp->Integral());
   Phi1_inp->Draw("SAME PLC");

   SmAndInp1->cd(2);
   Th1_sm->GetXaxis()->SetTitle("Polar angle #theta [deg]");
   Th1_sm->GetXaxis()->SetTitleSize(0.051);
   Th1_sm->GetXaxis()->SetLabelOffset(0.01);
   Th1_sm->GetYaxis()->SetTitle("Counts (normalized)");
   Th1_sm->GetYaxis()->SetTitleSize(0.051);
   Th1_sm->GetYaxis()->SetLabelOffset(0.01);
   Th1_sm->Scale(1./Th1_sm->Integral());
   Th1_sm->Draw("PLC");
   Th1_inp->Scale(1./Th1_inp->Integral());
   Th1_inp->Draw("SAME PLC");

   SmAndInp1->cd(3);
   En1_sm->GetXaxis()->SetTitle("Energy [MeV]");
   En1_sm->GetXaxis()->SetTitleSize(0.051);
   En1_sm->GetXaxis()->SetLabelOffset(0.01);
   En1_sm->GetYaxis()->SetTitle("Counts (normalized)");
   En1_sm->GetYaxis()->SetTitleSize(0.051);
   En1_sm->GetYaxis()->SetLabelOffset(0.01);
   En1_sm->Scale(1./En1_sm->Integral());
   En1_sm->Draw("PLC");
   En1_inp->Scale(1./En1_inp->Integral());
   En1_inp->Draw("SAME PLC");

   SmAndInp1->Write();
   SmAndInp1->SaveAs("SmAndInp1.png");
   SmAndInp1->SaveAs("SmAndInp1.root");

   SmAndInp2->cd(1);
   Phi2_sm->GetXaxis()->SetTitle("Azimuthal angle #phi [deg]");
   Phi2_sm->GetXaxis()->SetTitleSize(0.049);
   Phi2_sm->GetXaxis()->SetLabelOffset(0.01);
   Phi2_sm->GetYaxis()->SetTitle("Counts (normalized)");
   Phi2_sm->GetYaxis()->SetTitleSize(0.049);
   Phi2_sm->GetYaxis()->SetLabelOffset(0.01);
   Phi2_sm->GetXaxis()->SetNdivisions(9);
   Phi2_sm->Scale(1./Phi2_sm->Integral());
   Phi2_sm->Draw("PLC");
   Phi2_inp->Scale(1./Phi2_inp->Integral());
   Phi2_inp->Draw("SAME PLC");

   SmAndInp2->cd(2);
   Th2_sm->GetXaxis()->SetTitle("Polar angle #theta [deg]");
   Th2_sm->GetXaxis()->SetTitleSize(0.049);
   Th2_sm->GetXaxis()->SetLabelOffset(0.01);
   Th2_sm->GetYaxis()->SetTitle("Counts (normalized)");
   Th2_sm->GetYaxis()->SetTitleSize(0.049);
   Th2_sm->GetYaxis()->SetLabelOffset(0.01);
   Th2_sm->Scale(1./Th2_sm->Integral());
   Th2_sm->Draw("PLC");
   Th2_inp->Scale(1./Th2_inp->Integral());
   Th2_inp->Draw("SAME PLC");

   SmAndInp2->cd(3);
   En2_sm->GetXaxis()->SetTitle("Energy [MeV]");
   En2_sm->GetXaxis()->SetTitleSize(0.049);
   En2_sm->GetXaxis()->SetLabelOffset(0.01);
   En2_sm->GetYaxis()->SetTitle("Counts (normalized)");
   En2_sm->GetYaxis()->SetTitleSize(0.049);
   En2_sm->GetYaxis()->SetLabelOffset(0.01);
   En2_sm->Scale(1./En2_sm->Integral());
   En2_sm->Draw("PLC");
   En2_inp->Scale(1./En2_inp->Integral());
   En2_inp->Draw("SAME PLC");

   SmAndInp2->Write();
   SmAndInp2->SaveAs("SmAndInp2.png");
   SmAndInp2->SaveAs("SmAndInp2.root");

   fout->Close();
}
