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
   TFile *fout = new TFile("hist.root","RECREATE");

   TH1F* chi2 = new TH1F("chi2","chi2",200,0,20);
   TH1F* prob = new TH1F("prob","prob",200,0,1);


   TH1F* Res_spec_En1 = new TH1F("Res_spec_En1","Res_spec_En1",200,-8,8);
   TH1F* Res_spec_En2 = new TH1F("Res_spec_En2","Res_spec_En2",200,-8,8);
   // TH1F* Res_spec_En3 = new TH1F("Res_spec_En3","Res_spec_En3",200,-45,45);

   TH1F* QA_Res_En1 = new TH1F("QA_Res_En1","QA_Res_En1",200,-18,18);
   TH1F* QA_Res_En2 = new TH1F("QA_Res_En2","QA_Res_En2",200,-18,18);
   // TH1F* QA_Res_En3 = new TH1F("QA_Res_En3","QA_Res_En3",200,-120,120);

   TH1F* Res_spec_Th1 = new TH1F("Res_spec_Th1","Res_spec_Th1",200,-0.5,0.5);
   TH1F* Res_spec_Th2 = new TH1F("Res_spec_Th2","Res_spec_Th2",200,-0.5,0.5);
   // TH1F* Res_spec_Th3 = new TH1F("Res_spec_Th3","Res_spec_Th3",200,-15,15);

   TH1F* QA_Res_Th1 = new TH1F("QA_Res_Th1","QA_Res_Th1",200,-5,5);
   TH1F* QA_Res_Th2 = new TH1F("QA_Res_Th2","QA_Res_Th2",200,-5,5);
   // TH1F* QA_Res_Th3 = new TH1F("QA_Res_Th3","QA_Res_Th3",200,-50,50);

   TH1F* Res_spec_Ph1 = new TH1F("Res_spec_Ph1","Res_spec_Ph1",200,-0.25,0.25);
   TH1F* Res_spec_Ph2 = new TH1F("Res_spec_Ph2","Res_spec_Ph2",200,-0.25,0.25);
   // TH1F* Res_spec_Ph3 = new TH1F("Res_spec_Ph3","Res_spec_Ph3",200,-15,15);

   TH1F* QA_Res_Ph1 = new TH1F("QA_Res_Ph1","QA_Res_Ph1",200,-5,5);
   TH1F* QA_Res_Ph2 = new TH1F("QA_Res_Ph2","QA_Res_Ph2",200,-5,5);
   // TH1F* QA_Res_Ph3 = new TH1F("QA_Res_Ph3","QA_Res_Ph3",200,-40,40);
   
   TH1F* Pull_spec_En1 = new TH1F("Pull_spec_En1","Pull_spec_En1",200,-8,8);
   TH1F* Pull_spec_En2 = new TH1F("Pull_spec_En2","Pull_spec_En2",200,-8,8);
   // TH1F* Pull_spec_En3 = new TH1F("Pull_spec_En3","Pull_spec_En3",200,-30,30);

   TH1F* Pull_spec_Th1 = new TH1F("Pull_spec_Th1","Pull_spec_Th1",200,-5,5);
   TH1F* Pull_spec_Th2 = new TH1F("Pull_spec_Th2","Pull_spec_Th2",200,-5,5);
   // TH1F* Pull_spec_Th3 = new TH1F("Pull_spec_Th3","Pull_spec_Th3",200,-30,30);

   TH1F* Pull_spec_Ph1 = new TH1F("Pull_spec_Ph1","Pull_spec_Ph1",200,-5,5);
   TH1F* Pull_spec_Ph2 = new TH1F("Pull_spec_Ph2","Pull_spec_Ph2",200,-5,5);
   // TH1F* Pull_spec_Ph3 = new TH1F("Pull_spec_Ph3","Pull_spec_Ph3",200,-30,30);

   TH1F* MMProt2Neut = new TH1F("MMProt2Neut","MMProt2Neut",200,900,942);
   TH1F* MMProt2NeutPre = new TH1F("MMProt2NeutPre","MMProt2NeutPre",200,900,1000);

   TH1F* En1_inp = new TH1F("En1_inp","En1_inp",200,0,150);
   TH1F* En1_cond = new TH1F("En1_cond", "En1_cond",200,0,150);

   TH1F* Theta1=new TH1F("Theta1","Theta1",200,0,100);

   TH1F* hConverged = new TH1F("hConverged","hConverged",200,-10,10);

   TH2F* En1_smvsEn1_f = new TH2F("En1_smvsEn1_f","En1_smvsEn1_f",200,0,200,200,0,200);

   TH1F *Th1_sm=new TH1F("Th1_sm","Th1_sm",200,0,100);
   TH1F *Phi1_sm=new TH1F("Phi1_sm","Phi1_sm",200,-200,200);
   TH1F *En1_sm=new TH1F("En1_sm","En1_sm",200,-20,200);

   TH1F *Th1_fit=new TH1F("Th1_fit","Th1_fit",200,0,100);
   TH1F *Phi1_fit=new TH1F("Phi1_fit","Phi1_fit",200,-200,200);
   TH1F *En1_fit=new TH1F("En1_fit","En1_fit",200,0,200);

   TH1F *Th1_smC=new TH1F("Th1_smC","Th1_smC",200,0,100);
   TH1F *Phi1_smC=new TH1F("Phi1_smC","Phi1_smC",200,-200,200);
   TH1F *En1_smC=new TH1F("En1_smC","En1_smC",200,-20,200);

   TH1F *Th1_fitC=new TH1F("Th1_fitC","Th1_fitC",200,0,100);
   TH1F *Phi1_fitC=new TH1F("Phi1_fitC","Phi1_fitC",200,-200,200);
   TH1F *En1_fitC=new TH1F("En1_fitC","En1_fitC",200,0,200);

   TH1F *Th2_sm=new TH1F("Th2_sm","Th2_sm",200,0,100);
   TH1F *Phi2_sm=new TH1F("Phi2_sm","Phi2_sm",200,-200,200);
   TH1F *En2_sm=new TH1F("En2_sm","En2_sm",200,-20,200);

   TH1F *Th2_fit=new TH1F("Th2_fit","Th2_fit",200,0,100);
   TH1F *Phi2_fit=new TH1F("Phi2_fit","Phi2_fit",200,-200,200);
   TH1F *En2_fit=new TH1F("En2_fit","En2_fit",200,0,200);

   TH1F *Th2_smC=new TH1F("Th2_smC","Th2_smC",200,0,100);
   TH1F *Phi2_smC=new TH1F("Phi2_smC","Phi2_smC",200,-200,200);
   TH1F *En2_smC=new TH1F("En2_smC","En2_smC",200,-20,200);

   TH1F *Th2_fitC=new TH1F("Th2_fitC","Th2_fitC",200,0,100);
   TH1F *Phi2_fitC=new TH1F("Phi2_fitC","Phi2_fitC",200,-200,200);
   TH1F *En2_fitC=new TH1F("En2_fitC","En2_fitC",200,0,200);

   TH1F *Th3_sm=new TH1F("Th3_sm","Th3_sm",200,0,100);
   TH1F *Phi3_sm=new TH1F("Phi3_sm","Phi3_sm",200,-200,200);
   TH1F *En3_sm=new TH1F("En3_sm","En3_sm",200,0,200);

   TH1F *Th3_fit=new TH1F("Th3_fit","Th3_fit",200,0,100);
   TH1F *Phi3_fit=new TH1F("Phi3_fit","Phi3_fit",200,-200,200);
   TH1F *En3_fit=new TH1F("En3_fit","En3_fit",200,0,200);

   TH1F *Th3_smC=new TH1F("Th3_smC","Th3_smC",200,0,100);
   TH1F *Phi3_smC=new TH1F("Phi3_smC","Phi3_smC",200,-200,200);
   TH1F *En3_smC=new TH1F("En3_smC","En3_smC",200,0,200);

   TH1F *Th3_fitC=new TH1F("Th3_fitC","Th3_fitC",200,0,100);
   TH1F *Phi3_fitC=new TH1F("Phi3_fitC","Phi3_fitC",200,-200,200);
   TH1F *En3_fitC=new TH1F("En3_fitC","En3_fitC",200,0,200);

   TH1F *Pull_spec_En1prim = new TH1F("Pull_spec_En1prim","Pull_spec_En1prim",200,-10,10);
   TH1F *Pull_spec_Th1prim = new TH1F("Pull_spec_Th1prim","Pull_spec_Th1prim",200,-15,15);
   TH1F *Pull_spec_Ph1prim = new TH1F("Pull_spec_Ph1prim","Pull_spec_Ph1prim",200,-30,30);

   // Float_t sigma_En1, sigma_En2, sigma_En3;
   // Float_t sigma_Th1, sigma_Th2, sigma_Th3;
   // Float_t sigma_Ph1, sigma_Ph2, sigma_Ph3;

   Float_t sigma_En1, sigma_En2;
   Float_t sigma_Th1, sigma_Th2;
   Float_t sigma_Ph1, sigma_Ph2;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // cout<<Converged<<endl;
      // if (Converged==0) continue;
      
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if(Probability > 0.02)
      // {
      // cout<<Theta_1<<" "<<Energy_1<<" "<<Phi_1<<" "<<Converged<<endl;

      Res_spec_En1->Fill(Energy_1_sm-Energy_1);
      Res_spec_En2->Fill(Energy_2_sm-Energy_2);
      // Res_spec_En3->Fill(Energy_3_sm-Energy_3);

      // QA_Res_En1->Fill(Energy_1_sm-Energy_1_inp);
      // QA_Res_En2->Fill(Energy_2_sm-Energy_2_inp);
      // QA_Res_En3->Fill(Energy_3_sm-Energy_3_inp);

      QA_Res_En1->Fill(Energy_1_inp-Energy_1_sm);
      QA_Res_En2->Fill(Energy_2_inp-Energy_2_sm);
      // QA_Res_En3->Fill(Energy_3_inp-Energy_3_sm);
      
      Res_spec_Th1->Fill(Theta_1_sm-Theta_1);
      Res_spec_Th2->Fill(Theta_2_sm-Theta_2);
      // Res_spec_Th3->Fill(Theta_3_sm-Theta_3);

      QA_Res_Th1->Fill(Theta_1_sm-Theta_1_inp);
      QA_Res_Th2->Fill(Theta_2_sm-Theta_2_inp);
      // QA_Res_Th3->Fill(Theta_3_sm-Theta_3_inp);

      Res_spec_Ph1->Fill(Phi_1_sm-Phi_1);
      Res_spec_Ph2->Fill(Phi_2_sm-Phi_2);
      // Res_spec_Ph3->Fill(Phi_3_sm-Phi_3);

      QA_Res_Ph1->Fill(Phi_1_sm-Phi_1_inp);
      QA_Res_Ph2->Fill(Phi_2_sm-Phi_2_inp);
      // QA_Res_Ph3->Fill(Phi_3_sm-Phi_3_inp);

      MMProt2Neut->Fill(MissMassProt2Neut);
      MMProt2NeutPre->Fill(MissMassProt2NeutPre);

      En1_inp->Fill(Energy_1_inp);
      En1_cond->Fill(Energy_1);

      Theta1->Fill(Theta_1*Rad2Deg);
      hConverged->Fill(Converged);

      En1_smvsEn1_f->Fill(Energy_1_sm,Energy_1);

      Th1_smC->Fill(Theta_1_sm);
      Phi1_smC->Fill(Phi_1_sm);
      En1_smC->Fill(Energy_1_sm);

      Th1_fitC->Fill(Theta_1);
      Phi1_fitC->Fill(Phi_1);
      En1_fitC->Fill(Energy_1);

      Th2_smC->Fill(Theta_2_sm);
      Phi2_smC->Fill(Phi_2_sm);
      En2_smC->Fill(Energy_2_sm);

      Th2_fitC->Fill(Theta_2);
      Phi2_fitC->Fill(Phi_2);
      En2_fitC->Fill(Energy_2);

      Th3_smC->Fill(Theta_3_sm);
      Phi3_smC->Fill(Phi_3_sm);
      En3_smC->Fill(Energy_3_sm);

      Th3_fitC->Fill(Theta_3);
      Phi3_fitC->Fill(Phi_3);
      En3_fitC->Fill(Energy_3);
   // }
      Th1_sm->Fill(Theta_1_sm);
      Phi1_sm->Fill(Phi_1_sm);
      En1_sm->Fill(Energy_1_sm);

      Th1_fit->Fill(Theta_1);
      Phi1_fit->Fill(Phi_1);
      En1_fit->Fill(Energy_1);

      Th2_sm->Fill(Theta_2_sm);
      Phi2_sm->Fill(Phi_2_sm);
      En2_sm->Fill(Energy_2_sm);

      Th2_fit->Fill(Theta_2);
      Phi2_fit->Fill(Phi_2);
      En2_fit->Fill(Energy_2);

      Th3_sm->Fill(Theta_3_sm);
      Phi3_sm->Fill(Phi_3_sm);
      En3_sm->Fill(Energy_3_sm);

      Th3_fit->Fill(Theta_3);
      Phi3_fit->Fill(Phi_3);
      En3_fit->Fill(Energy_3);
}

   sigma_En1=Res_spec_En1->GetStdDev(1);
   sigma_En2=Res_spec_En2->GetStdDev(1);
   // sigma_En3=Res_spec_En3->GetStdDev(1);

   sigma_Th1=Res_spec_Th1->GetStdDev(1);
   sigma_Th2=Res_spec_Th2->GetStdDev(1);
   // sigma_Th3=Res_spec_Th3->GetStdDev(1);

   sigma_Ph1=Res_spec_Ph1->GetStdDev(1);
   sigma_Ph2=Res_spec_Ph2->GetStdDev(1);
   // sigma_Ph3=Res_spec_Ph3->GetStdDev(1);

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // if (Converged==0) continue;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // cout<<Theta_1<<" "<<Energy_1<<" "<<Phi_1<<endl;
      // if(Converged==1 && Probability>0.02)
      // if(Probability>0.02)
      // {      
      // Pull_spec_En1->Fill((Energy_1_sm*1000-Energy_1*1000)/sigma_En1);
      Pull_spec_En1->Fill((Energy_1_sm-Energy_1)/sigma_En1);
      Pull_spec_En2->Fill((Energy_2_sm-Energy_2)/sigma_En2);
      // Pull_spec_En3->Fill((Energy_3_sm*1000-Energy_3*1000)/sigma_En3);

      Pull_spec_Th1->Fill((Theta_1_sm-Theta_1)/sigma_Th1);
      Pull_spec_Th2->Fill((Theta_2_sm-Theta_2)/sigma_Th2);
      // Pull_spec_Th3->Fill((Theta_3_sm*Rad2Deg-Theta_3*Rad2Deg)/sigma_Th3);

      Pull_spec_Ph1->Fill((Phi_1_sm-Phi_1)/sigma_Ph1);
      Pull_spec_Ph2->Fill((Phi_2_sm-Phi_2)/sigma_Ph2);
      // Pull_spec_Ph3->Fill((Phi_3_sm*Rad2Deg-Phi_3*Rad2Deg)/sigma_Ph3);
      // }

      // if(Converged==1 && Probability >0.1)
      // if(Probability >0.2)
      // {
         Pull_spec_En1prim->Fill((Energy_1_sm-Energy_1)/sigma_En1);
         Pull_spec_Th1prim->Fill((Theta_1_sm-Theta_1)/sigma_Th1);
         Pull_spec_Ph1prim->Fill((Phi_1_sm-Phi_1)/sigma_Ph1);
      // }
   }

   fout->cd();

   // chi2->Write();
   // prob->Write();

   Res_spec_En1->Write();
   Res_spec_En2->Write();
   // Res_spec_En3->Write();

   QA_Res_En1->Write();
   QA_Res_En2->Write();
   // QA_Res_En3->Write();

   Res_spec_Th1->Write();
   Res_spec_Th2->Write();
   // Res_spec_Th3->Write();

   QA_Res_Th1->Write();
   QA_Res_Th2->Write();
   // QA_Res_Th3->Write();

   Res_spec_Ph1->Write();
   Res_spec_Ph2->Write();
   // Res_spec_Ph3->Write();

   QA_Res_Ph1->Write();
   QA_Res_Ph2->Write();
   // QA_Res_Ph3->Write();

   Pull_spec_En1->Write();
   Pull_spec_En2->Write();
   // Pull_spec_En3->Write();

   Pull_spec_Th1->Write();
   Pull_spec_Th2->Write();
   // Pull_spec_Th3->Write();

   Pull_spec_Ph1->Write();
   Pull_spec_Ph2->Write();
   // Pull_spec_Ph3->Write();

   // MMProt2Neut->Write();
   // MMProt2NeutPre->Write();

   En1_inp->Write();
   En1_cond->Write();
      
   Theta1->Write();

   hConverged->Write();

   En1_smvsEn1_f->Write();

   Th1_sm->Write();
   Phi1_sm->Write();
   En1_sm->Write();

   Th1_fit->Write();
   Phi1_fit->Write();
   En1_fit->Write();

   Pull_spec_En1prim->Write();
   Pull_spec_Th1prim->Write();
   Pull_spec_Ph1prim->Write();


   TCanvas *Residuals = new TCanvas("Residuals","Residuals",3000,3000);
   Residuals->Divide(3,2);
   TCanvas *QA_res = new TCanvas("QA_res","QA_res",3000,3000);
   QA_res->Divide(3,2);
   TCanvas *Pull = new TCanvas("Pull","Pull",3000,3000);
   Pull->Divide(3,2);
   TCanvas *Compare = new TCanvas("Compare","Compare",3000,3000);
   Compare->Divide(3,2);
   TCanvas *Compare2 = new TCanvas("Compare2","Compare2",3000,3000);
   Compare2->Divide(3,2);
   TCanvas *Pull_Prob = new TCanvas("Pull_Prob","Pull_Prob",3000,3000);
   Pull_Prob->Divide(3,1);
   TCanvas *QAandRes1 = new TCanvas("QAandRes1","QAandRes1",3000,3000);
   QAandRes1->Divide(3,2);
   TCanvas *QAandRes2 = new TCanvas("QAandRes2","QAandRes2",3000,3000);
   QAandRes2->Divide(3,2);

   Residuals->cd(1);
   Res_spec_En1->SetLineWidth(2);
   Res_spec_En1->Draw();
   Residuals->cd(4);
   Res_spec_En2->SetLineWidth(2);
   Res_spec_En2->Draw();
   // Residuals->cd(3);
   // Res_spec_En3->Draw();

   Residuals->cd(2);
   Res_spec_Th1->SetLineWidth(2);
   Res_spec_Th1->Draw();
   Residuals->cd(5);
   Res_spec_Th2->SetLineWidth(2);
   Res_spec_Th2->Draw();
   // Residuals->cd(6);
   // Res_spec_Th3->Draw();

   Residuals->cd(3);
   Res_spec_Ph1->SetLineWidth(2);
   Res_spec_Ph1->Draw();
   Residuals->cd(6);
   Res_spec_Ph2->SetLineWidth(2);
   Res_spec_Ph2->Draw();
   // Residuals->cd(9);
   // Res_spec_Ph3->Draw();

   Residuals->Write();
   Residuals->SaveAs("Residuals.png");

   QA_res->cd(1);
   QA_Res_En1->SetLineWidth(2);
   QA_Res_En1->Draw();
   QA_res->cd(4);
   QA_Res_En2->SetLineWidth(2);
   QA_Res_En2->Draw();
   // QA_res->cd(3);
   // QA_Res_En3->Draw();

   QA_res->cd(2);
   QA_Res_Th1->SetLineWidth(2);
   QA_Res_Th1->Draw();
   QA_res->cd(5);
   QA_Res_Th2->SetLineWidth(2);
   QA_Res_Th2->Draw();
   // QA_res->cd(6);
   // QA_Res_Th3->Draw();

   QA_res->cd(3);
   QA_Res_Ph1->SetLineWidth(2);
   QA_Res_Ph1->Draw();
   QA_res->cd(6);
   QA_Res_Ph2->SetLineWidth(2);
   QA_Res_Ph2->Draw();
   // QA_res->cd(9);
   // QA_Res_Ph3->Draw();

   QA_res->Write();
   QA_res->SaveAs("QA_res.png");

   Pull->cd(1);
   Pull_spec_En1->SetLineWidth(2);
   Pull_spec_En1->Draw();
   Pull->cd(4);
   Pull_spec_En2->SetLineWidth(2);
   Pull_spec_En2->Draw();
   // Pull->cd(3);
   // Pull_spec_En3->Draw();

   Pull->cd(2);
   Pull_spec_Th1->SetLineWidth(2);
   Pull_spec_Th1->Draw();
   Pull->cd(5);
   Pull_spec_Th2->SetLineWidth(2);
   Pull_spec_Th2->Draw();
   // Pull->cd(6);
   // Pull_spec_Th3->Draw();

   Pull->cd(3);
   Pull_spec_Ph1->SetLineWidth(2);
   Pull_spec_Ph1->Draw();
   Pull->cd(6);
   Pull_spec_Ph2->SetLineWidth(2);
   Pull_spec_Ph2->Draw();
   // Pull->cd(9);
   // Pull_spec_Ph3->Draw();

   Pull->Write();
   Pull->SaveAs("Pull.png");


   gStyle->SetPalette(1);

   Compare->cd(1);
   Th1_sm->SetLineWidth(2);
   Th1_sm->Scale(1./Th1_sm->Integral());
   Th1_sm->Draw("PLC");

   Th1_fit->SetLineWidth(2);
   Th1_fit->Scale(1./Th1_fit->Integral());
   Th1_fit->Draw("SAME PLC");//first is violet

   Compare->cd(2);
   Phi1_sm->SetLineWidth(2);
   Phi1_sm->Scale(1./Phi1_sm->Integral());
   Phi1_sm->Draw("PLC");

   Phi1_fit->SetLineWidth(2);
   Phi1_fit->Scale(1./Phi1_fit->Integral());
   Phi1_fit->Draw("SAME PLC");

   Compare->cd(3);
   En1_sm->SetLineWidth(2);
   En1_sm->Scale(1./En1_sm->Integral());
   En1_sm->Draw("PLC");

   En1_fit->SetLineWidth(2);
   En1_fit->Scale(1./En1_fit->Integral());
   En1_fit->Draw("SAME PLC");

   Compare->cd(4);
   Th1_smC->SetLineWidth(2);
   Th1_smC->Scale(1./Th1_smC->Integral());
   Th1_smC->Draw("PLC");

   Th1_fitC->SetLineWidth(2);
   Th1_fitC->Scale(1./Th1_fitC->Integral());
   Th1_fitC->Draw("SAME PLC");

   Compare->cd(5);
   Phi1_smC->SetLineWidth(2);
   Phi1_smC->Scale(1./Phi1_smC->Integral());
   Phi1_smC->Draw("PLC");

   Phi1_fitC->SetLineWidth(2);
   Phi1_fitC->Scale(1./Phi1_fitC->Integral());
   Phi1_fitC->Draw("SAME PLC");

   Compare->cd(6);
   En1_smC->SetLineWidth(2);
   En1_smC->Scale(1./En1_smC->Integral());
   En1_smC->Draw("PLC");

   En1_fitC->SetLineWidth(2);
   En1_fitC->Scale(1./En1_fitC->Integral());
   En1_fitC->Draw("SAME PLC");

   Compare->Write();
   Compare->SaveAs("Compare.png");

   gStyle->SetPalette(1);

   Compare2->cd(1);
   Th2_sm->SetLineWidth(2);
   Th2_sm->Scale(1./Th2_sm->Integral());
   Th2_sm->Draw("PLC");

   Th2_fit->SetLineWidth(2);
   Th2_fit->Scale(1./Th2_fit->Integral());
   Th2_fit->Draw("SAME PLC");//first is violet

   Compare2->cd(2);
   Phi2_sm->SetLineWidth(2);
   Phi2_sm->Scale(1./Phi2_sm->Integral());
   Phi2_sm->Draw("PLC");

   Phi2_fit->SetLineWidth(2);
   Phi2_fit->Scale(1./Phi2_fit->Integral());
   Phi2_fit->Draw("SAME PLC");

   Compare2->cd(3);
   En2_sm->SetLineWidth(2);
   En2_sm->Scale(1./En2_sm->Integral());
   En2_sm->Draw("PLC");

   En2_fit->SetLineWidth(2);
   En2_fit->Scale(1./En2_fit->Integral());
   En2_fit->Draw("SAME PLC");

   Compare2->cd(4);
   Th2_smC->SetLineWidth(2);
   Th2_smC->Scale(1./Th2_smC->Integral());
   Th2_smC->Draw("PLC");

   Th2_fitC->SetLineWidth(2);
   Th2_fitC->Scale(1./Th2_fitC->Integral());
   Th2_fitC->Draw("SAME PLC");

   Compare2->cd(5);
   Phi2_smC->SetLineWidth(2);
   Phi2_smC->Scale(1./Phi2_smC->Integral());
   Phi2_smC->Draw("PLC");

   Phi2_fitC->SetLineWidth(2);
   Phi2_fitC->Scale(1./Phi2_fitC->Integral());
   Phi2_fitC->Draw("SAME PLC");

   Compare2->cd(6);
   En2_smC->SetLineWidth(2);
   En2_smC->Scale(1./En2_smC->Integral());
   En2_smC->Draw("PLC");

   En2_fitC->SetLineWidth(2);
   En2_fitC->Scale(1./En2_fitC->Integral());
   En2_fitC->Draw("SAME PLC");

   Compare2->Write();
   Compare2->SaveAs("Compare2.png");

   gStyle->SetPalette(1);
   Pull_Prob->cd(1); //first is violet
   Pull_spec_En1->SetLineWidth(2);
   Pull_spec_En1->Draw("PLC");
   Pull_spec_En1prim->SetLineWidth(2);
   Pull_spec_En1prim->Draw("SAME PLC");
   Pull_Prob->cd(2);
   Pull_spec_Th1->SetLineWidth(2);
   Pull_spec_Th1->Draw("PLC");
   Pull_spec_Th1prim->SetLineWidth(2);
   Pull_spec_Th1prim->Draw("SAME PLC");
   Pull_Prob->cd(3);
   Pull_spec_Ph1->SetLineWidth(2);
   Pull_spec_Ph1->Draw("PLC");
   Pull_spec_Ph1prim->SetLineWidth(2);
   Pull_spec_Ph1prim->Draw("SAME PLC");

   Pull_Prob->Write();
   Pull_Prob->SaveAs("Pull_Prob.png");

   QAandRes1->cd(1);
   QA_Res_Th1->Draw();
   QAandRes1->cd(2);
   QA_Res_Ph1->Draw();
   QAandRes1->cd(3);
   QA_Res_En1->Draw();
   QAandRes1->cd(4);
   Res_spec_Th1->Draw();
   QAandRes1->cd(5);
   Res_spec_Ph1->Draw();
   QAandRes1->cd(6);
   Res_spec_En1->Draw();

   QAandRes1->Write();
   QAandRes1->SaveAs("QAandRes1.png");

   QAandRes2->cd(1);
   QA_Res_Th2->Draw();
   QAandRes2->cd(2);
   QA_Res_Ph2->Draw();
   QAandRes2->cd(3);
   QA_Res_En2->Draw();
   QAandRes2->cd(4);
   Res_spec_Th2->Draw();
   QAandRes2->cd(5);
   Res_spec_Ph2->Draw();
   QAandRes2->cd(6);
   Res_spec_En2->Draw();

   QAandRes2->Write();
   QAandRes2->SaveAs("QAandRes2.png");

   fout->Close();
}
