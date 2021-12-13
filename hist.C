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


   TH1F* Res_spec_En1 = new TH1F("Res_spec_En1","Res_spec_En1",200,-50,50);
   TH1F* Res_spec_En2 = new TH1F("Res_spec_En2","Res_spec_En2",200,-50,50);
   TH1F* Res_spec_En3 = new TH1F("Res_spec_En3","Res_spec_En3",200,-45,45);

   TH1F* QA_Res_En1 = new TH1F("QA_Res_En1","QA_Res_En1",200,-50,50);
   TH1F* QA_Res_En2 = new TH1F("QA_Res_En2","QA_Res_En2",200,-50,50);
   TH1F* QA_Res_En3 = new TH1F("QA_Res_En3","QA_Res_En3",200,-120,120);

   TH1F* Res_spec_Th1 = new TH1F("Res_spec_Th1","Res_spec_Th1",200,-20,20);
   TH1F* Res_spec_Th2 = new TH1F("Res_spec_Th2","Res_spec_Th2",200,-20,20);
   TH1F* Res_spec_Th3 = new TH1F("Res_spec_Th3","Res_spec_Th3",200,-15,15);

   TH1F* QA_Res_Th1 = new TH1F("QA_Res_Th1","QA_Res_Th1",200,-20,20);
   TH1F* QA_Res_Th2 = new TH1F("QA_Res_Th2","QA_Res_Th2",200,-20,20);
   TH1F* QA_Res_Th3 = new TH1F("QA_Res_Th3","QA_Res_Th3",200,-50,50);

   TH1F* Res_spec_Ph1 = new TH1F("Res_spec_Ph1","Res_spec_Ph1",200,-20,20);
   TH1F* Res_spec_Ph2 = new TH1F("Res_spec_Ph2","Res_spec_Ph2",200,-20,20);
   TH1F* Res_spec_Ph3 = new TH1F("Res_spec_Ph3","Res_spec_Ph3",200,-15,15);

   TH1F* QA_Res_Ph1 = new TH1F("QA_Res_Ph1","QA_Res_Ph1",200,-20,20);
   TH1F* QA_Res_Ph2 = new TH1F("QA_Res_Ph2","QA_Res_Ph2",200,-20,20);
   TH1F* QA_Res_Ph3 = new TH1F("QA_Res_Ph3","QA_Res_Ph3",200,-40,40);
   
   TH1F* Pull_spec_En1 = new TH1F("Pull_spec_En1","Pull_spec_En1",200,-30,30);
   TH1F* Pull_spec_En2 = new TH1F("Pull_spec_En2","Pull_spec_En2",200,-30,30);
   TH1F* Pull_spec_En3 = new TH1F("Pull_spec_En3","Pull_spec_En3",200,-30,30);

   TH1F* Pull_spec_Th1 = new TH1F("Pull_spec_Th1","Pull_spec_Th1",200,-30,30);
   TH1F* Pull_spec_Th2 = new TH1F("Pull_spec_Th2","Pull_spec_Th2",200,-30,30);
   TH1F* Pull_spec_Th3 = new TH1F("Pull_spec_Th3","Pull_spec_Th3",200,-30,30);

   TH1F* Pull_spec_Ph1 = new TH1F("Pull_spec_Ph1","Pull_spec_Ph1",200,-30,30);
   TH1F* Pull_spec_Ph2 = new TH1F("Pull_spec_Ph2","Pull_spec_Ph2",200,-30,30);
   TH1F* Pull_spec_Ph3 = new TH1F("Pull_spec_Ph3","Pull_spec_Ph3",200,-30,30);

   TH1F* MMProt2Neut = new TH1F("MMProt2Neut","MMProt2Neut",200,900,942);
   TH1F* MMProt2NeutPre = new TH1F("MMProt2NeutPre","MMProt2NeutPre",200,900,1000);

   TH1F* En1_inp = new TH1F("En1_inp","En1_inp",200,0,150);
   TH1F* En1_out = new TH1F("En1_out", "En1_out",200,0,150);

   TH1F* Theta1=new TH1F("Theta1","Theta1",200,0,100);

   TH1F* hConverged = new TH1F("hConverged","hConverged",200,-10,10);

   Float_t sigma_En1, sigma_En2, sigma_En3;
   Float_t sigma_Th1, sigma_Th2, sigma_Th3;
   Float_t sigma_Ph1, sigma_Ph2, sigma_Ph3;

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
      if (Converged==1){
      cout<<Theta_1<<" "<<Energy_1<<" "<<Phi_1<<" "<<Converged<<endl;
      // chi2->Fill(Chi2);
      // prob->Fill(Probability);
      // prob->Fill(TMath::Prob(Chi2,4));

      Res_spec_En1->Fill(Energy_1_sm-Energy_1);
      Res_spec_En2->Fill(Energy_2_sm-Energy_2);
      Res_spec_En3->Fill(Energy_3_sm-Energy_3);

      // QA_Res_En1->Fill(Energy_1_sm-Energy_1_inp);
      // QA_Res_En2->Fill(Energy_2_sm-Energy_2_inp);
      // QA_Res_En3->Fill(Energy_3_sm-Energy_3_inp);

      QA_Res_En1->Fill(Energy_1_inp-Energy_1_sm);
      QA_Res_En2->Fill(Energy_2_inp-Energy_2_sm);
      QA_Res_En3->Fill(Energy_3_inp-Energy_3_sm);
      
      Res_spec_Th1->Fill(Theta_1_sm-Theta_1);
      Res_spec_Th2->Fill(Theta_2_sm-Theta_2);
      Res_spec_Th3->Fill(Theta_3_sm-Theta_3);

      QA_Res_Th1->Fill(Theta_1_sm-Theta_1_inp);
      QA_Res_Th2->Fill(Theta_2_sm-Theta_2_inp);
      QA_Res_Th3->Fill(Theta_3_sm-Theta_3_inp);

      Res_spec_Ph1->Fill(Phi_1_sm-Phi_1);
      Res_spec_Ph2->Fill(Phi_2_sm-Phi_2);
      Res_spec_Ph3->Fill(Phi_3_sm-Phi_3);

      QA_Res_Ph1->Fill(Phi_1_sm-Phi_1_inp);
      QA_Res_Ph2->Fill(Phi_2_sm-Phi_2_inp);
      QA_Res_Ph3->Fill(Phi_3_sm-Phi_3_inp);

      MMProt2Neut->Fill(MissMassProt2Neut);
      MMProt2NeutPre->Fill(MissMassProt2NeutPre);

      En1_inp->Fill(Energy_1_inp);
      En1_out->Fill(Energy_1*1000);

      Theta1->Fill(Theta_1*Rad2Deg);
      hConverged->Fill(Converged);
   }
}

   sigma_En1=Res_spec_En1->GetStdDev(1);
   sigma_En2=Res_spec_En2->GetStdDev(1);
   sigma_En3=Res_spec_En3->GetStdDev(1);

   sigma_Th1=Res_spec_Th1->GetStdDev(1);
   sigma_Th2=Res_spec_Th2->GetStdDev(1);
   sigma_Th3=Res_spec_Th3->GetStdDev(1);

   sigma_Ph1=Res_spec_Ph1->GetStdDev(1);
   sigma_Ph2=Res_spec_Ph2->GetStdDev(1);
   sigma_Ph3=Res_spec_Ph3->GetStdDev(1);

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // if (Converged==0) continue;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // cout<<Theta_1<<" "<<Energy_1<<" "<<Phi_1<<endl;

      Pull_spec_En1->Fill((Energy_1_sm*1000-Energy_1*1000)/sigma_En1);
      Pull_spec_En2->Fill((Energy_2_sm*1000-Energy_2*1000)/sigma_En2);
      Pull_spec_En3->Fill((Energy_3_sm*1000-Energy_3*1000)/sigma_En3);

      Pull_spec_Th1->Fill((Theta_1_sm*Rad2Deg-Theta_1*Rad2Deg)/sigma_Th1);
      Pull_spec_Th2->Fill((Theta_2_sm*Rad2Deg-Theta_2*Rad2Deg)/sigma_Th2);
      Pull_spec_Th3->Fill((Theta_3_sm*Rad2Deg-Theta_3*Rad2Deg)/sigma_Th3);

      Pull_spec_Ph1->Fill((Phi_1_sm*Rad2Deg-Phi_1*Rad2Deg)/sigma_Ph1);
      Pull_spec_Ph2->Fill((Phi_2_sm*Rad2Deg-Phi_2*Rad2Deg)/sigma_Ph2);
      Pull_spec_Ph3->Fill((Phi_3_sm*Rad2Deg-Phi_3*Rad2Deg)/sigma_Ph3);

   }

   // Pull_spec_En1=Res_spec_En1;
   // Pull_spec_En1->Scale(sigma_En1);
   // Pull_spec_En1->SetTitle("Pull_spec_En1");
   // Pull_spec_En2=Res_spec_En2;
   // Pull_spec_En2->Scale(sigma_En2);
   // Pull_spec_En2->SetTitle("Pull_spec_En2");
   // Pull_spec_En3=Res_spec_En3;
   // Pull_spec_En3->Scale(sigma_En3);
   // Pull_spec_En3->SetTitle("Pull_spec_En3");

   // Pull_spec_Th1=Res_spec_Th1;
   // Pull_spec_Th1->Scale(sigma_Th1);
   // Pull_spec_Th1->SetTitle("Pull_spec_Th1");
   // Pull_spec_Th2=Res_spec_Th2;
   // Pull_spec_Th2->Scale(sigma_Th2);
   // Pull_spec_Th2->SetTitle("Pull_spec_Th2");
   // Pull_spec_Th3=Res_spec_Th3;
   // Pull_spec_Th3->Scale(sigma_Th3);
   // Pull_spec_Th3->SetTitle("Pull_spec_Th3");

   // Pull_spec_Ph1=Res_spec_Ph1;
   // Pull_spec_Ph1->Scale(sigma_Ph1);
   // Pull_spec_Ph2=Res_spec_Ph2;
   // Pull_spec_Ph2->Scale(sigma_Ph2);
   // Pull_spec_Ph3=Res_spec_Ph3;
   // Pull_spec_Ph3->Scale(sigma_Ph3);

   fout->cd();

   // chi2->Write();
   // prob->Write();

   Res_spec_En1->Write();
   Res_spec_En2->Write();
   Res_spec_En3->Write();

   QA_Res_En1->Write();
   QA_Res_En2->Write();
   QA_Res_En3->Write();

   Res_spec_Th1->Write();
   Res_spec_Th2->Write();
   Res_spec_Th3->Write();

   QA_Res_Th1->Write();
   QA_Res_Th2->Write();
   QA_Res_Th3->Write();

   Res_spec_Ph1->Write();
   Res_spec_Ph2->Write();
   Res_spec_Ph3->Write();

   QA_Res_Ph1->Write();
   QA_Res_Ph2->Write();
   QA_Res_Ph3->Write();

   Pull_spec_En1->Write();
   Pull_spec_En2->Write();
   Pull_spec_En3->Write();

   Pull_spec_Th1->Write();
   Pull_spec_Th2->Write();
   Pull_spec_Th3->Write();

   Pull_spec_Ph1->Write();
   Pull_spec_Ph2->Write();
   Pull_spec_Ph3->Write();

   // MMProt2Neut->Write();
   // MMProt2NeutPre->Write();

   En1_inp->Write();
   En1_out->Write();
      
   Theta1->Write();

   hConverged->Write();

   TCanvas *Residuals = new TCanvas("Residuals","Residuals",2000,2000);
   Residuals->Divide(3,3);
   TCanvas *QA_res = new TCanvas("QA_res","QA_res",2000,2000);
   QA_res->Divide(3,3);
   TCanvas *Pull = new TCanvas("Pull","Pull",2000,2000);
   Pull->Divide(3,3);
   // TCanvas *Various = new TCanvas("Various","Various",1000,1000);
   // Various->Divide(2,2);

   Residuals->cd(1);
   Res_spec_En1->Draw();
   Residuals->cd(2);
   Res_spec_En2->Draw();
   Residuals->cd(3);
   Res_spec_En3->Draw();

   Residuals->cd(4);
   Res_spec_Th1->Draw();
   Residuals->cd(5);
   Res_spec_Th2->Draw();
   Residuals->cd(6);
   Res_spec_Th3->Draw();

   Residuals->cd(7);
   Res_spec_Ph1->Draw();
   Residuals->cd(8);
   Res_spec_Ph2->Draw();
   Residuals->cd(9);
   Res_spec_Ph3->Draw();

   Residuals->Write();
   Residuals->SaveAs("Residuals.png");

   QA_res->cd(1);
   QA_Res_En1->Draw();
   QA_res->cd(2);
   QA_Res_En2->Draw();
   QA_res->cd(3);
   QA_Res_En3->Draw();

   QA_res->cd(4);
   QA_Res_Th1->Draw();
   QA_res->cd(5);
   QA_Res_Th2->Draw();
   QA_res->cd(6);
   QA_Res_Th3->Draw();

   QA_res->cd(7);
   QA_Res_Ph1->Draw();
   QA_res->cd(8);
   QA_Res_Ph2->Draw();
   QA_res->cd(9);
   QA_Res_Ph3->Draw();

   QA_res->Write();
   QA_res->SaveAs("QA_res.png");

   Pull->cd(1);
   Pull_spec_En1->Draw();
   Pull->cd(2);
   Pull_spec_En2->Draw();
   Pull->cd(3);
   Pull_spec_En3->Draw();

   Pull->cd(4);
   Pull_spec_Th1->Draw();
   Pull->cd(5);
   Pull_spec_Th2->Draw();
   Pull->cd(6);
   Pull_spec_Th3->Draw();

   Pull->cd(7);
   Pull_spec_Ph1->Draw();
   Pull->cd(8);
   Pull_spec_Ph2->Draw();
   Pull->cd(9);
   Pull_spec_Ph3->Draw();

   Pull->Write();
   Pull->SaveAs("Pull.png");

   // Various->cd(1);
   // chi2->Draw();
   // Various->cd(2);
   // prob->Draw();
   // Various->cd(3);
   // MMProt2Neut->Draw();
   // Various->cd(4);
   // MMProt2NeutPre->Draw();

   // Various->Write();
   // Various->SaveAs("Various.png");



   fout->Close();
}
