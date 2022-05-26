#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TFile.h"


void macro()
{
	const int n = 600;
	int i = 0;
	char name[n];
	char name1[n];
	float lr;
	int deltachi;

	TFile *f;

	double MMMean, MMStdDev, ProbMean, ProbStdDev;

	TCanvas *c1 = new TCanvas("c1","Total", 1000, 1000);
	TH1F *h = new TH1F("h","Total histogram", n, 0, n);
	h->SetMinimum(937.6);
	// h->Sumw2(kFALSE);

	for(lr=0.01;lr<=0.30;lr+=0.01)
	{
		for(deltachi=1;deltachi<=100;deltachi+=5)
		{
			sprintf(name,"out_%.2f_%d/Refit.root", lr, deltachi);
			// sprintf(name1, "out_%.2f_%d", lr, deltachi);
			f = new TFile(name, "READ");
			TH1F *hist;
			hist = (TH1F*)f->Get("hTotal");
			MMMean = hist->GetBinContent(1);
			MMStdDev = hist->GetBinContent(2);
			ProbMean = hist->GetBinContent(3);
			ProbStdDev = hist->GetBinContent(4);

			// cout<<"MMMean: "<<MMMean<<endl;
			h->SetBinContent(i + 1, MMStdDev);
		    // h->Fill(2, MMStdDev);
		    // h->Fill(3, ProbMean);
		    // h->Fill(4, ProbStdDev);
		    // if(i%10==0)
		    // {
		    // 	h->GetXaxis()->SetBinLabel(i + 1, name1);
		    // }
			
			i++;
		}

		sprintf(name1, "lr=%.2f", lr);
		h->GetXaxis()->SetBinLabel(i + 1, name1);
	}

	h->Draw("colz");
	c1->Write();
}





















// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 	for (int i=0;i<18;i++){
 //   for (int j=0;j<5;j++){
 //     int k=j*5+20;
 //     sprintf(name,"hresol1_pad%d_th%d_fast",i,k);
 //     hresol1_pad_th[i][j]= (TH1F*)file->Get(name);
 //   }
 // }
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!