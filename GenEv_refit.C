#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "time.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#define Pi_180 0.01745329251994
#define Rad2Deg 57.2957795
#define Deg2Rad 0.0174532925
int funny_sign_count=0;

TCanvas *Results = new TCanvas("Results","Results",1000,1000);
TCanvas *Various = new TCanvas("Various","Various",1000,1000);

// TMatrixD y = new TMatrixD(); // daughters * variables 
// TMatrixD V = new TMatrixD();
// y.ResizeTo(6,1);
// V.ResizeTo(6,6);
// y.Zero();
// V.Zero();

TMatrixD ytemp(6,1);
TMatrixD V(6,6);
Float_t mp = 938.272;
Float_t mn = 939.565;
Float_t md = 1875.612;
Int_t fNdf = 5; // 4+NumberOfConstraints
const TLorentzVector pinit(0.,0.,791.073,2973.88);
Float_t fChi2, fProb;
Int_t glob_event = 0.;
Int_t f_nan = 0;

void Counting(ULong64_t ent){
	if(ent%10000==0)
		{

	switch (funny_sign_count%4)
		{
	case 0:
		printf("  | \t");
		 break;
	case 1:
		printf("  / \t");
		break;
	case 2:
		printf("  --\t");
		break;
	case 3:
		printf("  \\ \t");
		break;
	}
		printf("Event_num=%lld\r",ent);
		fflush(stdout);
		funny_sign_count++;
	}
}


Float_t calc_MM(Float_t en1, Float_t en2, Float_t th1, Float_t th2, Float_t phi1, Float_t phi2, int cat){
	Float_t MMass=0;
  	Float_t Mp= 938.272;
  	Float_t Md=1875.612;
	Float_t Mn= 939.565;
	Float_t Ed=160.0;
	Float_t M1=0,M2=0,M3=0;
	if(cat==0){M1=Mp; M2=Mp; M3=Mn;} //Missing Mass dla pary proton-proton
	if(cat==1){M1=Mp, M2=Mn; M3=Mp;} //Missing Mass dla pary proton-neutron
	TVector3 p1,p2,p3,d;
	p1.SetMagThetaPhi(sqrt((en1*en1)+2*M1*en1), th1*Pi_180, phi1*Pi_180);
	p2.SetMagThetaPhi(sqrt((en2*en2)+2*M2*en2), th2*Pi_180, phi2*Pi_180);
	d.SetMagThetaPhi(sqrt((Ed*Ed)+2*Md*Ed), 0*Pi_180, 0*Pi_180);
	p3=d-(p1+p2);
	Float_t ER=160+Md-M3-en1-en2;
	MMass=sqrt(ER*ER - p3.Mag2());
	return MMass;
}

void show(TMatrixD matrix)
{
    for(int j=0;j<matrix.GetNcols();j++)
    {
    	for(int i=0;i<matrix.GetNrows();i++)
    	{
    		cout<<"("<<i<<", "<<j<<") = "<<matrix(i,j)<<endl;
    	}
    }

}

void SetValues(Float_t phi1, Float_t th1, Float_t en1, Float_t phi2, Float_t th2, Float_t en2)
{
	ytemp.Zero();
	ytemp(0,0)=phi1;
	ytemp(1,0)=th1;
	ytemp(2,0)=en1;
	ytemp(3,0)=phi2;
	ytemp(4,0)=th2;
	ytemp(5,0)=en2;

	// cout<<"Set Values "<<y(0,0)<<" "<<y(1,0)<<" "<<y(2,0)<<endl;
}

void SetCovariance(Float_t phi1_er, Float_t th1_er, Float_t en1_er, Float_t phi2_er, Float_t th2_er, Float_t en2_er)
{
	V.Zero();
	V(0,0) = phi1_er*phi1_er;
	V(1,1) = th1_er*th1_er;
	V(2,2) = en1_er*en1_er;
	V(3,3) = phi2_er*phi2_er;
	V(4,4) = th2_er*th2_er;
	V(5,5) = en2_er*en2_er;
	// cout<<"Set Covariance "<<V(0,0)<<" "<<V(1,1)<<" "<<V(2,2)<<endl;
}

TMatrixD f_eval(TMatrixD y)
{
	Int_t cov_dim = 3;
	TMatrixD d;
    d.ResizeTo(1, 1);
    Float_t P, Px, Py, Pz, E;
    Px = pinit.Px();
    Py = pinit.Py();
    Pz = pinit.Pz();
    E = pinit.E();
    // cout<<"pinit px: "<<pinit.Px()<<"pinit py: "<<pinit.Py()<<"pinit pz: "<<pinit.Pz()<<"pinit E: "<<pinit.E()<<endl;

    for (int q = 0; q < 2; q++) 
    {
        // E -= std::sqrt((1. / y(0 + q * cov_dim, 0)) *
        //                    (1. / y(0 + q * cov_dim, 0)) +
        //                mp * mp);
        // Px -= (1. / y(0 + q * cov_dim, 0)) *
        //       std::sin(y(1 + q * cov_dim, 0)) *
        //       std::cos(y(2 + q * cov_dim, 0));
        // Py -= (1. / y(0 + q * cov_dim, 0)) *
        //       std::sin(y(1 + q * cov_dim, 0)) *
        //       std::sin(y(2 + q * cov_dim, 0));
        // Pz -= (1. / y(0 + q * cov_dim, 0)) *
        //       std::cos(y(1 + q * cov_dim, 0));

       	E-=y(2+q*cov_dim,0)+mp;

        P=sqrt(pow(y(2+q*cov_dim,0)+mp,2)-mp*mp);
        Px-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*cos(y(0+q*cov_dim,0)*Deg2Rad);
		Py-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*sin(y(0+q*cov_dim,0)*Deg2Rad);
		Pz-=P*cos(y(1+q*cov_dim,0)*Deg2Rad);
		// cout<<"Ev: "<<glob_event<<" E: "<<E<<" P: "<<P<<" Px: "<<Px<<" Py: "<<Py<<" Pz: "<<Pz<<endl;
    }

    // cout<<"Ev: "<<glob_event<<" E: "<<E<<" P: "<<P<<" Px: "<<Px<<" Py: "<<Py<<" Pz: "<<Pz<<endl;

    d(0, 0) = std::pow(E, 2) - std::pow(Px, 2) - std::pow(Py, 2) - std::pow(Pz, 2) - mn * mn;
	return d;
}

TMatrixD Feta_eval(TMatrixD y)
{
	TMatrixD H;
    H.ResizeTo(1, 6);
    H.Zero();
    Int_t cov_dim=3;
    Float_t P, Px, Py, Pz, E;

    // Px += pinit.Px();
    // Py += pinit.Py();
    // Pz += pinit.Pz();
    // E += pinit.E();
    Px=0;
    Py=0;
    Pz=791.073;
    E=2973.88;


        for (int q = 0; q < 2; q++)
        {
            E-=y(2+q*cov_dim,0)+mp;

        	P=sqrt(pow(y(2+q*cov_dim,0)+mp,2)-mp*mp);
        	// P=sqrt(pow(E+mp,2)-mp*mp);
        	Px-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*cos(y(0+q*cov_dim,0)*Deg2Rad);
			Py-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*sin(y(0+q*cov_dim,0)*Deg2Rad);
			Pz-=P*cos(y(1+q*cov_dim,0)*Deg2Rad);
        }

        for (int q = 0; q < 2; q++)
        {
            double Pi = sqrt(pow(y(2+q*cov_dim,0)+mp,2)-mp*mp);
            double Ei = sqrt(Pi * Pi + mp * mp);

            // H(0, 2 + q * cov_dim) =
            //     2 * E * (Pi / Ei) -
            //     2 * std::sin(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::cos(y(0 + q * cov_dim, 0)*Deg2Rad) * Px -
            //     2 * std::sin(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::sin(y(0 + q * cov_dim, 0)*Deg2Rad) * Py +
            //     2 * std::cos(y(1 + q * cov_dim, 0)*Deg2Rad) * Pz;

            // H(0, 1 + q * cov_dim) =
            //     2 * Pi * std::cos(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::cos(y(0 + q * cov_dim, 0)*Deg2Rad) * Px +
            //     2 * Pi * std::cos(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::sin(y(0 + q * cov_dim, 0)*Deg2Rad) * Py -
            //     2 * Pi * std::sin(y(1 + q * cov_dim, 0)*Deg2Rad) * Pz;

            // H(0, 0 + q * cov_dim) =
            //     -2 * Pi * std::sin(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::sin(y(0 + q * cov_dim, 0)*Deg2Rad) * Px +
            //     2 * Pi * std::sin(y(1 + q * cov_dim, 0)*Deg2Rad) *
            //         std::cos(y(0 + q * cov_dim, 0)*Deg2Rad) * Py;

            H(0, 2 + q * cov_dim) = 
            	(-2 * Pi * E) / Ei -
            	2 * std::sin(y(1 + q * cov_dim, 0) * Deg2Rad) * std::cos(y(0 + 1 * cov_dim, 0)*Deg2Rad) * Px -
            	2 * std::sin(y(1 + q * cov_dim, 0) * Deg2Rad) * std::sin(y(0 + 1 * cov_dim, 0)*Deg2Rad) * Py +
            	2 * std::cos(y(1 + q * cov_dim, 0) * Deg2Rad) * Pz;

            H(0, 1 + q * cov_dim) =
            	0 -
            	Pi * std::cos(y(1 + q * cov_dim, 0) * Deg2Rad) * std::cos(y(0 + q * cov_dim, 0) * Deg2Rad) -
            	Pi * std::cos(y(1 + q * cov_dim, 0) * Deg2Rad) * std::sin(y(0 + q * cov_dim, 0) * Deg2Rad) -
            	Pi * std::sin(y(1 + q * cov_dim, 0) * Deg2Rad);

            H(0, 0 + q * cov_dim) =
            	0 +
            	Pi * std::sin(y(1 + q * cov_dim, 0) * Deg2Rad) * std::sin(y(0 + q * cov_dim, 0) * Deg2Rad) -
            	Pi * std::sin(y(1 + q * cov_dim, 0) * Deg2Rad) * std::cos(y(0 + q * cov_dim, 0) * Deg2Rad) -
            	0;

        }
    
    	// cout<<"Feta eval "<<H(0,0)<<" "<<H(0,1)<<" "<<H(0,2)<<endl;
 return H;

}


Int_t fit()
{
	double lr = 0.5;
    TMatrixD alpha0(fN * cov_dim, 1), alpha(fN * cov_dim, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    double chi2 = 1e6;
    TMatrixD D = Feta_eval(alpha);
    TMatrixD d = f_eval(alpha);

    for (int q = 0; q < 5; q++)
    {
        TMatrixD DT(D.GetNcols(), D.GetNrows());
        DT.Transpose(D);
        TMatrixD VD = D * V * DT;
        VD.Invert();

        TMatrixD delta_alpha = alpha - alpha0;
        TMatrixD lambda = VD * D * delta_alpha + VD * d;
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fN * cov_dim, 1);
        neu_alpha = alpha - lr * V * DT * lambda;

        double chisqrd = 0.;

        for (int p = 0; p < lambda.GetNrows(); p++)
        {
            chisqrd = lambdaT(0, p) * d(p, 0);
        }

        /* for checking convergence
        // three parameters are checked
        // 1. difference between measurements (for successive iterations) y
        // 2. difference between constraints (for successive iterations)  d
        // 3. difference between chi2 (for successive iterations)  chisqrd
        // check converge for 'y' measurements
        double sum0 = 0;
        for(int p=0; p<(fN*5); p++){
            sum0 += (neu_alpha(p,0)-alpha(p,0))*(neu_alpha(p,0)-alpha(p,0));
        }
        double d_const = fabs(d(0,0));
        if(fabs(chi2-chisqrd)<1e-3 && d_const<10 && sqrt(sum0)<1e-3){
            fIteration = q;
            fConverged = true;
            break;
        }
        */
        chi2 = chisqrd;
        alpha0 = alpha;
        alpha = neu_alpha;
        V = V - lr * V * DT * VD * D * V;
        D = Feta_eval(alpha);
        d = f_eval(alpha);
    }

    y = alpha;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    // cout<<"Ev: "<<glob_event<<" q: "<<q<<" Chi2 = "<<chi2<<" "<<" Prob = "<<TMath::Prob(chi2, fNdf)<<endl;

    // -----------------------------------------
    // Pull
    // -----------------------------------------
    // fPull.ResizeTo(fN * cov_dim, fN * cov_dim);
    // for (uint b = 0; b < (fN * cov_dim); b++)
    //     fPull(b, b) = -10000;

    // if (true)
    // {
    //     for (uint b = 0; b < (fN * cov_dim); b++)
    //     {
    //         double num = A0(b, 0) - alpha(b, 0);
    //         double dem = V0(b, b) - V(b, b);
    //         if (dem > 0) { fPull(b, b) = num / std::sqrt(dem); }
    //     }
    // }

    // updateDaughters();

    return fConverged; // for number of iterations greater than 1
    // return true; // for number of iterations equal to 1
}

TMatrixD GetFitVal(TMatrixD y)
{
	TMatrixD y_fit;
	TVector3 Convert;
	TLorentzVector *v;
	v = new TLorentzVector(0,0,0,0);
	y_fit.ResizeTo(6,1);
	Int_t cov_dim = 3;
	Int_t fN = 2;
	Float_t E = 0., P = 0., Px = 0., Py = 0., Pz = 0.;
	for(int val = 0; val < fN; val++)
	{
		E=y(2+val*cov_dim,0)+mp;

        P=sqrt(pow(y(2+val*cov_dim,0)+mp,2)-mp*mp);
        Px=P*sin(y(1+val*cov_dim,0)*Deg2Rad)*cos(y(0+val*cov_dim,0)*Deg2Rad);
		Py=P*sin(y(1+val*cov_dim,0)*Deg2Rad)*sin(y(0+val*cov_dim,0)*Deg2Rad);
		Pz=P*cos(y(1+val*cov_dim,0)*Deg2Rad);

    	// Convert.SetX(Px);
    	// Convert.SetY(Py);
    	// Convert.SetZ(Pz);

    	// y_fit(0 + val * cov_dim, 0) = Convert.Theta();
    	// y_fit(1 + val * cov_dim, 0) = Convert.Phi();
    	// y_fit(2 + val * cov_dim, 0) = Convert.Mag();

		v->SetPxPyPzE(Px,Py,Pz,E);

		y_fit(0 + val * cov_dim, 0) = v->Phi()*Rad2Deg;
		y_fit(1 + val * cov_dim, 0) = v->Theta()*Rad2Deg;
		y_fit(2 + val * cov_dim, 0) = v->E()-mp; //ODJAC MASE


	}

	return y_fit;
}

//   Float_t Mp= 938.272;
//   Float_t Md=1875.612;
//   Float_t Mn= 939.565;
//   Float_t Ed=160.0;
//   Float_t EpT=0.0;
//   TVector3 p1;
//   TVector3 n;
//   TVector3 d;
//   TVector3 pT;
//   TVector3 p2;
//   p1.SetMagThetaPhi(sqrt((Ep*Ep)+2*Mp*Ep), Th1[thp]*Pi_180, 0*Pi_180);
//   n.SetMagThetaPhi(sqrt((En*En)+2*Mn*En), Th2[thn]*Pi_180, Dphi[phi]*Pi_180);
//   d.SetMagThetaPhi(sqrt((Ed*Ed)+2*Md*Ed), 0, 0*Pi_180);
//   pT.SetMagThetaPhi(sqrt((EpT*EpT)+2*Mp*EpT), 0, 0*Pi_180);
//   p2=d-(p1+n);
//   Float_t ER=160+Md-Mp-En-Ep;
//   Float_t Mr=sqrt(ER*ER - p2.Mag2());
//
//
//   TLorentzVector pTl;
//   pTl.SetVectM(pT,Mp);
//   TLorentzVector pl;
//   pl.SetVectM(p2,Mp);
//   //dl*=0.5;
//
//     TLorentzVector p2l = pTl-pl;
// //cout<<"dl E="<<dl.E()<<" nl E="<<nl.E()<<" p2l E="<<p2l*p2l<<"\n";
//   Float_t Etr = -(p2l*p2l)/(2*Mp);
//
// //cout<<"pl E="<<pl.E()<<" pTl E="<<pTl.E()<<" p2l E="<<p2l*p2l<<" Mr="<<Mr<<" ER="<<ER<<" Etrp="<<Etr<<"\n";
//   return Etr;



int GenEv_refit(){
	clock_t fbegin;
	clock_t fend;
	fbegin = clock(); //Start clock

// TVector3 help;
// help.SetMagThetaPhi(sqrt((160*160)+2*1875.612*160), 0*Pi_180, 0*Pi_180);
// TLorentzVector pinit2(help, 1875.612+160.);
// TLorentzVector LVmp(0,0,0,mp);
// pinit2+=LVmp;
// TLorentzVector pinit(0,0,791.073,2973.88);

// cout<<"help px: "<<help.Px()<<"help py: "<<help.Py()<<"help pz: "<<help.Pz()<<endl;
// cout<<"pinit px: "<<pinit2.Px()<<"pinit py: "<<pinit2.Py()<<"pinit pz: "<<pinit2.Pz()<<"pinit E: "<<pinit2.E()<<endl;
// return 0;

// 	 TLorentzVector proj1, targ1, beam1;
//   TLorentzVector targC, beamC;
//   proj1.SetPxPyPzE(0,0,pion_mom,pion_en);
//   targ1.SetPxPyPzE(0,0,0,mp);
//   targC.SetPxPyPzE(0,0,0,mC);
//   beam1=proj1+targ1;
//   beamC=proj1+targC;
// 	 float  pion_mom=0.690; //GeV/c                                                                                               
//   float pion_en = sqrt( pion_mom*pion_mom + mpim*mpim );



	Float_t Th1,Th2,Th3,En1,En2,En3,Phi1,Phi2,Phi3,pTyp1,pTyp2,pTyp3; //Zmienne pierwotne z Generatora PLUTO
	Float_t Th1_sm,Th2_sm,Th3_sm,En1_sm,En2_sm,En3_sm,Phi1_sm,Phi2_sm,Phi3_sm; //Zmienne po rozmyciu
  Float_t sigma_Th1,sigma_Th2,sigma_Th3,sigma_En1,sigma_En2,sigma_En3,sigma_Phi1,sigma_Phi2,sigma_Phi3;	 //Sigma rozmyc
Float_t mean_Th1,mean_Th2,mean_Th3,mean_En1,mean_En2,mean_En3,mean_Phi1,mean_Phi2,mean_Phi3;	 //Sigma rozmyc

Float_t Th1_f,Th2_f,Th3_f,En1_f,En2_f,En3_f,Phi1_f,Phi2_f,Phi3_f; //Zmienne po rozmyciu
	// Float_t Th_p1,Phi_p1,En_p1,Th_p2,Phi_p2,En_p2,Th_n,Phi_n,En_n;
	// Float_t Th_p1_sm,Phi_p1_sm,En_p1_sm,Th_p2_sm,Phi_p2_sm,En_p2_sm,Th_n_sm,Phi_n_sm,En_n_sm;

	Int_t Converged; 



	TFile *f_old = new TFile("Output_GenEv.root","READ");

	TTree *tree_old;
	f_old->GetObject("T",tree_old);

		ULong64_t entries=(ULong64_t)tree_old->GetEntries();

	tree_old->SetBranchAddress("pTyp1",&pTyp1);
	tree_old->SetBranchAddress("pTyp2",&pTyp2);
	tree_old->SetBranchAddress("pTyp3",&pTyp3);


	tree_old->SetBranchAddress("Th1",&Th1);
	tree_old->SetBranchAddress("Phi1",&Phi1);
	tree_old->SetBranchAddress("En1",&En1);

	tree_old->SetBranchAddress("Th2",&Th2);
	tree_old->SetBranchAddress("Phi2",&Phi2);
	tree_old->SetBranchAddress("En2",&En2);

	tree_old->SetBranchAddress("Th3",&Th3);
	tree_old->SetBranchAddress("Phi3",&Phi3);
	tree_old->SetBranchAddress("En3",&En3);


	tree_old->SetBranchAddress("Th1_sm",&Th1_sm);
	tree_old->SetBranchAddress("Phi1_sm",&Phi1_sm);
	tree_old->SetBranchAddress("En1_sm",&En1_sm);

	tree_old->SetBranchAddress("Th2_sm",&Th2_sm);
	tree_old->SetBranchAddress("Phi2_sm",&Phi2_sm);
	tree_old->SetBranchAddress("En2_sm",&En2_sm);

	tree_old->SetBranchAddress("Th3_sm",&Th3_sm);
	tree_old->SetBranchAddress("Phi3_sm",&Phi3_sm);
	tree_old->SetBranchAddress("En3_sm",&En3_sm);

	tree_old->SetBranchAddress("sigma_Th1",&sigma_Th1);
	tree_old->SetBranchAddress("sigma_Phi1",&sigma_Phi1);
	tree_old->SetBranchAddress("sigma_En1",&sigma_En1);

	tree_old->SetBranchAddress("sigma_Th2",&sigma_Th2);
	tree_old->SetBranchAddress("sigma_Phi2",&sigma_Phi2);
	tree_old->SetBranchAddress("sigma_En2",&sigma_En2);

	tree_old->SetBranchAddress("sigma_Th3",&sigma_Th3);
	tree_old->SetBranchAddress("sigma_Phi3",&sigma_Phi3);
	tree_old->SetBranchAddress("sigma_En3",&sigma_En3);

	tree_old->SetBranchAddress("mean_Th1",&mean_Th1);
	tree_old->SetBranchAddress("mean_Phi1",&mean_Phi1);
	tree_old->SetBranchAddress("mean_En1",&mean_En1);

	tree_old->SetBranchAddress("mean_Th2",&mean_Th2);
	tree_old->SetBranchAddress("mean_Phi2",&mean_Phi2);
	tree_old->SetBranchAddress("mean_En2",&mean_En2);

	tree_old->SetBranchAddress("mean_Th3",&mean_Th3);
	tree_old->SetBranchAddress("mean_Phi3",&mean_Phi3);
	tree_old->SetBranchAddress("mean_En3",&mean_En3);

	TRandom3 *rnd = new TRandom3();
	mean_Th1=0;
		sigma_Th1=1;
	mean_Th2=0;
		sigma_Th2=1;
	mean_Th3=0;
		sigma_Th3=3;

	mean_Phi1=0;
		sigma_Phi1=1;
	mean_Phi2=0;
		sigma_Phi2=1;
	mean_Phi3=0;
		sigma_Phi3=3;

	mean_En1=0;
		sigma_En1=2;
	mean_En2=0;
		sigma_En2=2;
	mean_En3=0;
		sigma_En3=10;

	pTyp1=1;
	pTyp2=1;
	pTyp3=2;

	TFile *f = new TFile("./Output_GenEv_1.root","RECREATE");

	TTree *tree = new TTree("T","Smeared output");

	tree->Branch("pTyp1",&pTyp1,"pTyp1/F");
	tree->Branch("pTyp2",&pTyp2,"pTyp2/F");
	tree->Branch("pTyp3",&pTyp3,"pTyp3/F");

	tree->Branch("Th1",&Th1,"Th1/F");
	tree->Branch("Phi1",&Phi1,"Phi1/F");
	tree->Branch("En1",&En1,"En1/F");

	tree->Branch("Th2",&Th2,"Th2/F");
	tree->Branch("Phi2",&Phi2,"Phi2/F");
	tree->Branch("En2",&En2,"En2/F");

	tree->Branch("Th3",&Th3,"Th3/F");
	tree->Branch("Phi3",&Phi3,"Phi3/F");
	tree->Branch("En3",&En3,"En3/F");

	tree->Branch("Th1_sm",&Th1_sm,"Th1_sm/F");
	tree->Branch("Phi1_sm",&Phi1_sm,"Phi1_sm/F");
	tree->Branch("En1_sm",&En1_sm,"En1_sm/F");

	tree->Branch("Th2_sm",&Th2_sm,"Th2_sm/F");
	tree->Branch("Phi2_sm",&Phi2_sm,"Phi2_sm/F");
	tree->Branch("En2_sm",&En2_sm,"En2_sm/F");

	tree->Branch("Th3_sm",&Th3_sm,"Th3_sm/F");
	tree->Branch("Phi3_sm",&Phi3_sm,"Phi3_sm/F");
	tree->Branch("En3_sm",&En3_sm,"En3_sm/F");

	tree->Branch("sigma_Th1",&sigma_Th1,"sigma_Th1/F");
	tree->Branch("sigma_Phi1",&sigma_Phi1,"sigma_Phi1/F");
	tree->Branch("sigma_En1",&sigma_En1,"sigma_En1/F");

	tree->Branch("sigma_Th2",&sigma_Th2,"sigma_Th2/F");
	tree->Branch("sigma_Phi2",&sigma_Phi2,"sigma_Phi2/F");
	tree->Branch("sigma_En2",&sigma_En2,"sigma_En2/F");

	tree->Branch("sigma_Th3",&sigma_Th3,"sigma_Th3/F");
	tree->Branch("sigma_Phi3",&sigma_Phi3,"sigma_Phi3/F");
	tree->Branch("sigma_En3",&sigma_En3,"sigma_En3/F");

	tree->Branch("mean_Th1",&mean_Th1,"mean_Th1/F");
	tree->Branch("mean_Phi1",&mean_Phi1,"mean_Phi1/F");
	tree->Branch("mean_En1",&mean_En1,"mean_En1/F");

	tree->Branch("mean_Th2",&mean_Th2,"mean_Th2/F");
	tree->Branch("mean_Phi2",&mean_Phi2,"mean_Phi2/F");
	tree->Branch("mean_En2",&mean_En2,"mean_En2/F");

	tree->Branch("mean_Th3",&mean_Th3,"mean_Th3/F");
	tree->Branch("mean_Phi3",&mean_Phi3,"mean_Phi3/F");
	tree->Branch("mean_En3",&mean_En3,"mean_En3/F");

	TH1F *hMM1_sm= new TH1F("hMM1_sm","Smeared Missing Mass of 1st part.",100,900,1000);
	TH1F *hMM2_sm= new TH1F("hMM2_sm","Smeared Missing Mass of 2nd part.",100,900,1000);
	TH1F *hMM3_sm= new TH1F("hMM3_sm","Smeared Missing Mass of 3rd part.",100,900,1000);

	TH1F *hMM1= new TH1F("hMM1","Missing Mass of 1st part.",100,900,1000);
	TH1F *hMM2= new TH1F("hMM2","Missing Mass of 2nd part.",100,900,1000);
	TH1F *hMM3= new TH1F("hMM3","Missing Mass of 3rd part.",100,900,1000);




	TFile *ff = new TFile("./Refit.root","RECREATE");

	TTree *treef = new TTree("T","Fitted values");

	treef->Branch("Th1",&Th1,"Th1/F");
	treef->Branch("Phi1",&Phi1,"Phi1/F");
	treef->Branch("En1",&En1,"En1/F");

	treef->Branch("Th2",&Th2,"Th2/F");
	treef->Branch("Phi2",&Phi2,"Phi2/F");
	treef->Branch("En2",&En2,"En2/F");

	treef->Branch("Th1_sm",&Th1_sm,"Th1_sm/F");
	treef->Branch("Phi1_sm",&Phi1_sm,"Phi1_sm/F");
	treef->Branch("En1_sm",&En1_sm,"En1_sm/F");

	treef->Branch("Th2_sm",&Th2_sm,"Th2_sm/F");
	treef->Branch("Phi2_sm",&Phi2_sm,"Phi2_sm/F");
	treef->Branch("En2_sm",&En2_sm,"En2_sm/F");

	treef->Branch("Th1_f",&Th1_f,"Th1_f/F");
	treef->Branch("Phi1_f",&Phi1_f,"Phi1_f/F");
	treef->Branch("En1_f",&En1_f,"En1_f/F");

	treef->Branch("Th2_f",&Th2_f,"Th2_f/F");
	treef->Branch("Phi2_f",&Phi2_f,"Phi2_f/F");
	treef->Branch("En2_f",&En2_f,"En2_f/F");

	treef->Branch("Chi2",&fChi2,"Chi2/F");
	treef->Branch("Probability",&fProb,"Probability/F");

	treef->Branch("Converged",&Converged,"Converged/I");

	TH1F *hMM1_f= new TH1F("hMM1_f","Smeared Missing Mass of 1st part.",100,900,1000);
	TH1F *hMM2_f= new TH1F("hMM2_f","Smeared Missing Mass of 2nd part.",100,900,1000);
	TH1F *hMM3_f= new TH1F("hMM3_f","Smeared Missing Mass of 3rd part.",100,900,1000);

	TH1F *hChi2=new TH1F("Chi2","Chi2",100,0,30);
	TH1F *hProbability=new TH1F("Probability","Probability",100,0,2);

	TH1F *Th1_prev=new TH1F("Th1_prev","Th1_prev",200,0,100);
	TH1F *Phi1_prev=new TH1F("Phi1_prev","Phi1_prev",200,0,200);
	TH1F *Th2_prev=new TH1F("Th2_prev","Th2_prev",200,0,100);
	TH1F *Phi2_prev=new TH1F("Phi2_prev","Phi2_prev",200,0,200);
	TH1F *En1_prev=new TH1F("En1_prev","En1_prev",200,0,200);
	TH1F *En2_prev=new TH1F("En2_prev","En2_prev",200,0,200);

	TH1F *Th1_after=new TH1F("Th1_after","Th1_after",200,0,100);
	TH1F *Phi1_after=new TH1F("Phi1_after","Phi1_after",200,0,200);
	TH1F *Th2_after=new TH1F("Th2_after","Th2_after",200,0,100);
	TH1F *Phi2_after=new TH1F("Phi2_after","Phi2_after",200,0,200);
	TH1F *En1_after=new TH1F("En1_after","En1_after",200,0,200);
	TH1F *En2_after=new TH1F("En2_after","En2_after",200,0,200);


	// for(ULong64_t i=0;i<entries;i++){
	for(ULong64_t i=0;i<100000;i++){
		glob_event++;
		Counting(i);
		tree_old->GetEntry(i);
		fChi2=0;
		fProb=0;
		Th1_sm=Th1+rnd->Gaus(mean_Th1,sigma_Th1);
		Th2_sm=Th2+rnd->Gaus(mean_Th2,sigma_Th2);
		Th3_sm=Th3+rnd->Gaus(mean_Th3,sigma_Th3);

		Phi1_sm=Phi1+rnd->Gaus(mean_Phi1,sigma_Phi1);
		Phi2_sm=Phi2+rnd->Gaus(mean_Phi2,sigma_Phi2);
		Phi3_sm=Phi3+rnd->Gaus(mean_Phi3,sigma_Phi3);

		En1_sm=En1+rnd->Gaus(mean_En1,sigma_En1);
		En2_sm=En2+rnd->Gaus(mean_En2,sigma_En2);
		En3_sm=En3+rnd->Gaus(mean_En3,sigma_En3);

		//cyliczne warunki na Th (raczej hipotetyczne)

		if(Th1_sm<0) {Th1_sm=abs(Th1_sm); Phi1_sm+=180;};
    	if(Th2_sm<0) {Th2_sm=abs(Th2_sm); Phi2_sm+=180;};
		if(Th3_sm<0) {Th3_sm=abs(Th3_sm); Phi3_sm+=180;};

		if(Th1_sm>180) {Th1_sm=360-Th1_sm; Phi1_sm-=180;};
		if(Th2_sm>180) {Th2_sm=360-Th2_sm; Phi2_sm-=180;};
		if(Th3_sm>180) {Th3_sm=360-Th3_sm; Phi3_sm-=180;};

		//Cykliczne Warunki na Phi

		if(Phi1_sm>180) Phi1_sm-=360;
		if(Phi2_sm>180) Phi2_sm-=360;
		if(Phi3_sm>180) Phi3_sm-=360;

		if(Phi1_sm<-180) Phi1_sm+=360;
		if(Phi2_sm<-180) Phi2_sm+=360;
		if(Phi3_sm<-180) Phi3_sm+=360;


		Float_t MM1,MM2,MM3,MM1_sm,MM2_sm,MM3_sm,MM1_f,MM2_f,MM3_f;
		// MM1=calc_MM(En2,En3,Th2,Th3,Phi2,Phi3,1);
		// MM2=calc_MM(En1,En3,Th1,Th3,Phi1,Phi3,1);
		MM3=calc_MM(En1,En2,Th1,Th2,Phi1,Phi2,0);

		//cout<<"En1="<<En1<<" En2="<<En2<<" En3="<<En3<<endl;

		// MM1_sm=calc_MM(En2_sm,En3_sm,Th2_sm,Th3_sm,Phi2_sm,Phi3_sm,1);
		// MM2_sm=calc_MM(En1_sm,En3_sm,Th1_sm,Th3_sm,Phi1_sm,Phi3_sm,1);
		MM3_sm=calc_MM(En1_sm,En2_sm,Th1_sm,Th2_sm,Phi1_sm,Phi2_sm,0);

		// hMM1->Fill(MM1);
		// hMM2->Fill(MM2);
		hMM3->Fill(MM3);

		// hMM1_sm->Fill(MM1_sm);
		// hMM2_sm->Fill(MM2_sm);
		hMM3_sm->Fill(MM3_sm);

		tree->Fill();

		SetValues(Phi1_sm, Th1_sm, En1_sm, Phi2_sm, Th2_sm, En2_sm);
		SetCovariance(sigma_Phi1, sigma_Th1, sigma_En1, sigma_Phi2, sigma_Th2, sigma_En2);
		// SetCovariance(30.,30.,30.,30.,30.,30.);
		// f_eval();
		// Feta_eval();
		Converged = fit();
		// GetFitVal();

		TMatrixD FitVal = GetFitVal(ytemp);

		Phi1_f = FitVal(0,0);
		Th1_f = FitVal(1,0);
		En1_f = FitVal(2,0);

		Phi2_f = FitVal(3,0);
		Th2_f = FitVal(4,0);
		En2_f = FitVal(5,0);

		if(isnan(Th1_f)) f_nan++;

		// MM1_f=calc_MM(En2_f,En3_f,Th2_f,Th3_f,Phi2_f,Phi3_f,1);
		// MM2_f=calc_MM(En1_f,En3_f,Th1_f,Th3_f,Phi1_f,Phi3_f,1);
		MM3_f=calc_MM(En1_f,En2_f,Th1_f,Th2_f,Phi1_f,Phi2_f,0);

		// hMM1_f->Fill(MM1_f);
		// hMM2_f->Fill(MM2_f);
		hMM3_f->Fill(MM3_f);
		hChi2->Fill(fChi2);
		hProbability->Fill(fProb);

		Th1_prev->Fill(Th1_sm);
		Phi1_prev->Fill(Phi1_sm);
		Th2_prev->Fill(Th2_sm);
		Phi2_prev->Fill(Phi2_sm);
		En1_prev->Fill(En1_sm);
		En2_prev->Fill(En2_sm);

		Th1_after->Fill(Th1_f);
		Phi1_after->Fill(Phi1_f);
		Th2_after->Fill(Th2_f);
		Phi2_after->Fill(Phi2_f);
		En1_after->Fill(En1_f);
		En2_after->Fill(En2_f);


		// cout<<"Phi1: "<<Phi1_f<<" ";
		// cout<<"Phi2: "<<Phi2_f<<" ";
		// cout<<"Th1: "<<Th1_f<<" ";
		// cout<<"Th2: "<<Th2_f<<endl;


		treef->Fill();


	}
cout<<endl;
cout<<"ilosc nan: "<<f_nan<<endl;

f->cd();
tree->Write();
// hMM1->Write();
// hMM2->Write();
hMM3->Write();
// hMM1_sm->Write();
// hMM2_sm->Write();
hMM3_sm->Write();

Results->Divide(3,4);

Results->cd(1);
Th1_prev->Draw();
Results->cd(5);
Phi1_prev->Draw();
Results->cd(3);
Th2_prev->Draw();
Results->cd(7);
Phi2_prev->Draw();
Results->cd(2);
Th1_after->Draw();
Results->cd(6);
Phi1_after->Draw();
Results->cd(4);
Th2_after->Draw();
Results->cd(8);
Phi2_after->Draw();
Results->cd(9);
En1_prev->Draw();
Results->cd(10);
En1_after->Draw();
Results->cd(11);
En2_prev->Draw();
Results->cd(12);
En2_after->Draw();

Results->Write();
Results->SaveAs("Results.png");

ff->cd();
treef->Write();
// hMM1_f->Write();
// hMM2_f->Write();
hMM3_f->Write();
// hMM1_sm->Write();
// hMM2_sm->Write();
hMM3_sm->Write();
hMM3->Write();

Various->Divide(2,2); //przeniesc do hist.c

Various->cd(1);
hMM3_f->Draw();
Various->cd(2);
hMM3_sm->Draw();
Various->cd(3);
hChi2->Draw();
Various->cd(4);
hProbability->Draw();

Various->Write();
Various->SaveAs("Various.png");

f->Close();
ff->Close();
f_old->Close();

fend = clock();
printf("liczba eventow w pliku=%d\n",entries);
double time_spent;
time_spent = (double)(fend - fbegin) / CLOCKS_PER_SEC;
printf("Czas wykonywania programu =%f s\n",time_spent);
return 0;
}
