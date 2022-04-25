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

TCanvas *Results = new TCanvas("Results","Results",3000,3000);
TCanvas *Various = new TCanvas("Various","Various",3000,3000);

Int_t cov_dim = 3;
Int_t fN = 2;
TMatrixD ytemp(cov_dim * fN,1);
TMatrixD V(cov_dim * fN, cov_dim * fN);
TMatrixD m(fN,1);
float mp = 938.272;
float mn = 939.565;
float md = 1875.612;
Int_t fNdf = 1;
const TLorentzVector pinit(0.,0.,791.073,2973.88);
// const TLorentzVector pinit(0.,0.,791.073,3913.445);
float fChi2, fProb;
Int_t glob_event = 0.;
Int_t f_nan = 0;
float deltachi = 0;
Int_t iter = 0;
int nany = 0;

void Counting(ULong64_t ent)
{
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


float calc_MM(float en1, float en2, float th1, float th2, float phi1, float phi2, int cat)
{
    float MMass=0;
    float Mp= 938.272;
    float Md=1875.612;
    float Mn= 939.565;
    float Ed=160.0;
    float M1=0,M2=0,M3=0;
    if(cat==0){M1=Mp; M2=Mp; M3=Mn;} //Missing Mass dla pary proton-proton
    if(cat==1){M1=Mp, M2=Mn; M3=Mp;} //Missing Mass dla pary proton-neutron
    TVector3 p1,p2,p3,d;
    p1.SetMagThetaPhi(sqrt((en1*en1)+2*M1*en1), th1*Pi_180, phi1*Pi_180);
    p2.SetMagThetaPhi(sqrt((en2*en2)+2*M2*en2), th2*Pi_180, phi2*Pi_180);
    d.SetMagThetaPhi(sqrt((Ed*Ed)+2*Md*Ed), 0*Pi_180, 0*Pi_180);
    p3=d-(p1+p2);
    float ER=160+Md-M3-en1-en2;
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

// void SetValues(float phi1, float th1, float en1, float phi2, float th2, float en2, float phi3, float th3, float en3)
void SetValues(float phi1, float th1, float en1, float phi2, float th2, float en2)
{
    ytemp.Zero();
    ytemp(0,0)=phi1;
    ytemp(1,0)=th1;
    ytemp(2,0)=en1;
    ytemp(3,0)=phi2;
    ytemp(4,0)=th2;
    ytemp(5,0)=en2;
    // ytemp(6,0)=phi3;
    // ytemp(7,0)=th3;
    // ytemp(8,0)=en3;

    m.Zero();
    m(0,0)=mp;
    m(1,0)=mp;
    // m(2,0)=mn;
}

// void SetCovariance(float phi1_er, float th1_er, float en1_er, float phi2_er, float th2_er, float en2_er, float phi3_er, float th3_er, float en3_er)
void SetCovariance(float phi1_er, float th1_er, float en1_er, float phi2_er, float th2_er, float en2_er)

{
    V.Zero();
    V(0,0) = phi1_er*phi1_er;
    V(1,1) = th1_er*th1_er;
    V(2,2) = en1_er*en1_er;
    V(3,3) = phi2_er*phi2_er;
    V(4,4) = th2_er*th2_er;
    V(5,5) = en2_er*en2_er;
    // V(6,6) = phi3_er*phi3_er;
    // V(7,7) = th3_er*th3_er;
    // V(8,8) = en3_er*en3_er;
}

TMatrixD f_eval(TMatrixD y)
{
    TMatrixD d;
    d.ResizeTo(1, 1);
    float P, Px, Py, Pz, E;
    Px = pinit.Px();
    Py = pinit.Py();
    Pz = pinit.Pz();
    E = pinit.E();

    for (int q = 0; q < fN; q++) 
    {
        E-=y(2+q*cov_dim,0)+m(q,0);

        P=sqrt(pow(y(2+q*cov_dim,0)+m(q,0),2)-m(q,0)*m(q,0));
        Px-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*cos(y(0+q*cov_dim,0)*Deg2Rad);
        Py-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*sin(y(0+q*cov_dim,0)*Deg2Rad);
        Pz-=P*cos(y(1+q*cov_dim,0)*Deg2Rad);
        // cout<<"Ev: "<<glob_event<<" E: "<<E<<" P: "<<P<<" Px: "<<Px<<" Py: "<<Py<<" Pz: "<<Pz<<endl;

    }

    d(0, 0) = std::pow(E, 2) - std::pow(Px, 2) - std::pow(Py, 2) - std::pow(Pz, 2) - mn * mn;
    // cout<<"d: "<<d(0,0)<<endl;
    return d;
}

TMatrixD Feta_eval(TMatrixD y)
{
    TMatrixD H;
    H.ResizeTo(1, cov_dim * fN);
    H.Zero();
    float P, Px, Py, Pz, E;

    Px=0;
    Py=0;
    Pz=791.073;
    // E=3913.445;
    E=2973.88;

    for (int q = 0; q < fN; q++)
    {
        E-=y(2+q*cov_dim,0)+m(q,0);

        P=sqrt(pow(y(2+q*cov_dim,0)+m(q,0),2)-m(q,0)*m(q,0));
        Px-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*cos(y(0+q*cov_dim,0)*Deg2Rad);
        Py-=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*sin(y(0+q*cov_dim,0)*Deg2Rad);
        Pz-=P*cos(y(1+q*cov_dim,0)*Deg2Rad);
    }

    for (int q = 0; q < fN; q++)
    {
        float Pi = sqrt(pow(y(2+q*cov_dim,0)+m(q,0),2)-m(q,0)*m(q,0));
        float Ei = sqrt(Pi * Pi + m(q,0) * m(q,0));

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

    return H;
}


Int_t fit(double lr, int q, float deltachi)
{
    // double lr = 0.0006;
    // double lr = 0.01;
    TMatrixD alpha0(fN * cov_dim, 1), alpha(fN * cov_dim, 1);
    alpha0 = ytemp;
    alpha = alpha0;
    double chi2 = 1e6;
    TMatrixD D = Feta_eval(alpha);
    TMatrixD d = f_eval(alpha);
    // Int_t q = 0;
    Int_t fConverged = 0;
    iter=0;

    for (q = 0; q < 40; q++)
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
            // cout<<"lambdaT: "<<lambdaT(0,0)<<" d: "<<d(0,0)<<endl;
            chisqrd = lambdaT(0, p) * d(p, 0);
            // chisqrd = lambdaT(0, p) * d(p, 0)/10000;
            // cout<<"Ev: "<<glob_event<<" q: "<<q<<" chi2: "<<chi2<<" chisqrd: "<<chisqrd<<endl;
        }

        if(fabs(abs(chi2-chisqrd))<deltachi)
        {
            fConverged = 1;
            iter=q;
            break;
        }

        chi2 = chisqrd;
        alpha0 = alpha;
        alpha = neu_alpha;
        V = V - lr * V * DT * VD * D * V;
        D = Feta_eval(alpha);
        d = f_eval(alpha);
    }

    // fConverged = 1;
    ytemp = alpha;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    return fConverged;
}

TMatrixD GetFitVal(TMatrixD y)
{
    TMatrixD y_fit;
    TVector3 Convert;
    TLorentzVector *v;
    v = new TLorentzVector(0,0,0,0);
    y_fit.ResizeTo(cov_dim * fN,1);
    float E = 0., P = 0., Px = 0., Py = 0., Pz = 0.;
    for(int q = 0; q < fN; q++)
    {
        E=y(2+q*cov_dim,0)+m(q,0);

        P=sqrt(pow(y(2+q*cov_dim,0)+m(q,0),2)-m(q,0)*m(q,0));
        Px=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*cos(y(0+q*cov_dim,0)*Deg2Rad);
        Py=P*sin(y(1+q*cov_dim,0)*Deg2Rad)*sin(y(0+q*cov_dim,0)*Deg2Rad);
        Pz=P*cos(y(1+q*cov_dim,0)*Deg2Rad);

        v->SetPxPyPzE(Px,Py,Pz,E);

        y_fit(0 + q * cov_dim, 0) = v->Phi()*Rad2Deg;
        y_fit(1 + q * cov_dim, 0) = v->Theta()*Rad2Deg;
        y_fit(2 + q * cov_dim, 0) = v->E()-m(q,0);
    }

    return y_fit;
}

int GenEv_refit(double lr, int q, float deltachi)
{
    clock_t fbegin;
    clock_t fend;
    fbegin = clock();

    float Th1,Th2,Th3,En1,En2,En3,Phi1,Phi2,Phi3,pTyp1,pTyp2,pTyp3; //Zmienne pierwotne z Generatora PLUTO
    float Th1_sm,Th2_sm,Th3_sm,En1_sm,En2_sm,En3_sm,Phi1_sm,Phi2_sm,Phi3_sm; //Zmienne po rozmyciu
    float sigma_Th1,sigma_Th2,sigma_Th3,sigma_En1,sigma_En2,sigma_En3,sigma_Phi1,sigma_Phi2,sigma_Phi3;    //Sigma rozmyc
    float mean_Th1,mean_Th2,mean_Th3,mean_En1,mean_En2,mean_En3,mean_Phi1,mean_Phi2,mean_Phi3;   //Mean rozmyc
    float Th1_f,Th2_f,Th3_f,En1_f,En2_f,En3_f,Phi1_f,Phi2_f,Phi3_f; //Zmienne po rozmyciu

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

    TH1F *hMM1_f= new TH1F("hMM1_f","Smeared Missing Mass of 1st part.",200,920,955);
    TH1F *hMM2_f= new TH1F("hMM2_f","Smeared Missing Mass of 2nd part.",200,920,955);
    TH1F *hMM3_f= new TH1F("hMM3_f","Smeared Missing Mass of 3rd part.",200,900,1000);

    TH1F *hChi2=new TH1F("Chi2","Chi2",200,0,50);
    TH1F *hProbability=new TH1F("Probability","Probability",100,0,1);

    TH1F *Th1_prev=new TH1F("Th1_prev","Th1_prev",200,0,100);
    TH1F *Phi1_prev=new TH1F("Phi1_prev","Phi1_prev",200,-200,200);
    TH1F *Th2_prev=new TH1F("Th2_prev","Th2_prev",200,0,100);
    TH1F *Phi2_prev=new TH1F("Phi2_prev","Phi2_prev",200,-200,200);
    TH1F *En1_prev=new TH1F("En1_prev","En1_prev",200,-20,200);
    TH1F *En2_prev=new TH1F("En2_prev","En2_prev",200,-20,200);

    TH1F *Th1_after=new TH1F("Th1_after","Th1_after",200,0,100);
    TH1F *Phi1_after=new TH1F("Phi1_after","Phi1_after",200,-200,200);
    TH1F *Th2_after=new TH1F("Th2_after","Th2_after",200,0,100);
    TH1F *Phi2_after=new TH1F("Phi2_after","Phi2_after",200,-200,200);
    TH1F *En1_after=new TH1F("En1_after","En1_after",200,0,200);
    TH1F *En2_after=new TH1F("En2_after","En2_after",200,0,200);

    TH1F *hdeltachi=new TH1F("hdeltachi","hdeltachi",600,0,40);
    TH1F *hq=new TH1F("hq","hq",200,-1,100);
    TH1F *hnan=new TH1F("hnan","hnan",200,0,200);

    // for(ULong64_t i=0;i<entries;i++){
    for(ULong64_t i=0;i<15000;i++)
    {
        glob_event++;
        Counting(i);
        tree_old->GetEntry(i);

        float MM1,MM2,MM3,MM1_sm,MM2_sm,MM3_sm,MM1_f,MM2_f,MM3_f;

        // MM1=calc_MM(En2,En3,Th2,Th3,Phi2,Phi3,1);
        // MM2=calc_MM(En1,En3,Th1,Th3,Phi1,Phi3,1);
        MM3=calc_MM(En1,En2,Th1,Th2,Phi1,Phi2,0);

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

        // SetValues(Phi1_sm, Th1_sm, En1_sm, Phi2_sm, Th2_sm, En2_sm, Phi3_sm, Th3_sm, En3_sm);
        // SetCovariance(sigma_Phi1, sigma_Th1, sigma_En1, sigma_Phi2, sigma_Th2, sigma_En2, sigma_Phi3, sigma_Th3, sigma_En3);
        
        SetValues(Phi1_sm, Th1_sm, En1_sm, Phi2_sm, Th2_sm, En2_sm);
        SetCovariance(sigma_Phi1, sigma_Th1, sigma_En1, sigma_Phi2, sigma_Th2, sigma_En2);
        
        Converged = fit(lr, q, deltachi);
        
        TMatrixD FitVal = GetFitVal(ytemp);

        Phi1_f = FitVal(0,0);
        Th1_f = FitVal(1,0);
        En1_f = FitVal(2,0);

        Phi2_f = FitVal(3,0);
        Th2_f = FitVal(4,0);
        En2_f = FitVal(5,0);

        // Phi3_f = FitVal(6,0);
        // Th3_f = FitVal(7,0);
        // En3_f = FitVal(8,0);

        if(isnan(En1_f)) hnan->Fill(En1_sm);

        if(isnan(Phi1_f) || isnan(Th1_f) || isnan(En1_f)) f_nan++;

        // MM1_f=calc_MM(En2_f,En3_f,Th2_f,Th3_f,Phi2_f,Phi3_f,1);
        // MM2_f=calc_MM(En1_f,En3_f,Th1_f,Th3_f,Phi1_f,Phi3_f,1);
        MM3_f=calc_MM(En1_f,En2_f,Th1_f,Th2_f,Phi1_f,Phi2_f,0);

        // hMM1_f->Fill(MM1_f);
        // hMM2_f->Fill(MM2_f);
            
        hdeltachi->Fill(deltachi);
        hq->Fill(iter);

        
            Phi1_prev->Fill(Phi1_sm);
            Th1_prev->Fill(Th1_sm);
            En1_prev->Fill(En1_sm);
            Phi2_prev->Fill(Phi2_sm);
            Th2_prev->Fill(Th2_sm);
            En2_prev->Fill(En2_sm);

        //     if(Converged==1)
        // {
            hMM3_f->Fill(MM3_f);
            Phi1_after->Fill(Phi1_f);
            Th1_after->Fill(Th1_f);
            En1_after->Fill(En1_f);
            Phi2_after->Fill(Phi2_f);
            Th2_after->Fill(Th2_f);
            En2_after->Fill(En2_f);
            hProbability->Fill(fProb);
            hChi2->Fill(fChi2);
        // }

        treef->Fill();
    }

    cout<<endl;
    cout<<"ilosc nan: "<<f_nan<<endl;
    cout<<"nany: "<<nany<<endl;
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
    Th1_prev->SetLineWidth(2);
    Th1_prev->Draw();
    Results->cd(5);
    Phi1_prev->SetLineWidth(2);
    Phi1_prev->Draw();
    Results->cd(3);
    Th2_prev->SetLineWidth(2);
    Th2_prev->Draw();
    Results->cd(7);
    Phi2_prev->SetLineWidth(2);
    Phi2_prev->Draw();
    Results->cd(2);
    Th1_after->SetLineWidth(2);
    Th1_after->Draw();
    Results->cd(6);
    Phi1_after->SetLineWidth(2);
    Phi1_after->Draw();
    Results->cd(4);
    Th2_after->SetLineWidth(2);
    Th2_after->Draw();
    Results->cd(8);
    Phi2_after->SetLineWidth(2);
    Phi2_after->Draw();
    Results->cd(9);
    En1_prev->SetLineWidth(2);
    En1_prev->Draw();
    Results->cd(10);
    En1_after->SetLineWidth(2);
    En1_after->Draw();
    Results->cd(11);
    En2_prev->SetLineWidth(2);
    En2_prev->Draw();
    Results->cd(12);
    En2_after->SetLineWidth(2);
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
    hdeltachi->Write();
    hq->Write();
    hnan->Write();

    Various->Divide(2,2);

    Various->cd(1);
    hMM3_f->SetLineWidth(2);
    hMM3_f->Draw();
    Various->cd(2);
    hMM3_sm->SetLineWidth(2);
    hMM3_sm->Draw();
    Various->cd(3);
    hChi2->SetLineWidth(2);
    hChi2->Draw();
    Various->cd(4);
    hProbability->SetLineWidth(2);
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
