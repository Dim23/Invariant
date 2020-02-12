#define hTANA_cxx
#include "hTANA.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include <stdio.h>
#include "pid.C" 

int pn=14;TH1F *hmassW[14];TH1F *hmassE[14];TH1F *FitW[14];TH1F *FitE[14]; TH1F *invW;TH1F *invE;TH1F *invWE;TH1F *invALL;TH1F *normKW;TH1F *normKE;
//double A=0.5,B=3.5,d=(B-A)/(pn-1);
static const double bin_w[15]={0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.6,1.8,2.0,2.3};

//float *a=new float[pn];
//float *b=new float[pn];
double *binP=new double[pn];
double *sigW=new double[pn];
double *sigE=new double[pn];
double *sig2W=new double[pn];
double *sig3W=new double[pn];
//double *mean2W=new double[pn];
//double *mean3W=new double[pn];
double *sig2E=new double[pn];
double *sig3E=new double[pn];
//double *mean2E=new double[pn];
//double *mean3E=new double[pn];



int NN;
float pt,eta,m,beta,dtime,charge_mom;
const float CC=299792458,Mpi=0.019479835,Mk=0.24371698,Mpr=0.880354511,pi=3.141592654;
char strE[20],strW[20];


    TH1F *dtimW;
    TH1F *dtimE;

    TH1F *hm2E ;
    TH1F *hm2W  ;

    TH2F *hpidW;
    TH2F *hpidE;

    TH1F *htofW;
    TH1F *htofE;

    
    TFile *d_outfile;


void hTANA::Book_Hist(){
int NN=0;

for(int n=0;n<5;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",400,-0.5,1.5);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",400,-0.5,1.5);}

for(int n=5;n<10;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",200,-0.5,1.5);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",200,-0.5,1.5);}

for(int n=10;n<pn;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",100,-0.5,1.5);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",100,-0.5,1.5);}




dtimW=new TH1F("dtimeW","tof - t(exp for #pi),ns  TOF West Arm",200,-2,10);
dtimE=new TH1F("dtimeE","tof - t(exp for #pi),ns  TOF East Arm ",200,-2,10);

hm2E  = new TH1F("hm2E","mass^{2} for Tof.East",150,-0.2,1.8);
hm2W  = new TH1F("hm2W","mass^{2} for Tof.West",150,-0.2,1.8);
invW  = new TH1F("invW","inv West",1000,0,2);
invE  = new TH1F("invE","inv East",1000,0,2);
invWE  = new TH1F("invWE","inv West and East",500,0,7);
invALL  = new TH1F("invALL","inv West and East ALL",4000,0,8);

normKE  = new TH1F("normKE","NORM K East",50,-3,3);
normKW  = new TH1F("normKW","NORM K West",50,-3,3);

hpidW= new TH2F("hpidW","charge/momentum vs time of flight,  West Arm ",200,10,60,400,-2.2,2.2);
hpidE= new TH2F("hpidE","charge/momentum vs time of flight,  East Arm ",200,10,60,400,-2.2,2.2);

htofW= new TH1F("htofW","tof - t(exp for #pi),ns  TOF West Arm ",200,-2,3);
htofE= new TH1F("htofE","tof - t(exp for #pi),ns  TOF East Arm ",200,-3,3);

}


void hTANA::Loop()
{
//   In a ROOT session, you can do:
//      root> .L hTANA.C
//      root> hTANA t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
for(int n=0;n<pn;n++){ binP[n]=0.5*(bin_w[n]+bin_w[n+1]);}

Loop_imp();

/*if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
        
        Long64_t ientry = LoadTree(jentry);
        if(ientry%100000==0) cout << ientry<<endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

if(fabs(bbcz)<30 && cent>0 && cent<=80)
{
    for (int i=0;i<mh;i++)
    {
        pt=p[i]*sin(the0[i]);charge_mom=charge[i]/p[i];
    
        if (pt>1 && pt<1.5 && pltof[i]>0 && sigtof[i]<3) 
        {   
            m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1), 
            //m=p[i]*sqrt((pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1)),
            beta=p[i]/sqrt(pow(p[i],2)+pow(Mpi,2)),dtime=ttof[i]-pltof[i]*(1e+7)/(beta*CC);

            if (dcarm[i]==0 && etof[i]>0.002)
            {dtimE->Fill(dtime);htofE->Fill(dtime);hm2E->Fill(m);hpidE->Fill(ttof[i],charge_mom);}; 

            if (dcarm[i]==1 && etof[i]>60 && etof[i]<600)
            {dtimW->Fill(dtime);htofW->Fill(dtime);hm2W->Fill(m);hpidW->Fill(ttof[i],charge_mom);};
        }
    }
}}

*/
}

void hTANA::add_file(const char *file) {
  TFile *treefile = TFile::Open(file);
  TTree *tree = (TTree*)treefile->Get("mtree");
  if(tree == 0) {
    cout <<"htree is not found in "<<file<<endl;
    treefile->Close();
    return;
  }
  cout << file <<" is opened"<<endl;
  Init(tree);
  cout <<"one file processed"<<endl;
}

void hTANA::Loop_imp()
{
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast(),nbytes = 0, nb = 0;

for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if(ientry%100000==0) cout << ientry<< " all entri is "<< NN <<endl;
        nb = fChain->GetEntry(jentry);  nbytes += nb;


if(fabs(bbcz)<30 && cent>0 && cent<=80 )
{
    for (int i=0;i<mh;i++)
    {  
        pt=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);

for(int n=0;n<pn;n++){
 
        if ( pt>bin_w[n] && pt<bin_w[n+1] && sigtof[i]<3 && pltof[i]>0) 
        {   
            if (dcarm[i]==0 && etof[i]>0.002&&charge[i]>0)

            {hmassE[n]->Fill(m);NN+=1;}; 

            if (dcarm[i]==1 && etof[i]>60 && etof[i]<600)

            {hmassW[n]->Fill(m);NN+=1;};
        }}
    }
}}
cout << " Histograms for fitting was lopped"<<endl;
}

void hTANA::INV()
{
float a,b,c,m12,mk,pti,ptk;
int k;
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast(),nbytes = 0, nb = 0;

for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if(ientry%100000==0) cout << ientry<<endl;
        nb = fChain->GetEntry(jentry);  nbytes += nb;


    if(fabs(bbcz)<30 && cent>0 && cent<=80 )
        {
        for (int i=0;i<mh;i++)
            {  
            pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);

            if (dcarm[i]==0 && fabs(IsKaonE(m,pti))<3)

            {normKE->Fill(IsKaonE(m,pti));k=i;
            for (k;k<mh;k++)
            {
        mk=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);ptk=p[k]*sin(the0[k]);

            if (charge[i]!=charge[k] && fabs(IsKaonE(mk,ptk))<3 && dcarm[k]==0 && etof[i]>0.002){

            a=sqrt((p[i]*p[i]+Mk)*(p[k]*p[k]+Mk));b=cos(the0[i])*cos(the0[k]);
            m12=2*(Mk+a-(pti*ptk*cos(phi0[i]-phi0[k])+p[i]*p[k]*b));invE->Fill(sqrt(m12));invALL->Fill(sqrt(m12));}}}; 
        }   }


        if(fabs(bbcz)<30 && cent>0 && cent<=80 )
        {
        for (int i=0;i<mh;i++)
            {  
            pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);

            if (dcarm[i]==1 && fabs(IsKaonW(m,pti))<3)

            {normKW->Fill(IsKaonW(m,pti));k=i;
            for (k;k<mh;k++)
            {
        mk=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);ptk=p[k]*sin(the0[k]);

            if (charge[i]!=charge[k] && fabs(IsKaonW(mk,ptk))<3 && dcarm[k]==1 && etof[i]>60 && etof[i]<600){

            a=sqrt((p[i]*p[i]+Mk)*(p[k]*p[k]+Mk));b=cos(the0[i])*cos(the0[k]);
            m12=2*(Mk+a-(pti*ptk*cos(phi0[i]-phi0[k])+p[i]*p[k]*b));invW->Fill(sqrt(m12));invALL->Fill(sqrt(m12));}}
             }; 
        }}

if(fabs(bbcz)<30 && cent>0 && cent<=80 )
        {
        for (int i=0;i<mh;i++)
            {  
            pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);

            if (dcarm[i]==0&& fabs(IsKaonE(m,pti))<3)

            {k=i;
            for (k;k<mh;k++)
            {
        mk=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);ptk=p[k]*sin(the0[k]);

            if (dcarm[k]==1&&fabs(IsKaonW(mk,ptk))<3){

            a=sqrt((p[i]*p[i]+Mk)*(p[k]*p[k]+Mk));b=cos(the0[i])*cos(the0[k]);
            m12=2*(Mk+a-(pti*ptk*cos(phi0[i]-phi0[k])+p[i]*p[k]*b));invWE->Fill(sqrt(m12));invALL->Fill(sqrt(m12));}}
             }; 
        }}

    }
cout << " Histograms for fitting was lopped"<<endl;
TCanvas *MyW1 = new TCanvas ("CanSigma1","Test canvasqq1",1);MyW1->Divide(4,1);
MyW1 -> cd(1);invE->GetXaxis()->SetTitle("inv mass,Gev");invE->Draw();gStyle->SetOptStat(1111111);
c=invE->GetMaximum();TLine *line1=new TLine();
line1->DrawLine(1.02,0,1.02,c);cout <<gPad->GetUymax()<<endl;
MyW1 -> cd(2);invW->GetXaxis()->SetTitle("inv mass,Gev");invW->Draw();gStyle->SetOptStat(1111111);
c=invW->GetMaximum();TLine *line2=new TLine();
line2->DrawLine(1.02,0,1.02,c);cout <<MyW1->cd(2)->GetUymax()<<endl;
MyW1 -> cd(3);invWE->GetXaxis()->SetTitle("inv mass,Gev");invWE->Draw();
MyW1 -> cd(4);invALL->Draw();
TCanvas *MyW2 = new TCanvas ("CanSigma2","Test canvasqq2",1);MyW2->Divide(2,1);
MyW2 -> cd(1);normKE->Fit("gaus");normKE->Draw();MyW2 -> cd(2);normKW->Fit("gaus");normKW->Draw();
}


void hTANA::inv()
{
if (fChain == 0) return;
float a,b,c,m12,mk,pti,ptk;
int kp;
Long64_t nentries = fChain->GetEntriesFast(),nbytes = 0, nb = 0;

for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {    kp=0;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if(ientry%100000==0) cout << ientry<< " all entri is "<< kp <<endl;
        nb = fChain->GetEntry(jentry);  nbytes += nb;

        if(fabs(bbcz)<30 && cent>0 && cent<=80 )
        {
        for (int i=0;i<mh;i++)
            {   
                pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);
                if(dcarm[i]==0&& fabs(IsKaonE(m,pti))<3)
                {   kp+=1;
                    for(int k=i;k<mh;k++)
                    {
                        ptk=p[k]*sin(the0[k]);m=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);
                        if(charge[i]!=charge[k] &&dcarm[k]==0&& fabs(IsKaonE(m,ptk))<3)
                        {
                            a=sqrt((p[i]*p[i]+Mk)*(p[k]*p[k]+Mk));b=cos(the0[i])*cos(the0[k]);
                            m12=2*(Mk+a-(pti*ptk*cos(phi0[i]-phi0[k])+p[i]*p[k]*b));invE->Fill(sqrt(m12));
                        }
                    }                
                }
            cout << ientry<< " all entri is "<< kp <<endl;
            }
        }
    }
invE->Draw();
}

//Функция для фита;
double FitGaus(double *x, double *par) 
{
double g1,g2,g3,pdf;
g1=exp(-pow(x[0]-par[0],2)/(2*par[1]*par[1]))/(par[1]*pow(2*pi,0.5));
g2=exp(-pow(x[0]-par[2],2)/(2*par[3]*par[3]))/(par[3]*pow(2*pi,0.5));
g3=exp(-pow(x[0]-par[4],2)/(2*par[5]*par[5]))/(par[5]*pow(2*pi,0.5));
pdf=par[6]*(par[7]*g1+par[8]*g2+(1-par[7]-par[8])*g3);
    return pdf;
};

//Функция для фита;
double FitGaus2(double *x, double *par) 
{
double g1,g2,g3,pdf;
g1=exp(-pow(x[0]-par[0],2)/(2*par[1]*par[1]))/(par[1]*pow(2*pi,0.5));
g2=exp(-pow(x[0]-par[2],2)/(2*par[3]*par[3]))/(par[3]*pow(2*pi,0.5));
g3=exp(-pow(x[0]-par[4],2)/(2*par[5]*par[5]))/(par[5]*pow(2*pi,0.5));
pdf=par[6]*(par[8]*(par[7]*g1+(1-par[7])*g2)+(1-par[8])*g3);
    return pdf;
};
//Функция для опраксимации sigma и mean;
double Func(double x,double a0,double a1,double a2) 
{
    return a0*exp(a1*x+a2*x*x);
};

//Функция для опраксимации sigma и mean;
double Func2(double x,double a0,double a1,double a2,double a3,double a4) 
{
    return (a0+a1/x+a2/(x*x))*exp(a3*x+a4*x*x);
};

void hTANA::Fit_imp(){
TF1 *f1=new TF1("fc1",FitGaus,-0.5,1.5,9);
TF1 *f2=new TF1("fc2",FitGaus,-0.5,1.5,9);
double conW,meanm2W,meanm3W,sigmW,sigm2W,sigm3W;
double conE,meanm2E,meanm3E,sigmE,sigm2E,sigm3E;
for(int n=0;n<pn;n++){
conW=hmassW[n]->Integral(1, hmassW[n]->GetNbinsX(), "width");conE=hmassE[n]->Integral(1, hmassE[n]->GetNbinsX(), "width");


sigmE=Func2(binP[n],0.050268,-0.0463449,0.0117124,0.109743,0.136159);
sigm2E=Func2(binP[n],0.0485043,-0.04704,0.0138879,0.623618,-0.0250298);
sigm3E=Func2(binP[n],0.0515644,-0.0575289,0.03656,0.720548,-0.0436915);


f1->SetParameters(Mpi,sigmE,Mk,sigm2E,Mpr,sigm3E,conW,0.5,0.3);
//f1->SetParLimits(0,0,1.2*Mpi);f1->SetParLimits(2,0,1.2*Mk);f1->SetParLimits(0,0,1.2*Mpr);

f2->SetParameters(Mpi,sigmE,Mk,sigm2E,Mpr,sigm3E,conE,0.5,0.3);
//f2->SetParLimits(0,0,1.2*Mpi);f2->SetParLimits(2,0,1.2*Mk);f2->SetParLimits(0,0,1.2*Mpr);

FitW[n]=(TH1F*)hmassW[n]->Clone();FitW[n]->Fit(f1,"R");FitE[n]=(TH1F*)hmassE[n]->Clone();FitE[n]->Fit(f2,"R");

sigW[n]=f1->GetParameter(1);sigE[n]=f2->GetParameter(1);
sig2W[n]=f1->GetParameter(3);sig2E[n]=f2->GetParameter(3);
sig3W[n]=f1->GetParameter(5);sig3E[n]=f2->GetParameter(5);

    }
}


void hTANA::ana_end(const char *outfile) {

Fit_imp();
TGraph *grsig1W = new TGraph (pn,binP, sigW);grsig1W->SetName("gr_sigW"); grsig1W->GetYaxis()->SetTitle("<Sigma>");grsig1W->GetXaxis()->SetTitle("Pt,Gev/c");
TGraph *grsig1E = new TGraph (pn,binP, sigE);grsig1E->SetName("gr_sigE"); grsig1E->GetYaxis()->SetTitle("<Sigma>");grsig1E->GetXaxis()->SetTitle("Pt,Gev/c");
grsig1W->SetLineColor(1);grsig1E->SetLineColor(3);grsig1W->SetTitle("sigma for Pion, West");grsig1E->SetTitle("sigma for Pion, East");

TGraph *grsig2W = new TGraph (pn,binP, sig2W);grsig2W->SetName("gr_sig2W");// TGraph *grmean2W = new TGraph (pn,binP, mean2W);grmean2W->SetName("gr_mean2W");
TGraph *grsig3W = new TGraph (pn,binP, sig3W);grsig3W->SetName("gr_sig3W"); //TGraph *grmean3W = new TGraph (pn,binP, mean3W);grmean3W->SetName("gr_mean3W");
grsig2W->SetTitle("sigma for Kaon, West");grsig3W->SetTitle("sigma for Proton, West");
TGraph *grsig2E = new TGraph (pn,binP, sig2E);grsig2E->SetName("gr_sig2E"); //TGraph *grmean2E = new TGraph (pn,binP, mean2E);grmean2E->SetName("gr_mean2E");
TGraph *grsig3E = new TGraph (pn,binP, sig3E);grsig3E->SetName("gr_sig3E");// TGraph *grmean3E = new TGraph (pn,binP, mean3E);grmean3E->SetName("gr_mean3E");
grsig2E->SetTitle("sigma for Kaon, East");grsig3E->SetTitle("sigma for Proton, East");
TCanvas *MyW = new TCanvas ("CanSigma","Test canvasqq",1);
MyW -> Divide(3,2);
  MyW -> cd(1);grsig1W ->Draw();
  MyW -> cd(4);grsig1E -> Draw();
  MyW -> cd(2);grsig2W ->Draw();
  MyW -> cd(5);grsig2E -> Draw();
  MyW -> cd(3);grsig3W ->Draw();
  MyW -> cd(6);grsig3E -> Draw();

TCanvas *MyC = new TCanvas ("CanData","Test canvas",1);
  MyC -> Divide(3,2);
  MyC -> cd(1);hm2E ->Draw();
  MyC -> cd(4);hm2W -> Draw();
  MyC -> cd(2);hpidE->Draw("colz");
  MyC -> cd(5);hpidW->Draw("colz");
  MyC -> cd(3);htofE -> Draw();
  MyC -> cd(6);htofW -> Draw();

TCanvas *MyC2W = new TCanvas ("CanFitW","Test canvasW",1);
  MyC2W -> Divide(3,2);
  MyC2W -> cd(1);FitW[3] ->Draw();
  MyC2W -> cd(2);FitW[6] -> Draw();
  MyC2W -> cd(3);FitW[9]->Draw();
  MyC2W -> cd(4);FitW[11]->Draw();
  MyC2W -> cd(5);FitW[12] -> Draw();
  MyC2W -> cd(6);FitW[13] -> Draw();

TCanvas *MyC2E = new TCanvas ("CanFitE","Test canvasE",1);
  MyC2E -> Divide(3,2);
  MyC2E -> cd(1);FitE[3] ->Draw();
  MyC2E -> cd(2);FitE[6] -> Draw();
  MyC2E -> cd(3);FitE[9]->Draw();
  MyC2E -> cd(4);FitE[11]->Draw();
  MyC2E -> cd(5);FitE[12] -> Draw();
  MyC2E -> cd(6);FitE[13] -> Draw();
d_outfile = new TFile(outfile,"recreate");
d_outfile->cd();
MyW->Write();
//MyC->Write();
MyC2W->Write();MyC2E->Write();
grsig1W->Write();grsig1E->Write();
grsig2W->Write(); grsig3W->Write();
//grmean2W->Write(); grmean3W->Write();
grsig2E->Write(); grsig3E->Write();
//grmean2E->Write(); grmean3E->Write();

//hpidW->Write();hpidE->Write();
//htofW->Write();htofE->Write();
//dtimW->Write();dtimE->Write();
//hm2W->Write();hm2E->Write();
//hbbcz->Write();
for(int n=0;n<pn;n++){
hmassW[n]->Write();
hmassE[n]->Write();}
d_outfile->Close();}

