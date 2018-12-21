#include<cmath>
/*412V 413V 411V*/
void fit_residualplot_4(TString fname="../Input_rootfile/data_bhabha_l4.root") {
using namespace RooFit;
RooRealVar x("resolution", " residual in cm",-0.05, 0.05);
RooDataSet* data = new RooDataSet("data", "data",RooArgSet(x));


TString a="soumen_svdResolution";
TFile *f=new TFile(fname);
TString c=a+fname;
TDirectory* dir= (TDirectory*)f->Get("ntuple");
TNtuple *tre = (TNtuple*)f->Get(a);
int nentries=tre->GetEntries();
int eff1=0,eff2=0;

//file->cd();
float error,cls_error,intercept,cls,truehit;
float size,svd,cdc,pxd,layer,sensor,side;
tre->SetBranchAddress("size",&size);
tre->SetBranchAddress("PXD",&pxd);
tre->SetBranchAddress("SVD",&svd);
tre->SetBranchAddress("CDC",&cdc);
tre->SetBranchAddress("intercept_position",&intercept);
tre->SetBranchAddress("cluster_position",&cls);
tre->SetBranchAddress("layer",&layer);
tre->SetBranchAddress("sensor",&sensor);
tre->SetBranchAddress("side",&side);



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(TMath::Abs(intercept-cls)>0.040) continue;
  if(layer==4 && sensor==1 && svd>0 && pxd>0 && side==0){ 
     x.setVal(intercept-cls); 
     data->add(RooArgSet(x));
    }
    // cout<<x<<endl;  
}

//
// no of entries in the datafile
RooPlot *myframe = x.frame(Name("411V"),Title("411V residual distribution"),Bins(60));
data->plotOn(myframe);
RooRealVar mean1("#mu_{1}","mean_{1}",0.000,-0.002,0.002);
RooRealVar mean2("#mu_{2}","mean_{2}",0.001,-0.01,0.01);
RooRealVar sigma1("#sigma_{1}","sigma1",0.002,0,0.1);
RooRealVar sigma2("#sigma_{2}","sigma2",0.01,0,0.1);


RooGaussian signal1("signal1","signal1",x,mean1,sigma1);
RooGaussian signal2("signal2","signal2",x,mean2,sigma2);

RooRealVar sfrac("sfrac","Area fraction",0.6,0.0 ,01);
RooRealVar sfrac2("sfrac2","Area fraction2",0.9,0 ,01);

RooAddPdf signal("signal","model",RooArgList(signal1,signal2),sfrac);

signal.fitTo(*data);

//cout<<"sigma"<<sfrac2;//+(1-sfrac2)*sigma3*sigma3<<endl;

signal.paramOn(myframe);
TCanvas *c1 =new TCanvas("res","res",1000,1000);
c1->cd();

signal.plotOn(myframe);
cout << "chi^2 = " << myframe->chiSquare() << endl ;
signal.plotOn(myframe, Components("signal1"), LineColor(kRed),LineStyle(kDashed), LineWidth(2));
signal.plotOn(myframe, Components("signal2"), LineColor(kGreen),LineStyle(kDashed), LineWidth(2));
myframe->GetYaxis()->SetTitleOffset(1.5);
myframe->Draw();
}
