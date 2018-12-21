
void efficiency_calculation_6(TString fname="svdmumu_l6.root") {
TFile *file=new TFile("histograms/final_histograms_l6_data.root","RECREATE");

TString a="soumen_svdResolution";
TFile *f=new TFile(fname);
TNtuple *tre = (TNtuple*)f->Get(a);
int nentries=tre->GetEntries();

TString b="intercept";
TNtuple *tree = (TNtuple*)f->Get(b);
int nentries1=tree->GetEntries();
 
cout<<nentries<<" "<<nentries1<<endl;

file->cd();

int count_cluster[5][2]={0};
int count_intercept[5]={0};


Float_t L,S,PXD,SVD;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;

//   cout<<L_<<endl;   
   if(L_!=6) continue;
   
//     cout<<i<<endl;
   if(PXD_>0 && SVD_>0){
     count_intercept[S_-1]++;
     }
  }



float error,cls_error,intercept,cls,truehit;
float size,svd,cdc,pxd,layer,sensor,side;
tre->SetBranchAddress("size",&size);
tre->SetBranchAddress("PXD",&pxd);
tre->SetBranchAddress("SVD",&svd);
tre->SetBranchAddress("CDC",&cdc);
tre->SetBranchAddress("intercept_error",&error);
tre->SetBranchAddress("cluster_position_error",&cls_error);
tre->SetBranchAddress("intercept_position",&intercept);
tre->SetBranchAddress("cluster_position",&cls);
tre->SetBranchAddress("layer",&layer);
tre->SetBranchAddress("sensor",&sensor);
tre->SetBranchAddress("side",&side);
tre->SetBranchAddress("truehit_pos", &truehit);
//tre->SetBranchAddress("num", &num);
TH1F *residual[5][2];
for(int p=0;p<5;p++)
  for(int t=0;t<2;t++){
    TString s= "sensor";
    s+=(p+1);
    s+="_side_";
    s+=(t);
    s+="_size_";
    TString a,b;
    a=s;b=s;
    s+="cls-intercept",
    residual[p][t]= new TH1F(s,s,100,-0.02,0.02);
    residual[p][t]->GetXaxis()->SetTitle("residual(intercept-cluster) in cm");
    residual[p][t]->GetXaxis()->SetTitleOffset(0.5);
    residual[p][t]->GetXaxis()->SetTitleSize(0.07);
    residual[p][t]->SetTitle("Residual Distribution");
}



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(layer!=6) continue;
   if(svd<1|| pxd<1) continue;
   if(TMath::Abs(cls-intercept)<0.05) count_cluster[sensor_-1][side]++;
   residual[sensor_-1][side_]->Fill(cls-intercept);
  }

  gStyle->SetOptFit(kTRUE);

  cout<<"intercept:61= "<<count_intercept[0]<<endl;
  cout<<"cluster:61V= "<<count_cluster[0][0]<<"  cluster:61U= "<<count_cluster[0][1]<<endl;
  cout<<"intercept:62= "<<count_intercept[1]<<endl;
  cout<<"cluster:62V= "<<count_cluster[1][0]<<"  cluster:62U= "<<count_cluster[1][1]<<endl;
  cout<<"intercept:63= "<<count_intercept[2]<<endl;
  cout<<"cluster:63V= "<<count_cluster[2][0]<<"  cluster:63U= "<<count_cluster[2][1]<<endl;
  cout<<"intercept:64= "<<count_intercept[3]<<endl;
  cout<<"cluster:64V= "<<count_cluster[3][0]<<"  cluster:64U= "<<count_cluster[3][1]<<endl;
  cout<<"intercept:65= "<<count_intercept[4]<<endl;
  cout<<"cluster:65V= "<<count_cluster[4][0]<<"  cluster:65U= "<<count_cluster[4][1]<<endl;












file->Write();
}



