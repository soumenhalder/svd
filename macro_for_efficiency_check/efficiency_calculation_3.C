
void efficiency_calculation_3(TString fname="../Input_rootfile/data_bhabha_l3.root") {
TFile *file=new TFile("../histograms/final_histograms_l3_data_bhabha_for_subtracting.root","RECREATE");

TString a="soumen_svdResolution";
TFile *f=new TFile(fname);
TNtuple *tre = (TNtuple*)f->Get(a);
int nentries=tre->GetEntries();

TString b="intercept";
TNtuple *tree = (TNtuple*)f->Get(b);
int nentries1=tree->GetEntries();
 
cout<<nentries<<" "<<nentries1<<endl;

file->cd();

int count_cluster[2][2]={0};
int count_intercept[2]={0};


Float_t L,S,PXD,SVD,u_id,v_id,U,V;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
tree->SetBranchAddress("u_id",&u_id);
tree->SetBranchAddress("v_id",&v_id);
tree->SetBranchAddress("U",&U);
tree->SetBranchAddress("V",&V);

TH2F *uv=new TH2F("312uv","312uv",50,-4,4,50,-8,8);
for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;
   if(L_!=3) continue;
   if(S==2) uv->Fill(U,V); 
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
TH1F *residual[2][2];
//TH1F *true_residual[2][2];
for(int p=0;p<2;p++)
  for(int t=0;t<2;t++){
    TString s= "sensor";
    s+=(p+1);
    s+="_side_";
    s+=(t);
    s+="_size_";
    TString a,b;
    a=s;b=s;
    s+="cls-intercept";
    a+="cls-truehit";
    residual[p][t]= new TH1F(s,s,100,-1,1);
    residual[p][t]->GetXaxis()->SetTitle("residual(intercept-cluster) in cm");
    residual[p][t]->GetXaxis()->SetTitleOffset(0.5);
    residual[p][t]->GetXaxis()->SetTitleSize(0.07);
    residual[p][t]->SetTitle("Residual Distribution");

 //   true_residual[p][t]= new TH1F(a,a,100,-0.02,0.02);
 //   true_residual[p][t]->GetXaxis()->SetTitle("residual(truehit-cluster) in cm");
 //   true_residual[p][t]->GetXaxis()->SetTitleOffset(0.5);
 //   true_residual[p][t]->GetXaxis()->SetTitleSize(0.07);
 //   true_residual[p][t]->SetTitle("true Residual Distribution");
}



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(layer!=3) continue;
   if(svd>0 && pxd>0){
     if(TMath::Abs(cls-intercept)<0.05) count_cluster[sensor_-1][side]++;
     residual[sensor_-1][side_]->Fill(cls-intercept);
     //true_residual[sensor_-1][side_]->Fill(cls-truehit);
     }
  }

gStyle->SetOptFit(kTRUE);

 cout<<"intercept:31= "<<count_intercept[0]<<endl;
  cout<<"cluster:31V= "<<count_cluster[0][0]<<"  cluster:31U= "<<count_cluster[0][1]<<endl;
  cout<<"intercept:32= "<<count_intercept[1]<<endl;
  cout<<"cluster:32V= "<<count_cluster[1][0]<<"  cluster:32U= "<<count_cluster[1][1]<<endl;
file->Write();
}



