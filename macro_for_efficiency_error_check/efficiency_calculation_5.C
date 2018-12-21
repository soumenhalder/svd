
void efficiency_calculation_5(TString fname="../Input_rootfile/data_bhabha_l5.root") {
TFile *file=new TFile("../histograms/final_histograms_l5_data_bhabha_modified.root","RECREATE");

TString a="soumen_svdResolution";
TFile *f=new TFile(fname);
TNtuple *tre = (TNtuple*)f->Get(a);
int nentries=tre->GetEntries();

TString b="intercept";
TNtuple *tree = (TNtuple*)f->Get(b);
int nentries1=tree->GetEntries();
 
cout<<nentries<<" "<<nentries1<<endl;

file->cd();

int count_cluster[4][2]={0};
int count_intercept[4]={0};
int count_background_cluster[4][2]={0};


Float_t L,S,PXD,SVD,u_id,v_id,U,V,U_ID,V_ID;;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
tree->SetBranchAddress("U",&U);
tree->SetBranchAddress("V",&V);
tree->SetBranchAddress("u_id",&U_ID);
tree->SetBranchAddress("v_id",&V_ID);
for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;
    if(!(U_ID>20 && U_ID<748 && V_ID>20 && V_ID<492)) continue;
   if(U_ID>108 && U_ID<276) continue;
   if(L_!=5) continue;
   if(PXD_>0 && SVD_>0){
     count_intercept[S_-1]++;
     }
  }


float error,cls_error,intercept,cls,truehit;
float size,svd,cdc,pxd,layer,sensor,side,,u_id,v_id;
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
tre->SetBranchAddress("u_id", &u_id);
tre->SetBranchAddress("v_id", &v_id);

//tre->SetBranchAddress("num", &num);
TH1F *residual[4][2];
for(int p=0;p<4;p++)
  for(int t=0;t<2;t++){
    TString s= "sensor";
    s+=(p+1);
    s+="_side_";
    s+=(t);
    s+="_size_";
    TString a,b;
    a=s;b=s;
    s+="cls-intercept",
    residual[p][t]= new TH1F(s,s,700,-0.3,0.3);
    residual[p][t]->GetXaxis()->SetTitle("residual(intercept-cluster) in cm");
    residual[p][t]->GetXaxis()->SetTitleOffset(0.5);
    residual[p][t]->GetXaxis()->SetTitleSize(0.07);
    residual[p][t]->SetTitle("Residual Distribution");
}



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
//    cout<<u_id<<endl;
  if(u_id>108 && u_id<276) continue;
   if(!(u_id>20 && u_id<748 && v_id>20 && v_id<492)) continue;

   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(layer!=5) continue;
   if(svd>0 && pxd>0){
     if(TMath::Abs(cls-intercept)<0.05) count_cluster[sensor_-1][side]++;
     if(TMath::Abs(cls-intercept)>0.05 && TMath::Abs(cls-intercept)<0.1) count_background_cluster[sensor_-1][side]++;
     residual[sensor_-1][side_]->Fill(cls-intercept);
      }
  }

gStyle->SetOptFit(kTRUE);

 cout<<"intercept:51= "<<count_intercept[0]<<endl;
  cout<<"cluster:51V= "<<count_cluster[0][0]<<"  cluster:51U= "<<count_cluster[0][1]<<endl;
  cout<<"background cluster:51V= "<<count_background_cluster[0][0]<<"  cluster:51U= "<<count_background_cluster[0][1]<<endl; 
  cout<<"intercept:52= "<<count_intercept[1]<<endl;
  cout<<"cluster:52V= "<<count_cluster[1][0]<<"  cluster:52U= "<<count_cluster[1][1]<<endl;
  cout<<"background cluster:52V= "<<count_background_cluster[1][0]<<"  cluster:52U= "<<count_background_cluster[1][1]<<endl; 
  cout<<"intercept:53= "<<count_intercept[2]<<endl;
  cout<<"cluster:53V= "<<count_cluster[2][0]<<"  cluster:53U= "<<count_cluster[2][1]<<endl;
  cout<<"background cluster:53V= "<<count_background_cluster[2][0]<<"  cluster:53U= "<<count_background_cluster[2][1]<<endl; 
  cout<<"intercept:54= "<<count_intercept[3]<<endl;
  cout<<"cluster:54V= "<<count_cluster[3][0]<<"  cluster:54U= "<<count_cluster[3][1]<<endl;
  cout<<"background cluster:54V= "<<count_background_cluster[3][0]<<"  cluster:54U= "<<count_background_cluster[3][1]<<endl; 

   cout<<count_intercept[0]<<" "<<count_cluster[0][0]<<" "<<count_background_cluster[0][0]<<" "<<count_cluster[0][1]<<" "<<count_background_cluster[0][1]<<endl;
  cout<<count_intercept[1]<<" "<<count_cluster[1][0]<<" "<<count_background_cluster[1][0]<<" "<<count_cluster[1][1]<<" "<<count_background_cluster[1][1]<<endl;
  cout<<count_intercept[2]<<" "<<count_cluster[2][0]<<" "<<count_background_cluster[2][0]<<" "<<count_cluster[2][1]<<" "<<count_background_cluster[2][1]<<endl;
  cout<<count_intercept[3]<<" "<<count_cluster[3][0]<<" "<<count_background_cluster[3][0]<<" "<<count_cluster[3][1]<<" "<<count_background_cluster[3][1]<<endl;
  
file->Write();
}



