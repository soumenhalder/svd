
void efficiency_calculation_4(TString fname="../Input_rootfile/MC_l4_nocut.root") {
TFile *file=new TFile("histograms/final_histograms_l4_data.root","RECREATE");

TString a="soumen_svdResolution";
TFile *f=new TFile(fname);
TNtuple *tre = (TNtuple*)f->Get(a);
int nentries=tre->GetEntries();

TString b="intercept";
TNtuple *tree = (TNtuple*)f->Get(b);
int nentries1=tree->GetEntries();
 
cout<<nentries<<" "<<nentries1<<endl;

file->cd();

int count_cluster[3][2]={0};
int count_intercept[3]={0};
int count_clustersize[10][2]={0};

Float_t L,S,PXD,SVD,u_id,v_id,U,V;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
tree->SetBranchAddress("u_id",&u_id);
tree->SetBranchAddress("v_id",&v_id);
tree->SetBranchAddress("U",&U);
tree->SetBranchAddress("V",&V);


TH2F *UV[3][1];
for(int p=0;p<3;p++){
   TString s= "sensor";
    s+=(p+1);
    s+="uv";
    UV[p][0]= new TH2F(s,s,100,-5,5,100,-7,7);
    

}

for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;
   if(L_!=4) continue;
   UV[S-1][0]->Fill(U,V); 
   if(PXD_>0 && SVD_>0){
     count_intercept[S_-1]++;
     }
  }



float error,cls_error,intercept,cls,truehit,angle;
float size,svd,cdc,pxd,layer,sensor,side;
tre->SetBranchAddress("size",&size);
tre->SetBranchAddress("angle",&angle);
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

TH1F *residual[3][2][10];
TH1F *angleuv[3][2];
TH2F *residual_clspos[3][2];
//TH1F *true_residual[2][2];
for(int p=0;p<3;p++)
  for(int t=0;t<2;t++){
    for(int q=0;q<10;q++){
    TString s= "sensor";
    s+=(p+1);
    s+="_side_";
    s+=(t);
    s+="_size_";
    s+=(q+1);
    TString a,b;
    a=s;b=s;
    s+="cls-intercept";
    a+="residual_vs_clspos";
    residual[p][t][q]= new TH1F(s,s,100,-0.02,0.02);
    residual[p][t][q]->GetXaxis()->SetTitle("residual(intercept-cluster) in cm");
    residual[p][t][q]->GetXaxis()->SetTitleOffset(0.5);
    residual[p][t][q]->GetXaxis()->SetTitleSize(0.07);
    residual[p][t][q]->SetTitle("Residual Distribution");

   
   }
  TString s= "sensor";
    s+=(p+1);
    s+="_side_";
    s+=(t);
    TString a;
    a=s;
    s+="_angle_";
  angleuv[p][t]= new TH1F(s,s,60,0,180);
   angleuv[p][t]->GetXaxis()->SetTitle("angle in degree");
   angleuv[p][t]->GetXaxis()->SetTitleOffset(0.5);
   angleuv[p][t]->GetXaxis()->SetTitleSize(0.07);
   angleuv[p][t]->SetTitle("Incident angle"); 
   
   a+="residual vs cls pos"; 
   residual_clspos[p][t]= new TH2F(a,a,100,-0.02,0.02,100,-7,7);
//   residual_clspos[p][t]>GetXaxis()->SetTitle("residual(intercept-cluster) in cm");
//    residual_clspos[p][t]>GetXaxis()->SetTitle("cluster position in cm");
//   residual_clspos[p][t]>SetTitle("Residual vs cluster position correlation Distribution");

}



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(layer!=4) continue;
   if(size_>9) count_clustersize[9][side_]++;
   else count_clustersize[size_-1][side_]++;
   if(size_>9) size_=10;
   if(svd<1|| pxd<1) continue;
   if(TMath::Abs(cls-intercept)<0.05) count_cluster[sensor_-1][side]++;
   residual[sensor_-1][side_][size_-1]->Fill(cls-intercept);
   angleuv[sensor_-1][side_]->Fill(angle);
   residual_clspos[sensor_-1][side_]->Fill(cls-intercept,cls);
 //  true_residual[sensor_-1][side_]->Fill(cls-truehit);
    
  }

gStyle->SetOptFit(kTRUE);
   cout<<"cluster:51V= "<<count_cluster[0][0]<<"  cluster:51U= "<<count_cluster[0][1]<<endl;
   cout<<"intercept:52= "<<count_intercept[1]<<endl;
   cout<<"cluster:52V= "<<count_cluster[1][0]<<"  cluster:52U= "<<count_cluster[1][1]<<endl;
   cout<<"intercept:53= "<<count_intercept[2]<<endl;
   cout<<"cluster:53V= "<<count_cluster[2][0]<<"  cluster:53U= "<<count_cluster[2][1]<<endl;
   cout<<"intercept:54= "<<count_intercept[3]<<endl;
   cout<<"cluster:54V= "<<count_cluster[3][0]<<"  cluster:54U= "<<count_cluster[3][1]<<endl


  for(int i=0;i<10;i++) 
     cout<<i+1<<"  V="<<count_clustersize[i][0]<<"  :  U="<<count_clustersize[i][1]<<endl;
file->Write();
}



