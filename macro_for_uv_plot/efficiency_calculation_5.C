
void efficiency_calculation_5(TString fname="../Input_rootfile/data_bhabha_l5.root") {
TString b="intercept";
TFile *f=new TFile(fname);
TNtuple *tree = (TNtuple*)f->Get(b);
int nentries1=tree->GetEntries();


Float_t L,S,PXD,SVD,u_id,v_id,U,V;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
tree->SetBranchAddress("u_id",&u_id);
tree->SetBranchAddress("v_id",&v_id);
tree->SetBranchAddress("U",&U);
tree->SetBranchAddress("V",&V);
TH2F *uv = new TH2F("L5","L5",30,-4,4,150,-30,18);
for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;
   if(L_!=5) continue;
   if(S==1) V=V+11;
   if(S==3) V=V-12;
   if(S==4) V=V-24;
   if(!(PXD_>0 && SVD_>0)) continue;
   uv->Fill(U,V); 
  }

  TCanvas *c1= new TCanvas("c1","c1",300,1000);
  c1->cd();
  uv->SetStats(0);
  uv->Draw("colz");
}



