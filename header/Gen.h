#include <TChain.h>

class Gen{
public:
  double mu_Px[3];
  double mu_Py[3];
  double mu_Pz[3];
  double mu_E[3];
  double nu_Px;
  double nu_Py;
  double nu_Pz;
  double nu_E;
  double W_Px;
  double W_Py;
  double W_Pz;
  double W_E;
  double HN_Px;
  double HN_Py;
  double HN_Pz;
  double HN_E;
  //only for Wtollln
  //double gamma_star_found;
  //double gamma_star_mass;
  
  Gen(TChain *chain);
};

Gen::Gen(TChain *chain){
  TString varname[24] = {"mu_0_Px", "mu_0_Py", "mu_0_Pz", "mu_0_E", "mu_1_Px", "mu_1_Py", "mu_1_Pz", "mu_1_E", "mu_2_Px", "mu_2_Py", "mu_2_Pz", "mu_2_E", "nu_Px", "nu_Py", "nu_Pz", "nu_E", "W_Px", "W_Py", "W_Pz", "W_E", "HN_Px", "HN_Py", "HN_Pz", "HN_E"};
  
  for(int i=0; i<3; i++){
    chain->SetBranchAddress(varname[0+4*i], mu_Px+i);
    chain->SetBranchAddress(varname[1+4*i], mu_Py+i);
    chain->SetBranchAddress(varname[2+4*i], mu_Pz+i);
    chain->SetBranchAddress(varname[3+4*i], mu_E+i);
  }
  chain->SetBranchAddress(varname[12], &nu_Px);
  chain->SetBranchAddress(varname[13], &nu_Py);
  chain->SetBranchAddress(varname[14], &nu_Pz);
  chain->SetBranchAddress(varname[15], &nu_E);
  chain->SetBranchAddress(varname[16], &W_Px);
  chain->SetBranchAddress(varname[17], &W_Py);
  chain->SetBranchAddress(varname[18], &W_Pz);
  chain->SetBranchAddress(varname[19], &W_E);
  chain->SetBranchAddress(varname[20], &HN_Px);
  chain->SetBranchAddress(varname[21], &HN_Py);
  chain->SetBranchAddress(varname[22], &HN_Pz);
  chain->SetBranchAddress(varname[23], &HN_E);
  
  //only for Wtollln
  //chain->SetBranchAddress("gamma_star_found", &gamma_star_found);
  //chain->SetBranchAddress("gamma_star_mass", &gamma_star_mass);
  
}
