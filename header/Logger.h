#include <fstream>
#include <ctime>
#include <TString.h>
#include <TH1D.h>
#include <TMath.h>
#include <vector>

using namespace std;

class Logger{
public:
  TString folder_path, logfile_path;
  vector<double> N_entry_sec, Error_entry_sec;
  
  Logger(TString path);
  void start_log(int sample, int isPU, int ischi2, int additional_cut);
  void write_n_events(TString sample_name, double weight, int event, double PUweight);
  void save_n_entry_sec(double cont, double error);
  void write_n_entry_sec(vector<TString> sector);
  void write_onebin(TH1D *hist_sig, TH1D *hist_dat, TH1D *hist_bkg);
  void divide_log();
  TString GetLogfilePath();
};

Logger::Logger(TString path){
  folder_path = path;
}

void Logger::start_log(int sample, int isPU, int ischi2, int additional_cut){
  logfile_path = folder_path+TString::Itoa(sample, 10)+"_log.log";
  ofstream logfile(logfile_path, ios::trunc);
  const time_t ctt = time(0);
  logfile
  << "[TIME] " << asctime(localtime(&ctt)) << endl
  << "===============================================" << endl
  << "40, 50, 60, 100(no gen lvl info)" << endl
  << "=> " << sample << endl
  << "PileUP? (0 = no, 1 = yes)" << endl
  << "=> " << isPU << endl
  << "chi2 fit? (0 = no, 1 = yes)" << endl
  << "=> " << ischi2 << endl
  << "===============================================" << endl
  << "Additional Cut?" << endl
  << "  0 : No Additional Cut" << endl
  << "  1 : delta R > 0.5" << endl
  << "  2 : delta R > 0.5 & W_reco < 100 GeV" << endl
  << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
  << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda =0" << endl
  << "=> " << additional_cut << endl
  << "===============================================" << endl;
  logfile.close();
}

void Logger::write_n_events(TString sample_name, double weight, int event, double PUweight){
  ofstream logfile;
  logfile.open(logfile_path, ios::app);
  double Entry = weight*event;
  double Entry_PU = weight*event*PUweight;
  if(TMath::IsNaN(Entry)) Entry = 0;
  if(TMath::IsNaN(Entry_PU)) Entry_PU = 0;
  if(TMath::IsNaN(PUweight)) PUweight = 0;
  logfile
  << sample_name << "\t" << event << "\t" << weight << "\t" << PUweight << "\t" << Entry << "\t" << Entry_PU << endl;
  logfile.close();
}

void Logger::save_n_entry_sec(double cont, double error){
  N_entry_sec.push_back(cont);
  Error_entry_sec.push_back(error);
}

void Logger::write_n_entry_sec(vector<TString> sector){
  ofstream logfile;
  logfile.open(logfile_path, ios::app);
  if(sector.size() != N_entry_sec.size()){
    cout << "sector.size() != N_entry_sec.size()" << endl;
    logfile
    << "sector.size() != N_entry_sec.size()" << endl
    << "===============================================" << endl;
    logfile.close();
  }
  else{
    for(unsigned int i=0; i<N_entry_sec.size(); i++){
      logfile
      << sector.at(i) << "\t" << N_entry_sec.at(i) << "\t" << Error_entry_sec.at(i) << endl;
    }
    logfile
    << "===============================================" << endl;
    logfile.close();
  }
  
}

void Logger::divide_log(){
  ofstream logfile;
  logfile.open(logfile_path, ios::app);
  logfile
  << "===============================================" << endl;
  logfile.close();
}

void Logger::write_onebin(TH1D *hist_sig, TH1D *hist_dat, TH1D *hist_bkg){
  ofstream logfile;
  logfile.open(logfile_path, ios::app);
  logfile
  << "(# of events) for ..." << endl
  << '\t' << "-signal" << '\t' << hist_sig->GetBinContent(1) << endl
  << '\t' << "-data" << '\t' << hist_dat->GetBinContent(1) << endl
  << '\t' << "-bkg" << '\t' << hist_bkg->GetBinContent(1) << endl
  << "(Error) for ..." << endl
  << '\t' << "-signal" << '\t' << hist_sig->GetBinError(1) << endl
  << '\t' << "-data" << '\t' << hist_dat->GetBinError(1) << endl
  << '\t' << "-bkg" << '\t' << hist_bkg->GetBinError(1) << endl
  << "===============================================" << endl;
  logfile.close();
}

TString Logger::GetLogfilePath(){
  return logfile_path;
}