#include "PndTutAnaTask.h" 
class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

// **** some auxilliary functions in auxtut.C ****
// - FairRunAna* initrun(TString prefix, TString outfile, int min=-1, int max=-1) --> Init FairRunAna
// - initPDG()                     --> Initialize TDatabasePDG with EvtGen data; returns pointer to singleton
// - plotmyhistos()                --> Plots all histograms in current TDirectory on a autosized canvas
// - writemyhistos()               --> Writes all histos in current TFile
// - fillM(RhoCandList l, TH1* h)  --> Fill mass histogram h with masses of candidates in l
// - fillP(RhoCandList &l, TH1* h) --> Fill momentum histogram h with momenta of candidates in l
// - RemoveGeoManager()            --> Temporary fix for error on macro exit
// **** some auxilliary functions in auxtut.C ****
#include "auxtut.C"

// *** routine to only keep PID matched candidates in list
int SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
  int removed = 0;

  for (int ii = l.GetLength() - 1; ii >= 0; --ii) {
    if (!(ana->McTruthMatch(l[ii]))) {
      l.Remove(l[ii]);
      removed++;
    }
  }

  return removed;
}

void tut_ana(int nevts = 0, TString prefix = "D0_kpi")
{
  // *** some variables
  int i = 0, j = 0, k = 0, l = 0;
  gStyle->SetOptFit(1011);

  // *** the output file for FairRunAna
  TString OutFile = "out_dummy.root";

  // *** the files coming from the simulation
  TString inPidFile = prefix + "_pid.root"; // this file contains the PndPidCandidates and McTruth
  TString inParFile = prefix + "_par.root";

  // *** PID table with selection thresholds; can be modified by the user
  TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR")) + "/macro/params/all.par";

  // *** initialization
  // FairLogger::GetLogger()->SetLogToFile(kFALSE);
  FairRunAna *fRun = new FairRunAna();
  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  fRun->SetSource(new FairFileSource(inPidFile));

  // *** setup parameter database
  FairParRootFileIo *parIO = new FairParRootFileIo();
  parIO->open(inParFile);
  FairParAsciiFileIo *parIOPid = new FairParAsciiFileIo();
  parIOPid->open(pidParFile.Data(), "in");

  rtdb->setFirstInput(parIO);
  rtdb->setSecondInput(parIOPid);
  rtdb->setOutput(parIO);

  fRun->SetSink(new FairRootFileSink(OutFile));
  fRun->Init();
  
  // *** intialize TDatabasePDG with EvtGen names 
  TDatabasePDG::Instance()->ReadPDGTable(TString(gSystem->Getenv("VMCWORKDIR"))+"/input/pdg_table_evtgen.txt");
  // *** create shortcut for TDatabasePDG::Instance()
  TDatabasePDG *pdb = TDatabasePDG::Instance();
  

  // *** create an output file for all histograms
  TFile *out = TFile::Open(prefix + "_ana.root", "RECREATE");

  // *** create some histograms
  TH1F *hD0m_all = new TH1F("hD0m_all", "D0 mass (all);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_all = new TH1F("hD0barm_all", "D0bar mass (all);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_all = new TH1F("hpsim_all", "interaction point mass (all);m [GeV/c^{2}];entries", 200, 1.5, 7.5);

  TH1F *hD0m_lpid = new TH1F("hD0m_lpid", "D0 mass (loose pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_lpid = new TH1F("hD0barm_lpid", "D0bar mass (loose pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_lpid = new TH1F("hpsim_lpid", "interaction point mass (loose pid);m [GeV/c^{2}];entries", 200, 2.5, 7.5);

  TH1F *hD0m_tpid = new TH1F("hD0m_tpid", "D0 mass (tight pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_tpid = new TH1F("hD0barm_tpid", "D0bar mass (tight pid);m [GeV/c^{2}];entries", 200, 0, 4.5); 
  TH1F *hpsim_tpid = new TH1F("hpsim_tpid", "interaction point mass (tight pid);m [GeV/c^{2}];entries", 200, 2.5, 7.5);

  TH1F *hD0m_pid04 = new TH1F("hD0m_pid04", "D0 mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_pid04 = new TH1F("hD0barm_pid04", "D0bar mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_pid04 = new TH1F("hpsim_pid04", "interaction point mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_trpid = new TH1F("hD0m_trpid", "D0 mass (true pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_trpid = new TH1F("hD0barm_trpid", "D0bar mass (true pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_trpid = new TH1F("hpsim_trpid", "interaction point mass (true pid);m [GeV/c^{2}];entries", 200, 2.5, 7.5);

  TH1F *hD0m_ftm = new TH1F("hD0m_ftm", "D0 mass (full truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_ftm = new TH1F("hD0barm_ftm", "D0bar mass (full truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_ftm = new TH1F("hpsim_ftm", "interaction point mass (full truth match);m [GeV/c^{2}];entries", 200, 2.5, 7.5);

  TH1F *hD0m_nm = new TH1F("hD0m_nm", "D0 mass (no truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_nm = new TH1F("hD0barm_nm", "D0bar mass (no truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_nm = new TH1F("hpsim_nm", "interaction point mass (no truth match);m [GeV/c^{2}];entries", 200, 1.5, 7.5);

  TH1F *hD0m_diff = new TH1F("hD0m_diff", "D0 mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, -2, 2);
  TH1F *hD0barm_diff = new TH1F("hD0barm_diff", "D0bar mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, -2, 2);
  TH1F *hpsim_diff = new TH1F("hpsim_diff", "interaction point mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, 1.5, 7.5);

  TH1F *hD0m_vf = new TH1F("hD0m_vf", "D0 mass (vertex fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0m_4cf = new TH1F("hD0m_4cf", "D0 mass (4C fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0m_mcf = new TH1F("hD0m_mcf", "D0 mass (mass constraint fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0m_dc = new TH1F("hD0m_dc", "D0 mass (decay tree fit);m [GeV/c^{2}];entries", 200, 0, 4.5);

  TH1F *hD0barm_vf = new TH1F("hD0barm_vf", "D0bar mass (vertex fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_4cf = new TH1F("hD0barm_4cf", "D0bar mass (4C fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_mcf = new TH1F("hD0barm_mcf", "D0bar mass (mass constraint fit);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_dc = new TH1F("hD0barm_dc", "D0bar mass (decay tree fit);m [GeV/c^{2}];entries", 200, 0, 4.5);

  TH1F *hD0_chi2_vf = new TH1F("hD0_chi2_vf", "D0: #chi^{2} vertex fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hD0bar_chi2_vf = new TH1F("hD0bar_chi2_vf", "D0bar: #chi^{2} vertex fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hpsi_chi2_4c = new TH1F("hpsi_chi2_4c", "interaction point: #chi^{2} 4C fit;#chi^{2};entries", 100, 0, 250);
  TH1F *hD0_chi2_mf = new TH1F("hD0_chi2_mf", "D0: #chi^{2} mass fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hD0bar_chi2_mf = new TH1F("hD0bar_chi2_mf", "D0bar: #chi^{2} mass fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hpsi_chi2_dc = new TH1F("hpsi_chi2_dc", "interaction point: #chi^{2} decay tree fit;#chi^{2};entries", 100, 0, 100);

  TH1F *hD0_prob_vf = new TH1F("hD0_prob_vf", "D0: Prob vertex fit;P;entries", 100, 0, 1);
  TH1F *hD0bar_prob_vf = new TH1F("hD0bar_prob_vf", "D0bar: Prob vertex fit;P;entries", 100, 0, 1);
  TH1F *hpsi_prob_4c = new TH1F("hpsi_prob_4c", "interaction point: Prob 4C fit;P;entries", 100, 0, 1);
  TH1F *hD0_prob_mf = new TH1F("hD0_prob_mf", "D0: Prob mass fit;P;entries", 100, 0, 1);
  TH1F *hD0bar_prob_mf = new TH1F("hD0bar_prob_mf", "D0bar: Prob mass fit;P;entries", 100, 0, 1);
  TH1F *hpsi_prob_dc = new TH1F("hpsi_prob_dc", "interaction point: Prob decay tree fit;P;entries", 100, 0, 1);

  TH2F *hvpos = new TH2F("hvpos", "(x,y) projection of fitted decay vertex;x [cm];y [cm]", 100, -2, 2, 100, -2, 2);
 
  TH1F *hmkp = new TH1F("hmkp", "PdgCode of the Mother Particle of kplus", 800, -450, 450);
  TH1F *hmkp_ = new TH1F("hmkp_", "PdgCode of kplus", 800, -450, 450);
  TH1F *hmkm = new TH1F("hmkm", "PdgCode of the Mother Particle of kminus", 800, -450, 450);
  TH1F *hmkm_ = new TH1F("hmkm_", "PdgCode of kminus", 800, -450, 450);
  TH1F *hmpp = new TH1F("hmpp", "PdgCode of the Mother Particle of piplus", 800, -450, 450);
  TH1F *hmpp_ = new TH1F("hmpp_", "PdgCode of piplus", 800, -450, 450);
  TH1F *hmpm = new TH1F("hmpm", "PdgCode of the Mother Particle of piminus", 800, -450, 450);
  TH1F *hmpm_ = new TH1F("hmpm_", "PdgCode of piminus", 800, -450, 450);

  TH1F *hmkp_lo = new TH1F("hmkp_lo", "PdgCode of the Mother Particle of kplus_lo", 800, -450, 450);
  TH1F *hmkp_l = new TH1F("hmkp_l", "PdgCode of kplus_lo", 800, -450, 450);
  TH1F *hmkm_lo = new TH1F("hmkm_lo", "PdgCode of the Mother Particle of kminus_lo", 800, -450, 450);
  TH1F *hmkm_l = new TH1F("hmkm_l", "PdgCode of kminus_lo", 800, -450, 450);
  TH1F *hmpp_lo = new TH1F("hmpp_lo", "PdgCode of the Mother Particle of piplus_lo", 800, -450, 450);
  TH1F *hmpp_l = new TH1F("hmpp_l", "PdgCode of piplus_lo", 800, -450, 450);
  TH1F *hmpm_lo = new TH1F("hmpm_lo", "PdgCode of the Mother Particle of piminus_lo", 800, -450, 450);
  TH1F *hmpm_l = new TH1F("hmpm_l", "PdgCode of piminus_lo", 800, -450, 450);


  TH1F *hmkp_ti = new TH1F("hmkp_ti", "PdgCode of the Mother Particle of kplus_ti", 800, -450, 450);
  TH1F *hmkp_t = new TH1F("hmkp_t", "PdgCode of kplus_ti", 800, -450, 450);
  TH1F *hmkm_ti = new TH1F("hmkm_ti", "PdgCode of the Mother Particle of kminus_ti", 800, -450, 450);
  TH1F *hmkm_t = new TH1F("hmkm_t", "PdgCode of kminus_ti", 800, -450, 450);
  TH1F *hmpp_ti = new TH1F("hmpp_ti", "PdgCode of the Mother Particle of piplus_ti", 800, -450, 450);
  TH1F *hmpp_t = new TH1F("hmpp_t", "PdgCode of piplus_ti", 800, -450, 450);
  TH1F *hmpm_ti = new TH1F("hmpm_ti", "PdgCode of the Mother Particle of piminus_ti", 800, -450, 450);
  TH1F *hmpm_t = new TH1F("hmpm_t", "PdgCode of piminus_ti", 800, -450, 450);

 
  TH1F *hD0_momentum = new TH1F("hD0_momentum", "the momentum of D0", 200, 0, 15);
  TH1F *hD0_costheta = new TH1F("hD0_costheta", "the cos theta distribution of D0", 200, -1.5, 1.5);
  TH1F *hD0_phi = new TH1F("hD0_phi", "the phi distribution of D0", 200, -3.5, 3.5);
  TH1F *hD0bar_momentum = new TH1F("hD0bar_momentum", "the momentum of D0bar", 200, 0, 15);
  TH1F *hD0bar_costheta = new TH1F("hD0bar_costheta", "the cos theta distribution of D0bar", 200, -1.5, 1.5);
  TH1F *hD0bar_phi = new TH1F("hD0barm_phi", "the phi distribution of D0bar", 200, -3.5, 3.5);

  TH1F *hkp_momentum = new TH1F("hkp_momentum", "the momentum of kplus", 200, 0, 15);
  TH1F *hkp_costheta = new TH1F("hkp_costheta", "the cos theta distribution of kplus", 200, -1.5, 1.5);
  TH1F *hkp_phi = new TH1F("hkp_phi", "the phi distribution of kplus", 200, -3.5, 3.5);
  TH1F *hpm_momentum = new TH1F("hpm_momentum", "the momentum of piminus", 200, 0, 15);
  TH1F *hpm_costheta = new TH1F("hpm_costheta", "the cos theta distribution of piminus", 200, -1.5, 1.5);
  TH1F *hpm_phi = new TH1F("hpm_phi", "the phi distribution of piminus", 200, -3.5, 3.5);

  TH1F *hkm_momentum = new TH1F("hkm_momentum", "the momentum of kminus", 200, 0, 15);
  TH1F *hkm_costheta = new TH1F("hkm_costheta", "the cos theta distribution of kminus", 200, -1.5, 1.5);
  TH1F *hkm_phi = new TH1F("hkm_phi", "the phi distribution of kminus", 200, -3.5, 3.5);
  TH1F *hpp_momentum = new TH1F("hpp_momentum", "the momentum of piplus", 200, 0, 15);
  TH1F *hpp_costheta = new TH1F("hpp_costheta", "the cos theta distribution of piplus", 200, -1.5, 1.5);
  TH1F *hpp_phi = new TH1F("hpp_phi", "the phi distribution of piplus", 200, -3.5, 3.5);

  TH1F *hD0_m_resolution = new TH1F("hD0_m_resolution", "resolution of the mass of D0", 200, -2, 2);
  TH1F *hD0bar_m_resolution = new TH1F("hD0bar_m_resolution", "resolution of the mass of D0bar", 200, -2, 2); 
  TH1F *hD0_p_resolution = new TH1F("hD0_p_resolution", "resolution of the momentum of D0", 200, -2, 2);
  TH1F *hD0bar_p_resolution = new TH1F("hD0bar_p_resolution", "resolution of the momentum of D0bar", 200, -2, 2); 
  TH1F *hD0_cos_resolution = new TH1F("hD0_cos_resolution", "resolution of cos theta of D0", 200, -2, 2);
  TH1F *hD0bar_cos_resolution = new TH1F("hD0bar_cos_resolution", "resolution of cos theta of D0bar", 200, -2, 2); 
  TH1F *hD0_phi_resolution = new TH1F("hD0_phi_resolution", "resolution of phi of D0", 200, -2, 2);
  TH1F *hD0bar_phi_resolution = new TH1F("hD0bar_phi_resolution", "resolution of phi of D0bar", 200, -2, 2); 

  TH1F *hkp_m_resolution = new TH1F("hkp_m_resolution", "resolution of the mass of kplus", 200, -2, 2);
  TH1F *hpm_m_resolution = new TH1F("hpm_m_resolution", "resolution of the mass of piminus", 200, -2, 2); 
  TH1F *hkp_p_resolution = new TH1F("hkp_p_resolution", "resolution of the momentum of kplus", 200, -2, 2);
  TH1F *hpm_p_resolution = new TH1F("hpm_p_resolution", "resolution of the momentum of piminus", 200, -2, 2); 
  TH1F *hkp_cos_resolution = new TH1F("hkp_cos_resolution", "resolution of cos theta of kplus", 200, -2, 2);
  TH1F *hpm_cos_resolution = new TH1F("hpm_cos_resolution", "resolution of cos theta of piminus", 200, -2, 2); 
  TH1F *hkp_phi_resolution = new TH1F("hkp_phi_resolution", "resolution of phi of kplus", 200, -2, 2);
  TH1F *hpm_phi_resolution = new TH1F("hpm_phi_resolution", "resolution of phi of piminus", 200, -2, 2);

  TH1F *hkm_m_resolution = new TH1F("hkm_m_resolution", "resolution of the mass of kminus", 200, -2, 2);
  TH1F *hpp_m_resolution = new TH1F("hpp_m_resolution", "resolution of the mass of piplus", 200, -2, 2); 
  TH1F *hkm_p_resolution = new TH1F("hkm_p_resolution", "resolution of the momentum of kminus", 200, -2, 2);
  TH1F *hpp_p_resolution = new TH1F("hpp_p_resolution", "resolution of the momentum of piplus", 200, -2, 2); 
  TH1F *hkm_cos_resolution = new TH1F("hkm_cos_resolution", "resolution of cos theta of kminus", 200, -2, 2);
  TH1F *hpp_cos_resolution = new TH1F("hpp_cos_resolution", "resolution of cos theta of piplus", 200, -2, 2); 
  TH1F *hkm_phi_resolution = new TH1F("hkm_phi_resolution", "resolution of phi of kminus", 200, -2, 2);
  TH1F *hpp_phi_resolution = new TH1F("hpp_phi_resolution", "resolution of phi of piplus", 200, -2, 2);

  TH1F *hD0tmp_mass = new TH1F("hD0tmp_mass", "mass of D0tmp", 800, -450, 450);
  TH1F *hD0tmp_p = new TH1F("hD0tmp_p", "momentum of D0tmp", 200, 0, 15);
  TH1F *hD0tmp_cos = new TH1F("hD0tmp_cos", "costheta of D0tmp", 200, -1.5, 1.5);
  TH1F *hD0tmp_phi = new TH1F("hD0tmp_mass", "mass of D0tmp", 200, -2, 2);

  TH1F *hD0bartmp_mass = new TH1F("hD0bartmp_mass", "mass of D0bartmp", 800, -450, 450);
  TH1F *hD0bartmp_p = new TH1F("hD0bartmp_p", "momentum of D0bartmp", 200, 0, 15);
  TH1F *hD0bartmp_cos = new TH1F("hD0bartmp_cos", "costheta of D0bartmp", 200, -1.5, 1.5);
  TH1F *hD0bartmp_phi = new TH1F("hD0bartmp_mass", "mass of D0bartmp", 200, -2, 2);




  //
  // Now the analysis stuff comes...
  //

  // *** the data reader object
  PndAnalysis *theAnalysis = new PndAnalysis();
  if (nevts == 0)
    nevts = theAnalysis->GetEntries();

 

  // *** Mass selector for the D0 cands
  double m0_D0 = pdb->GetParticle("D0")->Mass(); // Get nominal PDG mass of the D0
  RhoMassParticleSelector *D0MassSel = new RhoMassParticleSelector("D0", m0_D0, 1.0);

  // *** the lorentz vector of the initial interaction point
  TLorentzVector ini(0, 0, 6.231552, 7.240065);

  // ***
  // the event loop
  // ***
  while (theAnalysis->GetEvent() && i++ < nevts) {
     // *** RhoCandLists for the analysis
  RhoCandList kplus_lo, kminus_lo, piplus_lo, piminus_lo, D0_lo, D0bar_lo, psi2s_lo, kplus_ti, kminus_ti, piplus_ti, piminus_ti, D0_ti, D0bar_ti, psi2s_ti, kplus, kminus, piplus, piminus, D0, D0bar, psi2s;
    if ((i % 100) == 0)
      cout << "evt " << i << endl;

    // *** Select with no PID info ('All'); type and mass are set
    theAnalysis->FillList(kplus, "KaonAllPlus");
    theAnalysis->FillList(kminus, "KaonAllMinus");
    theAnalysis->FillList(piplus, "PionAllPlus");
    theAnalysis->FillList(piminus, "PionAllMinus");

    //kplus loop
    for(j = 0; j < kplus.GetLength(); ++j){
      RhoCandidate *kplusMcTruth = kplus[j]->GetMcTruth();
      if(kplusMcTruth){
        hmkp_->Fill(kplus[j]->PdgCode());
        RhoCandidate *kplusMother = kplusMcTruth->TheMother();
        if(kplusMother){
          hmkp->Fill(kplusMother->PdgCode());
        }
      }
    }

    //kminus loop
    for(j = 0; j < kminus.GetLength(); ++j){
      RhoCandidate *kminusMcTruth = kminus[j]->GetMcTruth();
      if(kminusMcTruth){
        hmkm_->Fill(kminus[j]->PdgCode());
        RhoCandidate *kminusMother = kminusMcTruth->TheMother();
        if(kminusMother){
          hmkm->Fill(kminusMother->PdgCode());
        }
      }
    }

    //piplus loop
    for(j = 0; j < piplus.GetLength(); ++j){
      RhoCandidate *piplusMcTruth = piplus[j]->GetMcTruth();
      if(piplusMcTruth){
        hmpp_->Fill(piplus[j]->PdgCode());
        RhoCandidate *piplusMother = piplusMcTruth->TheMother();
        if(piplusMother){
          hmpp->Fill(piplusMother->PdgCode());
        }
      }
    }

    //piminus loop
    for(j = 0; j < piminus.GetLength(); ++j){
      RhoCandidate *piminusMcTruth = piminus[j]->GetMcTruth();
      if(piminusMcTruth){
        hmpm_->Fill(piminus[j]->PdgCode());
        RhoCandidate *piminusMother = piminusMcTruth->TheMother();
        if(piminusMother){
           hmpm->Fill(piminusMother->PdgCode());
        }
      }
    }


    // *** combinatorics for D0 -> k+ pi-, and set type to D0
    D0.Combine(kplus, piminus, pdb->GetParticle("D0"));

    // *** combinatorics for D0bar -> k- pi+, and set type to D0bar
    D0bar.Combine(kminus, piplus, pdb->GetParticle("D0bar"));

    // ***
    // *** do the TRUTH MATCH for D0
    // ***

    for (j = 0; j < D0.GetLength(); ++j) {
      hD0m_all->Fill(D0[j]->M());

      if (theAnalysis->McTruthMatch(D0[j])) {
        hD0m_ftm->Fill(D0[j]->M());
        hD0m_diff->Fill(D0[j]->GetMcTruth()->M() - D0[j]->M());
      } else
        hD0m_nm->Fill(D0[j]->M());
    }

    // ***
    // *** do VERTEX FIT (D0)
    // ***
    for (j = 0; j < D0.GetLength(); ++j) {
     RhoKinVtxFitter vtxfitter(D0[j]); // instantiate a vertex fitter
      vtxfitter.Fit();

      double chi2_vtx = vtxfitter.GetChi2(); // access chi2 of fit
      double prob_vtx = vtxfitter.GetProb(); // access probability of fit
      hD0_chi2_vf->Fill(chi2_vtx);
      hD0_prob_vf->Fill(prob_vtx);

      if (prob_vtx > 0.01) // when good enough, fill some histos
      {
        RhoCandidate *jfit = D0[j]->GetFit(); // access the fitted cand
        TVector3 jVtx = jfit->Pos();            // and the decay vertex position

        hD0m_vf->Fill(jfit->M());
        hvpos->Fill(jVtx.X(), jVtx.Y());
      }
    }

    // *** some rough mass selection
    D0.Select(D0MassSel);

    // ***
    // *** do the TRUTH MATCH for D0bar
    // ***

    for (j = 0; j < D0bar.GetLength(); ++j) {
      hD0barm_all->Fill(D0bar[j]->M());

      if (theAnalysis->McTruthMatch(D0bar[j])) {
        hD0barm_ftm->Fill(D0bar[j]->M());
        hD0barm_diff->Fill(D0bar[j]->GetMcTruth()->M() - D0bar[j]->M());
      } else
        hD0barm_nm->Fill(D0bar[j]->M());
    }

    // ***
    // *** do VERTEX FIT (D0bar)
    // ***
    for (j = 0; j < D0bar.GetLength(); ++j) {
      RhoKinVtxFitter vtxfitter(D0bar[j]); // instantiate a vertex fitter
      vtxfitter.Fit();

      double chi2_vtx = vtxfitter.GetChi2(); // access chi2 of fit
      double prob_vtx = vtxfitter.GetProb(); // access probability of fit
      hD0bar_chi2_vf->Fill(chi2_vtx);
      hD0bar_prob_vf->Fill(prob_vtx);

      if (prob_vtx > 0.01) // when good enough, fill some histos
      {
        RhoCandidate *jfit = D0bar[j]->GetFit(); // access the fitted cand
        TVector3 jVtx = jfit->Pos();            // and the decay vertex position

        hD0barm_vf->Fill(jfit->M());
        hvpos->Fill(jVtx.X(), jVtx.Y());
      }
    }

    // *** some rough mass selection

    // *** combinatorics for psi(2S) -> D0 D0bar; set type to pbarpSystem (since it was generated like this)
    psi2s.Combine(D0, D0bar, pdb->GetParticle("pbarpSystem") );


    // ***
    // *** do the TRUTH MATCH for psi(2S)
    // ***
    /*
    for (j = 0; j < psi2s.GetLength(); ++j) {
      hpsim_all->Fill(psi2s[j]->M());

      if (theAnalysis->McTruthMatch(psi2s[j])) {
        hpsim_ftm->Fill(psi2s[j]->M());
        hpsim_diff->Fill(psi2s[j]->GetMcTruth()->M() - psi2s[j]->M());
      } else
        hpsim_nm->Fill(psi2s[j]->M());
    } */




    // ***
    // *** do 4C FIT (initial psi(2S) system)
    // ***
    for (j = 0; j < psi2s.GetLength(); ++j) {
      RhoKinFitter fitter(psi2s[j]); // instantiate the kin fitter in psi(2S)
      fitter.Add4MomConstraint(ini); // set 4 constraint
      fitter.Fit();                  // do fit

      double chi2_4c = fitter.GetChi2(); // get chi2 of fit
      double prob_4c = fitter.GetProb(); // access probability of fit
      hpsi_chi2_4c->Fill(chi2_4c);
      hpsi_prob_4c->Fill(prob_4c);

      if (prob_4c > 0.01) // when good enough, fill some histo
      {
        RhoCandidate *jfit = psi2s[j]->Daughter(0)->GetFit(); // get fitted J/psi

        hD0m_4cf->Fill(jfit->M());
        hD0barm_4cf->Fill(jfit->M());
      }
    }

    // ***
    // *** do DECAY TREE FIT (initial psi(2S) system)
    // ***
    for (j = 0; j < psi2s.GetLength(); ++j) {
      // *** decay tree fitter
      RhoDecayTreeFitter fittree(psi2s[j], ini);
      fittree.Fit();

      double chi2_m = fittree.GetChi2(); // get chi2 of fit
      double prob_m = fittree.GetProb(); // access probability of fit
      hpsi_chi2_dc->Fill(chi2_m);
      hpsi_prob_dc->Fill(prob_m);

      if (prob_m > 0.00001) // when good enough, fill some histo
      {
        RhoCandidate *jfit = psi2s[j]->Daughter(0)->GetFit(); // access the fitted cand
        hD0m_dc->Fill(jfit->M());
        hD0barm_dc->Fill(jfit->M());
      }
    }


    // ***
    // *** do MASS CONSTRAINT FIT (D0)
    // ***
    for (j = 0; j < D0.GetLength(); ++j) {
      RhoKinFitter mfitter(D0[j]);      // instantiate the RhoKinFitter in psi(2S)
      mfitter.AddMassConstraint(m0_D0); // add the mass constraint
      mfitter.Fit();                      // do fit

      double chi2_m = mfitter.GetChi2(); // get chi2 of fit
      double prob_m = mfitter.GetProb(); // access probability of fit
      hD0_chi2_mf->Fill(chi2_m);
      hD0_prob_mf->Fill(prob_m);

      if (prob_m > 0.01) // when good enough, fill some histo
      {
        RhoCandidate *jfit = D0[j]->GetFit(); // access the fitted cand
        hD0m_mcf->Fill(jfit->M());
        hD0barm_mcf->Fill(jfit->M());
      }
    }

    // ***
    // *** TRUE PID combinatorics
    // ***

    // *** do MC truth match for PID type
    SelectTruePid(theAnalysis, kplus);
    SelectTruePid(theAnalysis, kminus);
    SelectTruePid(theAnalysis, piplus);
    SelectTruePid(theAnalysis, piminus);

    // *** all combinatorics again with true PID
    D0.Combine(kplus, piminus);
    for (j = 0; j < D0.GetLength(); ++j)
      hD0m_trpid->Fill(D0[j]->M());
    D0.Select(D0MassSel);

    D0bar.Combine(kminus, piplus);
    for (j = 0; j < D0bar.GetLength(); ++j)
      hD0barm_trpid->Fill(D0bar[j]->M());
    D0bar.Select(D0MassSel);

    psi2s.Combine(D0, D0bar);
    for (j = 0; j < psi2s.GetLength(); ++j)
      hpsim_trpid->Fill(psi2s[j]->M());


    // ***
    // *** LOOSE PID combinatorics
    // ***

    // *** and again with PidAlgoMvd;PidAlgoStt;PidAlgoDrc and loose selection
    theAnalysis->FillList(kplus_lo, "KaonLoosePlus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(kminus_lo, "KaonLooseMinus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piplus_lo, "PionLoosePlus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piminus_lo, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
   
    //kplus_lo loop, pdg code of the partice, pdg code of the mother particle of kplus 
    for(j = 0; j < kplus_lo.GetLength(); ++j){
      RhoCandidate *kplusMcTruth = kplus_lo[j]->GetMcTruth();
      hkp_momentum->Fill(kplus_lo[j]->P());
      hkp_costheta->Fill(kplus_lo[j]->GetPosition().CosTheta());
      hkp_phi->Fill(kplus_lo[j]->GetPosition().Phi());
      if(kplusMcTruth){
         hmkp_l->Fill(kplus_lo[j]->PdgCode());
         RhoCandidate *kplusMother = kplusMcTruth->TheMother();
         if(kplusMother){
           hmkp_lo->Fill(kplusMother->PdgCode());
        }
      }
    }

    //kminus_lo loop
    for(j = 0; j < kminus_lo.GetLength(); ++j){
      RhoCandidate *kminusMcTruth = kminus_lo[j]->GetMcTruth();
      hkm_momentum->Fill(kminus_lo[j]->P());
      hkm_costheta->Fill(kminus_lo[j]->GetPosition().CosTheta());
      hkm_phi->Fill(kminus_lo[j]->GetPosition().Phi());
      if(kminusMcTruth){
        hmkm_l->Fill(kminus_lo[j]->PdgCode());
        RhoCandidate *kminusMother = kminusMcTruth->TheMother();
        if(kminusMother){
          hmkm_lo->Fill(kminusMother->PdgCode());
        }
      }
    }

    //piplus_lo loop
    for(j = 0; j < piplus_lo.GetLength(); ++j){
      RhoCandidate *piplusMcTruth = piplus_lo[j]->GetMcTruth();
      hpp_momentum->Fill(piplus_lo[j]->P());
      hpp_costheta->Fill(piplus_lo[j]->GetPosition().CosTheta());
      hpp_phi->Fill(piplus_lo[j]->GetPosition().Phi());

      if(piplusMcTruth){
        hmpp_l->Fill(piplus_lo[j]->PdgCode());
        RhoCandidate *piplusMother = piplusMcTruth->TheMother();
        if(piplusMother){
          hmpp_lo->Fill(piplusMother->PdgCode());
        }
      }
    }

    //piminus_lo loop
    for(j = 0; j < piminus_lo.GetLength(); ++j){
      RhoCandidate *piminusMcTruth = piminus_lo[j]->GetMcTruth();
      hpm_momentum->Fill(piminus_lo[j]->P());
      hpm_costheta->Fill(piminus_lo[j]->GetPosition().CosTheta());
      hpm_phi->Fill(piminus_lo[j]->GetPosition().Phi());
      if(piminusMcTruth){
        hmpm_l->Fill(piminus_lo[j]->PdgCode());
        RhoCandidate *piminusMother = piminusMcTruth->TheMother();
        if(piminusMother){
           hmpm_lo->Fill(piminusMother->PdgCode());
        }
      }
    }



    

    
  


    D0_lo.Combine(kplus_lo, piminus_lo);
    for (j = 0; j < D0_lo.GetLength(); ++j){
      hD0m_lpid->Fill(D0_lo[j]->M());
      hD0_momentum->Fill(D0_lo[j]-> P());
      hD0_costheta->Fill(D0_lo[j]->GetPosition().CosTheta());
      hD0_phi->Fill(D0_lo[j]->GetPosition().Phi());
    D0_lo.Select(D0MassSel);
    }

    D0bar_lo.Combine(kminus_lo, piplus_lo);
    for (j = 0; j < D0bar_lo.GetLength(); ++j){
      hD0barm_lpid->Fill(D0bar_lo[j]->M());
      hD0bar_momentum->Fill(D0bar_lo[j]->P());
      hD0bar_costheta->Fill(D0bar_lo[j]->GetPosition().CosTheta());
      hD0bar_phi->Fill(D0bar_lo[j]->GetPosition().Phi());
    D0bar_lo.Select(D0MassSel);
    }

    psi2s.Combine(D0_lo, D0bar_lo);
    for (j = 0; j < psi2s_lo.GetLength(); ++j){
      hpsim_lpid->Fill(psi2s_lo[j]->M());
// get daughters
      RhoCandidate *D0_tmp = psi2s_lo[j]->Daughter(0);
      RhoCandidate *D0bar_tmp = psi2s_lo[j]->Daughter(1);
      RhoCandidate *kp = D0_tmp->Daughter(0);
      RhoCandidate *pim = D0_tmp->Daughter(1);
      RhoCandidate *pip = D0bar_tmp->Daughter(1);
      RhoCandidate *km = D0bar_tmp->Daughter(0);


  const int kaon_plus_id = 321;
  const int kaon_minus_id = -321;
  const int pion_plus_id = 211;
  const int pion_minus_id = -211;
  const int D0_id = 421;
  const int D0bar_id = -421;
  RhoCandidate *kplusMcTruth = kp->GetMcTruth();
  RhoCandidate *kminusMCTruth = km->GetMcTruth();
  RhoCandidate *piplusMcTruth = pip->GetMcTruth();
  RhoCandidate *piminusMcTruth = pim->GetMcTruth();
 
 if(kplusMcTruth && kminusMCTruth && piplusMcTruth && piminusMcTruth)
  if(kaon_plus_id == kplusMcTruth->PdgCode()){
          if(kaon_minus_id == kminusMCTruth->PdgCode()){
                  if(pion_plus_id == piplusMcTruth->PdgCode()){
                          if(pion_minus_id == piminusMcTruth->PdgCode()){ 
                                  RhoCandidate *kp_mother_MC = kplusMcTruth->TheMother();               
                                  if(kp_mother_MC->PdgCode() == D0_id){
                                         RhoCandidate *pim_mother_MC = piminusMcTruth->TheMother();                                        
                                          if(pim_mother_MC->PdgCode() == D0_id){                                              
                                                  double D0_mass_truth_matched = D0_tmp->M();
                                                  double D0_mass_MC = kp_mother_MC->M();
                                                  double D0_momentum_truth_matched = D0_tmp->P();
                                                  double D0_momentum_MC = kp_mother_MC->P();
                                                  double D0_cos_truth_matched = D0_tmp->GetPosition().CosTheta();
                                                  double D0_cos_MC = kp_mother_MC->GetPosition().CosTheta();
                                                  double D0_phi_truth_matched = D0_tmp->GetPosition().Phi();
                                                  double D0_phi_MC = kp_mother_MC->GetPosition().Phi(); 

                                                  hD0tmp_mass->Fill(D0_mass_truth_matched);
                                                  hD0tmp_p->Fill(D0_momentum_truth_matched);
                                                  hD0tmp_cos->Fill(D0_cos_truth_matched);
                                                  hD0tmp_phi->Fill(D0_phi_truth_matched);

                                                  hD0_m_resolution->Fill(D0_mass_truth_matched - D0_mass_MC); 
                                                  hD0_p_resolution->Fill(D0_momentum_truth_matched - D0_momentum_MC); 
                                                  hD0_cos_resolution->Fill(D0_cos_truth_matched - D0_cos_MC);
                                                  hD0_phi_resolution->Fill(D0_phi_truth_matched - D0_phi_MC);
                                          }
                                  }
                                  RhoCandidate *km_mother_MC = kminusMCTruth->TheMother();
                                  if(km_mother_MC->PdgCode() == D0bar_id){
                                         RhoCandidate *pip_mother_MC = piplusMcTruth->TheMother();
                                          if(pip_mother_MC->PdgCode() == D0_id){
                                                  double D0bar_mass_truth_matched = D0bar_tmp->M();
                                                  double D0bar_mass_MC = km_mother_MC->M();
                                                  double D0bar_momentum_truth_matched = D0bar_tmp->P();
                                                  double D0bar_momentum_MC = km_mother_MC->P();
                                                  double D0bar_cos_truth_matched = D0bar_tmp->GetPosition().CosTheta();
                                                  double D0bar_cos_MC = km_mother_MC->GetPosition().CosTheta();
                                                  double D0bar_phi_truth_matched = D0bar_tmp->GetPosition().Phi();
                                                  double D0bar_phi_MC = km_mother_MC->GetPosition().Phi(); 

                                                  hD0bartmp_mass->Fill(D0bar_mass_truth_matched);
                                                  hD0bartmp_p->Fill(D0bar_momentum_truth_matched);
                                                  hD0bartmp_cos->Fill(D0bar_cos_truth_matched);
                                                  hD0bartmp_phi->Fill(D0bar_phi_truth_matched);

                                                  hD0bar_m_resolution->Fill(D0bar_mass_truth_matched - D0bar_mass_MC);
                                                  hD0bar_p_resolution->Fill(D0bar_momentum_truth_matched - D0bar_momentum_MC); 
                                                  hD0bar_cos_resolution->Fill(D0bar_cos_truth_matched - D0bar_cos_MC);
                                                  hD0bar_phi_resolution->Fill(D0bar_phi_truth_matched - D0bar_phi_MC);
                                          }
                                  } 
                          }
                  }
          }
  }
}








    // ***
    // *** TIGHT PID combinatorics
    // ***

    // *** and again with PidAlgoMvd;PidAlgoStt and tight selection
    theAnalysis->FillList(kplus_ti, "KaonTightPlus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(kminus_ti, "KaonTightMinus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piplus_ti, "PionLoosePlus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piminus_ti, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    
    //kplus_ti loop
    for(j = 0; j < kplus_ti.GetLength(); ++j){
      RhoCandidate *kplusMcTruth = kplus_ti[j]->GetMcTruth();
      if(kplusMcTruth){
        hmkp_t->Fill(kplus_ti[j]->PdgCode());
        RhoCandidate *kplusMother = kplusMcTruth->TheMother();
        if(kplusMother){
          hmkp_ti->Fill(kplusMother->PdgCode());
        }
      }
    }

    //kminus_ti loop
    for(j = 0; j < kminus_ti.GetLength(); ++j){
      RhoCandidate *kminusMcTruth = kminus_ti[j]->GetMcTruth();
      if(kminusMcTruth){
        hmkm_t->Fill(kminus_ti[j]->PdgCode());
        RhoCandidate *kminusMother = kminusMcTruth->TheMother();
        if(kminusMother){
          hmkm_ti->Fill(kminusMother->PdgCode());
        }
      }
    }

    //piplus_ti loop
    for(j = 0; j < piplus_ti.GetLength(); ++j){
      RhoCandidate *piplusMcTruth = piplus_ti[j]->GetMcTruth();
      if(piplusMcTruth){
        hmpp_t->Fill(piplus_ti[j]->PdgCode());
        RhoCandidate *piplusMother = piplusMcTruth->TheMother();
        if(piplusMother){
          hmpp_ti->Fill(piplusMother->PdgCode());
        }
      }
    }

    //piminus_ti loop
    for(j = 0; j < piminus_ti.GetLength(); ++j){
      RhoCandidate *piminusMcTruth = piminus_ti[j]->GetMcTruth();
      if(piminusMcTruth){
        hmpm_t->Fill(piminus_ti[j]->PdgCode());
        RhoCandidate *piminusMother = piminusMcTruth->TheMother();
        if(piminusMother){
           hmpm_ti->Fill(piminusMother->PdgCode());
        }
      }
    }

    D0_ti.Combine(kplus_ti, piminus_ti);
    for (j = 0; j < D0_ti.GetLength(); ++j)
      hD0m_tpid->Fill(D0_ti[j]->M());
    D0_ti.Select(D0MassSel);

    D0bar_ti.Combine(kminus_ti, piplus_ti);
    for (j = 0; j < D0bar_ti.GetLength(); ++j)
      hD0barm_tpid->Fill(D0bar_ti[j]->M());
    D0bar_ti.Select(D0MassSel);

    psi2s_ti.Combine(D0_ti, D0bar_ti);
    for (j = 0; j < psi2s_ti.GetLength(); ++j)
      hpsim_tpid->Fill(psi2s_ti[j]->M());
      
    // ***
    // *** P>0.4 PID combinatorics
    // ***

    // *** and again with PidAlgoMvd;PidAlgoStt and tight selection
    theAnalysis->FillList(kplus, "KaonPlus0.4", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(kminus, "KaonMinus0.4", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piplus, "PionLoosePlus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");

    D0.Combine(kplus, piminus);
    for (j = 0; j < D0.GetLength(); ++j)
      hD0m_pid04->Fill(D0[j]->M());
    D0.Select(D0MassSel);

    D0bar.Combine(kminus, piplus);
    for (j = 0; j < D0bar.GetLength(); ++j)
      hD0barm_pid04->Fill(D0bar[j]->M());
    D0bar.Select(D0MassSel);

    psi2s.Combine(D0, D0bar);
    for (j = 0; j < psi2s.GetLength(); ++j)
      hpsim_pid04->Fill(psi2s[j]->M());
  }

  // *** write out all the histos
  out->cd();

  hD0m_all->Write();
  hD0barm_all->Write();
  hpsim_all->Write();
  hD0m_lpid->Write();
  hD0barm_lpid->Write();
  hpsim_lpid->Write();
  hD0m_tpid->Write();
  hD0barm_tpid->Write();
  hpsim_tpid->Write();
  hD0m_trpid->Write();
  hD0barm_trpid->Write();
  hpsim_trpid->Write();

  hD0m_ftm->Write();
  hD0barm_ftm->Write();
  hpsim_ftm->Write();
  hD0m_nm->Write();
  hD0barm_nm->Write();
  hpsim_nm->Write();

  hpsim_diff->Write();
  hD0barm_diff->Write();
  hD0m_diff->Write();

  hD0m_vf->Write();
  hD0m_4cf->Write();
  hD0m_mcf->Write();
  hD0m_dc->Write();

  hD0barm_vf->Write();
  hD0barm_4cf->Write();
  hD0barm_mcf->Write();
  hD0barm_dc->Write();

  hD0_chi2_vf->Write();
  hD0bar_chi2_vf->Write(); 
  hpsi_chi2_4c->Write();
  hD0_chi2_mf->Write();
  hD0bar_chi2_mf->Write();
  hpsi_chi2_dc->Write();
 
  hD0_prob_vf->Write();
  hD0bar_prob_vf->Write();
  hpsi_prob_4c->Write();
  hD0_prob_mf->Write();
  hD0bar_prob_mf->Write(); 
  hpsi_prob_dc->Write();

  hvpos->Write();
  hmkp->Write();
  hmkm->Write();
  hmpp->Write();
  hmpm->Write();

  hmkp_lo->Write();
  hmkm_lo->Write();
  hmpp_lo->Write();
  hmpm_lo->Write();

  hmkp_l->Write();
  hmkm_l->Write();
  hmpp_l->Write();
  hmpm_l->Write();


  hmkp_ti->Write();
  hmkm_ti->Write();
  hmpp_ti->Write();
  hmpm_ti->Write();

  hmkp_t->Write();
  hmkm_t->Write();
  hmpp_t->Write();
  hmpm_t->Write();

  hD0_momentum->Write();
  hD0_costheta->Write();
  hD0_phi->Write();
  hD0bar_momentum->Write();
  hD0bar_costheta->Write();
  hD0bar_phi->Write();

  hkp_momentum->Write();
  hkp_costheta->Write();
  hkp_phi->Write();
  hpm_momentum->Write();
  hpm_costheta->Write();
  hpm_phi->Write();

  hpp_momentum->Write();
  hpp_costheta->Write();
  hpp_phi->Write();
  hkm_momentum->Write();
  hkm_costheta->Write();
  hkm_phi->Write();

  hD0_m_resolution->Write();
  hD0_p_resolution->Write();
  hD0_cos_resolution->Write();
  hD0_phi_resolution->Write();

  hD0bar_m_resolution->Write();
  hD0bar_p_resolution->Write();
  hD0bar_cos_resolution->Write();
  hD0bar_phi_resolution->Write();

  hkp_m_resolution->Write();
  hkp_p_resolution->Write();
  hkp_cos_resolution->Write();
  hkp_phi_resolution->Write();

  hkm_m_resolution->Write();
  hkm_p_resolution->Write();
  hkm_cos_resolution->Write();
  hkm_phi_resolution->Write();

  hpp_m_resolution->Write();
  hpp_p_resolution->Write();
  hpp_cos_resolution->Write();
  hpp_phi_resolution->Write();

  hpm_m_resolution->Write();
  hpm_p_resolution->Write();
  hpm_cos_resolution->Write();
  hpm_phi_resolution->Write();

  out->Save();

  // *** temporaty fix to avoid error on macro exit
  RemoveGeoManager();
}
