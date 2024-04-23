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

void tut_ana(int nevts = 0, TString prefix = "signal")
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
  TH1F *hpsim_all = new TH1F("hpsim_all", "#psi(2S) mass (all);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_lpid = new TH1F("hD0m_lpid", "D0 mass (loose pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_lpid = new TH1F("hD0barm_lpid", "D0bar mass (loose pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_lpid = new TH1F("hpsim_lpid", "#psi(2S) mass (loose pid);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_tpid = new TH1F("hD0m_tpid", "D0 mass (tight pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_tpid = new TH1F("hD0barm_tpid", "D0bar mass (tight pid);m [GeV/c^{2}];entries", 200, 0, 4.5); 
  TH1F *hpsim_tpid = new TH1F("hpsim_tpid", "#psi(2S) mass (tight pid);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_pid04 = new TH1F("hD0m_pid04", "D0 mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_pid04 = new TH1F("hD0barm_pid04", "D0bar mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_pid04 = new TH1F("hpsim_pid04", "#psi(2S) mass (P_{PID} > 0.4);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_trpid = new TH1F("hD0m_trpid", "D0 mass (true pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_trpid = new TH1F("hD0barm_trpid", "D0bar mass (true pid);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_trpid = new TH1F("hpsim_trpid", "#psi(2S) mass (true pid);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_ftm = new TH1F("hD0m_ftm", "D0 mass (full truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_ftm = new TH1F("hD0barm_ftm", "D0bar mass (full truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_ftm = new TH1F("hpsim_ftm", "#psi(2S) mass (full truth match);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_nm = new TH1F("hD0m_nm", "D0 mass (no truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hD0barm_nm = new TH1F("hD0barm_nm", "D0bar mass (no truth match);m [GeV/c^{2}];entries", 200, 0, 4.5);
  TH1F *hpsim_nm = new TH1F("hpsim_nm", "#psi(2S) mass (no truth match);m [GeV/c^{2}];entries", 200, 0, 5);

  TH1F *hD0m_diff = new TH1F("hD0m_diff", "D0 mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, -2, 2);
  TH1F *hD0barm_diff = new TH1F("hD0barm_diff", "D0bar mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, -2, 2);
  TH1F *hpsim_diff = new TH1F("hpsim_diff", "#psi(2S) mass diff to truth;#Deltam [GeV/c^{2}];entries", 100, -2, 2);

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
  TH1F *hpsi_chi2_4c = new TH1F("hpsi_chi2_4c", "#psi(2S): #chi^{2} 4C fit;#chi^{2};entries", 100, 0, 250);
  TH1F *hD0_chi2_mf = new TH1F("hD0_chi2_mf", "D0: #chi^{2} mass fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hD0bar_chi2_mf = new TH1F("hD0bar_chi2_mf", "D0bar: #chi^{2} mass fit;#chi^{2};entries", 100, 0, 10);
  TH1F *hpsi_chi2_dc = new TH1F("hpsi_chi2_dc", "#psi(2S): #chi^{2} decay tree fit;#chi^{2};entries", 100, 0, 100);

  TH1F *hD0_prob_vf = new TH1F("hD0_prob_vf", "D0: Prob vertex fit;P;entries", 100, 0, 1);
  TH1F *hD0bar_prob_vf = new TH1F("hD0bar_prob_vf", "D0bar: Prob vertex fit;P;entries", 100, 0, 1);
  TH1F *hpsi_prob_4c = new TH1F("hpsi_prob_4c", "#psi(2S): Prob 4C fit;P;entries", 100, 0, 1);
  TH1F *hD0_prob_mf = new TH1F("hD0_prob_mf", "D0: Prob mass fit;P;entries", 100, 0, 1);
  TH1F *hD0bar_prob_mf = new TH1F("hD0bar_prob_mf", "D0bar: Prob mass fit;P;entries", 100, 0, 1);
  TH1F *hpsi_prob_dc = new TH1F("hpsi_prob_dc", "#psi(2S): Prob decay tree fit;P;entries", 100, 0, 1);

  TH2F *hvpos = new TH2F("hvpos", "(x,y) projection of fitted decay vertex;x [cm];y [cm]", 100, -2, 2, 100, -2, 2);

  //
  // Now the analysis stuff comes...
  //

  // *** the data reader object
  PndAnalysis *theAnalysis = new PndAnalysis();
  if (nevts == 0)
    nevts = theAnalysis->GetEntries();

  // *** RhoCandLists for the analysis
  RhoCandList kplus, kminus, piplus, piminus, D0, D0bar, psi2s;

  // *** Mass selector for the D0 cands
  double m0_D0 = pdb->GetParticle("D0")->Mass(); // Get nominal PDG mass of the D0
  RhoMassParticleSelector *D0MassSel = new RhoMassParticleSelector("D0", m0_D0, 1.0);

  // *** the lorentz vector of the initial psi(2S)
  TLorentzVector ini(0, 0, 6.231552, 7.240065);

  // ***
  // the event loop
  // ***
  while (theAnalysis->GetEvent() && i++ < nevts) {
    if ((i % 100) == 0)
      cout << "evt " << i << endl;

    // *** Select with no PID info ('All'); type and mass are set
    theAnalysis->FillList(kplus, "KaonAllPlus");
    theAnalysis->FillList(kminus, "KaonAllMinus");
    theAnalysis->FillList(piplus, "PionAllPlus");
    theAnalysis->FillList(piminus, "PionAllMinus");

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

    for (j = 0; j < psi2s.GetLength(); ++j) {
      hpsim_all->Fill(psi2s[j]->M());

      if (theAnalysis->McTruthMatch(psi2s[j])) {
        hpsim_ftm->Fill(psi2s[j]->M());
        hpsim_diff->Fill(psi2s[j]->GetMcTruth()->M() - psi2s[j]->M());
      } else
        hpsim_nm->Fill(psi2s[j]->M());
    }

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
    theAnalysis->FillList(kplus, "KaonLoosePlus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(kminus, "KaonLooseMinus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piplus, "PionLoosePlus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");

    D0.Combine(kplus, piminus);
    for (j = 0; j < D0.GetLength(); ++j)
      hD0m_lpid->Fill(D0[j]->M());
    D0.Select(D0MassSel);

    D0bar.Combine(kminus, piplus);
    for (j = 0; j < D0bar.GetLength(); ++j)
      hD0barm_lpid->Fill(D0bar[j]->M());
    D0bar.Select(D0MassSel);

    psi2s.Combine(D0, D0bar);
    for (j = 0; j < psi2s.GetLength(); ++j)
      hpsim_lpid->Fill(psi2s[j]->M());

    // ***
    // *** TIGHT PID combinatorics
    // ***

    // *** and again with PidAlgoMvd;PidAlgoStt and tight selection
    theAnalysis->FillList(kplus, "KaonTightPlus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(kminus, "KaonTightMinus", "PidAlgoMdtHardCuts;PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piplus, "PionLoosePlus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
    theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");

    D0.Combine(kplus, piminus);
    for (j = 0; j < D0.GetLength(); ++j)
      hD0m_tpid->Fill(D0[j]->M());
    D0.Select(D0MassSel);

    D0bar.Combine(kminus, piplus);
    for (j = 0; j < D0bar.GetLength(); ++j)
      hD0barm_tpid->Fill(D0bar[j]->M());
    D0bar.Select(D0MassSel);

    psi2s.Combine(D0, D0bar);
    for (j = 0; j < psi2s.GetLength(); ++j)
      hpsim_tpid->Fill(psi2s[j]->M());
      
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

  out->Save();

  // *** temporaty fix to avoid error on macro exit
  RemoveGeoManager();
}
