from code_generation.quantity import Quantity

lumi = Quantity("lumi")
puweight = Quantity("puweight")
puweight_up = Quantity("puweight_up")
puweight_down = Quantity("puweight_down")
prefireweight = Quantity("prefiring_wgt")

base_taus_mask = Quantity("base_taus_mask")
good_taus_mask = Quantity("good_taus_mask")
base_muons_mask = Quantity("base_muons_mask")
good_muons_mask = Quantity("good_muons_mask")
veto_muons_mask = Quantity("veto_muons_mask")
veto_muons_mask_2 = Quantity("veto_muons_mask_2")
muon_veto_flag = Quantity("extramuon_veto")
base_electrons_mask = Quantity("base_electrons_mask")
good_electrons_mask = Quantity("good_electrons_mask")
veto_electrons_mask = Quantity("veto_electrons_mask")
veto_electrons_mask_2 = Quantity("veto_electrons_mask_2")
electron_veto_flag = Quantity("extraelec_veto")
base_photons_mask = Quantity("base_photons_mask")
base_pfcands_mask = Quantity("base_pfcands_mask")
jet_id_mask = Quantity("jet_id_mask")
jet_puid_mask = Quantity("jet_puid_mask")
jet_overlap_veto_mask = Quantity("jet_overlap_veto_mask")
jet_electron_overlap_veto_mask = Quantity("jet_electron_overlap_veto_mask")
jet_muon_overlap_veto_mask = Quantity("jet_muon_overlap_veto_mask")
jet_photon_overlap_veto_mask = Quantity("jet_photon_overlap_veto_mask")
pfcand_electron_overlap_veto_mask = Quantity("pfcand_electron_overlap_veto_mask")
pfcand_muon_overlap_veto_mask = Quantity("pfcand_muon_overlap_veto_mask")
pfcand_photon_overlap_veto_mask = Quantity("pfcand_photon_overlap_veto_mask")
good_jets_mask = Quantity("good_jets_mask")
good_fatjets_mask = Quantity("good_fatjets_mask")
good_bjets_mask = Quantity("good_bjets_mask")
Tau_pt_ele_corrected = Quantity("Tau_pt_ele_corrected")
Tau_pt_ele_mu_corrected = Quantity("Tau_pt_mu_corrected")
Tau_pt_corrected = Quantity("Tau_pt_corrected")
Tau_mass_corrected = Quantity("Tau_mass_corrected")
Jet_pt_corrected = Quantity("Jet_pt_corrected")
Jet_mass_corrected = Quantity("Jet_mass_corrected")
dileptonpair = Quantity("dileptonpair")
dielectronpair = Quantity("dielectronpair")
gen_dileptonpair = Quantity("gen_dileptonpair")
truegenpair = Quantity("truegenpair")
good_jet_collection = Quantity("good_jet_collection")
good_fatjet_collection = Quantity("good_fatjet_collection")
good_bjet_collection = Quantity("good_bjet_collection")

base_zs_mask = Quantity("base_zs_mask")
z_p4 = Quantity("z_p4")
z_pt = Quantity("z_pt")

nelectrons = Quantity("nelectrons")
nmuons = Quantity("nmuons")
ntaus = Quantity("ntaus")
p4_1 = Quantity("p4_1")
p4_1_uncorrected = Quantity("p4_1_uncorrected")
pt_1 = Quantity("pt_1")
eta_1 = Quantity("eta_1")
phi_1 = Quantity("phi_1")
p4_2 = Quantity("p4_2")
p4_2_uncorrected = Quantity("p4_2_uncorrected")
pt_2 = Quantity("pt_2")
eta_2 = Quantity("eta_2")
phi_2 = Quantity("phi_2")
mass_1 = Quantity("mass_1")
mass_2 = Quantity("mass_2")
dxy_1 = Quantity("dxy_1")
dxy_2 = Quantity("dxy_2")
dz_1 = Quantity("dz_1")
dz_2 = Quantity("dz_2")
q_1 = Quantity("q_1")
q_2 = Quantity("q_2")
iso_1 = Quantity("iso_1")
iso_2 = Quantity("iso_2")
decaymode_1 = Quantity("decaymode_1")
decaymode_2 = Quantity("decaymode_2")
gen_match_1 = Quantity("gen_match_1")
gen_match_2 = Quantity("gen_match_2")
gen_tau_pt_1 = Quantity("gen_tau_pt_1")
gen_tau_pt_2 = Quantity("gen_tau_pt_2")
gen_tau_eta_1 = Quantity("gen_tau_eta_1")
gen_tau_eta_2 = Quantity("gen_tau_eta_2")
gen_tau_phi_1 = Quantity("gen_tau_phi_1")
gen_tau_phi_2 = Quantity("gen_tau_phi_2")
taujet_pt_1 = Quantity("taujet_pt_1")
taujet_pt_2 = Quantity("taujet_pt_2")
gen_taujet_pt_1 = Quantity("gen_taujet_pt_1")
gen_taujet_pt_2 = Quantity("gen_taujet_pt_2")


# Combined event quantities
m_vis = Quantity("m_vis")
pt_vis = Quantity("pt_vis")
eta_vis = Quantity("eta_vis")
phi_vis = Quantity("phi_vis")

pzetamissvis = Quantity("pzetamissvis")
mTdilepton = Quantity("mTdileptonMET")
mt_1 = Quantity("mt_1")
mt_2 = Quantity("mt_2")
pt_tt = Quantity("pt_tt")
pt_ttjj = Quantity("pt_ttjj")
mt_tot = Quantity("mt_tot")

good_pfcands_mask = Quantity("good_pfcands_mask")
npfcands = Quantity("npfcands")

higgsdaughters = Quantity("higgsdaughters")
ChargedPFCands = Quantity("ChargedPFCands")
NeutralPFCands = Quantity("NeutralPFCands")
PhotonPFCands = Quantity("PhotonPFCands")
fromPV = Quantity("fromPV")

ps1Pair = Quantity("ps1Pair")
ps2Pair = Quantity("ps2Pair")
ps_1_d_1_p4 = Quantity("ps_1_d_1_p4")
ps_1_d_2_p4 = Quantity("ps_1_d_2_p4")
ps_2_d_1_p4 = Quantity("ps_2_d_1_p4")
ps_2_d_2_p4 = Quantity("ps_2_d_2_p4")

ps_1_mass = Quantity("ps_1_mass")
ps_2_mass = Quantity("ps_2_mass")

pfcands_iso = Quantity("pfcands_iso")
pfcands_pdgid = Quantity("pfcands_pdgid")
pfcands_pt = Quantity("pfcands_pt")
pfcands_eta = Quantity("pfcands_eta")
pfcands_phi = Quantity("pfcands_phi")
pfcands_mass = Quantity("pfcands_mass")

daughter_isos = Quantity("daughter_isos")

d1_p4 = Quantity("d1_p4")
d1_pt = Quantity("d1_pt")
d1_eta = Quantity("d1_eta")
d1_phi = Quantity("d1_phi")
d1_mass = Quantity("d1_mass")
d1_prompt = Quantity("d1_prompt")
d1_iso = Quantity("d1_iso")
d2_p4 = Quantity("d2_p4")
d2_pt = Quantity("d2_pt")
d2_eta = Quantity("d2_eta")
d2_phi = Quantity("d2_phi")
d2_mass = Quantity("d2_mass")
d2_prompt = Quantity("d2_prompt")
d2_iso = Quantity("d2_iso")
d3_p4 = Quantity("d3_p4")
d3_pt = Quantity("d3_pt")
d3_eta = Quantity("d3_eta")
d3_phi = Quantity("d3_phi")
d3_mass = Quantity("d3_mass")
d3_prompt = Quantity("d3_prompt")
d3_iso = Quantity("d3_iso")
d4_p4 = Quantity("d4_p4")
d4_pt = Quantity("d4_pt")
d4_eta = Quantity("d4_eta")
d4_phi = Quantity("d4_phi")
d4_mass = Quantity("d4_mass")
d4_prompt = Quantity("d4_prompt")
d4_iso = Quantity("d4_iso")

truth_d1_p4 = Quantity("truth_d1_p4")
truth_d1_pt = Quantity("truth_d1_pt")
truth_d1_eta = Quantity("truth_d1_eta")
truth_d1_phi = Quantity("truth_d1_phi")

truth_d2_p4 = Quantity("truth_d2_p4")
truth_d2_pt = Quantity("truth_d2_pt")
truth_d2_eta = Quantity("truth_d2_eta")
truth_d2_phi = Quantity("truth_d2_phi")

truth_d3_p4 = Quantity("truth_d3_p4")
truth_d3_pt = Quantity("truth_d3_pt")
truth_d3_eta = Quantity("truth_d3_eta")
truth_d3_phi = Quantity("truth_d3_phi")

truth_d4_p4 = Quantity("truth_d4_p4")
truth_d4_pt = Quantity("truth_d4_pt")
truth_d4_eta = Quantity("truth_d4_eta")
truth_d4_phi = Quantity("truth_d4_phi")

truth_H_p4 = Quantity("truth_H_p4")
truth_H_pt = Quantity("truth_H_pt")
truth_H_eta = Quantity("truth_H_eta")
truth_H_phi = Quantity("truth_H_phi")
truth_H_mass = Quantity("truth_H_mass")

H_p4 = Quantity("H_p4")
H_pt = Quantity("H_pt")
H_eta = Quantity("H_eta")
H_phi = Quantity("H_phi")
H_mass = Quantity("H_mass")

ps_mass_12 = Quantity("ps_mass_12")
ps_mass_14 = Quantity("ps_mass_14")
ps_mass_23 = Quantity("ps_mass_23")
ps_mass_34 = Quantity("ps_mass_34")

njets = Quantity("njets")
nbtag = Quantity("nbtag")
jet_p4_1 = Quantity("jet_p4_1")
jpt_1 = Quantity("jpt_1")
jeta_1 = Quantity("jeta_1")
jphi_1 = Quantity("jphi_1")
jmass_1 = Quantity("jmass_1")
jbtagDeepB_1 = Quantity("jbtagDeepB_1")
jbtagDeepFlavQG_1 = Quantity("jbtagDeepFlavQG_1")
jchEmEF_1 = Quantity("jchEmEF_1")
jchHEF_1 = Quantity("jchHEF_1")
jneEmEF_1 = Quantity("jneEmEF_1")
jneHEF_1 = Quantity("jneHEF_1")
jet_p4_2 = Quantity("jet_p4_2")
jpt_2 = Quantity("jpt_2")
jeta_2 = Quantity("jeta_2")
jphi_2 = Quantity("jphi_2")
jmass_2 = Quantity("jmass_2")
jbtagDeepB_2 = Quantity("jbtagDeepB_2")
jbtagDeepFlavQG_2 = Quantity("jbtagDeepFlavQG_2")
jchEmEF_2 = Quantity("jchEmEF_2")
jchHEF_2 = Quantity("jchHEF_2")
jneEmEF_2 = Quantity("jneEmEF_2")
jneHEF_2 = Quantity("jneHEF_2")
jet_p4_3 = Quantity("jet_p4_3")
jpt_3 = Quantity("jpt_3")
jeta_3 = Quantity("jeta_3")
jphi_3 = Quantity("jphi_3")
jmass_3 = Quantity("jmass_3")
jbtagDeepB_3 = Quantity("jbtagDeepB_3")
jbtagDeepFlavQG_3 = Quantity("jbtagDeepFlavQG_3")
jchEmEF_3 = Quantity("jchEmEF_3")
jchHEF_3 = Quantity("jchHEF_3")
jneEmEF_3 = Quantity("jneEmEF_3")
jneHEF_3 = Quantity("jneHEF_3")
jet_p4_4 = Quantity("jet_p4_4")
jpt_4 = Quantity("jpt_4")
jeta_4 = Quantity("jeta_4")
jphi_4 = Quantity("jphi_4")
jmass_4 = Quantity("jmass_4")
jbtagDeepB_4 = Quantity("jbtagDeepB_4")
jbtagDeepFlavQG_4 = Quantity("jbtagDeepFlavQG_4")
jchEmEF_4 = Quantity("jchEmEF_4")
jchHEF_4 = Quantity("jchHEF_4")
jneEmEF_4 = Quantity("jneEmEF_4")
jneHEF_4 = Quantity("jneHEF_4")
jet_p4_5 = Quantity("jet_p4_5")
jpt_5 = Quantity("jpt_5")
jeta_5 = Quantity("jeta_5")
jphi_5 = Quantity("jphi_5")
jmass_5 = Quantity("jmass_5")
jbtagDeepB_5 = Quantity("jbtagDeepB_5")
jbtagDeepFlavQG_5 = Quantity("jbtagDeepFlavQG_5")
jchEmEF_5 = Quantity("jchEmEF_5")
jchHEF_5 = Quantity("jchHEF_5")
jneEmEF_5 = Quantity("jneEmEF_5")
jneHEF_5 = Quantity("jneHEF_5")

nfatjets = Quantity("nfatjets")
fatjet_p4_1 = Quantity("fatjet_p4_1")
fjpt_1 = Quantity("fjpt_1")
fjeta_1 = Quantity("fjeta_1")
fjphi_1 = Quantity("fjphi_1")
fjmass_1 = Quantity("fjmass_1")
fjmsoftdrop_1 = Quantity("fjmsoftdrop_1")
fjn2b1_1 = Quantity("fjn2b1_1")
fjn3b1_1 = Quantity("fjn3b1_1")
fjparticleNet_H4qvsQCD_1 = Quantity("fjparticleNet_H4qvsQCD_1")
fjparticleNet_QCD_1 = Quantity("fjparticleNet_QCD_1")
fjparticleNet_mass_1 = Quantity("fjparticleNet_mass_1")
fatjet_p4_2 = Quantity("fatjet_p4_2")
fjpt_2 = Quantity("fjpt_2")
fjeta_2 = Quantity("fjeta_2")
fjphi_2 = Quantity("fjphi_2")
fjmass_2 = Quantity("fjmass_2")
fjmsoftdrop_2 = Quantity("fjmsoftdrop_2")
fjn2b1_2 = Quantity("fjn2b1_2")
fjn3b1_2 = Quantity("fjn3b1_2")
fjparticleNet_H4qvsQCD_2 = Quantity("fjparticleNet_H4qvsQCD_2")
fjparticleNet_QCD_2 = Quantity("fjparticleNet_QCD_2")
fjparticleNet_mass_2 = Quantity("fjparticleNet_mass_2")
fatjet_p4_3 = Quantity("fatjet_p4_3")
fjpt_3 = Quantity("fjpt_3")
fjeta_3 = Quantity("fjeta_3")
fjphi_3 = Quantity("fjphi_3")
fjmass_3 = Quantity("fjmass_3")
fjmsoftdrop_3 = Quantity("fjmsoftdrop_3")
fjn2b1_3 = Quantity("fjn2b1_3")
fjn3b1_3 = Quantity("fjn3b1_3")
fjparticleNet_H4qvsQCD_3 = Quantity("fjparticleNet_H4qvsQCD_3")
fjparticleNet_QCD_3 = Quantity("fjparticleNet_QCD_3")
fjparticleNet_mass_3 = Quantity("fjparticleNet_mass_3")

mjj = Quantity("mjj")
bjet_p4_1 = Quantity("bjet_p4_1")
bpt_1 = Quantity("bpt_1")
beta_1 = Quantity("beta_1")
bphi_1 = Quantity("bphi_1")
btag_value_1 = Quantity("btag_value_1")
bjet_p4_2 = Quantity("bjet_p4_2")
bpt_2 = Quantity("bpt_2")
beta_2 = Quantity("beta_2")
bphi_2 = Quantity("bphi_2")
btag_value_2 = Quantity("btag_value_2")

dielectron_veto = Quantity("dielectron_veto")
dimuon_veto = Quantity("dimuon_veto")
dilepton_veto = Quantity("dilepton_veto")

## Gen Quantities
gen_p4_1 = Quantity("gen_p4_1")
gen_pt_1 = Quantity("gen_pt_1")
gen_eta_1 = Quantity("gen_eta_1")
gen_phi_1 = Quantity("gen_phi_1")
gen_mass_1 = Quantity("gen_mass_1")
gen_pdgid_1 = Quantity("gen_pdgid_1")

gen_p4_2 = Quantity("gen_p4_2")
gen_pt_2 = Quantity("gen_pt_2")
gen_eta_2 = Quantity("gen_eta_2")
gen_phi_2 = Quantity("gen_phi_2")
gen_mass_2 = Quantity("gen_mass_2")
gen_pdgid_2 = Quantity("gen_pdgid_2")
gen_m_vis = Quantity("gen_m_vis")
gen_pt_vis = Quantity("gen_pt_vis")
gen_eta_vis = Quantity("gen_eta_vis")

genpart_pdgid = Quantity("genpart_pdgid")
genpart_status = Quantity("genpart_status")
genpart_eta = Quantity("genpart_eta")
genpart_phi = Quantity("genpart_phi")

topPtReweightWeight = Quantity("topPtReweightWeight")
ZPtMassReweightWeight = Quantity("ZPtMassReweightWeight")

## HTXS quantities
ggh_NNLO_weight = Quantity("ggh_NNLO_weight")
THU_ggH_Mu = Quantity("THU_ggH_Mu")
THU_ggH_Res = Quantity("THU_ggH_Res")
THU_ggH_Mig01 = Quantity("THU_ggH_Mig01")
THU_ggH_Mig12 = Quantity("THU_ggH_Mig12")
THU_ggH_VBF2j = Quantity("THU_ggH_VBF2j")
THU_ggH_VBF3j = Quantity("THU_ggH_VBF3j")
THU_ggH_PT60 = Quantity("THU_ggH_PT60")
THU_ggH_PT120 = Quantity("THU_ggH_PT120")
THU_ggH_qmtop = Quantity("THU_ggH_qmtop")
THU_qqH_TOT = Quantity("THU_qqH_TOT")
THU_qqH_PTH200 = Quantity("THU_qqH_PTH200")
THU_qqH_Mjj60 = Quantity("THU_qqH_Mjj60")
THU_qqH_Mjj120 = Quantity("THU_qqH_Mjj120")
THU_qqH_Mjj350 = Quantity("THU_qqH_Mjj350")
THU_qqH_Mjj700 = Quantity("THU_qqH_Mjj700")
THU_qqH_Mjj1000 = Quantity("THU_qqH_Mjj1000")
THU_qqH_Mjj1500 = Quantity("THU_qqH_Mjj1500")
THU_qqH_25 = Quantity("THU_qqH_25")
THU_qqH_JET01 = Quantity("THU_qqH_JET01")

## MET quantities
met_p4 = Quantity("met_p4")
recoil_genboson_p4_vec = Quantity("recoil_genboson_p4_vec")
genbosonmass = Quantity("genbosonmass")
npartons = Quantity("npartons")
met_p4_leptoncorrected = Quantity("met_p4_leptoncorrected")
met_p4_jetcorrected = Quantity("met_p4_jetcorrected")
met_p4_recoilcorrected = Quantity("met_p4_recoilcorrected")
met = Quantity("met")
metphi = Quantity("metphi")
metSumEt = Quantity("metSumEt")
metcov00 = Quantity("metcov00")
metcov01 = Quantity("metcov01")
metcov10 = Quantity("metcov10")
metcov11 = Quantity("metcov11")
met_uncorrected = Quantity("met_uncorrected")
metphi_uncorrected = Quantity("metphi_uncorrected")

## PFMET quantities
pfmet = Quantity("pfmet")
pfmet_p4 = Quantity("pfmet_p4")
pfmetphi = Quantity("pfmetphi")
pfmet_uncorrected = Quantity("pfmet_uncorrected")
pfmetphi_uncorrected = Quantity("pfmetphi_uncorrected")
pfmet_p4_leptoncorrected = Quantity("pfmet_p4_leptoncorrected")
pfmet_p4_jetcorrected = Quantity("pfmet_p4_jetcorrected")
pfmet_p4_recoilcorrected = Quantity("pfmet_p4_recoilcorrected")

## embedding quantities
emb_genweight = Quantity("emb_genweight")
emb_initialMETEt = Quantity("emb_initialMETEt")
emb_initialMETphi = Quantity("emb_initialMETphi")
emb_initialPuppiMETEt = Quantity("emb_initialPuppiMETEt")
emb_initialPuppiMETphi = Quantity("emb_initialPuppiMETphi")
emb_isMediumLeadingMuon = Quantity("emb_isMediumLeadingMuon")
emb_isMediumTrailingMuon = Quantity("emb_isMediumTrailingMuon")
emb_isTightLeadingMuon = Quantity("emb_isTightLeadingMuon")
emb_isTightTrailingMuon = Quantity("emb_isTightTrailingMuon")
emb_InitialPairCandidates = Quantity("emb_InitialPairCandidates")
emb_SelectionOldMass = Quantity("emb_SelectionOldMass")
emb_SelectionNewMass = Quantity("emb_SelectionNewMass")

# safe nano PFCand vars
nano_pfcands_eta = Quantity("nano_pfcands_eta")
nano_pfcands_mass = Quantity("nano_pfcands_mass")
nano_pfcands_phi = Quantity("nano_pfcands_phi")
nano_pfcands_pt = Quantity("nano_pfcands_pt")
nano_pfcands_charge = Quantity("nano_pfcands_charge")
nano_pfcands_pdgId = Quantity("nano_pfcands_pdgId")
nano_pfcands_d0 = Quantity("nano_pfcands_d0")
nano_pfcands_d0Err = Quantity("nano_pfcands_d0Err")
nano_pfcands_dz = Quantity("nano_pfcands_dz")
nano_pfcands_dzErr = Quantity("nano_pfcands_dzErr")
nano_pfcands_trkQuality = Quantity("nano_pfcands_trkQuality")

# unsafe nano PFCand vars

nano_pfcands_PuppuWeight = Quantity("nano_pfcands_PuppuWeight")
nano_pfcands_PuppuWeightNoLep = Quantity("nano_pfcands_PuppuWeightNoLep")
nano_pfcands_trkChi2 = Quantity("nano_pfcands_trkChi2")
nano_pfcands_trkEta = Quantity("nano_pfcands_trkEta")
nano_pfcands_trkPt = Quantity("nano_pfcands_trkPt")
nano_pfcands_vtxChi2 = Quantity("nano_pfcands_vtxChi2")

nano_pfcands_lostInnerHits = Quantity("nano_pfcands_lostInnerHits")

nano_pfcands_pvAssocQuality = Quantity("nano_pfcands_pvAssocQuality")

# safe nano Jet vars
nano_jet_eta = Quantity("nano_jet_eta")
nano_jet_mass = Quantity("nano_jet_mass")
nano_jet_phi = Quantity("nano_jet_phi")
nano_jet_pt = Quantity("nano_jet_pt")
nano_jet_btagDeepB = Quantity("nano_jet_btagDeepB")
nano_jet_btagDeepFlavQG = Quantity("nano_jet_btagDeepFlavQG")
nano_jet_chEmEF = Quantity("nano_jet_chEmEF")
nano_jet_chHEF = Quantity("nano_jet_chHEF")
nano_jet_neEmEF = Quantity("nano_jet_neEmEF")
nano_jet_neHEF = Quantity("nano_jet_neHEF")
nano_jet_jetId = Quantity("nano_jet_jetId")

# Fat Jet vars
nano_fatjet_eta = Quantity("nano_fatjet_eta")
nano_fatjet_mass = Quantity("nano_fatjet_mass")
nano_fatjet_phi = Quantity("nano_fatjet_phi")
nano_fatjet_pt = Quantity("nano_fatjet_pt")
nano_fatjet_msoftdrop = Quantity("nano_fatjet_msoftdrop")
nano_fatjet_n2b1 = Quantity("nano_fatjet_n2b1")
nano_fatjet_n3b1 = Quantity("nano_fatjet_n3b1")
nano_fatjet_particleNet_H4qvsQCD = Quantity("nano_fatjet_particleNet_H4qvsQCD")
nano_fatjet_particleNet_QCD = Quantity("nano_fatjet_particleNet_QCD")
nano_fatjet_particleNet_mass = Quantity("nano_fatjet_particleNet_mass")

# sample flags
is_data = Quantity("is_data")
is_embedding = Quantity("is_embedding")
is_ttbar = Quantity("is_ttbar")
is_dyjets = Quantity("is_dyjets")
is_wjets = Quantity("is_wjets")
is_ggh_htautau = Quantity("is_ggh_htautau")
is_vbf_htautau = Quantity("is_vbf_htautau")
is_diboson = Quantity("is_diboson")

# Electron Weights
id_wgt_ele_wp90nonIso_1 = Quantity("id_wgt_ele_wp90nonIso_1")
id_wgt_ele_wp90nonIso_2 = Quantity("id_wgt_ele_wp90nonIso_2")
id_wgt_ele_wp80nonIso_1 = Quantity("id_wgt_ele_wp80nonIso_1")
id_wgt_ele_wp80nonIso_2 = Quantity("id_wgt_ele_wp80nonIso_2")
id_wgt_ele_1 = Quantity("id_wgt_ele_1")
id_wgt_ele_2 = Quantity("id_wgt_ele_2")
trigger_wgt_ele_1 = Quantity("trigger_wgt_ele_1")
trigger_wgt_ele_2 = Quantity("trigger_wgt_ele_2")
trigger_wgt_ele_1_up = Quantity("trigger_wgt_ele_1_up")
trigger_wgt_ele_2_up = Quantity("trigger_wgt_ele_2_up")
trigger_wgt_ele_1_down = Quantity("trigger_wgt_ele_1_down")
trigger_wgt_ele_2_down = Quantity("trigger_wgt_ele_2_down")
# Muon weights
id_wgt_mu_1 = Quantity("id_wgt_mu_1")
id_wgt_mu_2 = Quantity("id_wgt_mu_2")
iso_wgt_mu_1 = Quantity("iso_wgt_mu_1")
iso_wgt_mu_2 = Quantity("iso_wgt_mu_2")
trigger_wgt_mu_1 = Quantity("trigger_wgt_mu_1")
trigger_wgt_mu_2 = Quantity("trigger_wgt_mu_2")
