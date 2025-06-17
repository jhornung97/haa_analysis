#ifndef GUARD_HAA_H
#define GUARD_HAA_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

namespace haa {

ROOT::RDF::RNode ClosestToHiggsMassAlgo(ROOT::RDF::RNode df, 
					const std::string &str_pfcand_pt,
					const std::string &str_pfcand_eta,
					const std::string &str_pfcand_phi,
					const std::string &str_pfcand_mass,
					const std::string &str_pfcand_charge,
					const std::string &str_pfcand_mask,
					const std::string &str_daughteridx);

ROOT::RDF::RNode FourHardestPFCandsAlgo(ROOT::RDF::RNode df,
					const std::string &str_pfcand_pt,
					const std::string &str_pfcand_mask,
					const std::string &str_daughteridx);

ROOT::RDF::RNode ChargePairsAlgo(ROOT::RDF::RNode df,
				 const std::string &str_pfcand_pt,
				 const std::string &str_pfcand_eta,
				 const std::string &str_pfcand_phi,
				 const std::string &str_pfcand_mass,
				 const std::string &str_pfcand_charge,
				 const std::string &str_pfcand_mask,
				 const std::string &str_daughteridx);
  
ROOT::RDF::RNode GetHiggsP4(ROOT::RDF::RNode df, 
			    const std::string &str_d1_p4, 
			    const std::string &str_d2_p4,
			    const std::string &str_d3_p4, 
			    const std::string &str_d4_p4, 
			    const std::string &str_H_p4);
		
ROOT::RDF::RNode GetTrueDaughterP4s(ROOT::RDF::RNode df, 
				    const std::string &str_ngenpart, 
				    const std::string &str_genpart_pdgid,
				    const std::string &str_genpart_genpartidxmother, 
				    const std::string &str_genpart_pt,
				    const std::string &str_genpart_eta, 
				    const std::string &str_genpart_phi,
				    const std::string &str_genpart_mass, 
				    const std::string &str_truth_d1_p4,
				    const std::string &str_truth_d2_p4,
				    const std::string &str_truth_d3_p4, 
				    const std::string &str_truth_d4_p4,
				    const std::string &str_truth_h_p4); 

ROOT::RDF::RNode GetPseudoScalars(ROOT::RDF::RNode df, 
					const std::string &daughterIdx, 
					const std::string &str_pfcand_pt, 
					const std::string &str_pfcand_eta, 
					const std::string &str_pfcand_phi, 
                    const std::string &str_pfcand_mass, 
					const std::string &str_pfcand_charge, 
					const std::string &str_pfcand_mask,
					const std::string &str_highPtPair, 
					const std::string &str_lowPtPair);

ROOT::RDF::RNode GetMinMassDiff(ROOT::RDF::RNode df, 
					const std::string &daughterIdx, 
					const std::string &str_pfcand_pt, 
					const std::string &str_pfcand_eta, 
					const std::string &str_pfcand_phi, 
                    const std::string &str_pfcand_mass, 
					const std::string &str_pfcand_charge, 
					const std::string &str_pfcand_mask, 
					const std::string &str_ps1Pair, 
                    const std::string &str_ps2Pair);

ROOT::RDF::RNode pfCandIso(ROOT::RDF::RNode df, 
					const std::string &str_pf_cand_iso, 
					const std::string &str_pfcand_pt, 
					const std::string &str_pfcand_eta, 
					const std::string &str_pfcand_phi, 
					const std::string &str_pfcand_mass, 
					const std::string &str_higgsdaughters, 
					const std::string &str_pfcand_charged_hadron_mask, 
					const std::string &str_pfcand_neutral_hadron_mask, 
					const std::string &str_pfcand_photon_mask, 
					const std::string &str_pfcand_fromPV_mask); 

ROOT::RDF::RNode getGenPt(ROOT::RDF::RNode df, const std::string &str_genpart_mask, const std::string &str_genpart_pt, const std::string &str_gen_pt);

}

#endif
