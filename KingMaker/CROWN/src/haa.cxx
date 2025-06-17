#ifndef GUARD_HAA_H
#define GUARD_HAA_H

#include "../include/haa.hxx"
#include "../include/RoccoR.hxx"
#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TMinuit.h"
#include "TVector2.h"
#include "correction.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <vector>
#include <iostream>
#include <algorithm>

namespace haa {

ROOT::RDF::RNode GetTrueDaughterP4s(ROOT::RDF::RNode df, const std::string &str_ngenpart, const std::string &str_genpart_pdgid, 
				    const std::string &str_genpart_genpartidxmother, const std::string &str_genpart_pt,
				    const std::string &str_genpart_eta, const std::string &str_genpart_phi,  
				    const std::string &str_genpart_mass, const std::string &str_truth_d1_p4,
				    const std::string &str_truth_d2_p4, const std::string &str_truth_d3_p4, const std::string &str_truth_d4_p4, const std::string &str_truth_h_p4) {
  auto truth_didx = [](const unsigned int ngenpart, const ROOT::RVec<int> genpart_pdgid, const ROOT::RVec<int> genpart_genpartidxmother, const ROOT::RVec<float> genpart_pt) {
    
    ROOT::RVec<int> truedaughteridxs;
    std::vector<int> truedaughteridxs_pos;
    std::vector<int> truedaughteridxs_neg;
    
    int ndaughters = 0;

    for (int i = 0; i < ngenpart; i++) {
      if (std::abs(genpart_pdgid[genpart_genpartidxmother[i]]) == 9000006 && genpart_pdgid[i] > 0){
      	truedaughteridxs_pos.push_back(i);
        ndaughters++;
      }
      if (std::abs(genpart_pdgid[genpart_genpartidxmother[i]]) == 9000006 && genpart_pdgid[i] < 0) {
      	truedaughteridxs_neg.push_back(i);
        ndaughters++; 
      }
    }

    std::sort(truedaughteridxs_pos.begin(),truedaughteridxs_pos.end(), [&genpart_pt](int a, int b) {
	    return genpart_pt[a] > genpart_pt[b];
      }
      );

    std::sort(truedaughteridxs_neg.begin(),truedaughteridxs_neg.end(), [&genpart_pt](int a, int b) {
      return genpart_pt[a] > genpart_pt[b];
      }
      );

    truedaughteridxs = {truedaughteridxs_pos[0], truedaughteridxs_neg[0], truedaughteridxs_pos[1], truedaughteridxs_neg[1]};

    return truedaughteridxs;
  };

  auto truth_d1_p4 = [](const ROOT::RVec<int> truedaughteridxs, const ROOT::RVec<float> genpart_pt, const ROOT::RVec<float> genpart_eta, 
			const ROOT::RVec<float> genpart_phi, const ROOT::RVec<float> genpart_mass) {
    return ROOT::Math::PtEtaPhiMVector(genpart_pt[truedaughteridxs[0]], genpart_eta[truedaughteridxs[0]], genpart_phi[truedaughteridxs[0]], genpart_mass[truedaughteridxs[0]]); 
  };
  
  auto truth_d2_p4 = [](const ROOT::RVec<int> truedaughteridxs, const ROOT::RVec<float> genpart_pt, const ROOT::RVec<float> genpart_eta,
			const ROOT::RVec<float> genpart_phi, const ROOT::RVec<float> genpart_mass) {
    return ROOT::Math::PtEtaPhiMVector(genpart_pt[truedaughteridxs[1]], genpart_eta[truedaughteridxs[1]], genpart_phi[truedaughteridxs[1]], genpart_mass[truedaughteridxs[1]]);
  };
  
  auto truth_d3_p4 = [](const ROOT::RVec<int> truedaughteridxs, const ROOT::RVec<float> genpart_pt, const ROOT::RVec<float> genpart_eta,
			const ROOT::RVec<float> genpart_phi, const ROOT::RVec<float> genpart_mass) {
    return ROOT::Math::PtEtaPhiMVector(genpart_pt[truedaughteridxs[2]], genpart_eta[truedaughteridxs[2]], genpart_phi[truedaughteridxs[2]], genpart_mass[truedaughteridxs[2]]);
  };
  
  auto truth_d4_p4 = [](const ROOT::RVec<int> truedaughteridxs, const ROOT::RVec<float> genpart_pt, const ROOT::RVec<float> genpart_eta,
			const ROOT::RVec<float> genpart_phi, const ROOT::RVec<float> genpart_mass) {
    return ROOT::Math::PtEtaPhiMVector(genpart_pt[truedaughteridxs[3]], genpart_eta[truedaughteridxs[3]], genpart_phi[truedaughteridxs[3]], genpart_mass[truedaughteridxs[3]]);
  };

  auto truth_h_p4 = [](const unsigned int ngenpart, const ROOT::RVec<int> genpart_pdgid, const ROOT::RVec<int> genpart_genpartidxmother, 
		       const ROOT::RVec<float> genpart_pt, const ROOT::RVec<float> genpart_eta, const ROOT::RVec<float> genpart_phi, 
		       const ROOT::RVec<float> genpart_mass) {
    int thidx = 0;
    for (int i = 0; i < ngenpart; i++) {
      if (genpart_pdgid[genpart_genpartidxmother[i]] == 25) {
	      thidx = genpart_genpartidxmother[i];
      }
    }
    return ROOT::Math::PtEtaPhiMVector(genpart_pt[thidx], genpart_eta[thidx], genpart_phi[thidx], genpart_mass[thidx]);
  };

  auto df1 = df.Define("truedaughteridxs", truth_didx, {str_ngenpart, str_genpart_pdgid, str_genpart_genpartidxmother, str_genpart_pt});
  auto df2 = df1.Define(str_truth_d1_p4, truth_d1_p4, {"truedaughteridxs", str_genpart_pt, str_genpart_eta, str_genpart_phi, str_genpart_mass});
  auto df3 = df2.Define(str_truth_d2_p4, truth_d2_p4, {"truedaughteridxs", str_genpart_pt, str_genpart_eta, str_genpart_phi, str_genpart_mass});
  auto df4 = df3.Define(str_truth_d3_p4, truth_d3_p4, {"truedaughteridxs", str_genpart_pt, str_genpart_eta, str_genpart_phi, str_genpart_mass});
  auto df5 = df4.Define(str_truth_d4_p4, truth_d4_p4, {"truedaughteridxs", str_genpart_pt, str_genpart_eta, str_genpart_phi, str_genpart_mass});
  auto df6 = df5.Define(str_truth_h_p4, truth_h_p4, {str_ngenpart, str_genpart_pdgid, str_genpart_genpartidxmother, str_genpart_pt, str_genpart_eta, str_genpart_phi, str_genpart_mass});
  
  return df6;
  
}

ROOT::RDF::RNode ClosestToHiggsMassAlgo(ROOT::RDF::RNode df, const std::string &str_pfcand_pt, const std::string &str_pfcand_eta, const std::string &str_pfcand_phi, 
					const std::string &str_pfcand_mass, const std::string &str_pfcand_charge, const std::string &str_pfcand_mask, const std::string &str_daughteridxs) {
  Logger::get("HiggsSelection")->debug("Setting up algorithm");
  auto hidx = [](const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<float> pfcand_eta, const ROOT::RVec<float> pfcand_phi, const ROOT::RVec<float> pfcand_mass, const ROOT::RVec<int> pfcand_charge,
		 const ROOT::RVec<int> pfcand_mask) {

    Logger::get("HiggsSelectionAlgo")
    ->debug("Running algorithm on all hadrons (ID = 211 or ID == 1)");
    ROOT::RVec<int> selected_hadrons = {-1, -1, -1, -1};
    auto original_pfcand_indices = ROOT::VecOps::Nonzero(pfcand_mask);
    
    const auto good_pts = ROOT::VecOps::Take(pfcand_pt, original_pfcand_indices);
    const auto good_etas = ROOT::VecOps::Take(pfcand_eta, original_pfcand_indices);
    const auto good_phis = ROOT::VecOps::Take(pfcand_phi, original_pfcand_indices);
    const auto good_masses = ROOT::VecOps::Take(pfcand_mass, original_pfcand_indices);
    const auto good_charges = ROOT::VecOps::Take(pfcand_charge, original_pfcand_indices);

    std::sort(original_pfcand_indices.begin(), original_pfcand_indices.end(), [&pfcand_pt](int a, int b) {
	return pfcand_pt[a] > pfcand_pt[b];
      }
      );

    if (original_pfcand_indices.size() < 4) {
      return selected_hadrons; 
    }
    
    auto fourVecs = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
	good_pts, good_etas, good_phis, good_masses);

    float hmass_difference = -1.0;
                                                                                                     //float psmass_difference = 0.5;
    float hmass_candidate = -1.0;
    auto selected_hadron_indices = std::vector<int>{-1, -1, -1, -1};

    auto combinations = ROOT::VecOps::Combinations(original_pfcand_indices, 4);
    
    for (unsigned long n = 0; n < combinations[0].size(); n++) {
      auto hadron_1 = fourVecs[combinations[0][n]]; 
      auto hadron_2 = fourVecs[combinations[1][n]];
      auto hadron_3 = fourVecs[combinations[2][n]];
      auto hadron_4 = fourVecs[combinations[3][n]];
      int count = 0;
      for (int i = 0; i < 4; i++) {
	for (int j = i+1; j < 4; j++) {
	  float DeltaR = ROOT::Math::VectorUtil::DeltaR(fourVecs[combinations[i][n]], fourVecs[combinations[j][n]]);
	  int charge_1 = good_charges[combinations[i][n]];
	  int charge_2 = good_charges[combinations[j][n]];
	                                                                                             //float psmass_candidate = (fourVecs[combinations[i][n]] + fourVecs[combinations[j][n]]).M();
	  if (DeltaR < 0.1 && charge_1 != charge_2) {                                                //if (DeltaR < 0.1 &&(std::abs(1.5 - psmass_candidate) < psmass_difference) && charge_1 != charge_2) {
	    count++;
	  }
	}
      }
      hmass_candidate = (hadron_1 + hadron_2 + hadron_3 + hadron_4).M();
      if (count >= 2) {
	//if (std::abs(125.2 - hmass_candidate) < hmass_difference || hmass_difference < 0) {
	//hmass_difference = std::abs(125.2 - hmass_candidate);
	  selected_hadron_indices[0] = original_pfcand_indices[combinations[0][n]];
	  selected_hadron_indices[1] = original_pfcand_indices[combinations[1][n]];
	  selected_hadron_indices[2] = original_pfcand_indices[combinations[2][n]];
	  selected_hadron_indices[3] = original_pfcand_indices[combinations[3][n]];
	  //	}
      }
    }

    selected_hadrons = {static_cast<int>(selected_hadron_indices[0]), 
			static_cast<int>(selected_hadron_indices[1]), 
			static_cast<int>(selected_hadron_indices[2]), 
			static_cast<int>(selected_hadron_indices[3])};

    std::sort(selected_hadrons.begin(), selected_hadrons.end(), [&pfcand_pt](int a, int b) {
        return pfcand_pt[a] > pfcand_pt[b];
      }
      );

    return selected_hadrons;
  };
  
  auto df1 = df.Define(str_daughteridxs, hidx, {str_pfcand_pt, str_pfcand_eta, str_pfcand_phi, str_pfcand_mass, str_pfcand_charge, str_pfcand_mask});
  
  return df1;

}

ROOT::RDF::RNode FourHardestPFCandsAlgo(ROOT::RDF::RNode df, const std::string &str_pfcand_pt, const std::string &str_pfcand_mask, const std::string &str_daughteridxs) {
  auto hidx = [](const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<int> pfcand_mask) {

    auto original_pfcand_indices = ROOT::VecOps::Nonzero(pfcand_mask);
    ROOT::RVec<int> selected_hadrons = {-1, -1, -1, -1}; 
    
    std::sort(original_pfcand_indices.begin(), original_pfcand_indices.end(), [&pfcand_pt](int a, int b) {
	return pfcand_pt[a] > pfcand_pt[b];
      }
      );
    selected_hadrons = {static_cast<int>(original_pfcand_indices[0]),
                        static_cast<int>(original_pfcand_indices[1]),
                        static_cast<int>(original_pfcand_indices[2]),
                        static_cast<int>(original_pfcand_indices[3])};

    return selected_hadrons;
  };
  auto df1 = df.Define(str_daughteridxs, hidx, {str_pfcand_pt, str_pfcand_mask});
  return df1;
}

ROOT::RDF::RNode ChargePairsAlgo(ROOT::RDF::RNode df, const std::string &str_pfcand_pt, const std::string &str_pfcand_eta, const std::string &str_pfcand_phi, const std::string &str_pfcand_mass, const std::string &str_pfcand_charge, 
				 const std::string &str_pfcand_mask, const std::string &str_daughteridxs) {
/*
  auto iso = [](const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<float> pfcand_eta, const ROOT::RVec<float> pfcand_phi, ROOT::RVec<float> pfcand_mass){
    ROOT::RVec<float> sums_pt;
    auto pfcand = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pfcand_pt, pfcand_eta, pfcand_phi, pfcand_mass);
    for (int idx = 0; idx < pfcand_pt.size(); idx++) {
      float sum_pt = 0;
      for (int i = 0; i < pfcand_pt.size(); i++) {
        if (i == idx) continue;
          float dR = ROOT::Math::VectorUtil::DeltaR(pfcand[idx], pfcand[i]);
        if (dR < 0.4) {
          sum_pt += pfcand_pt[i];
        }
      }
      sums_pt.push_back(sum_pt);
    }
    return sums_pt;
  };
*/
  auto hidx = [](const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<int> pfcand_charge, const ROOT::RVec<int> pfcand_mask) {

    auto original_pfcand_indices = ROOT::VecOps::Nonzero(pfcand_mask);
    ROOT::RVec<int> selected_hadrons = {-1, -1, -1, -1};

    if (original_pfcand_indices.size() < 4) {
      return selected_hadrons;
    }
    
    std::sort(original_pfcand_indices.begin(), original_pfcand_indices.end(), [&pfcand_pt](int a, int b) {
	return pfcand_pt[a] > pfcand_pt[b];
      }
      );
    int countpos = 0;
    int countneg = 0;
    ROOT::RVec<int> two_hardest_pos = {-1, -1};
    ROOT::RVec<int> two_hardest_neg = {-1, -1};
    for (int i = 0; i < original_pfcand_indices.size(); i++) {
      if (pfcand_charge[original_pfcand_indices[i]] == 1 && countpos < 2) {
	two_hardest_pos[countpos] = original_pfcand_indices[i];
	countpos++;
      }
      if (pfcand_charge[original_pfcand_indices[i]] == -1 && countneg < 2) {
	two_hardest_neg[countneg] = original_pfcand_indices[i];
	countneg++;
      }
    }

    if (countpos < 2 || countneg < 2) {
      return selected_hadrons;
    }

    selected_hadrons = {static_cast<int>(two_hardest_pos[0]), 
			static_cast<int>(two_hardest_neg[0]),
			static_cast<int>(two_hardest_pos[1]),
			static_cast<int>(two_hardest_neg[1])};
    return selected_hadrons;
  };
  //auto df1 = df.Define("PFCands_iso", iso, {str_pfcand_pt, str_pfcand_eta, str_pfcand_phi, str_pfcand_mass});
  auto df1 = df.Define(str_daughteridxs, hidx, {str_pfcand_pt, str_pfcand_charge, str_pfcand_mask});
  return df1;
  }


ROOT::RDF::RNode GetHiggsP4(ROOT::RDF::RNode df, const std::string &str_d1_p4, const std::string &str_d2_p4, 
			    const std::string &str_d3_p4, const std::string &str_d4_p4, const std::string &str_H_p4) {
  auto H_p4 = [](const ROOT::Math::PtEtaPhiMVector d1_p4, const ROOT::Math::PtEtaPhiMVector d2_p4, const ROOT::Math::PtEtaPhiMVector d3_p4,                                                                
		 const ROOT::Math::PtEtaPhiMVector d4_p4) {
    auto H_p4 = d1_p4 + d2_p4 + d3_p4 + d4_p4;
    return H_p4;
  }; 
  auto df1 = df.Define(str_H_p4, H_p4, {str_d1_p4, str_d2_p4, str_d3_p4, str_d4_p4});
  return df1;
}

ROOT::RDF::RNode GetPseudoScalars(ROOT::RDF::RNode df, const std::string &daughterIdx, const std::string &str_pfcand_pt, const std::string &str_pfcand_eta, const std::string &str_pfcand_phi, 
                                  const std::string &str_pfcand_mass, const std::string &str_pfcand_charge, const std::string &str_pfcand_mask, const std::string &str_highPtPair, 
                                  const std::string &str_lowPtPair) {
  auto highPtPair = [](const ROOT::RVec<int> daughterIdx, const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<float> pfcand_eta, const ROOT::RVec<float> pfcand_phi, const ROOT::RVec<float> pfcand_mass, const ROOT::RVec<int> pfcand_charge, const ROOT::RVec<int> pfcand_mask) {
    
    const auto good_pts = ROOT::VecOps::Take(pfcand_pt, daughterIdx);
    const auto good_etas = ROOT::VecOps::Take(pfcand_eta, daughterIdx);
    const auto good_phis = ROOT::VecOps::Take(pfcand_phi, daughterIdx);
    const auto good_masses = ROOT::VecOps::Take(pfcand_mass, daughterIdx);
    const auto good_charges = ROOT::VecOps::Take(pfcand_charge, daughterIdx);

    ROOT::RVec<int> highPtPairIdxs = {-1, -1};  
    float highPt = 0;

    auto combinations = ROOT::VecOps::Combinations(daughterIdx, 2);
    for (int n = 0; n < combinations[0].size(); n++) {
      auto daughter_1 = ROOT::Math::PtEtaPhiMVector(good_pts[combinations[0][n]], good_etas[combinations[0][n]], good_phis[combinations[0][n]], good_masses[combinations[0][n]]);
      auto daughter_2 = ROOT::Math::PtEtaPhiMVector(good_pts[combinations[1][n]], good_etas[combinations[1][n]], good_phis[combinations[1][n]], good_masses[combinations[1][n]]);
      int charge_1 = good_charges[combinations[0][n]];
      int charge_2 = good_charges[combinations[1][n]];
      if (charge_1 != charge_2) {
        auto ps = daughter_1 + daughter_2;
        if (ps.Pt() > highPt) {
          highPt = ps.Pt();
          highPtPairIdxs = {static_cast<int>(daughterIdx[combinations[0][n]]), static_cast<int>(daughterIdx[combinations[1][n]])};
        }
      }
    }
    return highPtPairIdxs;
  };

  auto lowPtPair = [](const ROOT::RVec<int> daughterIdx, const ROOT::RVec<int> highPtPairIdxs) {
    ROOT::RVec<int> lowPtPairIdxs = {-1, -1};
    for (auto idx : daughterIdx) {
      if (idx != highPtPairIdxs[0] && idx != highPtPairIdxs[1]) {
        if (lowPtPairIdxs[0] == -1) {
          lowPtPairIdxs[0] = idx;
        } else {
          lowPtPairIdxs[1] = idx;
        }
      }
    }

    return lowPtPairIdxs;

  };

  auto df1 = df.Define(str_highPtPair, highPtPair, {daughterIdx, str_pfcand_pt, str_pfcand_eta, str_pfcand_phi, str_pfcand_mass, str_pfcand_charge, str_pfcand_mask});
  auto df2 = df1.Define(str_lowPtPair, lowPtPair, {daughterIdx, str_highPtPair});

  return df2;
    
}

ROOT::RDF::RNode GetMinMassDiff(ROOT::RDF::RNode df, const std::string &daughterIdx, const std::string &str_pfcand_pt, const std::string &str_pfcand_eta, const std::string &str_pfcand_phi, 
                               const std::string &str_pfcand_mass, const std::string &str_pfcand_charge, const std::string &str_pfcand_mask, const std::string &str_ps1Pair, 
                               const std::string &str_ps2Pair) {
  auto minMassDiff = [](const ROOT::RVec<int> daughterIdx, const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<float> pfcand_eta, 
                        const ROOT::RVec<float> pfcand_phi, const ROOT::RVec<float> pfcand_mass, const ROOT::RVec<int> pfcand_charge, const ROOT::RVec<int> pfcand_mask) {
                        
    //const auto good_pts = ROOT::VecOps::Take(pfcand_pt, daughterIdx);
    //const auto good_etas = ROOT::VecOps::Take(pfcand_eta, daughterIdx);
    //const auto good_phis = ROOT::VecOps::Take(pfcand_phi, daughterIdx);
    //const auto good_masses = ROOT::VecOps::Take(pfcand_mass, daughterIdx);
    //const auto good_charges = ROOT::VecOps::Take(pfcand_charge, daughterIdx);

    ROOT::RVec<int> minMassDiffIdxs = {-1, -1};
    ROOT::RVec<int> daughters = daughterIdx;
    std::sort(daughters.begin(), daughters.end());
    float minMassDiff = 9999.0;
    do {
      if (pfcand_charge[daughters[0]] != pfcand_charge[daughters[1]] && pfcand_charge[daughters[2]] != pfcand_charge[daughters[3]]) {
        auto daughter_1 = ROOT::Math::PtEtaPhiMVector(pfcand_pt[daughters[0]], pfcand_eta[daughters[0]], pfcand_phi[daughters[0]], pfcand_mass[daughters[0]]);
        auto daughter_2 = ROOT::Math::PtEtaPhiMVector(pfcand_pt[daughters[1]], pfcand_eta[daughters[1]], pfcand_phi[daughters[1]], pfcand_mass[daughters[1]]);
        auto daughter_3 = ROOT::Math::PtEtaPhiMVector(pfcand_pt[daughters[2]], pfcand_eta[daughters[2]], pfcand_phi[daughters[2]], pfcand_mass[daughters[2]]);
        auto daughter_4 = ROOT::Math::PtEtaPhiMVector(pfcand_pt[daughters[3]], pfcand_eta[daughters[3]], pfcand_phi[daughters[3]], pfcand_mass[daughters[3]]);
        auto ps1 = daughter_1 + daughter_2;
        auto ps2 = daughter_3 + daughter_4;
        float massDiff = std::abs(ps1.M() - ps2.M());
        if (massDiff < minMassDiff) {
          minMassDiff = massDiff;
          minMassDiffIdxs = {daughters[0], daughters[1]};
        }
      }
    }
    while (std::next_permutation(daughters.begin(), daughters.end()));

    return minMassDiffIdxs;
  };

  auto partners = [](const ROOT::RVec<int> daughterIdx, const ROOT::RVec<int> minMassDiffIdxs) {
    ROOT::RVec<int> partners = {-1, -1};
    for (auto idx : daughterIdx) {
      if (idx != minMassDiffIdxs[0] && idx != minMassDiffIdxs[1]) {
        if (partners[0] == -1) {
          partners[0] = idx;
        } else {
          partners[1] = idx;
        }
      }
    }
    return partners;
  };

  auto df1 = df.Define(str_ps1Pair, minMassDiff, {daughterIdx, str_pfcand_pt, str_pfcand_eta, str_pfcand_phi, str_pfcand_mass, str_pfcand_charge, str_pfcand_mask});
  auto df2 = df1.Define(str_ps2Pair, partners, {daughterIdx, str_ps1Pair});

  return df2;

}

ROOT::RDF::RNode pfCandIso(ROOT::RDF::RNode df, const std::string &str_pf_cand_iso, const std::string &str_pfcand_pt, const std::string &str_pfcand_eta, const std::string &str_pfcand_phi, 
                           const std::string &str_pfcand_mass, const std::string &str_higgsdaughters, const std::string &str_pfcand_charged_hadron_mask, 
                           const std::string &str_pfcand_neutral_hadron_mask, const std::string &str_pfcand_photon_mask, const std::string &str_pfcand_fromPV_mask) {
  auto iso = [](const ROOT::RVec<float> pfcand_pt, const ROOT::RVec<float> pfcand_eta, const ROOT::RVec<float> pfcand_phi, const ROOT::RVec<float> pfcand_mass, 
                const ROOT::RVec<int> higgsdaughters, const ROOT::RVec<int> pfcand_charged_hadron_mask, const ROOT::RVec<int> pfcand_neutral_hadron_mask,
                const ROOT::RVec<int> pfcand_photon_mask, const ROOT::RVec<int> pfcand_fromPV_mask) {
    ROOT::RVec<float> Rel_isos;

    auto daughter_pts =  ROOT::VecOps::Take(pfcand_pt, higgsdaughters);
    auto daughter_etas = ROOT::VecOps::Take(pfcand_eta, higgsdaughters);
    auto daughter_phis = ROOT::VecOps::Take(pfcand_phi, higgsdaughters);
    ROOT::RVec<float> daughter_masses = {0.493677, 0.493677, 0.493677, 0.493677};

    auto daughter_p4s = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(daughter_pts, daughter_etas, daughter_phis, daughter_masses);

    auto non_PU_charged_hadron_indices = ROOT::VecOps::Nonzero(pfcand_charged_hadron_mask && pfcand_fromPV_mask);
    auto neutral_hadron_indices = ROOT::VecOps::Nonzero(pfcand_neutral_hadron_mask);
    auto photon_indices = ROOT::VecOps::Nonzero(pfcand_photon_mask);
    auto PU_charged_hadron_indices = ROOT::VecOps::Nonzero(pfcand_charged_hadron_mask && !pfcand_fromPV_mask);

    auto charged_hadron_pts_from_pv = ROOT::VecOps::Take(pfcand_pt, non_PU_charged_hadron_indices);
    auto charged_hadron_etas_from_pv = ROOT::VecOps::Take(pfcand_eta, non_PU_charged_hadron_indices);
    auto charged_hadron_phis_from_pv = ROOT::VecOps::Take(pfcand_phi, non_PU_charged_hadron_indices);
    auto charged_hadron_masses_from_pv = ROOT::VecOps::Take(pfcand_mass, non_PU_charged_hadron_indices);

    auto charged_hadron_p4s_from_pv = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(charged_hadron_pts_from_pv, charged_hadron_etas_from_pv, charged_hadron_phis_from_pv, charged_hadron_masses_from_pv);

    auto neutral_hadron_pts = ROOT::VecOps::Take(pfcand_pt, neutral_hadron_indices);
    auto neutral_hadron_etas = ROOT::VecOps::Take(pfcand_eta, neutral_hadron_indices);
    auto neutral_hadron_phis = ROOT::VecOps::Take(pfcand_phi, neutral_hadron_indices);
    auto neutral_hadron_masses = ROOT::VecOps::Take(pfcand_mass, neutral_hadron_indices);

    auto neutral_hadron_p4s = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(neutral_hadron_pts, neutral_hadron_etas, neutral_hadron_phis, neutral_hadron_masses);

    auto photon_pts = ROOT::VecOps::Take(pfcand_pt, photon_indices);
    auto photon_etas = ROOT::VecOps::Take(pfcand_eta, photon_indices);
    auto photon_phis = ROOT::VecOps::Take(pfcand_phi, photon_indices);
    auto photon_masses = ROOT::VecOps::Take(pfcand_mass, photon_indices);

    auto photon_p4s = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(photon_pts, photon_etas, photon_phis, photon_masses);

    auto charged_hadron_pts_not_from_pv = ROOT::VecOps::Take(pfcand_pt, PU_charged_hadron_indices);
    auto charged_hadron_etas_not_from_pv = ROOT::VecOps::Take(pfcand_eta, PU_charged_hadron_indices);
    auto charged_hadron_phis_not_from_pv = ROOT::VecOps::Take(pfcand_phi, PU_charged_hadron_indices);
    auto charged_hadron_masses_not_from_pv = ROOT::VecOps::Take(pfcand_mass, PU_charged_hadron_indices);

    auto charged_hadron_p4s_not_from_pv = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(charged_hadron_pts_not_from_pv, charged_hadron_etas_not_from_pv, charged_hadron_phis_not_from_pv, charged_hadron_masses_not_from_pv);
    
    for (int i = 0; i < daughter_p4s.size(); i++) {
      float sum_charged_hadron_pt_from_pv = 0;
      float sum_neutral_hadron_pt = 0;
      float sum_photon_pt = 0;
      float sum_charged_hadron_pt_not_from_pv = 0;
      
      for (int j = 0; j < charged_hadron_p4s_from_pv.size(); j++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(daughter_p4s[i], charged_hadron_p4s_from_pv[j]);
        if (dR < 0.4 && std::find(higgsdaughters.begin(), higgsdaughters.end(), non_PU_charged_hadron_indices[j]) == higgsdaughters.end()) {
          sum_charged_hadron_pt_from_pv += charged_hadron_pts_from_pv[j];
        }
      }
      
      for (int j = 0; j < neutral_hadron_p4s.size(); j++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(daughter_p4s[i], neutral_hadron_p4s[j]);
        if (dR < 0.4) {
          sum_neutral_hadron_pt += neutral_hadron_pts[j];
        }
      }
      for (int j = 0; j < photon_p4s.size(); j++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(daughter_p4s[i], photon_p4s[j]);
        if (dR < 0.4) {
          sum_photon_pt += photon_pts[j];
        }
      }
      for (int j = 0; j < charged_hadron_p4s_not_from_pv.size(); j++) {
        float dR = ROOT::Math::VectorUtil::DeltaR(daughter_p4s[i], charged_hadron_p4s_not_from_pv[j]);
        if (dR < 0.4 && std::find(higgsdaughters.begin(), higgsdaughters.end(), PU_charged_hadron_indices[j]) == higgsdaughters.end()) {
          sum_charged_hadron_pt_not_from_pv += charged_hadron_pts_not_from_pv[j];
        }
      }
      float Rel_iso = (sum_charged_hadron_pt_from_pv + std::max(0.0, sum_neutral_hadron_pt + sum_photon_pt - 0.5*sum_charged_hadron_pt_not_from_pv))/daughter_p4s[i].Pt();
      Rel_isos.push_back(Rel_iso);
    }

    return Rel_isos;
  };


  auto df1 = df.Define(str_pf_cand_iso, iso, {str_pfcand_pt, str_pfcand_eta, str_pfcand_phi, str_pfcand_mass, str_higgsdaughters, str_pfcand_charged_hadron_mask, str_pfcand_neutral_hadron_mask, str_pfcand_photon_mask, str_pfcand_fromPV_mask});
  return df1;
}

ROOT::RDF::RNode getGenPt(ROOT::RDF::RNode df, const std::string &str_genpart_mask, const std::string &str_genpart_pt, const std::string &str_gen_pt){ 
  auto df1 = df.Define(
		       str_gen_pt,
		       [](const ROOT::RVec<int> genpart_mask, const ROOT::RVec<float> genpart_pt){
			 auto original_genpart_indices = ROOT::VecOps::Nonzero(genpart_mask);
			 const auto good_pts = ROOT::VecOps::Take(genpart_pt, original_genpart_indices);
			 return good_pts[0];
		       },
		       {str_genpart_mask, str_genpart_pt}
		       );
  return df1;
}

}

#endif /* GUARD_HAA_H */
