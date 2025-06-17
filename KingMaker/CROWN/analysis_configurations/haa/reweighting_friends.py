from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List, Union
import os
#from .producers import muon_sf_friends as muon_sf_friends
from .producers import pfcands as pfcands
from .producers import event as event
from .quantities import output as q
from .quantities import nanoAOD as nanoAOD
from code_generation.configuration import Configuration
from code_generation.modifiers import EraModifier
from code_generation.rules import RemoveProducer
from code_generation.systematics import SystematicShift


def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
    quantities_map: Union[str, None] = None,
):
    # for the test, we provide a quantities map
    '''
    quantities_map = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dyjets_shift_quantities_map.json",
    )
    '''
    configuration = Configuration(
        era,
        sample,
        scopes,
        shifts,
        available_sample_types,
        available_eras,
        available_scopes,
        #quantities_map,
    )

    configuration.add_config_parameters(
        ["mm"],
        {
            "PU_reweighting_file": EraModifier(
                {
                    "2016preVFP": "data/jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json.gz",
                    "2016postVFP": "data/jsonpog-integration/POG/LUM/2016postVFP_UL/puWeights.json.gz",
                    "2017": "data/jsonpog-integration/POG/LUM/2017_UL/puWeights.json.gz",
                    "2018": "data/jsonpog-integration/POG/LUM/2018_UL/puWeights.json.gz",
                }
            ),
            "PU_reweighting_era": EraModifier(
                {
                    "2016preVFP": "Collisions16_UltraLegacy_goldenJSON",
                    "2016postVFP": "Collisions16_UltraLegacy_goldenJSON",
                    "2017": "Collisions17_UltraLegacy_goldenJSON",
                    "2018": "Collisions18_UltraLegacy_goldenJSON",
                }
            ),
        },
    )

    configuration.add_producers(
        ["global"],
        [
            #event.PUweights,
            event.Lumi,
        ],
    )

    configuration.add_producers(
        ["mm"],
        [
            event.PUweights,
        ],
    )

    configuration.add_outputs(
        ["global"],
        [
            nanoAOD.PV_npvsGood,
            nanoAOD.Pileup_nTrueInt,
            #q.puweight,
            q.lumi,
            nanoAOD.event,
            nanoAOD.run,
            nanoAOD.luminosityBlock,
            nanoAOD.Pileup_nPU,
        ],
    )
    configuration.add_outputs(
        ["mm"],
        [
            q.puweight,
        ],
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
