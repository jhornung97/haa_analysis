ifeq ($(strip $(FWCore/Version)),)
FWCoreVersion := self/FWCore/Version
FWCore/Version := FWCoreVersion
FWCoreVersion_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
FWCoreVersion_EX_USE := $(foreach d, self cmssw ,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
FWCoreVersion_EX_LIB   := FWCoreVersion
ALL_EXTERNAL_PRODS += FWCoreVersion
FWCoreVersion_CLASS := LIBRARY
FWCore/Version_relbigobj+=FWCoreVersion
endif
ifeq ($(strip $(GeneratorInterfaceCore_plugins)),)
GeneratorInterfaceCore_plugins := self/src/GeneratorInterface/Core/plugins
GeneratorInterfaceCore_plugins_LOC_USE := self cmssw FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger FWCore/Utilities SimDataFormats/GeneratorProducts GeneratorInterface/Core FWCore/SharedMemory
ALL_EXTERNAL_PLUGIN_PRODS += GeneratorInterfaceCore_plugins
GeneratorInterface/Core_relbigobj+=GeneratorInterfaceCore_plugins
endif
ifeq ($(strip $(GeneratorInterface/Core)),)
GeneratorInterfaceCore := self/GeneratorInterface/Core
GeneratorInterface/Core := GeneratorInterfaceCore
GeneratorInterfaceCore_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
GeneratorInterfaceCore_EX_USE := $(foreach d, self cmssw FWCore/Concurrency FWCore/ServiceRegistry FWCore/Utilities FWCore/Framework SimDataFormats/GeneratorProducts GeneratorInterface/LHEInterface heppdt boost clhep lhapdf f77compiler root,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
GeneratorInterfaceCore_EX_LIB   := GeneratorInterfaceCore
ALL_EXTERNAL_PRODS += GeneratorInterfaceCore
GeneratorInterfaceCore_CLASS := LIBRARY
GeneratorInterface/Core_relbigobj+=GeneratorInterfaceCore
endif
ifeq ($(strip $(GeneratorInterfaceGenFiltersPlugins)),)
GeneratorInterfaceGenFiltersPlugins := self/src/GeneratorInterface/GenFilters/plugins
GeneratorInterfaceGenFiltersPlugins_LOC_USE := self cmssw FWCore/ParameterSet FWCore/Framework FWCore/Utilities SimDataFormats/GeneratorProducts DataFormats/HepMCCandidate DataFormats/JetReco
ALL_EXTERNAL_PLUGIN_PRODS += GeneratorInterfaceGenFiltersPlugins
GeneratorInterface/GenFilters_relbigobj+=GeneratorInterfaceGenFiltersPlugins
endif
ifeq ($(strip $(GeneratorInterface/GenFilters)),)
GeneratorInterfaceGenFilters := self/GeneratorInterface/GenFilters
GeneratorInterface/GenFilters := GeneratorInterfaceGenFilters
GeneratorInterfaceGenFilters_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
GeneratorInterfaceGenFilters_EX_USE := $(foreach d, self cmssw fastjet boost FWCore/PluginManager FWCore/ParameterSet FWCore/Framework FWCore/Utilities SimDataFormats/GeneratorProducts heppdt clhep root GeneratorInterface/Pythia6Interface GeneratorInterface/Core pythia8 SimGeneral/HepPDTRecord DataFormats/GeometryVector DataFormats/GeometrySurface TrackPropagation/SteppingHelixPropagator MagneticField/Records TrackingTools/TrajectoryState TrackingTools/TrajectoryParametrization TrackingTools/Records CommonTools/UtilAlgos FWCore/ServiceRegistry CommonTools/BaseParticlePropagator TrackingTools/GeomPropagators DataFormats/HepMCCandidate DataFormats/JetReco DataFormats/EgammaReco DataFormats/Math,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += GeneratorInterfaceGenFilters
GeneratorInterfaceGenFilters_CLASS := LIBRARY
GeneratorInterface/GenFilters_relbigobj+=GeneratorInterfaceGenFilters
endif
ifeq ($(strip $(SimDataFormats/ValidationFormats)),)
SimDataFormatsValidationFormats := self/SimDataFormats/ValidationFormats
SimDataFormats/ValidationFormats := SimDataFormatsValidationFormats
SimDataFormatsValidationFormats_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimDataFormatsValidationFormats_EX_USE := $(foreach d, self cmssw DataFormats/Common DataFormats/GeometryVector clhep geant4core expat,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimDataFormatsValidationFormats_LCGDICTS  := x 
SimDataFormatsValidationFormats_EX_LIB   := SimDataFormatsValidationFormats
ALL_EXTERNAL_PRODS += SimDataFormatsValidationFormats
SimDataFormatsValidationFormats_CLASS := LIBRARY
SimDataFormats/ValidationFormats_relbigobj+=SimDataFormatsValidationFormats
endif
ifeq ($(strip $(SimG4CMSCaloPlugins)),)
SimG4CMSCaloPlugins := self/src/SimG4CMS/Calo/plugins
SimG4CMSCaloPlugins_LOC_USE := self cmssw FWCore/PluginManager SimG4Core/Watcher SimG4CMS/Calo SimDataFormats/GeneratorProducts
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CMSCaloPlugins
SimG4CMS/Calo_relbigobj+=SimG4CMSCaloPlugins
endif
ifeq ($(strip $(SimG4CMS/Calo)),)
SimG4CMSCalo := self/SimG4CMS/Calo
SimG4CMS/Calo := SimG4CMSCalo
SimG4CMSCalo_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CMSCalo_EX_USE := $(foreach d, self cmssw FWCore/PluginManager SimGeneral/GFlash SimG4Core/SensitiveDetector SimG4Core/Notification DataFormats/EcalDetId DataFormats/HcalDetId DataFormats/ForwardDetId DataFormats/Math CondFormats/HcalObjects CondFormats/DataRecord SimDataFormats/SimHitMaker SimDataFormats/CaloHit SimDataFormats/CaloTest Geometry/HcalCommonData Geometry/EcalCommonData Geometry/HGCalCommonData DetectorDescription/Core FWCore/ParameterSet FWCore/MessageLogger FWCore/ServiceRegistry CommonTools/UtilAlgos boost clhep geant4core hepmc root rootmath,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CMSCalo_EX_LIB   := SimG4CMSCalo
ALL_EXTERNAL_PRODS += SimG4CMSCalo
SimG4CMSCalo_CLASS := LIBRARY
SimG4CMS/Calo_relbigobj+=SimG4CMSCalo
endif
ifeq ($(strip $(SimG4CMSFP420Plugins)),)
SimG4CMSFP420Plugins := self/src/SimG4CMS/FP420/plugins
SimG4CMSFP420Plugins_LOC_FLAGS_REM_CXXFLAGS   := -Werror=unused-variable
SimG4CMSFP420Plugins_LOC_USE := self cmssw SimG4CMS/FP420
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CMSFP420Plugins
SimG4CMS/FP420_relbigobj+=SimG4CMSFP420Plugins
endif
ifeq ($(strip $(SimG4CMS/FP420)),)
SimG4CMSFP420 := self/SimG4CMS/FP420
SimG4CMS/FP420 := SimG4CMSFP420
SimG4CMSFP420_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CMSFP420_EX_USE := $(foreach d, self cmssw FWCore/PluginManager SimG4Core/Watcher SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Application SimG4CMS/Calo SimDataFormats/SimHitMaker DataFormats/GeometryVector SimDataFormats/CaloHit DetectorDescription/Core FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger SimGeneral/HepPDTRecord boost geant4core root,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CMSFP420_EX_LIB   := SimG4CMSFP420
ALL_EXTERNAL_PRODS += SimG4CMSFP420
SimG4CMSFP420_CLASS := LIBRARY
SimG4CMS/FP420_relbigobj+=SimG4CMSFP420
endif
ifeq ($(strip $(SimG4CMSForwardPlugins)),)
SimG4CMSForwardPlugins := self/src/SimG4CMS/Forward/plugins
SimG4CMSForwardPlugins_LOC_USE := self cmssw FWCore/PluginManager SimG4Core/Watcher SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Physics SimG4CMS/Calo SimG4CMS/Forward DataFormats/Math SimDataFormats/SimHitMaker SimDataFormats/CaloHit SimDataFormats/CaloTest SimDataFormats/Forward DetectorDescription/Core FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger geant4core root rootmath
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CMSForwardPlugins
SimG4CMS/Forward_relbigobj+=SimG4CMSForwardPlugins
endif
ifeq ($(strip $(SimG4CMS/Forward)),)
SimG4CMSForward := self/SimG4CMS/Forward
SimG4CMS/Forward := SimG4CMSForward
SimG4CMSForward_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CMSForward_EX_USE := $(foreach d, self cmssw FWCore/PluginManager FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger SimG4Core/Watcher SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Physics SimG4CMS/Calo SimG4CMS/Tracker DataFormats/ForwardDetId DataFormats/Math SimDataFormats/SimHitMaker SimDataFormats/CaloHit SimDataFormats/Forward DetectorDescription/Core Geometry/HGCalCommonData Geometry/MTDCommonData boost clhep geant4core root rootmath,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CMSForward_EX_LIB   := SimG4CMSForward
ALL_EXTERNAL_PRODS += SimG4CMSForward
SimG4CMSForward_CLASS := LIBRARY
SimG4CMS/Forward_relbigobj+=SimG4CMSForward
endif
ifeq ($(strip $(SimG4CMS/Muon)),)
SimG4CMSMuon := self/SimG4CMS/Muon
SimG4CMS/Muon := SimG4CMSMuon
SimG4CMSMuon_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CMSMuon_EX_USE := $(foreach d, self cmssw SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Physics SimDataFormats/SimHitMaker Geometry/MuonNumbering boost geant4core hepmc,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += SimG4CMSMuon
SimG4CMSMuon_CLASS := LIBRARY
SimG4CMS/Muon_relbigobj+=SimG4CMSMuon
endif
ifeq ($(strip $(SimG4CMSTrackerPlugins)),)
SimG4CMSTrackerPlugins := self/src/SimG4CMS/Tracker/plugins
SimG4CMSTrackerPlugins_LOC_USE := self cmssw SimG4CMS/Tracker FWCore/ParameterSet SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Physics Geometry/TrackerNumberingBuilder geant4core SimDataFormats/SimHitMaker
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CMSTrackerPlugins
SimG4CMS/Tracker_relbigobj+=SimG4CMSTrackerPlugins
endif
ifeq ($(strip $(SimG4CMS/Tracker)),)
SimG4CMSTracker := self/SimG4CMS/Tracker
SimG4CMS/Tracker := SimG4CMSTracker
SimG4CMSTracker_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CMSTracker_EX_USE := $(foreach d, self cmssw FWCore/ParameterSet SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Physics Geometry/TrackerNumberingBuilder geant4core SimDataFormats/SimHitMaker,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CMSTracker_EX_LIB   := SimG4CMSTracker
ALL_EXTERNAL_PRODS += SimG4CMSTracker
SimG4CMSTracker_CLASS := LIBRARY
SimG4CMS/Tracker_relbigobj+=SimG4CMSTracker
endif
ifeq ($(strip $(SimG4CoreApplicationPlugins)),)
SimG4CoreApplicationPlugins := self/src/SimG4Core/Application/plugins
SimG4CoreApplicationPlugins_LOC_USE := self cmssw FWCore/Concurrency FWCore/Framework FWCore/ParameterSet SimDataFormats/GeneratorProducts SimDataFormats/CaloHit SimDataFormats/Track SimDataFormats/TrackingHit SimDataFormats/Vertex geant4core hepmc SimG4Core/Application
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CoreApplicationPlugins
SimG4Core/Application_relbigobj+=SimG4CoreApplicationPlugins
endif
ifeq ($(strip $(SimG4Core/Application)),)
SimG4CoreApplication := self/SimG4Core/Application
SimG4Core/Application := SimG4CoreApplication
SimG4CoreApplication_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreApplication_EX_USE := $(foreach d, self cmssw DataFormats/Common DataFormats/Math SimDataFormats/GeneratorProducts SimDataFormats/Forward SimDataFormats/Track SimDataFormats/Vertex SimG4Core/Generators SimG4Core/Geometry SimG4Core/MagneticField SimG4Core/Notification SimG4Core/CustomPhysics SimG4Core/Physics SimG4Core/SensitiveDetector SimG4Core/Watcher SimGeneral/HepPDTRecord SimGeneral/GFlash FWCore/ParameterSet FWCore/PluginManager FWCore/Framework FWCore/Utilities MagneticField/Engine MagneticField/Records clhep xerces-c geant4core hepmc heppdt,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreApplication_EX_LIB   := SimG4CoreApplication
ALL_EXTERNAL_PRODS += SimG4CoreApplication
SimG4CoreApplication_CLASS := LIBRARY
SimG4Core/Application_relbigobj+=SimG4CoreApplication
endif
ifeq ($(strip $(SimG4CoreCustomPhysicsPlugins)),)
SimG4CoreCustomPhysicsPlugins := self/src/SimG4Core/CustomPhysics/plugins
SimG4CoreCustomPhysicsPlugins_LOC_USE := self cmssw SimG4Core/Watcher SimG4Core/Notification FWCore/ParameterSet FWCore/ServiceRegistry geant4core clhep boost SimG4Core/CustomPhysics
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CoreCustomPhysicsPlugins
SimG4Core/CustomPhysics_relbigobj+=SimG4CoreCustomPhysicsPlugins
endif
ifeq ($(strip $(SimG4Core/CustomPhysics)),)
SimG4CoreCustomPhysics := self/SimG4Core/CustomPhysics
SimG4Core/CustomPhysics := SimG4CoreCustomPhysics
SimG4CoreCustomPhysics_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreCustomPhysics_EX_USE := $(foreach d, self cmssw FWCore/Framework FWCore/PluginManager FWCore/MessageLogger SimG4Core/Physics SimG4Core/PhysicsLists geant4core clhep boost,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreCustomPhysics_EX_LIB   := SimG4CoreCustomPhysics
ALL_EXTERNAL_PRODS += SimG4CoreCustomPhysics
SimG4CoreCustomPhysics_CLASS := LIBRARY
SimG4Core/CustomPhysics_relbigobj+=SimG4CoreCustomPhysics
endif
ifeq ($(strip $(SimG4Core/Generators)),)
SimG4CoreGenerators := self/SimG4Core/Generators
SimG4Core/Generators := SimG4CoreGenerators
SimG4CoreGenerators_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreGenerators_EX_USE := $(foreach d, self cmssw FWCore/ParameterSet FWCore/MessageLogger SimDataFormats/GeneratorProducts boost heppdt geant4core rootmath expat,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreGenerators_EX_LIB   := SimG4CoreGenerators
ALL_EXTERNAL_PRODS += SimG4CoreGenerators
SimG4CoreGenerators_CLASS := LIBRARY
SimG4Core/Generators_relbigobj+=SimG4CoreGenerators
endif
ifeq ($(strip $(SimG4Core/Geometry)),)
SimG4CoreGeometry := self/SimG4Core/Geometry
SimG4Core/Geometry := SimG4CoreGeometry
SimG4CoreGeometry_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreGeometry_EX_USE := $(foreach d, self cmssw DetectorDescription/Core geant4core Geometry/Records FWCore/MessageLogger FWCore/Utilities xerces-c expat,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreGeometry_EX_LIB   := SimG4CoreGeometry
ALL_EXTERNAL_PRODS += SimG4CoreGeometry
SimG4CoreGeometry_CLASS := LIBRARY
SimG4Core/Geometry_relbigobj+=SimG4CoreGeometry
endif
ifeq ($(strip $(SimG4Core/HelpfulWatchers)),)
SimG4CoreHelpfulWatchers := self/SimG4Core/HelpfulWatchers
SimG4Core/HelpfulWatchers := SimG4CoreHelpfulWatchers
SimG4CoreHelpfulWatchers_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreHelpfulWatchers_EX_USE := $(foreach d, self cmssw FWCore/Framework FWCore/ParameterSet SimG4Core/Watcher SimG4Core/Notification boost geant4core CommonTools/UtilAlgos MagneticField/Engine MagneticField/Records SimG4Core/Physics,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += SimG4CoreHelpfulWatchers
SimG4CoreHelpfulWatchers_CLASS := LIBRARY
SimG4Core/HelpfulWatchers_relbigobj+=SimG4CoreHelpfulWatchers
endif
ifeq ($(strip $(SimG4Core/MagneticField)),)
SimG4CoreMagneticField := self/SimG4Core/MagneticField
SimG4Core/MagneticField := SimG4CoreMagneticField
SimG4CoreMagneticField_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreMagneticField_EX_USE := $(foreach d, self cmssw FWCore/PluginManager FWCore/ParameterSet boost geant4core expat,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreMagneticField_EX_LIB   := SimG4CoreMagneticField
ALL_EXTERNAL_PRODS += SimG4CoreMagneticField
SimG4CoreMagneticField_CLASS := LIBRARY
SimG4Core/MagneticField_relbigobj+=SimG4CoreMagneticField
endif
ifeq ($(strip $(SimG4Core/Notification)),)
SimG4CoreNotification := self/SimG4Core/Notification
SimG4Core/Notification := SimG4CoreNotification
SimG4CoreNotification_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreNotification_EX_USE := $(foreach d, self cmssw SimDataFormats/Forward SimDataFormats/Track SimDataFormats/Vertex DataFormats/Math geant4core FWCore/MessageLogger rootmath expat hepmc,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreNotification_EX_LIB   := SimG4CoreNotification
ALL_EXTERNAL_PRODS += SimG4CoreNotification
SimG4CoreNotification_CLASS := LIBRARY
SimG4Core/Notification_relbigobj+=SimG4CoreNotification
endif
ifeq ($(strip $(SimG4Core/Physics)),)
SimG4CorePhysics := self/SimG4Core/Physics
SimG4Core/Physics := SimG4CorePhysics
SimG4CorePhysics_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CorePhysics_EX_USE := $(foreach d, self cmssw DetectorDescription/Core FWCore/ParameterSet geant4core heppdt boost sigcpp expat,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CorePhysics_EX_LIB   := SimG4CorePhysics
ALL_EXTERNAL_PRODS += SimG4CorePhysics
SimG4CorePhysics_CLASS := LIBRARY
SimG4Core/Physics_relbigobj+=SimG4CorePhysics
endif
ifeq ($(strip $(SimG4CorePhysicsListsPlugins)),)
SimG4CorePhysicsListsPlugins := self/src/SimG4Core/PhysicsLists/plugins
SimG4CorePhysicsListsPlugins_LOC_USE := self cmssw FWCore/ParameterSet FWCore/MessageLogger FWCore/PluginManager SimG4Core/Physics geant4core clhep heppdt boost SimG4Core/PhysicsLists
ALL_EXTERNAL_PLUGIN_PRODS += SimG4CorePhysicsListsPlugins
SimG4Core/PhysicsLists_relbigobj+=SimG4CorePhysicsListsPlugins
endif
ifeq ($(strip $(SimG4Core/PhysicsLists)),)
SimG4CorePhysicsLists := self/SimG4Core/PhysicsLists
SimG4Core/PhysicsLists := SimG4CorePhysicsLists
SimG4CorePhysicsLists_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CorePhysicsLists_EX_USE := $(foreach d, self cmssw FWCore/MessageLogger SimG4Core/Physics SimG4Core/MagneticField geant4core clhep heppdt boost,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CorePhysicsLists_EX_LIB   := SimG4CorePhysicsLists
ALL_EXTERNAL_PRODS += SimG4CorePhysicsLists
SimG4CorePhysicsLists_CLASS := LIBRARY
SimG4Core/PhysicsLists_relbigobj+=SimG4CorePhysicsLists
endif
ifeq ($(strip $(SimG4Core/PrintGeomInfo)),)
SimG4CorePrintGeomInfo := self/SimG4Core/PrintGeomInfo
SimG4Core/PrintGeomInfo := SimG4CorePrintGeomInfo
SimG4CorePrintGeomInfo_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CorePrintGeomInfo_EX_USE := $(foreach d, self cmssw SimG4Core/Notification SimG4Core/Watcher SimG4Core/Geometry FWCore/ParameterSet geant4core boost,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += SimG4CorePrintGeomInfo
SimG4CorePrintGeomInfo_CLASS := LIBRARY
SimG4Core/PrintGeomInfo_relbigobj+=SimG4CorePrintGeomInfo
endif
ifeq ($(strip $(SimG4Core/SensitiveDetector)),)
SimG4CoreSensitiveDetector := self/SimG4Core/SensitiveDetector
SimG4Core/SensitiveDetector := SimG4CoreSensitiveDetector
SimG4CoreSensitiveDetector_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimG4CoreSensitiveDetector_EX_USE := $(foreach d, self cmssw FWCore/MessageLogger FWCore/PluginManager SimG4Core/Geometry SimG4Core/Notification geant4core,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimG4CoreSensitiveDetector_EX_LIB   := SimG4CoreSensitiveDetector
ALL_EXTERNAL_PRODS += SimG4CoreSensitiveDetector
SimG4CoreSensitiveDetector_CLASS := LIBRARY
SimG4Core/SensitiveDetector_relbigobj+=SimG4CoreSensitiveDetector
endif
ifeq ($(strip $(SimRomanPotSimFP420Plugins)),)
SimRomanPotSimFP420Plugins := self/src/SimRomanPot/SimFP420/plugins
SimRomanPotSimFP420Plugins_LOC_USE := self cmssw SimRomanPot/SimFP420
ALL_EXTERNAL_PLUGIN_PRODS += SimRomanPotSimFP420Plugins
SimRomanPot/SimFP420_relbigobj+=SimRomanPotSimFP420Plugins
endif
ifeq ($(strip $(SimRomanPot/SimFP420)),)
SimRomanPotSimFP420 := self/SimRomanPot/SimFP420
SimRomanPot/SimFP420 := SimRomanPotSimFP420
SimRomanPotSimFP420_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
SimRomanPotSimFP420_EX_USE := $(foreach d, self cmssw FWCore/PluginManager SimDataFormats/TrackingHit DataFormats/Math SimDataFormats/TrackerDigiSimLink SimDataFormats/CrossingFrame DataFormats/FP420Digi DataFormats/Common SimGeneral/HepPDTRecord SimG4Core/Watcher SimG4Core/SensitiveDetector SimG4Core/Notification SimG4Core/Application SimG4CMS/Calo SimG4CMS/FP420 SimDataFormats/SimHitMaker SimDataFormats/CaloHit Mixing/Base DetectorDescription/Core FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger boost geant4core root clhep hepmc gsl,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
SimRomanPotSimFP420_EX_LIB   := SimRomanPotSimFP420
ALL_EXTERNAL_PRODS += SimRomanPotSimFP420
SimRomanPotSimFP420_CLASS := LIBRARY
SimRomanPot/SimFP420_relbigobj+=SimRomanPotSimFP420
endif
ifeq ($(strip $(TauAnalysis/MCEmbeddingTools)),)
src_TauAnalysis_MCEmbeddingTools := self/TauAnalysis/MCEmbeddingTools
TauAnalysis/MCEmbeddingTools  := src_TauAnalysis_MCEmbeddingTools
src_TauAnalysis_MCEmbeddingTools_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
src_TauAnalysis_MCEmbeddingTools_EX_USE := $(foreach d, DataFormats/TrackReco DataFormats/Common self cmssw,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += src_TauAnalysis_MCEmbeddingTools
endif

ifeq ($(strip $(TauAnalysisMCEmbeddingToolsAuto)),)
TauAnalysisMCEmbeddingToolsAuto := self/src/TauAnalysis/MCEmbeddingTools/plugins
TauAnalysisMCEmbeddingToolsAuto_LOC_USE := self cmssw FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/PatCandidates pythia8 GeneratorInterface/Pythia8Interface GeneratorInterface/ExternalDecays GeneratorInterface/Herwig6Interface DataFormats/RPCRecHit DataFormats/CaloRecHit DataFormats/EcalRecHit DataFormats/HcalRecHit TrackingTools/TrackAssociator TrackingTools/TransientTrack DataFormats/TrackReco CondFormats/BeamSpotObjects
ALL_EXTERNAL_PLUGIN_PRODS += TauAnalysisMCEmbeddingToolsAuto
TauAnalysis/MCEmbeddingTools_relbigobj+=TauAnalysisMCEmbeddingToolsAuto
endif
ifeq ($(strip $(Validation/EcalHits)),)
ValidationEcalHits := self/Validation/EcalHits
Validation/EcalHits := ValidationEcalHits
ValidationEcalHits_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
ValidationEcalHits_EX_USE := $(foreach d, self cmssw DataFormats/Common root geant4core boost FWCore/ParameterSet SimG4Core/Watcher SimG4CMS/Calo DataFormats/EcalDetId SimDataFormats/GeneratorProducts SimDataFormats/ValidationFormats DQMServices/Core rootmath DataFormats/Math,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += ValidationEcalHits
ValidationEcalHits_CLASS := LIBRARY
Validation/EcalHits_relbigobj+=ValidationEcalHits
endif
ifeq ($(strip $(Validation/HcalHits)),)
ValidationHcalHits := self/Validation/HcalHits
Validation/HcalHits := ValidationHcalHits
ValidationHcalHits_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
ValidationHcalHits_EX_USE := $(foreach d, self cmssw FWCore/Framework FWCore/ParameterSet FWCore/MessageLogger geant4core boost SimG4Core/Notification SimG4Core/Watcher SimG4CMS/Calo Geometry/HcalCommonData DataFormats/HcalDetId SimDataFormats/CaloHit SimDataFormats/ValidationFormats SimDataFormats/GeneratorProducts hepmc DataFormats/Math rootmath DQMServices/Core DataFormats/HepMCCandidate,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += ValidationHcalHits
ValidationHcalHits_CLASS := LIBRARY
Validation/HcalHits_relbigobj+=ValidationHcalHits
endif
