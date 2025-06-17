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
ifeq ($(strip $(L1Trigger/RPCTechnicalTrigger)),)
L1TriggerRPCTechnicalTrigger := self/L1Trigger/RPCTechnicalTrigger
L1Trigger/RPCTechnicalTrigger := L1TriggerRPCTechnicalTrigger
L1TriggerRPCTechnicalTrigger_BuildFile    := $(RELEASETOP)/.SCRAM/$(SCRAM_ARCH)/MakeData/DirCache.mk
L1TriggerRPCTechnicalTrigger_EX_USE := $(foreach d, self cmssw FWCore/Framework FWCore/MessageLogger FWCore/ParameterSet FWCore/ServiceRegistry DataFormats/RPCDigi DataFormats/L1GlobalTrigger DataFormats/L1GlobalMuonTrigger DataFormats/Common CondFormats/RPCObjects CondFormats/DataRecord Geometry/Records Geometry/RPCGeometry root,$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
ALL_EXTERNAL_PRODS += L1TriggerRPCTechnicalTrigger
L1TriggerRPCTechnicalTrigger_CLASS := LIBRARY
L1Trigger/RPCTechnicalTrigger_relbigobj+=L1TriggerRPCTechnicalTrigger
endif
