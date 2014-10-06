#void forest2yskim_jetSkim_forestV3(TString inputFile_="forestFiles",
#                                   std::string MinbiasFname = "for"
#                                   float cutphotonPt  = 35, 
#                                   sampleType colli=kHIDATA,
#                                   TString jetAlgo="akPu3PF",
#                                   bool doMix = true,
#                                   TString triggerSelection="",
#                                   bool doJetResCorrection = 0, 
#                                   int smearingCentBin = -1,
#                                   bool useGenJetColl = 0
#                                   )
#

#enum sampleType { kHIDATA, kHIMC, kPPDATA, kPPMC, kPADATA, kPAMC};
#                       0    1      2       3       4        5



# pbpb data


kHIDATA=0
kHIMC=1
kPPDATA=2

for jetAlgo in "akVs3Calo"
do
    input="forestFiles/dataPbPb/forest_pbpb_photon30_isoUpdated.root"
    inputSkim="forestFiles/dataPbPb/skim_collId_kHIDATA_jetAlgo_akVs3Calo_HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root"
    root -l -q -b 'forest2yskim_jetSkim.C++("'$input'", "'$inputSkim'", 30,   '$kHIDATA', "'$jetAlgo'", 1)'
done



