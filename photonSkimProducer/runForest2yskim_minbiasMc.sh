

#void forest2yskim_minbias_tracks(TString inputFile_="forestFiles/HiForest4/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root",
#                                   sampleType colli=kHIDATA,
#                                   int maxEvent = -1
#                                   )


####void forest2yskim_minbias_forestV3(TString inputFile_="forestFiles/HiForest4/HiForest_PbPb_MinBias_Track8_Jet5_GR_R_53_LV2B_merged_forest_0.root",
    ####                               sampleType colli=kHIDATA,
        #                           int maxEvent = -1,
         #                          bool useGenJetColl = 0,
          #                          TString jetAlgo="akPu3PF"
            #                        )
                                                                                                               
#enum sampleType { kHIDATA, kHIMC, kPPDATA, kPPMC, kPADATA, kPAMC};                                                          
#                   0 (X)    1 (x)  2       3       4        5                                                  


# PbPb mc
#root -l -q -b 'forest2yskim_minbias_tracks.C++("forestFiles/HiForest4/HiForest_HYDJET_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root",1,-1)'
root -l -q -b 'forest2yskim_minbias_forestV3.C++("forestFiles/HiForest4/HiForest_HYDJET_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root",  "forestFiles/HiForest4/trackSkim_collId_kHIMC_HiForest_HYDJET_Track8_Jet21_STARTHI53_LV1_merged_forest_0.root",1, -1)'
