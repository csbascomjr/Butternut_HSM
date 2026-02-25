library(dplyr)
library(RSQLite)
DataBase=dbConnect(SQLite(), "../shared/FIADB_ENTIRE/SQLite_FIADB_ENTIRE.db")
#DataBase=dbConnect(SQLite(), "../../FIADB_ENTIRE/SQLite_FIADB_ENTIRE.db")

AllSpecies_2020=dbGetQuery(DataBase,"SELECT 
state,
county,
EVAL_GRP, 
SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) osc_tns, 
CASE WHEN SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) <> 0 
     THEN ABS(SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) / 
          SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) * 100) ELSE 0 END AS SE_OF_ESTIMATE_PCT, 
SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) SE_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE) VAR_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOTAL_PLOTS) TOTAL_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOT_POP_AREA_ACRES) TOT_POP_AC 
FROM (SELECT POP_EVAL_GRP_CN, 
             EVAL_GRP, 
             PHASE_SUMMARY.state,
             PHASE_SUMMARY.county,
             SUM(COALESCE(YSUM_HD, 0) * PHASE_1_SUMMARY.EXPNS) ESTIMATE, 
             PHASE_1_SUMMARY.N TOTAL_PLOTS, 
             SUM(PHASE_SUMMARY.NUMBER_PLOTS_IN_DOMAIN) DOMAIN_PLOTS, 
             SUM(PHASE_SUMMARY.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
             TOTAL_AREA * TOTAL_AREA / PHASE_1_SUMMARY.N * 
              ((SUM(W_H * PHASE_1_SUMMARY.N_H * (((COALESCE(YSUM_HD_SQR, 0) / 
                PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H) * 
               (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / (PHASE_1_SUMMARY.N_H - 1)))) + 1 
              / PHASE_1_SUMMARY.N * (SUM((1 - W_H) * PHASE_1_SUMMARY.N_H * 
             (((COALESCE(YSUM_HD_SQR, 0) / PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) 
              / PHASE_1_SUMMARY.N_H) * (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / 
               (PHASE_1_SUMMARY.N_H - 1))))) VAR_OF_ESTIMATE, TOTAL_AREA TOT_POP_AREA_ACRES 
  FROM (SELECT PEV.CN EVAL_CN, 
               PEG.EVAL_GRP, 
               PEG.CN POP_EVAL_GRP_CN, 
               POP_STRATUM.ESTN_UNIT_CN, 
               POP_STRATUM.EXPNS, 
               POP_STRATUM.CN POP_STRATUM_CN, 
               P1POINTCNT / (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) W_H, 
              (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N_PRIME, 
               P1POINTCNT N_PRIME_H, 
              (SELECT SUM(EU_S.AREA_USED) FROM POP_ESTN_UNIT EU_S WHERE EU_S.CN = POP_STRATUM.ESTN_UNIT_CN) TOTAL_AREA, 
              (SELECT SUM(STR.P2POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N, 
               POP_STRATUM.P2POINTCNT N_H FROM  POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
               JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
               JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
               JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
               WHERE PEG.EVAL_GRP IN (12020,52020,92020,102020,122020,132020,172020,
                                               182020,192020,212020,222020,232020,242020,
                                               252020,262020,272020,282020,292020,332020,
                                               342020,362020,372020,392020,422020,442020,
                                               452020,472020,502020,512020,542020,552020) AND PET.EVAL_TYP = 'EXPCURR') PHASE_1_SUMMARY 
LEFT OUTER JOIN (SELECT 
                        state,
                        county,
                        POP_STRATUM_CN,
                        ESTN_UNIT_CN,
                        EVAL_CN, 
                        SUM(Y_HID_ADJUSTED) YSUM_HD, 
                        SUM(Y_HID_ADJUSTED * Y_HID_ADJUSTED) YSUM_HD_SQR, 
                        COUNT(*) NUMBER_PLOTS_IN_DOMAIN, 
                        SUM(CASE WHEN Y_HID_ADJUSTED IS NULL THEN 0 WHEN Y_HID_ADJUSTED = 0 THEN 0 ELSE 1 END) NON_ZERO_PLOTS 
                  FROM (SELECT 
                        SUM(COND.CONDPROP_UNADJ * COND.CARBON_SOIL_ORG * CASE COND.PROP_BASIS WHEN 'MACR' 
                            THEN POP_STRATUM.ADJ_FACTOR_MACR ELSE POP_STRATUM.ADJ_FACTOR_SUBP END) AS Y_HID_ADJUSTED,
                        PEU.CN ESTN_UNIT_CN, 
                        PEV.CN EVAL_CN, 
                        POP_STRATUM.CN POP_STRATUM_CN, 
                        PLOT.CN PLT_CN,
                        PLOT.STATECD state,
                        PLOT.COUNTYCD county
                  FROM POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
                                        JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
                                        JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
                                        JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
                                        JOIN POP_PLOT_STRATUM_ASSGN POP_PLOT_STRATUM_ASSGN ON (POP_PLOT_STRATUM_ASSGN.STRATUM_CN = POP_STRATUM.CN) 
                                        JOIN PLOT ON (POP_PLOT_STRATUM_ASSGN.PLT_CN = PLOT.CN) 
                                        JOIN PLOTGEOM ON (PLOT.CN = PLOTGEOM.CN) 
                                        JOIN COND ON (COND.PLT_CN = PLOT.CN) 
                  WHERE 
                      COND.COND_STATUS_CD = 1 
                  AND COND.CONDPROP_UNADJ IS NOT NULL 
                  AND COND.CARBON_SOIL_ORG IS NOT NULL 
                  AND PET.EVAL_TYP = 'EXPCURR' 
                  AND PEG.EVAL_GRP IN (12020,52020,92020,102020,122020,132020,172020,
                                               182020,192020,212020,222020,232020,242020,
                                               252020,262020,272020,282020,292020,332020,
                                               342020,362020,372020,392020,422020,442020,
                                               452020,472020,502020,512020,542020,552020)  
                  AND  1 = 1 
          GROUP BY PEU.CN, 
                   PEV.CN, 
                   POP_STRATUM.CN, 
                   PLOT.CN,
                   state,
                   county) PLOT_SUMMARY 
          GROUP BY POP_STRATUM_CN, 
                   ESTN_UNIT_CN, 
                   EVAL_CN,
                   state,
                   county) PHASE_SUMMARY 
          ON (PHASE_1_SUMMARY.POP_STRATUM_CN = PHASE_SUMMARY.POP_STRATUM_CN 
          AND PHASE_1_SUMMARY.EVAL_CN = PHASE_SUMMARY.EVAL_CN 
          AND PHASE_1_SUMMARY.ESTN_UNIT_CN = PHASE_SUMMARY.ESTN_UNIT_CN) 
          GROUP BY PHASE_1_SUMMARY.POP_EVAL_GRP_CN, 
                   PHASE_1_SUMMARY.EVAL_GRP,
                   PHASE_1_SUMMARY.ESTN_UNIT_CN, 
                   PHASE_1_SUMMARY.TOTAL_AREA, 
                   PHASE_1_SUMMARY.N,
                   PHASE_SUMMARY.state,
                   PHASE_SUMMARY.county) ESTIMATE_BY_ESTN_UNIT 
          WHERE NON_ZERO_PLOTS IS NOT NULL 
          GROUP BY POP_EVAL_GRP_CN, 
                   EVAL_GRP,
                   state,
                   county;")
AllSpecies_2020
beep(5)

write.csv(AllSpecies_2020, "./OSC_CountyLvl_2020.csv", row.names=FALSE)


##From GET ATTRIBUTE function
dbGetQuery(DataBase,"SELECT 
state,
county,
EVAL_GRP, 
SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) osc_tns, 
CASE WHEN SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) <> 0 
     THEN ABS(SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) / 
          SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) * 100) ELSE 0 END AS SE_OF_ESTIMATE_PCT, 
SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) SE_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE) VAR_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOTAL_PLOTS) TOTAL_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOT_POP_AREA_ACRES) TOT_POP_AC 
FROM (SELECT POP_EVAL_GRP_CN, 
             EVAL_GRP, 
             PHASE_SUMMARY.state,
             PHASE_SUMMARY.county,
             SUM(COALESCE(YSUM_HD, 0) * PHASE_1_SUMMARY.EXPNS) ESTIMATE, 
             PHASE_1_SUMMARY.N TOTAL_PLOTS, 
             SUM(PHASE_SUMMARY.NUMBER_PLOTS_IN_DOMAIN) DOMAIN_PLOTS, 
             SUM(PHASE_SUMMARY.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
             TOTAL_AREA * TOTAL_AREA / PHASE_1_SUMMARY.N * 
              ((SUM(W_H * PHASE_1_SUMMARY.N_H * (((COALESCE(YSUM_HD_SQR, 0) / 
                PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H) * 
               (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / (PHASE_1_SUMMARY.N_H - 1)))) + 1 
              / PHASE_1_SUMMARY.N * (SUM((1 - W_H) * PHASE_1_SUMMARY.N_H * 
             (((COALESCE(YSUM_HD_SQR, 0) / PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) 
              / PHASE_1_SUMMARY.N_H) * (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / 
               (PHASE_1_SUMMARY.N_H - 1))))) VAR_OF_ESTIMATE, TOTAL_AREA TOT_POP_AREA_ACRES 
  FROM (SELECT PEV.CN EVAL_CN, 
               PEG.EVAL_GRP, 
               PEG.CN POP_EVAL_GRP_CN, 
               POP_STRATUM.ESTN_UNIT_CN, 
               POP_STRATUM.EXPNS, 
               POP_STRATUM.CN POP_STRATUM_CN, 
               P1POINTCNT / (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) W_H, 
              (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N_PRIME, 
               P1POINTCNT N_PRIME_H, 
              (SELECT SUM(EU_S.AREA_USED) FROM POP_ESTN_UNIT EU_S WHERE EU_S.CN = POP_STRATUM.ESTN_UNIT_CN) TOTAL_AREA, 
              (SELECT SUM(STR.P2POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N, 
               POP_STRATUM.P2POINTCNT N_H FROM  POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
               JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
               JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
               JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
               WHERE PEG.EVAL_GRP IN (442020) AND PET.EVAL_TYP = 'EXPCURR') PHASE_1_SUMMARY 
LEFT OUTER JOIN (SELECT 
                        state,
                        county,
                        POP_STRATUM_CN,
                        ESTN_UNIT_CN,
                        EVAL_CN, 
                        SUM(Y_HID_ADJUSTED) YSUM_HD, 
                        SUM(Y_HID_ADJUSTED * Y_HID_ADJUSTED) YSUM_HD_SQR, 
                        COUNT(*) NUMBER_PLOTS_IN_DOMAIN, 
                        SUM(CASE WHEN Y_HID_ADJUSTED IS NULL THEN 0 WHEN Y_HID_ADJUSTED = 0 THEN 0 ELSE 1 END) NON_ZERO_PLOTS 
                  FROM (SELECT 
                        SUM(COND.CONDPROP_UNADJ * COND.CARBON_SOIL_ORG * CASE COND.PROP_BASIS WHEN 'MACR' 
                            THEN POP_STRATUM.ADJ_FACTOR_MACR ELSE POP_STRATUM.ADJ_FACTOR_SUBP END) AS Y_HID_ADJUSTED,
                        PEU.CN ESTN_UNIT_CN, 
                        PEV.CN EVAL_CN, 
                        POP_STRATUM.CN POP_STRATUM_CN, 
                        PLOT.CN PLT_CN,
                        PLOT.STATECD state,
                        PLOT.COUNTYCD county
                  FROM POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
                                        JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
                                        JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
                                        JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
                                        JOIN POP_PLOT_STRATUM_ASSGN POP_PLOT_STRATUM_ASSGN ON (POP_PLOT_STRATUM_ASSGN.STRATUM_CN = POP_STRATUM.CN) 
                                        JOIN PLOT ON (POP_PLOT_STRATUM_ASSGN.PLT_CN = PLOT.CN) 
                                        JOIN PLOTGEOM ON (PLOT.CN = PLOTGEOM.CN) 
                                        JOIN COND ON (COND.PLT_CN = PLOT.CN) 
                  WHERE 
                      COND.COND_STATUS_CD = 1 
                  AND COND.CONDPROP_UNADJ IS NOT NULL 
                  AND COND.CARBON_SOIL_ORG IS NOT NULL 
                  AND PET.EVAL_TYP = 'EXPCURR' 
                  AND PEG.EVAL_GRP IN (442020)  
                  AND  1 = 1 
          GROUP BY PEU.CN, 
                   PEV.CN, 
                   POP_STRATUM.CN, 
                   PLOT.CN,
                   state,
                   county) PLOT_SUMMARY 
          GROUP BY POP_STRATUM_CN, 
                   ESTN_UNIT_CN, 
                   EVAL_CN,
                   state,
                   county) PHASE_SUMMARY 
          ON (PHASE_1_SUMMARY.POP_STRATUM_CN = PHASE_SUMMARY.POP_STRATUM_CN 
          AND PHASE_1_SUMMARY.EVAL_CN = PHASE_SUMMARY.EVAL_CN 
          AND PHASE_1_SUMMARY.ESTN_UNIT_CN = PHASE_SUMMARY.ESTN_UNIT_CN) 
          GROUP BY PHASE_1_SUMMARY.POP_EVAL_GRP_CN, 
                   PHASE_1_SUMMARY.EVAL_GRP,
                   PHASE_1_SUMMARY.ESTN_UNIT_CN, 
                   PHASE_1_SUMMARY.TOTAL_AREA, 
                   PHASE_1_SUMMARY.N,
                   PHASE_SUMMARY.state,
                   PHASE_SUMMARY.county) ESTIMATE_BY_ESTN_UNIT 
          WHERE NON_ZERO_PLOTS IS NOT NULL 
          GROUP BY POP_EVAL_GRP_CN, 
                   EVAL_GRP,
                   state,
                   county;")

###From the above code, 5.2.2025 @9:16AM
state county EVAL_GRP    osc_tns SE_OF_ESTIMATE_PCT SE_OF_ESTIMATE VAR_OF_ESTIMATE TOTAL_PLOTS NON_ZERO_PLOTS TOT_POP_AC
1    44      1   442020   579825.8           0.000000            0.0               0         204              3   670107.7
2    44      3   442020  3318184.6           0.000000            0.0               0         204             17   670107.7
3    44      5   442020  1068748.9           0.000000            0.0               0         172              6   574656.7
4    44      7   442020 11985251.5           2.586569       310006.8     96104232682         231             61   781971.6
5    44      9   442020  7982286.2           0.000000            0.0               0         204             43   670107.7

##From EVALIDATOR
#fips   osc         non_zero_plots
#44001  591881       3
#44003  3322676     17
#44005  1073922     6
#44007  12004524    61
#44009  8027472     43


####Test a couple states
dbGetQuery(DataBase,"SELECT 
state,
county,
EVAL_GRP, 
SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) osc_tns, 
CASE WHEN SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) <> 0 
     THEN ABS(SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) / 
          SUM(ESTIMATE_BY_ESTN_UNIT.ESTIMATE) * 100) ELSE 0 END AS SE_OF_ESTIMATE_PCT, 
SQRT(SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE)) SE_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.VAR_OF_ESTIMATE) VAR_OF_ESTIMATE, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOTAL_PLOTS) TOTAL_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
     SUM(ESTIMATE_BY_ESTN_UNIT.TOT_POP_AREA_ACRES) TOT_POP_AC 
FROM (SELECT POP_EVAL_GRP_CN, 
             EVAL_GRP, 
             PHASE_SUMMARY.state,
             PHASE_SUMMARY.county,
             SUM(COALESCE(YSUM_HD, 0) * PHASE_1_SUMMARY.EXPNS) ESTIMATE, 
             PHASE_1_SUMMARY.N TOTAL_PLOTS, 
             SUM(PHASE_SUMMARY.NUMBER_PLOTS_IN_DOMAIN) DOMAIN_PLOTS, 
             SUM(PHASE_SUMMARY.NON_ZERO_PLOTS) NON_ZERO_PLOTS, 
             TOTAL_AREA * TOTAL_AREA / PHASE_1_SUMMARY.N * 
              ((SUM(W_H * PHASE_1_SUMMARY.N_H * (((COALESCE(YSUM_HD_SQR, 0) / 
                PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H) * 
               (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / (PHASE_1_SUMMARY.N_H - 1)))) + 1 
              / PHASE_1_SUMMARY.N * (SUM((1 - W_H) * PHASE_1_SUMMARY.N_H * 
             (((COALESCE(YSUM_HD_SQR, 0) / PHASE_1_SUMMARY.N_H) - ((COALESCE(YSUM_HD, 0) 
              / PHASE_1_SUMMARY.N_H) * (COALESCE(YSUM_HD, 0) / PHASE_1_SUMMARY.N_H))) / 
               (PHASE_1_SUMMARY.N_H - 1))))) VAR_OF_ESTIMATE, TOTAL_AREA TOT_POP_AREA_ACRES 
  FROM (SELECT PEV.CN EVAL_CN, 
               PEG.EVAL_GRP, 
               PEG.CN POP_EVAL_GRP_CN, 
               POP_STRATUM.ESTN_UNIT_CN, 
               POP_STRATUM.EXPNS, 
               POP_STRATUM.CN POP_STRATUM_CN, 
               P1POINTCNT / (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) W_H, 
              (SELECT SUM(STR.P1POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N_PRIME, 
               P1POINTCNT N_PRIME_H, 
              (SELECT SUM(EU_S.AREA_USED) FROM POP_ESTN_UNIT EU_S WHERE EU_S.CN = POP_STRATUM.ESTN_UNIT_CN) TOTAL_AREA, 
              (SELECT SUM(STR.P2POINTCNT) FROM POP_STRATUM STR WHERE STR.ESTN_UNIT_CN = POP_STRATUM.ESTN_UNIT_CN) N, 
               POP_STRATUM.P2POINTCNT N_H FROM  POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
               JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
               JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
               JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
               WHERE PEG.EVAL_GRP IN (542020) AND PET.EVAL_TYP = 'EXPCURR') PHASE_1_SUMMARY 
LEFT OUTER JOIN (SELECT 
                        state,
                        county,
                        POP_STRATUM_CN,
                        ESTN_UNIT_CN,
                        EVAL_CN, 
                        SUM(Y_HID_ADJUSTED) YSUM_HD, 
                        SUM(Y_HID_ADJUSTED * Y_HID_ADJUSTED) YSUM_HD_SQR, 
                        COUNT(*) NUMBER_PLOTS_IN_DOMAIN, 
                        SUM(CASE WHEN Y_HID_ADJUSTED IS NULL THEN 0 WHEN Y_HID_ADJUSTED = 0 THEN 0 ELSE 1 END) NON_ZERO_PLOTS 
                  FROM (SELECT 
                        SUM(COND.CONDPROP_UNADJ * COND.CARBON_SOIL_ORG * CASE COND.PROP_BASIS WHEN 'MACR' 
                            THEN POP_STRATUM.ADJ_FACTOR_MACR ELSE POP_STRATUM.ADJ_FACTOR_SUBP END) AS Y_HID_ADJUSTED,
                        PEU.CN ESTN_UNIT_CN, 
                        PEV.CN EVAL_CN, 
                        POP_STRATUM.CN POP_STRATUM_CN, 
                        PLOT.CN PLT_CN,
                        PLOT.STATECD state,
                        PLOT.COUNTYCD county
                  FROM POP_EVAL_GRP PEG JOIN POP_EVAL_TYP PET ON (PET.EVAL_GRP_CN = PEG.CN) 
                                        JOIN POP_EVAL PEV ON (PEV.CN = PET.EVAL_CN) 
                                        JOIN POP_ESTN_UNIT PEU ON (PEV.CN = PEU.EVAL_CN) 
                                        JOIN POP_STRATUM POP_STRATUM ON (PEU.CN = POP_STRATUM.ESTN_UNIT_CN) 
                                        JOIN POP_PLOT_STRATUM_ASSGN POP_PLOT_STRATUM_ASSGN ON (POP_PLOT_STRATUM_ASSGN.STRATUM_CN = POP_STRATUM.CN) 
                                        JOIN PLOT ON (POP_PLOT_STRATUM_ASSGN.PLT_CN = PLOT.CN) 
                                        JOIN PLOTGEOM ON (PLOT.CN = PLOTGEOM.CN) 
                                        JOIN COND ON (COND.PLT_CN = PLOT.CN) 
                  WHERE 
                      COND.COND_STATUS_CD = 1 
                  AND COND.CONDPROP_UNADJ IS NOT NULL 
                  AND COND.CARBON_SOIL_ORG IS NOT NULL 
                  AND PET.EVAL_TYP = 'EXPCURR' 
                  AND PEG.EVAL_GRP IN (542020)  
                  AND  1 = 1 
          GROUP BY PEU.CN, 
                   PEV.CN, 
                   POP_STRATUM.CN, 
                   PLOT.CN,
                   state,
                   county) PLOT_SUMMARY 
          GROUP BY POP_STRATUM_CN, 
                   ESTN_UNIT_CN, 
                   EVAL_CN,
                   state,
                   county) PHASE_SUMMARY 
          ON (PHASE_1_SUMMARY.POP_STRATUM_CN = PHASE_SUMMARY.POP_STRATUM_CN 
          AND PHASE_1_SUMMARY.EVAL_CN = PHASE_SUMMARY.EVAL_CN 
          AND PHASE_1_SUMMARY.ESTN_UNIT_CN = PHASE_SUMMARY.ESTN_UNIT_CN) 
          GROUP BY PHASE_1_SUMMARY.POP_EVAL_GRP_CN, 
                   PHASE_1_SUMMARY.EVAL_GRP,
                   PHASE_1_SUMMARY.ESTN_UNIT_CN, 
                   PHASE_1_SUMMARY.TOTAL_AREA, 
                   PHASE_1_SUMMARY.N,
                   PHASE_SUMMARY.state,
                   PHASE_SUMMARY.county) ESTIMATE_BY_ESTN_UNIT 
          WHERE NON_ZERO_PLOTS IS NOT NULL 
          GROUP BY POP_EVAL_GRP_CN, 
                   EVAL_GRP,
                   state,
                   county;")


