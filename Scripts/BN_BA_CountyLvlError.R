library(dplyr)
library(RSQLite)

DataBase=dbConnect(SQLite(), "../shared/FIADB_ENTIRE/SQLite_FIADB_ENTIRE.db")
#setwd("/Volumes/Onni_Drive/*****_Projects/Butternut_Landscape_Genomics")
#DataBase=dbConnect(SQLite(), "../../FIADB_ENTIRE/SQLite_FIADB_ENTIRE.db")

County_Lvl_Error = dbGetQuery(DataBase, "SELECT 
      state,
      county,
      long,
      lat,
      SUM(estimate_by_estn_unit.estimate) AS basal_area_ftsqrd,
      CASE
        WHEN SUM(estimate_by_estn_unit.estimate) <> 0 THEN
          abs(sqrt(sum(estimate_by_estn_unit.var_of_estimate)) / sum(estimate_by_estn_unit.estimate) * 100)
        ELSE
          0
      END AS CI_pct_of_estimate_68pct,
      sqrt(sum(estimate_by_estn_unit.var_of_estimate)) AS se_of_estimate_68pct,
      SUM(estimate_by_estn_unit.var_of_estimate) AS var_of_estimate,
      SUM(estimate_by_estn_unit.total_plots) AS total_plots, 
      SUM(estimate_by_estn_unit.non_zero_plots) AS non_zero_plots, 
      SUM(estimate_by_estn_unit.total_population_area_acres) AS total_population_acres
    FROM (SELECT pop_eval_grp_cn,
    phase_2_summary.state state,
    phase_2_summary.county county,
    phase_2_summary.plot plot,
    phase_2_summary.long long,
    phase_2_summary.lat lat,
    estn_unit_cn,
    SUM(total_area*(coalesce(ysum_hd, 0) / phase_1_summary.n_h) * w_h) estimate,
     SUM(phase_1_summary.n_h) total_plots,
     SUM(phase_2_summary.number_plots_in_domain) domain_plots,
     SUM(phase_2_summary.non_zero_plots) non_zero_plots,
     total_area*total_area/SUM(phase_1_summary.n_h) *
       ((SUM(w_h * phase_1_summary.n_h *
               (((coalesce(ysum_hd_sqr, 0)/phase_1_summary.n_h)-
                   ((coalesce(ysum_hd, 0)/phase_1_summary.n_h)*
                      (coalesce(ysum_hd, 0)/phase_1_summary.n_h)))/
                  (phase_1_summary.n_h-1)))) +
          1/SUM(phase_1_summary.n_h) *
          (SUM((1-w_h) * phase_1_summary.n_h *
                 (((COALESCE(ysum_hd_sqr, 0)/phase_1_summary.n_h) -
                     ((COALESCE(ysum_hd, 0)/phase_1_summary.n_h) *
                        (COALESCE(ysum_hd, 0)/phase_1_summary.n_h)))/
                    (phase_1_summary.n_h-1))))) var_of_estimate,
     total_area total_population_area_acres
     
--Phase 1 Summary			
     FROM(SELECT pop_eval_grp.eval_grp,
          pop_eval_grp.cn pop_eval_grp_cn,
          pop_stratum.estn_unit_cn,
          pop_stratum.cn pop_stratum_cn,
          CAST(p1pointcnt AS REAL)/
            (SELECT SUM(str.p1pointcnt)
             FROM pop_stratum str
             WHERE str.estn_unit_cn=pop_stratum.estn_unit_cn) w_h,
          (SELECT SUM(str.p1pointcnt)
           FROM pop_stratum str
           WHERE str.estn_unit_cn=pop_stratum.estn_unit_cn) n_prime,
          p1pointcnt n_prime_h,
          (SELECT SUM(eu_s.area_used)
           FROM pop_estn_unit eu_s
           WHERE eu_s.cn=pop_stratum.estn_unit_cn) total_area,
          pop_stratum.p2pointcnt n_h
          FROM pop_eval_grp,
          pop_eval,
          pop_eval_typ,
          pop_estn_unit,
          pop_stratum
          WHERE pop_eval_typ.eval_grp_cn=pop_eval_grp.cn
          AND pop_eval.cn=pop_eval_typ.eval_cn
          AND pop_estn_unit.cn=pop_stratum.estn_unit_cn
          AND CAST(pop_eval_grp.eval_grp IN (12020,
	                         52020,
	                         92020,
	                         102020,
	                         112020,
	                         122020,
	                         132020,
	                         172020,
	                         182020,
	                         192020,
	                         212020,
	                         222020,
	                         232020,
	                         242020,
	                         252020,
	                         262020,
	                         272020,
	                         282020,
	                         292020,
	                         332020,
	                         342020,
	                         362020,
	                         372020,
	                         392020,
	                         422020,
	                         442020,
	                         452020,
	                         472020,
	                         502020,
	                         512020,
	                         542020,
	                         552020) AS INTEGER) --statecd/year to update
          AND pop_eval_typ.eval_typ='EXPVOL') phase_1_summary,
     
--Phase 2 Summary
(SELECT pop_stratum_cn,
  state,
  county,
  plot,
  long,
  lat,
  SUM(y_hid_adjusted) ysum_hd,
  SUM(y_hid_adjusted*y_hid_adjusted) ysum_hd_sqr,
  COUNT(*) number_plots_in_domain,
  SUM(CASE y_hid_adjusted
      WHEN 0 THEN 
      0 
      WHEN NULL THEN
      0
      ELSE 
      1
      END) non_zero_plots
  FROM(SELECT pop_stratum.cn pop_stratum_cn,
       p.cn plt_cn,
       p.plot plot,
       p.lon long,
       p.lat lat,
       t.statecd state,
       t.countycd county,
       SUM(COALESCE(T.DIA*T.DIA*0.005454*T.TPA_UNADJ * 
      CASE WHEN T.DIA IS NULL THEN POP_STRATUM.ADJ_FACTOR_SUBP 
      ELSE 
      CASE MIN(T.DIA, 5 - 0.001) WHEN T.DIA THEN POP_STRATUM.ADJ_FACTOR_MICR 
      ELSE 
      CASE MIN(T.DIA, COALESCE(P.MACRO_BREAKPOINT_DIA, 9999) - 0.001) WHEN T.DIA THEN POP_STRATUM.ADJ_FACTOR_SUBP 
      ELSE 
      POP_STRATUM.ADJ_FACTOR_MACR END END END, 0)) y_hid_adjusted
       FROM pop_eval_grp, 
       pop_eval_typ, 
       pop_eval, 
       pop_estn_unit, 
       pop_stratum, 
       pop_plot_stratum_assgn, 
       plot p, 
       plotgeom, 
       cond c, 
       tree t
       WHERE pop_eval_typ.eval_grp_cn=pop_eval_grp.cn 
       AND pop_eval.cn=pop_eval_typ.eval_cn 
       AND pop_eval.cn=pop_estn_unit.eval_cn 
       AND pop_estn_unit.cn=pop_stratum.estn_unit_cn 
       AND pop_plot_stratum_assgn.stratum_cn=pop_stratum.cn 
       AND pop_plot_stratum_assgn.plt_cn=p.cn 
       AND p.cn=plotgeom.cn 
       AND c.plt_cn=p.cn 
       AND t.plt_cn=c.plt_cn 
       AND t.condid=c.condid
       AND c.cond_status_cd = 1
       AND t.spcd IN (601)
       AND t.statuscd=1
       AND pop_eval_typ.eval_typ='EXPVOL'
       AND CAST(pop_eval_grp.eval_grp IN (12020,
	                         52020,
	                         92020,
	                         102020,
	                         112020,
	                         122020,
	                         132020,
	                         172020,
	                         182020,
	                         192020,
	                         212020,
	                         222020,
	                         232020,
	                         242020,
	                         252020,
	                         262020,
	                         272020,
	                         282020,
	                         292020,
	                         332020,
	                         342020,
	                         362020,
	                         372020,
	                         392020,
	                         422020,
	                         442020,
	                         452020,
	                         472020,
	                         502020,
	                         512020,
	                         542020,
	                         552020) AS INTEGER)
       GROUP BY pop_stratum.cn,
       p.cn,
       t.statecd,
       t.countycd)
  GROUP BY pop_stratum_cn,
            state,
            county) phase_2_summary
WHERE phase_1_summary.pop_stratum_cn = phase_2_summary.pop_stratum_cn
GROUP BY pop_eval_grp_cn,
         estn_unit_cn,
         phase_1_summary.total_area,
         phase_2_summary.state,
         phase_2_summary.county) estimate_by_estn_unit 
GROUP BY pop_eval_grp_cn,
         state,
         county")  
County_Lvl_Error <- County_Lvl_Error %>% distinct(.keep_all = TRUE)
County_Lvl_Error
County_Lvl_Error$basal_area_msqrd <- County_Lvl_Error$basal_area_ftsqrd/10.764
County_Lvl_Error$se_basal_area_msqrd <- County_Lvl_Error$se_of_estimate_68pct/10.764
County_Lvl_Error$ba_msqrd_per_ha <- County_Lvl_Error$basal_area_msqrd/(County_Lvl_Error$total_population_acres/2.471)
County_Lvl_Error$se_ba_msqrd_per_ha <- County_Lvl_Error$se_basal_area_msqrd/(County_Lvl_Error$total_population_acres/2.471)

write.csv(County_Lvl_Error,"./Results/20250207_PlotLvlError/County_Level_Error_BasalArea_2020_v2.csv", row.names=FALSE)



