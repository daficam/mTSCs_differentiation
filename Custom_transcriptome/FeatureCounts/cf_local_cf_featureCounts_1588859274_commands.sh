#!/bin/bash

# Cluster Flow local script
# Pipeline: featureCounts
# Created at 14:47, 07-05-2020

/home/rsh46/.clusterflow/modules/featureCounts.cfmod.pl --run_fn Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_featureCounts_235 --prev_job_id start_000 --cores 4 --mem 5G  2>&1 | sed s/^/###CF_cf_featureCounts_1588859274_featureCounts_235:/ >> Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts_clusterFlow.txt


/storage/Software/packages/clusterflow/modules/cf_run_finished.cfmod.pl --run_fn Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_email_run_complete_686 --prev_job_id null --cores 1 --mem 4G --param hide_log_header=true --param outfn=Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts_clusterFlow.txt >> Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts_clusterFlow.txt


/home/rsh46/.clusterflow/modules/featureCounts.cfmod.pl --run_fn lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_featureCounts_339 --prev_job_id start_000 --cores 4 --mem 5G  2>&1 | sed s/^/###CF_cf_featureCounts_1588859274_featureCounts_339:/ >> lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt


/storage/Software/packages/clusterflow/modules/cf_run_finished.cfmod.pl --run_fn lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_email_run_complete_049 --prev_job_id null --cores 1 --mem 4G --param hide_log_header=true --param outfn=lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt >> lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt


/home/rsh46/.clusterflow/modules/featureCounts.cfmod.pl --run_fn lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_featureCounts_499 --prev_job_id start_000 --cores 4 --mem 5G  2>&1 | sed s/^/###CF_cf_featureCounts_1588859274_featureCounts_499:/ >> lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt


/storage/Software/packages/clusterflow/modules/cf_run_finished.cfmod.pl --run_fn lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_email_run_complete_960 --prev_job_id null --cores 1 --mem 4G --param hide_log_header=true --param outfn=lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt >> lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt


/home/rsh46/.clusterflow/modules/featureCounts.cfmod.pl --run_fn lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_featureCounts_547 --prev_job_id start_000 --cores 4 --mem 5G  2>&1 | sed s/^/###CF_cf_featureCounts_1588859274_featureCounts_547:/ >> lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts_clusterFlow.txt


/storage/Software/packages/clusterflow/modules/cf_run_finished.cfmod.pl --run_fn lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_email_run_complete_875 --prev_job_id null --cores 1 --mem 4G --param hide_log_header=true --param outfn=lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts_clusterFlow.txt >> lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts_clusterFlow.txt


/storage/Software/packages/clusterflow/modules/cf_runs_all_finished.cfmod.pl --run_fn lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts.run --job_id cf_featureCounts_1588859274_email_pipeline_complete_050 --prev_job_id cf_featureCounts_1588859274_email_run_complete_875 --cores 1 --mem 4G --param summary_module=true --param hide_log_header=true --param outfn_0=lane3273_CGATGT_20_DM_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt --param outfn_1=Sample1009_lane7_TS_GFP_ACTTGA_L007_R1_trimmed.bam_featureCounts_clusterFlow.txt --param outfn_2=lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts_clusterFlow.txt --param outfn_3=lane3_scr-1_ATCACGAT_L003_R1_trimmed.bam_featureCounts_clusterFlow.txt >> lane8_24-DES1_ACAGTGAT_L008_R1_trimmed.bam_featureCounts_clusterFlow.txt


