Collection of scripts for paper analysis. Run after pipeline. 


Some scripts should be run before others:
0-1 EDTA.sh
0-2 te_anno.sh

1-0 qc_ichip2.R
1-1a merge_bam.sh
1-1b chromhmm_makeRef.sh
1-2a bedmap.sh
1-2b log2r.sh
1-2c macs2.sh
1-2d chromhmm_binarize.sh
1-3a plot_genome.sh (may need some iChIP_prep.R lines to be run to get list of 5kb genes)
1-3b chromhmm_learnModel.sh

2-1 iChIP_prep.R
2-2a iChIP_plots.R
2-2b iChIP_plots_ac.R


P.S. I tried to tidy up the scripts but it is far from a "press and play" pipeline.
If something is confusing or unclear it's entirely my fault, and please contact me (Sean) for help.
I'd be happy to assist where I can!