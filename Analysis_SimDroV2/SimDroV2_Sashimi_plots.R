# Created 26 Nov 2014

#### metadata


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata





#############################
# sashimi plots
#############################

### test
# cd /usr/local/lib/python2.7/dist-packages/misopy-0.4.9-py2.7-linux-x86_64.egg/misopy/sashimi_plot/
#
# python ../index_gff.py --index test-data/events.gff /home/gosia/MISO/test-data/event-data/
#
#   python plot.py --plot-event "chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-" /home/gosia/MISO/test-data/event-data/ ../settings/sashimi_plot_settings.txt --output-dir /home/gosia/MISO/test-plot



miso_path <- "/usr/local/lib/python2.7/dist-packages/misopy-0.4.9-py2.7-linux-x86_64.egg/misopy/"
gff3 <- "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gff3"
miso_index <- "/home/Shared/data/annotation/Drosophila/Ensembl70/miso_0.4.9_index_byGosia"


####### create miso index - MUST CREATE AGAIN - problem with paths
# cmd <- paste0(miso_path, "index_gff.py --index ", gff3, " ", miso_index)
# cat(cmd)
# system(cmd)



######### entries for sashimi_plot_settings.txt

# ## Number of mapped reads in each sample
# cmd <- paste0("samtools view -c -F 4 ", bam_prefix, bam_files, "\n")
# #cmd <- paste0("samtools view -c -f 1 -F 12 sam/", bam.files, "\n")
# cat(cmd)
#
# for(i in 1:length(cmd))
#   system(cmd[i])

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plot_sashimi.R")

coverages <- c(40504622,40798950,40241213,42380898,40202220,41438659)
outPath <- "/home/gosia/Multinomial_project/Simulations_drosophila_V2/Sashimi_plots/"
bam_prefix <- "/home/Shared/tmp/Simulation/simulation_drosophila_V2/1_reads/bam/"
bam_files <- paste0("sim_samp", 1:6, "_s.bam")

colors <- c(rep("#3366FF", 3), rep("#00CC00", 3))
logged <- "False"
fig_width <- 20
fig_height <- 10
font_size <- 5
nyticks <- 2
nxticks <- 5



plot_sashimi(miso_path = miso_path, miso_index = miso_index, genes = "FBgn0000639", outPath = outPath ,bam_prefix = bam_prefix, bam_files = bam_files, coverages = coverages, colors = colors, ymax = 40, logged = logged, fig_width = fig_width, fig_height = fig_height, font_size = font_size, nyticks = nyticks, nxticks = nxticks)


plot_sashimi(miso_path = miso_path, miso_index = miso_index, genes = "FBgn0000173", outPath = outPath ,bam_prefix = bam_prefix, bam_files = bam_files, coverages = coverages, colors = colors, ymax = 20, logged = logged, fig_width = fig_width, fig_height = fig_height, font_size = font_size, nyticks = nyticks, nxticks = nxticks)







































































