### Created 26 Nov 2014 

MHmakeRandomString <- function(n=1, lenght=12)
{
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                 lenght, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}



create_sashimi_plot_settings <- function(outPath, bam_prefix, bam_files, coverages, colors, ymax = 300, logged = "False", fig_width = 20, fig_height = 10, font_size = 8, nyticks = 3, nxticks = 4)
{
	

	sashimi_plot_settings <- c("[data]",
	paste0("bam_prefix = ", bam_prefix),
	"miso_prefix =",

	paste0("bam_files = [", paste0(paste0("\"",bam_files, "\""), collapse=",") ,"]"),

	"miso_files = []",

	"[plotting]",
	"# Dimensions of figure to be plotted (in inches)",
	paste0("fig_width = ",fig_width),
	paste0("fig_height = ",fig_height),
	"# Factor to scale down introns and exons by",
	"intron_scale = 30",
	"exon_scale = 4",
	"# Whether to use a log scale or not when plotting",
	paste0("logged = ",logged),
	paste0("font_size = ", font_size),

	"bar_posteriors = False",

	"# Max y-axis",
	paste0("ymax = ",ymax),

	"# Axis tick marks",
	paste0("nyticks = ",nyticks),
	paste0("nxticks = ",nxticks),

	"# Whether to show axis labels",
	"show_ylabel = False",
	"show_xlabel = False",

	"# Whether to plot posterior distributions inferred by MISO",
	"show_posteriors = False ",

	"# Whether to plot the number of reads in each junction",
	"number_junctions = True",

	"resolution = .5",
	"posterior_bins = 40",
	"gene_posterior_ratio = 5",

	"# List of colors for read denisites of each sample",
	paste0("colors = [", paste0(paste0("\"",colors, "\""), collapse=",") ,"]"),

	"# Number of mapped reads in each sample",
	"# (Used to normalize the read density for RPKM calculation)",
	paste0("coverages = [", paste0(coverages, collapse = ","), "]"),

	"# Bar color for Bayes factor distribution",
	"# plots (--plot-bf-dist)",
	"# Paint them blue",
	"bar_color = \"b\"",

	"# Bayes factors thresholds to use for --plot-bf-dist",
	"bf_thresholds = [0, 1, 2, 5, 10, 20]", sep=" ")
	
	
	file <- paste0(outPath, "/sashimi_plot_settings_", MHmakeRandomString(), ".txt")
	
	write.table(sashimi_plot_settings, file=file, quote=FALSE, col.names=FALSE, row.names=FALSE)
	
	return(file)
	
}


# create_sashimi_plot_settings(outPath, bam_prefix, bam_files, coverages, colors, ymax = 300, logged = "False", fig_width = 20, fig_height = 10, font_size = 8, nyticks = 3, nxticks = 4)




plot_sashimi <- function(miso_path, miso_index, genes, outPath,bam_prefix, bam_files, coverages, colors, ymax = 300, logged = "False", fig_width = 20, fig_height = 10, font_size = 8, nyticks = 3, nxticks = 4)
{
	### settings 
	settings <- create_sashimi_plot_settings(outPath = outPath, bam_prefix = bam_prefix, bam_files = bam_files, coverages = coverages, colors = colors, ymax = ymax, logged = logged, fig_width = fig_width, fig_height = fig_height, font_size = font_size, nyticks = nyticks, nxticks = nxticks)
	
	### plotting
	cmd <- paste0(miso_path, "sashimi_plot/plot.py --plot-event ", genes, " ", miso_index, " ", settings, " --output-dir ", outPath)

	for(i in 1:length(cmd))
	  system(cmd[i])
	
	file.remove(settings)
	
}


