library("optparse")

option_list = list(make_option("--profile_file", type = "character", default = NULL, help  = "One profile file to be plotted."),
				   make_option("--sample_sheet", type = "character", default = NULL, help  = "A samplesheet to describe multiple profile files."),
				   make_option("--outTag", type = "character", default = NULL, help = "Prefix for output files."))		   
args <- parse_args(OptionParser(option_list=option_list))

# check parameters
if(is.null(args$profile_file))
	stop("\nArgument profile_file is required.\n")
if(is.null(args$outTag))
	stop("Argument outTag is required.\n")

library('ggplot2')

# args = list()
# args$profile_file = 'ff_profile.tsv'
# args$outTag = 'ff'

# a function to order windows
orderWindows = function(window_names){
	chunks = strsplit(window_names, '_')
	group = sapply(chunks, function(z) z[1])
	grank = as.numeric(sapply(chunks, function(z) z[2]))
	names(group) = names(grank) = window_names

	# order each group
	ugroup = c('EU', 'PU', 'P', 'PD', 'B', 'ED')
	res = c()
	for(i in 1:6){
		idx = group == ugroup[i]
		tmp = names(sort(grank[idx], decreasing=(i < 4)))
		res = c(res, tmp)
	}
	out = data.frame(name = res, group = group[res], rank = 1:length(res))
	return(out)
}

tx = read.table(file = args$profile_file, row.names = 1, header = TRUE, sep = '\t')
tmp = orderWindows(rownames(tx))
tx = tx[tmp$name, ]
out = data.frame(tx, tmp)

# plot
pdf(paste0(args$outTag, '_plot.pdf'), width = 8, height = 4)
	p <- ggplot(out, aes(x = rank, y = Score_All, color = group))
	p <- p + theme_bw() + labs(x = '', y = '', title = '') + ylim(0,1)
	p <- p + theme(plot.title = element_text(hjust = 0.5))
	p <- p + geom_line(linetype="dashed", color="blue", linewidth=0.5) + geom_point(size=0.8)
	p
dev.off()
