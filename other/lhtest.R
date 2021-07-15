library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
folder <- dirname(args[1])
df <- read.table(args[1])
df[,3] <- 2*(df[,1]-df[,2])
names(df) <- c('global', 'multi', 'diff')

global.lnk <- as.numeric(args[2]) # -162576.606204
multi.lnk <- as.numeric(args[3]) # -149055.330013

print(df)

h <- ggplot(df, aes(x=diff)) + geom_histogram() + xlim(-12, 2) + xlab("2*(global_lnL-multi_lnL)") + ylab("Count") +
    theme_bw() +
    theme(# first 5 to make background white and draw axes lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          # now making stuff transparent
          plot.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          ## axis.ticks = element_line(),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.background=element_rect(fill = "transparent"),
          legend.title=element_text(size=12),
          legend.text=element_text(size=12)
    )

pdf(paste0(folder, "/", "lk_null.pdf"), width=4, height=4)
h + annotate("text", -8, 10, label=paste0("p-value = \n(Counts < ", 2*(global.lnk-multi.lnk), ")/100"), size=3)
dev.off()

