# this script compared the correlation based analysis with the kinetic modeling result in the reference study
library(xlsx)

info = read.xlsx('allosteric_annotation.xlsx',sheetName = 'info')
# skip the unknowns(incomplete GPR)
info = info[info$consistency !='unknown',]
info$enzyme_contribution = info$enzyme_contribution * 100 # percentage

# PLOT
dev.off()
pdf('figures/correlation_vs_simmer.pdf',width = 4,height = 4.3)
boxplot(enzyme_contribution~correlated, data = info, outline= F,ylim = c(0,80),col = "white",
        ylab = 'metabolic leverage from enzyme (%)')
stripchart(enzyme_contribution~correlated, data = info,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  
dev.off()
