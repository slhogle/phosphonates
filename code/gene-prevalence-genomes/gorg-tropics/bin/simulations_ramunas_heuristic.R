library(tidyverse)

mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

df <- NULL
df <- data.frame()

for (i in 2:100) {
  for (j in c(1000, 1500, 2000, 2500)){
    for (k in c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8)){
      #genomes with gene
      y <- i
      #total genomes observed
      n <- 100
      #apparent frequency
      apparent_freq <- y/n
      # distribution for length of gene drawn from normal distribution
      p <- rnorm(y, mean=j, sd=150)
      # distribution for length of assemblies drawn from normal distribution
      g <- rnorm(n, mean=4e5, sd=1e4)
      # distribution for estimated completeness drawn from normal distribution
      c <- mysamp(n=n, m=k, s=0.1, lwr=0, upr=1, nnorm=1000)
      
      # observed fraction bp in gene
      #fobs <- sum(p)/sum(g)
      # true fraction bp assuming
      #ftrue_bp <- fobs/mean(c)
      # convert to genes/genome units
      #ftrue_g <- ftrue_bp * mean(g)/mean(p)
      #ftrue_g1 <- (y/n)/mean(c)
      
      # observed fraction bp in gene
      fobs <- sum(p)/sum(g)
      tlengenome <- sum(g/c)
      tlengene <- fobs*tlengenome
      # convert to genes/genome units
      ftrue_g <- tlengene/mean(p)/n
      
      # Ramunas' heuristic
      ftrue_g1 <- (y/n)/mean(c)
      
      df1 <- data.frame(
        genomes_w_gene = i,
        gene_len = j,
        completeness = k,
        estimated_freq = ftrue_g*100, 
        apparent_freq = apparent_freq*100,
        estimated_freq_heuristic = ftrue_g1*100
      )
      
      df <- bind_rows(df, df1)
    }
  }
}

t <- df %>%
  pivot_longer(c(estimated_freq, estimated_freq_heuristic)) %>%
  filter(value < 100)

p1 <- ggplot(t) + 
  geom_point(aes(x=apparent_freq, y=value, color=name, shape=factor(gene_len)), size=0.75) +
  facet_wrap(~ completeness, nrow=3, scales="free") +
  labs(x="Observed frequency [%]", y="Corrected frequency [%]", 
       shape="Gene length", color="Measure") 

p1 

ggsave(filename="/Users/shane/Desktop/test.png",
       device="png", plot = p1, width = 7, height = 3.5,  units = "in", dpi=300)