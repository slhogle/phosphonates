Formatting PEPM abundances
================

  - [READ DATA](#read-data)
  - [FORMAT CORE AND PEPM RPKM](#format-core-and-pepm-rpkm)
      - [SAR11](#sar11)
      - [PROCHLOROCOCCUS](#prochlorococcus)
      - [BACTERIOPLANKTON](#bacterioplankton)
  - [ADD ENVIRONMENTAL DATA](#add-environmental-data)
      - [READ](#read)
      - [FORMAT](#format)
      - [SAVE](#save)

``` r
library(here)
library(tidyverse)
```

# READ DATA

``` r
pepm <- read_tsv(here::here("input", "pepm.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   sampleID = col_character(),
    ##   prochlorococcus = col_double(),
    ##   sar11 = col_double(),
    ##   bactarc = col_double()
    ## )

``` r
bactarc <- read_tsv(here::here("input", "core_fams_bacterioplankton.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   sampleID = col_character(),
    ##   COG0012 = col_double(),
    ##   COG0016 = col_double(),
    ##   COG0018 = col_double(),
    ##   COG0172 = col_double(),
    ##   COG0215 = col_double(),
    ##   COG0495 = col_double(),
    ##   COG0525 = col_double(),
    ##   COG0533 = col_double(),
    ##   COG0541 = col_double(),
    ##   COG0552 = col_double()
    ## )

``` r
pro <- read_tsv(here::here("input", "core_fams_prochlorococcus.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   sampleID = col_character(),
    ##   RpsJ = col_double(),
    ##   PsbJ = col_double(),
    ##   ComplexI_49KD = col_double(),
    ##   RpsG = col_double(),
    ##   DUF3172 = col_double(),
    ##   AtpD = col_double(),
    ##   AtpA = col_double(),
    ##   PetB = col_double(),
    ##   LepA = col_double(),
    ##   CbbP = col_double()
    ## )

``` r
sar11 <- read_tsv(here::here("input", "core_fams_sar11.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   sampleID = col_character(),
    ##   trpS = col_double(),
    ##   DUF3470 = col_double(),
    ##   ACP = col_double(),
    ##   RpsL = col_double(),
    ##   NQO = col_double(),
    ##   CarD = col_double(),
    ##   UBQ = col_double(),
    ##   TPR_2 = col_double(),
    ##   RelA = col_double(),
    ##   CJE = col_double()
    ## )

# FORMAT CORE AND PEPM RPKM

## SAR11

TPR\_2 (sar11) is an outlier and always severely undercounted in
GEOTRACES - maybe not core? To stringent of cutoff? It’s not such a big
problem with Tara because I think the sequencing depth is so much
greater than the GEOTRACES thus SAR11 coverage is generally much higher.
Excluding it here to be safe

``` r
sar11cmb <- sar11 %>%
  pivot_longer(trpS:CJE) %>%
  filter(name != "TPR_2") %>% # exclude TPR_2
  group_by(sampleID) %>%
  summarize(coreRPKM=median(value)) %>% 
  left_join(., dplyr::select(pepm, sampleID, pepmRPKM=sar11)) %>%
  filter(coreRPKM > 75) %>%
  mutate(tax="sar11")
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## Joining, by = "sampleID"

## PROCHLOROCOCCUS

``` r
procmb <- pro %>%
  pivot_longer(RpsJ:CbbP) %>%
  group_by(sampleID) %>%
  summarize(coreRPKM=median(value)) %>% 
  left_join(., dplyr::select(pepm, sampleID, pepmRPKM=prochlorococcus)) %>%
  filter(coreRPKM > 75) %>%
  mutate(tax="prochlorococcus")
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## Joining, by = "sampleID"

## BACTERIOPLANKTON

``` r
bactarccmb <- bactarc %>%
  pivot_longer(COG0012:COG0552) %>%
  group_by(sampleID) %>%
  summarize(coreRPKM=median(value)) %>% 
  left_join(., dplyr::select(pepm, sampleID, pepmRPKM=bactarc)) %>%
  filter(coreRPKM > 75) %>%
  mutate(tax="bactarc")
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## Joining, by = "sampleID"

Time series metagenome samples (HOT/BATS)

``` r
ts <- c("S0582", "S0583", "S0584", "S0585", "S0586", "S0587", "S0567", 
"S0568", "S0569", "S0606", "S0607", "S0608", "S0588", "S0589", 
"S0590", "S0591", "S0592", "S0593", "S0570", "S0571", "S0572", 
"S0609", "S0610", "S0611", "S0612", "S0613", "S0614", "S0573", 
"S0574", "S0575", "S0624", "S0625", "S0626", "S0576", "S0577", 
"S0578", "S0594", "S0595", "S0596", "S0615", "S0616", "S0617", 
"S0618", "S0619", "S0620", "S0621", "S0622", "S0623", "S0597", 
"S0598", "S0599", "S0600", "S0601", "S0602", "S0603", "S0604", 
"S0605", "S0579", "S0580", "S0581", "S0629", "S0630", "S0528", 
"S0529", "S0530", "S0531", "S0532", "S0533", "S0534", "S0535", 
"S0536", "S0537", "S0538", "S0539", "S0540", "S0541", "S0542", 
"S0543", "S0544", "S0545", "S0546", "S0547", "S0548", "S0549", 
"S0550", "S0551", "S0175", "S0174", "S0173", "S0558", "S0559", 
"S0560", "S0552", "S0553", "S0554", "S0555", "S0556", "S0557", 
"S0561", "S0562", "S0563", "S0501", "S0502", "S0503", "S0564", 
"S0565", "S0566", "S0504", "S0505", "S0506", "S0507", "S0508", 
"S0509", "S0510", "S0511", "S0512", "S0513", "S0514", "S0515", 
"S0516", "S0517", "S0518", "S0519", "S0520", "S0521", "S0522", 
"S0523", "S0524", "S0525", "S0526", "S0527", "S0627", "S0628"
)
```

Separate out timeseries from global

``` r
pepm.global <- bind_rows(procmb, sar11cmb, bactarccmb) %>%
  filter(!(sampleID %in% ts))

pepm.timeseries <- bind_rows(procmb, sar11cmb, bactarccmb) %>%
  filter(sampleID %in% ts)
```

# ADD ENVIRONMENTAL DATA

So the GEOTRACES data license is not GPL3 and you must request access
and be approved before being allowed to download. So I am not including
the GEOTRACES data here since it would be in violation of their license
agreement. Unfortunately, this means that this is not truly reproducible
as is. You can fill in the missing data for yourself by following the
code available here - <https://doi.org/10.5281/zenodo.3689249>.

``` r
env.global <- readRDS(here::here("input-old", "geotraces_tara_combined_imputed_new"))

as.na <- function(x) {x[] <- NA; x}

env.global.tara <- env.global %>%
  filter(str_detect(sampleID, "tara_")) %>%
  mutate(across(all_of(dput(colnames(env.global.tara)[61:98])), as.na))

env.global.simons <- env.global %>%
  filter(!(str_detect(sampleID, "tara_"))) %>%
  mutate(across(all_of(dput(colnames(env.global.tara)[61:98])), as.na)) %>%
  mutate(across(all_of(dput(colnames(env.global.tara)[6:13])), as.na))

env.global <- bind_rows(env.global.simons, env.global.tara)

env.global1 <- env.global[colSums(!is.na(env.global)) > 0]

write_tsv(env.global1, here::here("input", "geotraces_tara_env_data.tsv"))
```

## READ

``` r
env.global <- read_tsv(here::here("input", "geotraces_tara_env_data.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   sampleID = col_character(),
    ##   section = col_character(),
    ##   replicate = col_character(),
    ##   LongCode = col_character(),
    ##   LongWind = col_character(),
    ##   LongDesc = col_character(),
    ##   ocean = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
env.timeseries <- read_tsv(here::here("input", "timeseries_env_data.tsv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   sampleID = col_character(),
    ##   sample_category = col_character(),
    ##   section = col_character(),
    ##   date = col_date(format = ""),
    ##   ocean = col_character(),
    ##   environment = col_character(),
    ##   LongCode = col_character(),
    ##   LongWind = col_character(),
    ##   LongDesc = col_character()
    ## )
    ## See spec(...) for full column specifications.

## FORMAT

Format a global dataset. The W/M notation is for beta binomial
regresion. The RA (relative abundance) is just for beta distribution. We
will use the beta distribution with a random forest in the Ranger
package.

``` r
pepm.global <- bind_rows(procmb, sar11cmb, bactarccmb) %>%
  select(sampleID, tax, W=pepmRPKM, M=coreRPKM) %>%
  mutate(RA=W/M) %>%
  right_join(., env.global) %>%
  drop_na(tax) %>%
  filter(depth <= 250)
```

    ## Joining, by = "sampleID"

``` r
pepm.ts <- bind_rows(procmb, sar11cmb, bactarccmb) %>%
  select(sampleID, tax, W=pepmRPKM, M=coreRPKM) %>%
  mutate(RA=W/M) %>%
  right_join(., env.timeseries) %>%
  drop_na(tax) %>%
  filter(depth <= 250)
```

    ## Joining, by = "sampleID"

Make a sumarized version of the timeseries at HOT and BATS.

``` r
pepm.ts.summary <- pepm.ts %>%
  group_by(tax, section, environment) %>%
  summarize_at(.vars=c("W","M", "RA"),
               .funs=median) %>%
  ungroup() %>%
  mutate(sampleID=case_when(section=="BATS" & environment == "surface" ~ "S0582",
                            section=="BATS" & environment == "dcm" ~ "S0583",
                            section=="BATS" & environment == "mixed.layer" ~ "S0584",
                            section=="HOT" & environment == "surface" ~ "S0528",
                            section=="HOT" & environment == "dcm" ~ "S0529",
                            section=="HOT" & environment == "mixed.layer" ~ "S0530")) %>%
  mutate(dcm.layer=ifelse(environment=="dcm", 1, 0),
         dcm.max=ifelse(environment=="dcm", 1, 0)) %>%
  select(-section, -environment) %>%
  left_join(., env.timeseries) %>%
  select(-cruise_id, -environment, -sample_category)
```

    ## Joining, by = "sampleID"

Make the final combined dataset. Here we only select the variables we
want to use to model TBDT relative abundance, exclude the southern ocean
samples and any sample that is not in the euphotic zone (here deeper
than 250 meters), remove any samples that have NA values, and add 1e-7
pseudovalue to RA zero values to be compatible with beta distribution
split rule in `ranger::ranger`.

For macronutrients and nitrite replace 0 values with with random value
centered on the lowest nonzero value. This reduces weird zero-inflation
issues.

Also exclude very high relative abundances \> 0.5. These all happen
pretty much in Prochlorococcus when the marker gene RPKM goes way low
and we still have pepm RPKM (that may or may not be noise)

``` r
pepm.global.final <- bind_rows(pepm.ts.summary, pepm.global) %>%
  filter(depth < 250) %>%
  filter(ocean != "southern_ocean") %>%
  mutate(ocean=factor(ocean, levels=c("atlantic_north", "atlantic_south", "pacific_north",
                                      "pacific_south", "med_red", "indian"))) %>%
  mutate(dcm.layer=as.factor(ceiling(dcm.layer)),
         dcm.max=as.factor(ceiling(dcm.max))) %>%
  mutate(nitrate_dissolved_umol.kg=ifelse(nitrate_dissolved_umol.kg==0, 
                                          abs(rnorm(n(), mean=0.001, sd=0.0005)), 
                                          nitrate_dissolved_umol.kg),
         nitrite_dissolved_umol.kg=ifelse(nitrite_dissolved_umol.kg==0, 
                                          abs(rnorm(n(), mean=0.0005, sd=0.00008)),
                                          nitrite_dissolved_umol.kg),
         phosphate_dissolved_umol.kg=ifelse(phosphate_dissolved_umol.kg==0, 
                                            abs(rnorm(n(), mean=0.001, sd=0.0005)),
                                            phosphate_dissolved_umol.kg),
         silicate_dissolved_umol.kg=ifelse(silicate_dissolved_umol.kg==0, 
                                           abs(rnorm(n(), mean=0.05, sd=0.005)), silicate_dissolved_umol.kg),
         chla_norm=ifelse(chla_norm==0, rnorm(n(), mean=0.08, sd=0.001), chla_norm)) %>%
  mutate(RA=ifelse(RA>0.5, abs(rnorm(n(), mean=0.5, sd=0.05)), RA))
```

## SAVE

``` r
saveRDS(pepm.global.final, file=here::here("output", "pepm.global.final"))
saveRDS(pepm.ts, file=here::here("output", "pepm.timeseries.final"))
```
