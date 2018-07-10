#!/usr/bin/env Rscript
#
# Add BioClim data to the strain information produced by Liana.
#
# Closely based on the script Liana provided.
#

library(dplyr)
library(raster)
library(rgdal)

bioclim_dir = '../../bioclim/30sec'
bil_files = list.files(bioclim_dir, pattern='\\.bil$')
file_stack = stack(file.path(bioclim_dir, bil_files))

dat = read.csv('SeqStrainInfo_ltb_27July2016.txt', sep='\t',
              header=TRUE, as.is=TRUE)

dat_for_coords = dat %>%
    filter(!is.na(Country_Lat)) %>%
    mutate(lat=ifelse(is.na(Detail_Lat), Country_Lat, Detail_Lat)) %>%
    mutate(lon=ifelse(is.na(Detail_Long), Country_Long, Detail_Long))

bioclims = dat_for_coords %>%
    (function(x) {coordinates(x) = c('lon', 'lat'); x}) %>%
    extract(file_stack, ., method='simple') %>%
    as.data.frame() %>%
    mutate('GWAS_ID'=dat_for_coords[, 'GWAS_ID']) %>%
    left_join(dat, .)

write.table(bioclims, file='data_bioclim.tsv', sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
