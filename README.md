# lst_gap_fill
This module is part of an unsuccessful attempt of filling the gaps in MODIS LST using ERA5 Skin temperature.
The project was discarded because the results were not as good as anticipated.

Description of the files

1.LST_ERA_data_match
  Reads ERA5 skin temperature data and MODIS LST data
  Then selects the ERA5 data closest to MODIS data
  Builds netcdf file for LST as well as chosen ERA5 data
  
 2. lst_skt_fill_module
    It consists of functions used for gap filling
  
 3. LST_ERA_gapfill_20210228
    It contains the actual process of gapfilling used
