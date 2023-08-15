# Montgomery2023_urban_no2_hotspot

Repository for "Intraurban NO2 Hotspot Detection across Multiple Air Quality Products"

Includes:

  preprocessing/
    
    s1-clean-eclipse-append-cmaq.ipynb    # formatting eclipse to grid, eclipse data is outdated but need this first step
    
    s2-format_input_data_on_cmaq.ipynb    # add demographics & land use to grid
    
    s3-eclipse_interpolated_daily.r       # create daily interpolated maps from Eclipse
    
    s4-download_TropOMI_directly.ipynb    # download TropOMI data locally
    
    s5-regrid_tropomi_simple.py           # regrid TropOMI to cmaq grid
  
  calculations/
    
    cluster.py                            # create clusters, check significant attributes in clusters
    
    check.py                              # compare products against eachother and EPA
    
    check_distance_and_height.py          # analyze highway hotspots
  
  data/
    
    eclipse_daily_qaqc.csv                # hourly eclipse data
    
    monthly.json                          # monthly shapefile with gridded eclipse, tropomi, & wrf-cmaq

Other data (not generated from this work):

demographic data: ACS from https://www.nhgis.org/

land-use and vehicle data: streets, industrial zoning, and bus speeds from https://data.cityofchicago.org/

wrf-cmaq model: available from EPA/CMAS
