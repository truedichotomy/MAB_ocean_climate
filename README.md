# MAB_ocean_climate

Citation: 

Wallace, E. J., Looney, L. B., & Gong, D. (2018). Multi-decadal trends and variability in temperature and salinity in the Mid-Atlantic Bight, Georges Bank, and Gulf of Maine. Journal of Marine Research, 76(5-6), 163-215.

Abstract:

Increasing attention is being placed on the regional impact of climate change. This study focuses on the decadal scale variabilities of temperature and salinity in the Mid-Atlantic Bight (MAB), Georges Bank (GB), and Gulf of Maine (GOM) from 1977 to 2016 using hydrographic survey data from the National Oceanic and Atmospheric Administration (NOAA) Northeast Fisheries Science Center. The MAB (as defined by the shelf regions from Cape Hatteras to Cape Cod) experienced warming rates of 0.57 ◦C per decade during the Winter/Spring season (Jan–Apr) and 0.47 ◦C per decade during the Fall/Winter season (Sep–Dec). The GOM and GB, on the other hand, warmed at approximately half the rate of the MAB over the same time span (1977–2016). We found that rates of warming vary on decadal time scales. From 1977 to 1999, significant temperature increases (>0.6 ◦C/decade) were found in the southern regions of the MAB during the Winter/Spring season. During the same period, significant freshening (stronger than −0.2/decade) was found in GB and the northern regions of the MAB during the Winter/Spring and Summer seasons. From 1999 to 2016, on the other hand, we found no significant trends in temperature and few significant trends in salinity with the exceptions of some northern MAB regions showing significant salting. Interannual variability in shelf salinity can in part be attributed to river discharge variability in the Hudson River and Chesapeake Bay. However, decadal scale change in shelf salinity cannot be attributed to changes in river discharge as there were no significant decadal scale changes in river outflow. Variability in along-shelf freshwater transport and saline intrusions from offshore were the likely drivers of long-term changes in MAB shelf-salinity.

Data Analysis:

Software required:
- MATLAB
- MATLAB Mapping Toolbox
- Gibbs SeaWater (GSW) Oceanographic Toolbox (http://www.teos-10.org/software.htm)

Aanalysis of CTD data:
Run 'dg_grid_regions2D.m' first to load CTD data & calculate grids.

Analysis of river discharge data:
Run 'dg_load_hudson_chesapeake_discharge.m' first to load river discharge data.
