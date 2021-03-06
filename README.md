General Notes 
==============

This document summarizes the paleo-lake model created by Julian Cross. This document was compiled by Bryce Glenn and Julian Cross. There are three separate components to fully run the model. For specific details on how to run each model see the separate README files created by Julian Cross. The model relies on three different scenarios for the Ross Ice Sheet (RIS) details are described below.

Model Components
==============
1.  ICEMELT energy and mass balance model

    -   icemelt\_cross\_v05.f95

    -   FORTRAN code

    -   This is run on the PSU cluster

2.  Meltwater model

    -   A MATLAB script

3.  Lake-level model

    -   A MATLAB script

Terrain and Domain Definition Grids
==============

There are three scenarios for this model, no RIS (modern topography), minimum incursion RIS, and maximum incursion RIS. Each scenario has its own set of terrain and domain definition grids (ASCII file format) that are read by the ICEMELT model.

-   Domain definition grids indicates which cells to run in the ICEMELT model. Glaciers in the model domain (Taylor Valley) given a unique drainage basin code (10 to 90), stored in this file.

-   Mean annual temperature grids are used to shift initial ice column temperatures to a local value.

-   Terrain grids are 250-m DEMs describing the topography of Taylor Valley under each scenario.

**Table 1.** Differences between terrain and domain definition grids modern model versus the paleo model.

|                              |             Modern Model            |                          Paleo Model                          |
|------------------------------|:-----------------------------------:|:-------------------------------------------------------------:|
|       Domain Definition Grid | `tv\_basins\_surface\_jmc\_ekh.txt` | `tv\_basins\_surface\_min.txt` `tv\_basins\_surface\_max.txt` |
|      Annual Temperature Grid | `T\_avg\_all.txt`                   | `T\_avg\_min.txt` `T\_avg\_max.txt`                           |
|                 Terrain Grid | `tv\_dem250.txt`                    | `tv\_mindem250.txt` `tv\_maxdem250.txt`                       |

For the minimum intrusion RIS scenario, the extent of the RIS lobe in Taylor Valley was digitized as a polygon feature in ArcGIS based on Hall et al. 2000. This extent polygon feature was converted to points (10 m interval?). Elevation values from the 30-m Taylor Valley DEM were extracted to these points. Other points were added within the RIS lobe extent feature and given elevation values, corresponding to the 100, 200, 400, 750 and 1000 meter contours.

For the maximum intrusion RIS scenario, the extent of the RIS lobe in Taylor Valley was mapped in ArcGIS, reaching as far west as defile between the Nussbaum Riegel and Suess Glacier. The outside extent followed valley-wall contours from an initial height at the Suess Glacier of 300 meters ?? to 1500 meters at the mouth of the valley. Modern elevations were extracted to points along this extent polygon feature following the same method described above. Additional elevation points were added up-glacier corresponding to the 400, 750 and 1000 meter contours.

The elevation points described above for each scenario were interpolated to a 30-m grid using the ?? interpolation method.

The new RIS DEMs were merged with the exiting 250-m DEM...

For paleo- scenarios, cells are added to the modern model domain definition grid based on the mapped RIS extents described above. These new zones of RIS ice are given the drainage basin code of 90.

Using the mapped RIS extents described above, mean annual temperatures were calculated for each new RIS cell in MATLAB using the 13 years of MICROMET data.

Climate Variables 
==============

The ICEMELT model relies on several climate variables. Meteorological forcing data are stored in binary format MICROMET input and albedo data are stored in ASCII format input files.

**Table 2.** Differences between MICROMET and albedo for the modern model versus the paleo model.

|          	|                                                                                                     Modern Model                                                                                                    	|                                                                                               Paleo Model                                                                                               	|
|:--------:	|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:	|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:	|
| MICROMET 	| - Modified temperature and vapor pressure lapse rates are used in the Fryxell Basin (basin codes 50-82).                                                                                                            	| - Spatially uniform temperature and vapor pressure lapse rates.<br>- Meteorological conditions based on modified terrain grids (e.g. `tv\_mindem250.txt`).                                              	|
| Albedo   	| - A combination of MODIS and in-situ measured albedo values from Taylor, Canada, and Commonwealth glaciers are used as representative of albedo for glaciers in Bonney, Hoare and Fryxell watersheds, respectively. 	| - Initially Commonwealth Glacier albedo is used for the Minimum RIS.<br>- Subsequent adjustments to albedo are made using a fixed value or offset/multiplier applied to the temporally variable albedo. 	|

Hypsometry and Watershed Boundary Delineation
==============

Lake Model Time Step and Modeled Period
==============

The modern model is run with a daily time-step, while the paleo- model uses an annual time-step. Test runs over the 16 historical years showed that with daily and annual versions of the model, water balances converged at the end of each water year. The paleo-model is setup to use repeating sequences of the 16 historical years across the selected modeled period (e.g. a model period of 10,000 years would the same 16-year sequence repeated 625 times).
