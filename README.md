### Physically Based Adaptive empirical local patches for CaMa-Flood hydrodynamic model
==========================================================
### Local patached area derived using spatial auto-correlation
==========================================================
1. s01-CaMa_sim.sh - run CaMa-Flood
2. s02-conv_bin2nc.sh - convert binary files to netCDF4
3. s03-remove_trend.sh - remove trend lines
4. s04-remove_season.sh - remove seasonality
5. s05-standardaize.sh - standardize the data
6. s06-semivarince.sh - calculate the experimental semivariances
7. s07-weightage.sh - calcualte the spatial dependency weightages
8. s08-localiz_para.sh - write localization parameters to text files
9. s09-localizMS_para.sh - localization parameters only for main stem
==========================================================
### Reference:
1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S.,A Physically Based Empirical Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829, 2019
2. Revel, M., Yamazaki, D., & Kanae, S.,Model Based Observation Localization Weighting Function for Amazon Mainstream, Journal of Japan Society of Civil Engineers, Ser. B1 (Hydraulic Engineering), 74(5), I_157-I_162,https://doi.org/10.2208/jscejhe.74.5_I_157,2018 