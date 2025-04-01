# Supplemental-ERPV-Hotspot


This repository contains the supplemental materials and BRDF model fitting code for the paper:

**"An enhanced RPV model to better capture hotspot signatures in vegetation canopy reflectance observed by the geostationary meteorological satellite Himawari-8."**

---

##  Overview

This work focuses on improving the fitting of vegetation canopy reflectance using an enhanced RPV model with hotspot sensitivity. The repository includes:

📁 A curated dataset of hotspot observations, including BRFs and fitted results from multiple BRDF models

📊  Fitting code for various BRDF models described in the study

---

##  CSV Column Descriptions

Each `.csv` file corresponds to a specific location and observation period. The following fields are included:

| Column Name         | Description |
|---------------------|-------------|
| `TIME`              | Observation timestamp (Local) |
| `RED`               | Red band reflectance|
| `NIR`               | Near-infrared band reflectance |
| `SAZ`               | Satellite zenith angle in degrees |
| `SOZ`               | Solar zenith angle in degrees |
| `SAA`               | Satellite azimuth angle in degrees |
| `SOA`               | Solar azimuth angle in degrees |
| `SA`                | Scattering angle in degrees|
| `ERPV_RED_FIT`      | Fitted RED reflectance using **Enhanced RPV model** |
| `ERPV_NIR_FIT`      | Fitted NIR reflectance using **Enhanced RPV model** |
| `RPV_RED_FIT`       | Fitted RED reflectance using **Original RPV model** |
| `RPV_NIR_FIT`       | Fitted NIR reflectance using **Original RPV model** |
| `JIAO2016_RED_FIT`  | Fitted RED reflectance using **Jiao2016 model** |
| `JIAO2016_NIR_FIT`  | Fitted NIR reflectance using **Jiao2016 model** |
| `RTLSR_RED_FIT`     | Fitted RED reflectance using **RTLSR model** |
| `RTLSR_NIR_FIT`     | Fitted NIR reflectance using **RTLSR model** |
| `MAIGNAN2004_RED_FIT` | Fitted RED reflectance using **Maignan2004 model)** |
| `MAIGNAN2004_NIR_FIT` | Fitted NIR reflectance using **Maignan2004 model)** |
| `NDVI`              | Normalized Difference Vegetation Index computed as (NIR - RED) / (NIR + RED) |

---
📌 Note:
The name of each .csv file contains only the site name and vegetation type (e.g., 0,0_AUC_EBF_2018_0329.csv).
To obtain detailed information such as geographic coordinates and the full site name, please refer to the supplementary file siteinfo.csv.


##  BRDF Model Descriptions

This repository provides code implementations for multiple BRDF models used in the study, including:

Enhanced RPV model (ERPV): An Enhanced RPV model developed in this study to better capture hotspot effects.

Other commonly used BRDF models: Original RPV, RTLSR, Jiao2016, and Maignan2004.

📚 Citation Notice

If you use the ERPV model or the dataset in this repository, please cite our paper:
Yang Wei, Qiao Zhi, Li Wei, Ma Xuanlong, and Ichii Kazuhito. (2025).
An Enhanced RPV Model to Better Capture Hotspot Signatures in Vegetation Canopy Reflectance Observed by the Geostationary Meteorological Satellite Himawari-8.
Science of Remote Sensing.

If you use or extend other BRDF models (e.g., Original RPV, RTLSR, Jiao2016, Maignan2004), please cite the original publications. Brief descriptions and corresponding citations for these models are provided in our paper.

---

## ⚙️ How to Use

⚙️ How to Use
1. Prepare the data
Place your input .csv files in a local directory.
Each file should be named according to site and observation metadata (e.g., 0,0_AUC_EBF_2018_0329.csv).

Each file must contain at least the following columns:
TIME, SOZ, SAZ, SOA, SAA, RED, NIR
where:

SOZ, SAZ, SOA, and SAA represent solar/satellite zenith and azimuth angles.

RED and NIR represent red and near-infrared reflectance.

2. Place the model scripts
Ensure the following custom Python modules are placed in the same directory as the main script:

Func_LSR_Kgeo.py

Func_RossThick_Kvol.py

Func_Maignan2004_Kvol.py

These contain kernel functions required by the RTLSR, Jiao2016, and Maignan2004 models.

3. Edit the input file path
In the main script (e.g., main_fitting_script.py), set the path to the .csv file you wish to process:


file_path = r'F:/data/your_csv_file.csv'
Note: Use raw string format (r'...') or double backslashes for Windows paths.

4. Run the script
Execute the script in your Python environment (e.g., Spyder, VSCode, or terminal):

5. Output
The script prints fitting metrics (R², RMSE, BIAS) and model parameters for each band and model.

Results are stored in a dictionary named fitting_results.

You may extend the code to export fitted values or plots.

📜 License

This project is released under the MIT license. See the LICENSE file for details.

🙌 Acknowledgements

Created by Qiao Zhi.
If you have any questions, feedback, or suggestions, feel free to open an Issue on this repository.
We appreciate your interest and welcome any contributions that support open and reproducible research!
