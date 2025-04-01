# -*- coding: utf-8 -*-
"""
Main script for fitting BRDF models to directional reflectance data.

This script reads observation data (angles + reflectance) from a CSV file,
fits multiple BRDF models to both RED and NIR bands, and prints the results.


@author: QIAO ZHI
"""



import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression
import Func_LSR_Kgeo
import Func_RossThick_Kvol
import Func_Maignan2004_Kvol

# Define Improved RPV model
def Improved_RPV_func(xdata, rho, kei, theta, C1, C2):
    sza, vza, saa, vaa = np.radians(xdata[:, 0]), np.radians(xdata[:, 1]), np.radians(xdata[:, 2]), np.radians(xdata[:, 3])
    raa = vaa - saa
    G = np.sqrt((np.tan(sza))**2 + (np.tan(vza))**2 - 2 * np.tan(sza) * np.tan(vza) * np.cos(raa))
    cos_g = np.cos(vza) * np.cos(sza) + np.sin(vza) * np.sin(sza) * np.cos(raa)
    g = np.degrees(np.arccos(cos_g))
    H = 1 + (1 - rho) / (1 + G) * (1 + C1 * np.exp(-g / C2))
    F_HG = (1 - theta**2) / (1 + 2 * theta * cos_g + theta**2)**1.5
    var1 = np.cos(sza)**(kei - 1)
    var2 = np.cos(vza)**(kei - 1)
    var3 = (np.cos(sza) + np.cos(vza))**(kei - 1)
    M_I = var1 * var2 * var3
    rho_sfc = rho * M_I * F_HG * H
    return rho_sfc

# Define Original RPV model
def Original_RPV_func(xdata, rho, kei, theta):
    sza, vza, saa, vaa = np.radians(xdata[:, 0]), np.radians(xdata[:, 1]), np.radians(xdata[:, 2]), np.radians(xdata[:, 3])
    raa = vaa - saa
    G = np.sqrt((np.tan(sza))**2 + (np.tan(vza))**2 - 2 * np.tan(sza) * np.tan(vza) * np.cos(raa))
    cos_g = np.cos(vza) * np.cos(sza) + np.sin(vza) * np.sin(sza) * np.cos(raa)
    H = 1 + (1 - rho) / (1 + G)
    F_HG = (1 - theta**2) / (1 + 2 * theta * cos_g + theta**2)**1.5
    var1 = np.cos(sza)**(kei - 1)
    var2 = np.cos(vza)**(kei - 1)
    var3 = (np.cos(sza) + np.cos(vza))**(kei - 1)
    M_I = var1 * var2 * var3
    rho_sfc = rho * M_I * F_HG * H
    return rho_sfc

# Define RTLSR model
def RTLSR_func(xdata, ydata):
    sza, vza, saa, vaa = xdata[:, 0], xdata[:, 1], xdata[:, 2], xdata[:, 3]
    K_geo = Func_LSR_Kgeo.LiSparseReciproal_Kgeo(sza, vza, saa, vaa)
    K_vol = Func_RossThick_Kvol.RossThick_Kvol(sza, vza, saa, vaa)
    X = np.column_stack([K_vol, K_geo])
    model = LinearRegression().fit(X, ydata)
    y_fit = model.predict(X)
    return model, y_fit

# Define Jiao2016 model
def Jiao2016_func(xdata, f0, f1, f2, C1, C2):
    sza, vza, saa, vaa = np.radians(xdata[:, 0]), np.radians(xdata[:, 1]), np.radians(xdata[:, 2]), np.radians(xdata[:, 3])
    raa = vaa - saa
    cos_kai = np.cos(sza) * np.cos(vza) + np.sin(sza) * np.sin(vza) * np.cos(raa)
    kai = np.arccos(cos_kai)
    numerator = (0.5 * np.pi - kai) * cos_kai + np.sin(kai)
    denominator = np.cos(sza) + np.cos(vza)
    Hotspot_factor = 1 + C1 * np.exp(-np.degrees(kai) / C2)
    K_vol = Hotspot_factor * numerator / denominator - 0.5 * np.pi
    K_geo = Func_LSR_Kgeo.LiSparseReciproal_Kgeo(sza, vza, saa, vaa)
    rho_sfc = f0 + f1 * K_geo + f2 * K_vol
    return rho_sfc

# Define Maignan2004 model
def Maignan2004_func(xdata, ydata):
    sza, vza, saa, vaa = xdata[:, 0], xdata[:, 1], xdata[:, 2], xdata[:, 3]
    K_geo = Func_LSR_Kgeo.LiSparseReciproal_Kgeo(sza, vza, saa, vaa)
    K_vol = Func_Maignan2004_Kvol.Maignan2004_Kvol(sza, vza, saa, vaa)
    X = np.column_stack([K_vol, K_geo])
    model = LinearRegression().fit(X, ydata)
    y_fit = model.predict(X)
    return model, y_fit

# === Set input CSV file path ===
file_path = r'sample/sample.csv'

# Read data
data_df = pd.read_csv(file_path)


# Extract necessary data
Solar_zenith = np.array(data_df['SOZ'])
View_zenith = np.array(data_df['SAZ'])
Solar_azimuth = np.array(data_df['SOA'])
View_azimuth = np.array(data_df['SAA'])
RED = np.array(data_df['RED'])   
NIR = np.array(data_df['NIR']) 

# Filter data

xdata = np.column_stack((Solar_zenith, View_zenith, Solar_azimuth, View_azimuth))

# Models and bands
models = {
    'Improved_RPV': Improved_RPV_func,
    'Original_RPV': Original_RPV_func,
    'RTLSR': RTLSR_func,
    'Jiao2016': Jiao2016_func,
    'Maignan2004': Maignan2004_func
}
bands = {'RED': RED, 'NIR': NIR}

# Fit models
fitting_results = {}
for band_name, ydata in bands.items():
    fitting_results[band_name] = {}
    for model_name, model_func in models.items():
        try:
            if model_name in ['Improved_RPV', 'Original_RPV', 'Jiao2016']:
                # Nonlinear curve fitting
                popt, _ = curve_fit(model_func, xdata, ydata, maxfev=5000)
                y_fit = model_func(xdata, *popt)
                r2 = r2_score(ydata, y_fit)
                bias = np.mean(y_fit - ydata)
                rmse = np.sqrt(mean_squared_error(ydata, y_fit))
                fitting_results[band_name][model_name] = {'r2': r2, 'bias': bias, 'rmse': rmse, 'parameters': popt, 'y_fit': y_fit}
                print(f"File: {file_path}, Band: {band_name}, Model: {model_name}")
                print(f"R²: {r2:.4f}, BIAS: {bias:.4f}, RMSE: {rmse:.4f}")
                print(f"Fitted Parameters: {popt}")
            elif model_name in ['RTLSR', 'Maignan2004']:
                # Linear fitting for RTLSR and Maignan2004
                model, y_fit = model_func(xdata, ydata)
                r2 = r2_score(ydata, y_fit)
                bias = np.mean(y_fit - ydata)
                rmse = np.sqrt(mean_squared_error(ydata, y_fit))
                fitting_results[band_name][model_name] = {'r2': r2, 'bias': bias, 'rmse': rmse, 'coefficients': model.coef_, 'y_fit': y_fit}
                print(f"File: {file_path}, Band: {band_name}, Model: {model_name}")
                print(f"R²: {r2:.4f}, BIAS: {bias:.4f}, RMSE: {rmse:.4f}")
                print(f"Coefficients: {model.coef_}")
        except Exception as e:
            print(f"Model {model_name} fitting failed for {band_name} in {file_path}: {e}")
