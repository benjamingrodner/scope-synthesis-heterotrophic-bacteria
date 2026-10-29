# Framework generated with Gemini 3 free tier 28 Oct 2026
# prompt:
# Please write a python script to plot a color map of a region of the globe defined by a window of latitude and longitude
# Edited by Ben Grodner

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def load_copernicus_netcdf(file_path):
    """
    Loads a NetCDF4 file 

    Returns:
    - xarray.Dataset or xarray.DataArray subsetted to the region
    """
    print(f"Opening NetCDF file: {file_path}")

    # Open dataset using xarray (requires engine='netcdf4')
    ds = xr.open_dataset(file_path, engine="netcdf4")

    return ds

def average_over_time(ds, val_name):
    return ds[val_name].mean(axis=0)

def subset_latlon(df, step):
    return df.isel(
        {
            "latitude": slice(None, None, step),
            "longitude": slice(None, None, step),
        }
    )

# def sub_df_metat(df, cruise, date_re):
#     bool_dp = (df["depth"] > depths[0]) & (df["depth"] < depths[1])
#     bool_dt = df['Gradients'] == cruise
#     return df.loc[bool_dp & bool_dt, ['lat','lon']]

def prep_dfs_metat(fn, dict_cruise_xyshift):
    dfs = {}
    df = pd.read_csv(fn)
    bool_dp = df["depth"] < 200
    bool_ns = ['ns' in s for s in df['SampleID']]
    # dfs.append(df[bool_ns & bool_dp])
    for cruise, shft in dict_cruise_xyshift.items():
        bool_c = df['Gradients'] == cruise
        df_sub = df[bool_dp & bool_ns & bool_c].copy()
        df_sub["lon_old"] = df_sub["lon"]
        df_sub["lat_old"] = df_sub["lat"]
        for ls in shft:
            bool_lgt = (df_sub['lon_old'] > ls['lonlims'][0]) 
            bool_llt = (df_sub['lon_old'] < ls['lonlims'][1])
            df_sub.loc[bool_lgt & bool_llt, 'lon'] = (
                df_sub.loc[bool_lgt & bool_llt, 'lon_old'] + ls['shift'][0]
            )
            df_sub.loc[bool_lgt & bool_llt, "lat"] = (
                df_sub.loc[bool_lgt & bool_llt, "lat_old"] + ls["shift"][1]
            )
        dfs['metaT_' + cruise] = df_sub
    return dfs

def prep_dfs_asv_g123(fmt, dict_cruise_xyshift):
    dfs = {}
    for crs, lonshift in dict_cruise_xyshift.items():
        fn = glob.glob(fmt.format(cruise=crs))[0]
        df = pd.read_csv(fn)
        df = df.loc[df['depth'] < 200, ['lat','lon','sampleID']].drop_duplicates()
        df['lon_old'] = df['lon']
        df['lon'] = df['lon_old'] + lonshift
        dfs[f'ASVs_{crs}'] = df
    return dfs


def prep_dfs_asv_5(fn, dict_cruise_xyshift):
    dfs = {}
    df = pd.read_csv(fn)
    # dfs.append(df[bool_ns & bool_dp])
    for cruise, shft in dict_cruise_xyshift.items():
        df_sub = df.copy()
        df_sub["lon_old"] = df_sub["lon"]
        df_sub["lat_old"] = df_sub["lat"]
        for ls in shft:
            bool_lgt = df_sub["lon_old"] > ls["lonlims"][0]
            bool_llt = df_sub["lon_old"] < ls["lonlims"][1]
            df_sub.loc[bool_lgt & bool_llt, "lon"] = (
                df_sub.loc[bool_lgt & bool_llt, "lon_old"] + ls["shift"][0]
            )
            df_sub.loc[bool_lgt & bool_llt, "lat"] = (
                df_sub.loc[bool_lgt & bool_llt, "lat_old"] + ls["shift"][1]
            )
        dfs["metaG_" + cruise] = df_sub
    return dfs


def prep_dfs_metag(FN_METAG, dict_cruise_xyshift):
    dfs = {}
    df = pd.read_csv(FN_METAG)
    for crs, lonshift in dict_cruise_xyshift.items():
        df_sub = df.loc[df["cruise"] == crs]
        df_sub["lon_old"] = df_sub["lon"]
        df_sub["lon"] = df_sub["lon_old"] + lonshift
        dfs[f"metaG_{crs}"] = df_sub
    return dfs

def prep_dfs_influx(fmt, dict_cruise_xyshift):
    dfs = {}
    for crs, shft in dict_cruise_xyshift.items():
        fn = glob.glob(fmt.format(cruise=crs))[0]
        df = pd.read_csv(fn)
        if 'count_bacteria' in df.columns:
            nacol = 'count_bacteria'
            st = 'file'
        elif 'cell_abundance_bacteria':
            nacol = "cell_abundance_bacteria"
            st = 'station'
        df_sub = df.loc[(df["depth"] < 200) & (~df[nacol].isna()), ["lat", "lon", 'depth', st]].drop_duplicates()
        df_sub["lon_old"] = df_sub["lon"]
        df_sub["lat_old"] = df_sub["lat"]
        for ls in shft:
            bool_lgt = df_sub["lon_old"] > ls["lonlims"][0]
            bool_llt = df_sub["lon_old"] < ls["lonlims"][1]
            df_sub.loc[bool_lgt & bool_llt, "lon"] = (
                df_sub.loc[bool_lgt & bool_llt, "lon_old"] + ls["shift"][0]
            )
            df_sub.loc[bool_lgt & bool_llt, "lat"] = (
                df_sub.loc[bool_lgt & bool_llt, "lat_old"] + ls["shift"][1]
            )
        dfs[f"influx_{crs}"] = df_sub
    return dfs


def plot_geographic_colormap(lat_min, lat_max, lon_min, lon_max, 
                             color_grid=None, modality_locs={}, 
                             modality_params={}):
    """
    Plots a color map over a specific latitude/longitude window of the globe.

    Parameters:
    - lat_min, lat_max: Latitude bounds in degrees (-90 to 90)
    - lon_min, lon_max: Longitude bounds in degrees (-180 to 180)
    - resolution: Number of grid points along each axis for the data grid
    """
    if color_grid is None:
        # 1. Create a spatial grid for latitude and longitude
        resolution = 100
        lats = np.linspace(lat_min, lat_max, resolution)
        lons = np.linspace(lon_min, lon_max, resolution)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # 2. Generate sample data (replace this with your actual data grid, e.g., temperature, elevation, etc.)
        # Example function generating a 2D wave surface over the region
        data = np.sin(np.radians(lat_grid) * 4) + np.cos(np.radians(lon_grid) * 4)
    else:
        lon_grid, lat_grid = np.meshgrid(color_grid["longitude"], color_grid["latitude"])
        data = color_grid.values

    # 3. Initialize matplotlib figure with a PlateCarree map projection
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    fig_cbar, ax_cbar = plt.subplots(figsize=(1.5, 0.15))  # Tall and narrow layout    

    # 4. Set map extent [longitude_min, longitude_max, latitude_min, latitude_max]
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # 5. Add geographic features for context
    ax.add_feature(cfeature.LAND, facecolor='lightgrey', alpha=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor='black')
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.8)

    # Add gridlines with latitude/longitude labels
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(),
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 6, "color": "black"}
    gl.ylabel_style = {"size": 6, "color": "black"}
    print('Plotting...')
    # 6. Plot the 2D color map (pcolormesh or contourf)
    mesh = ax.pcolormesh(
        lon_grid, 
        lat_grid, 
        data, 
        transform=ccrs.PlateCarree(),
        cmap='Greens',      # Change colormap as needed (e.g., 'coolwarm', 'terrain', 'plasma')
        shading='auto',
        alpha=0.75,
        vmin=0,
        vmax=0.35
    )

    # 7. Add a colorbar
    cbar = fig_cbar.colorbar(mesh, cax=ax_cbar, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label('Chlorophyll (mg/m^3)', fontsize=7)
    
    # Plot cruise track
    ax.plot([-158, -158], [22.5, 43], "k")
    ax.plot([-155.5, -140], [18.5, -0.5], "k")
    ax.plot([-140, -140], [-4, 18.5], "k")
    ax.plot([-140, -119], [18.5, 31.5], "k")

    # Plot modality locs
    for m, locs in modality_locs.items():
        col = modality_params[m]['color']
        ax.scatter(locs['lon'],locs['lat'], 
                   label=m, alpha=0.3, color=col,
                   s=20
                   )

    ax.legend(fontsize=6)
    # ax.tick_params(axis='both', which='major', labelsize=6)
    # plt.title(f'Color Map ({lat_min}°N to {lat_max}°N, {lon_min}°E to {lon_max}°E)', fontsize=14, pad=15)
    plt.tight_layout()

    # Show plot (or save using plt.savefig('map.png', dpi=300))
    # plt.show()
    # plt.savefig("map.png", dpi=300)
    return fig, ax, fig_cbar


if __name__ == '__main__':
    # Define your bounding box window
    DIR_OUT = 'results/summary_figure'

    LAT_MIN, LAT_MAX = -10.0, 45.0    # e.g., Mediterranean / Southern Europe region
    LON_MIN, LON_MAX = -180.0, -108.0

    NETCDF_FILE = "data/copernicus/180-110W_10S-45N_2026_04_20-27_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_1793256551732.nc"
    STEP = 3
    DPI = 300

    FN_METAT = 'data/cmap/metaT/tblNPac_metat_samples_data.csv'
    FMT_ASV_123 = "data/cmap/ASVs/{cruise}_ASVs/*.csv"
    FN_METAG = "data/cmap/metaG/metaG_latlon.csv"
    FN_ASV_5 = "data/cmap/ASVs/G5_ASVs/G5_ASVs_latlon.csv"
    FMT_INFLUX = "data/cmap/influx/CMAP-Bulk-Download/Influx_Stations_Gradients_{cruise}*/*.csv"

    ds = load_copernicus_netcdf(NETCDF_FILE)
    df_latlon = average_over_time(ds, val_name='CHL')
    df_latlon_sub = subset_latlon(df_latlon, step=STEP) if STEP > 1 else df_latlon

    # MetaT
    s1 = 3
    s2 = 10.5
    s3 = 18
    s4 = 3
    s5 = 7
    s6 = 9
    dict_cruise_xyshift = {
        "Gradients1": [{"lonlims": [-160, -150], "shift": [-s3, 0]}],
        "Gradients2": [
            {
                "lonlims": [
                    -160,
                    -150,
                ],
                "shift": [-s2, 0],
            }
        ],
        "Gradients3": [{"lonlims": [-160, -150], "shift": [-s1, 0]}],
        "Gradients4": [
            {
                "lonlims": [-160, -141],
                "shift": [-np.sin([np.pi / 3.8]) * s5, -np.cos([np.pi / 3.8]) * s5],
            },
            {"lonlims": [-141.5, -139.1], "shift": [s4, 0]},
            {
                "lonlims": [-139, -110],
                "shift": [-np.cos([np.pi / 3.5]) * s6, np.sin([np.pi / 3.5]) * s6],
            },
        ],
        "Gradients5": [
            {
                "lonlims": [-160, -141],
                "shift": [
                    -np.sin([np.pi / 3.8]) * s4,
                    -np.cos([np.pi / 3.8]) * s4,
                ],
            },
            {"lonlims": [-141.5, -139.1], "shift": [s6, 0]},
            {
                "lonlims": [-139, -110],
                "shift": [-np.cos([np.pi / 3.5]) * s4, np.sin([np.pi / 3.5]) * s4],
            },
        ],
    }
    dict_metat_locs = prep_dfs_metat(FN_METAT, dict_cruise_xyshift)

    # ASVs 123
    s1 = 6
    s2 = 13.5
    s3 = 21
    dict_cruise_xyshift = {
        "G1": -s3,
        "G2": -s2,
        "G3": -s1,
    }
    dict_asv123_locs = prep_dfs_asv_g123(FMT_ASV_123, dict_cruise_xyshift)

    # ASV 5
    s1 = 4.5
    s2 = 7.5
    dict_cruise_xyshift = {
        "G5": [
            {
                "lonlims": [-160, -141],
                "shift": [
                    -np.sin([np.pi / 3.8]) * s1,
                    -np.cos([np.pi / 3.8]) * s1,
                ],
            },
            {"lonlims": [-141.5, -139.2], "shift": [s2, 0]},
            {
                "lonlims": [-139.1, -110],
                "shift": [-np.cos([np.pi / 3.5]) * s1, np.sin([np.pi / 3.5]) * s1],
            },
        ],
    }
    dict_asv5_locs = prep_dfs_asv_5(FN_ASV_5, dict_cruise_xyshift)

    # metaG
    s1 = 4.5
    s2 = 12
    s3 = 19.5
    dict_cruise_xyshift = {
        "G1": -s3,
        "G2": -s2,
        "G3": -s1,
    }
    dict_metag_locs = prep_dfs_metag(FN_METAG, dict_cruise_xyshift)

    # Influx
    s1 = 1.5
    s2 = 9
    s3 = 16.5
    s4 = 1.5
    s5 = 5.5
    s6 = 4.5
    s7 = 10.5
    s8 = 7.5
    dict_cruise_xyshift = {
        "2016": [{"lonlims": [-160, -150], "shift": [-s3, 0]}],
        "2017": [
            {
                "lonlims": [
                    -160,
                    -150,
                ],
                "shift": [-s2, 0],
            }
        ],
        "2019": [{"lonlims": [-160, -150], "shift": [-s1, 0]}],
        "2021": [
            {
                "lonlims": [-160, -141],
                "shift": [-np.sin([np.pi / 3.8]) * s5, -np.cos([np.pi / 3.8]) * s5],
            },
            {"lonlims": [-141.5, -139.1], "shift": [s6, 0]},
            {
                "lonlims": [-139, -110],
                "shift": [-np.cos([np.pi / 3.5]) * s8, np.sin([np.pi / 3.5]) * s8],
            },
        ],
        "2023": [
            {
                "lonlims": [-160, -141],
                "shift": [
                    -np.sin([np.pi / 3.8]) * s4,
                    -np.cos([np.pi / 3.8]) * s4,
                ],
            },
            {"lonlims": [-141.5, -139.1], "shift": [s7, 0]},
            {
                "lonlims": [-139, -110],
                "shift": [-np.cos([np.pi / 3.5]) * s4, np.sin([np.pi / 3.5]) * s4],
            },
        ],
    }
    dict_influx_locs = prep_dfs_influx(FMT_INFLUX, dict_cruise_xyshift)

    # Merge dicts
    dict_modality_locs = dict_metat_locs | dict_metag_locs | dict_asv123_locs | dict_asv5_locs | dict_influx_locs

    # Plot
    dict_modality_params = {
        "metaT_Gradients1": {"color": "#ffbb78"},
        "metaT_Gradients2": {"color": "#ffbb78"},
        "metaT_Gradients3": {"color": "#ffbb78"},
        "metaT_Gradients4": {"color": "#ffbb78"},
        "metaT_Gradients5": {"color": "#ffbb78"},
        "ASVs_G1": {"color": "#1f77b4"},
        "ASVs_G2": {"color": "#1f77b4"},
        "ASVs_G3": {"color": "#1f77b4"},
        "metaG_G1": {"color": "#9467bd"},
        "metaG_G2": {"color": "#9467bd"},
        "metaG_G3": {"color": "#9467bd"},
        "metaG_G5": {"color": "#9467bd"},
        "influx_2016": {"color": "#8c564b"},
        "influx_2017": {"color": "#8c564b"},
        "influx_2019": {"color": "#8c564b"},
        "influx_2021": {"color": "#8c564b"},
        "influx_2023": {"color": "#8c564b"},
    }

    fig, ax, fig_cbar = plot_geographic_colormap(
        LAT_MIN, LAT_MAX, LON_MIN, LON_MAX, 
        color_grid=df_latlon_sub,
        modality_locs=dict_modality_locs,
        modality_params=dict_modality_params,
    )
    print('Saving...')
    bn_out = f'{DIR_OUT}/sample_map_multimodal'
    os.makedirs(DIR_OUT, exist_ok=True)
    fig.savefig(bn_out + '.pdf', dpi=DPI)

    # Save the separate colorbar figure
    fig_cbar.savefig(f'{bn_out}_cbar.pdf', dpi=DPI, bbox_inches='tight')
    plt.close(fig_cbar)

    print(f'Saved to {bn_out}')
