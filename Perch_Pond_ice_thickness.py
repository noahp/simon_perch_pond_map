from dataclasses import dataclass

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry import Point


def read_excel_file(excel_file_path):
    # Step 1: Read the Excel file containing observations
    observations_df = pd.read_excel(excel_file_path)

    @dataclass
    class Observation:
        latitude: float
        longitude: float
        thickness: float

    # Python dictionary, keys are the date of observations, values are lists of Observations
    obs_bydate = {}

    for index, row in observations_df.iterrows():
        date = row["Date"]
        obs = Observation(row["y"], row["x"], row["Thickness"])
        if date not in obs_bydate:
            obs_bydate[date] = []
        obs_bydate[date].append(obs)

    # Debug, print out the dict
    # for date in obs_bydate:
    #     print(date)
    #     for obs in obs_bydate[date]:
    #         print(f"  {obs.latitude}, {obs.longitude}, {obs.thickness}")

    return obs_bydate


def read_shapefile(shapefile_path):
    # Step 2: Read the shapefile containing the boundary polygon
    boundary_gdf = gpd.read_file(shapefile_path)

    # filter out desired polygon
    for index, shape in boundary_gdf.iterrows():
        if shape["GNIS_ID"] == "00872485":
            print(shape)
            return shape

    raise ValueError("No shape found with GNIS_ID 00872485")


def filter_observations(observations_df, boundary_gdf):
    # Step 3: Filter observations that fall within the boundary polygon
    boundary_polygon = boundary_gdf.geometry.iloc[0]
    within_boundary = observations_df.apply(
        lambda row: boundary_polygon.contains(Point(row["longitude"], row["latitude"])),
        axis=1,
    )
    observations_within_boundary = observations_df[within_boundary]

    return observations_within_boundary


def interpolate_thickness(observations_within_boundary, boundary_polygon):
    # Step 4: Interpolate thickness values within the boundary polygon
    x = observations_within_boundary["longitude"]
    y = observations_within_boundary["latitude"]
    thickness = observations_within_boundary["thickness"]
    interpolator = LinearNDInterpolator(list(zip(x, y)), thickness)

    # Step 5: Create an interpolated map of thickness values within the polygon
    # Create a grid of points within the boundary polygon
    min_x, min_y, max_x, max_y = boundary_polygon.bounds
    x_grid, y_grid = np.meshgrid(
        np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100)
    )
    points = np.vstack([x_grid.ravel(), y_grid.ravel()]).T

    # Interpolate thickness values at grid points
    thickness_grid = interpolator(points).reshape(x_grid.shape)

    return thickness_grid


def plot_thickness_map(boundary_gdf, thickness_grid, min_x, max_x, min_y, max_y):
    # Plot interpolated thickness map
    fig, ax = plt.subplots(figsize=(10, 8))
    boundary_gdf.plot(ax=ax, color="lightgray", edgecolor="black")
    im = ax.imshow(
        thickness_grid,
        extent=(min_x, max_x, min_y, max_y),
        origin="lower",
        cmap="viridis",
        alpha=0.5,
    )
    plt.colorbar(im, label="Thickness")
    ax.set_title("Interpolated Thickness Map within Boundary Polygon")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.show()


# Python will run this if the script is the main file, i.e. `python <this file>.py`
if __name__ == "__main__":
    excel_file_path = "S123_11f180d520ad4d4a8a82480bab5d5754_EXCEL.xlsx"

    # load the observations, bucketed by date of observation
    obs_by_date = read_excel_file(excel_file_path)

    # confirm the shapefile is sane
    shapefile_path = "New_Hampshire_Hydrography_Dataset_(Waterbody)/New_Hampshire_Hydrography_Dataset_(Waterbody).shp"
    pond_shape = read_shapefile(shapefile_path)

    # filter_observations(excel_file_path, shapefile_path)

    # interpolate_thickness(excel_file_path, "Perch_Pond.shp")

    # plot_thickness_map(excel_file_path, "Perch_Pond.shp")
