#!/usr/bin/env python3
"""
Script to download OpenStreetMap residential landuse data, save it to files, and visualize it.
"""

import os
import osmnx as ox
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
from matplotlib_scalebar.scalebar import ScaleBar
import argparse

# Set the OSMnx configuration for newer versions of OSMnx
ox.settings.use_cache = True
ox.settings.log_console = True

def download_osm_data(location=None, bbox=None, tags={'landuse': 'residential'}, save_to_disk=True):
    """
    Download OpenStreetMap data for a specific location or bounding box.
    
    Args:
        location (str): Name of the location
        bbox (tuple): Bounding box as (north, south, east, west)
        tags (dict): Tags to filter OSM features (default: residential landuse)
        save_to_disk (bool): Whether to save the data to disk
        
    Returns:
        urban_areas (GeoDataFrame): GeoDataFrame containing urban area polygons
    """
    if location is None and bbox is None:
        raise ValueError("Either location or bbox must be provided")
    
    if location:
        print(f"Downloading OpenStreetMap data for {location}...")
        try:
            urban_areas = ox.features.features_from_place(location, tags=tags)
            print(f"Downloaded {len(urban_areas)} residential landuse polygons")
        except Exception as e:
            print(f"Error downloading residential areas: {e}")
            urban_areas = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    else:
        print(f"Downloading OpenStreetMap data for bounding box {bbox}...")
        north, south, east, west = bbox
        try:
            urban_areas = ox.features.features_from_bbox(north, south, east, west, tags=tags)
            print(f"Downloaded {len(urban_areas)} residential landuse polygons")
        except Exception as e:
            print(f"Error downloading residential areas: {e}")
            urban_areas = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    
    if save_to_disk and not urban_areas.empty:
        # Create output directory if it doesn't exist
        os.makedirs('output', exist_ok=True)
        
        # Create filename based on location or bbox
        name_part = location.replace(" ", "_").replace(",", "").lower() if location else f"bbox_{north}_{south}_{east}_{west}"
        
        # Save as GeoJSON for compatibility
        output_geojson = f'output/urban_areas_{name_part}.geojson'
        urban_areas.to_file(output_geojson, driver='GeoJSON')
        print(f"Saved GeoJSON to {output_geojson}")
        
        # Save as parquet for efficient storage and querying
        output_parquet = f'output/urban_areas_{name_part}.parquet'
        urban_areas.to_parquet(output_parquet)
        print(f"Saved Parquet to {output_parquet}")
    
    return urban_areas

def visualize_urban_areas(urban_areas, title=None, output_path=None, figsize=(15, 15)):
    """
    Create a visualization of urban areas (residential landuse) on a map.
    
    Args:
        urban_areas: GeoDataFrame with urban area polygons
        title: Title for the map (optional)
        output_path: Path to save the output image (optional)
        figsize: Size of the figure as tuple (width, height)
    """
    if urban_areas.empty:
        print("No data to visualize")
        return
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot urban areas
    urban_areas.plot(ax=ax, facecolor='#FF9999', alpha=0.5, edgecolor='red', linewidth=1, label='Residential Areas')
    
    # Add basemap
    try:
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)
    except Exception as e:
        print(f"Could not add basemap: {e}")
    
    # Add scale bar
    ax.add_artist(ScaleBar(1))
    
    # Remove axis
    ax.set_axis_off()
    
    # Set title and legend
    if title:
        plt.title(title, fontsize=15)
    else:
        plt.title('OpenStreetMap Urban Areas', fontsize=15)
    plt.legend(loc='upper right')
    
    # Save figure if output path is provided
    if output_path:
        # Fix: Ensure output_path has a directory component
        directory = os.path.dirname(output_path)
        if directory:  # If there is a directory component
            os.makedirs(directory, exist_ok=True)
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Map saved to {output_path}")
    
    # Show figure
    plt.tight_layout()
    plt.show()

def visualize_urban_areas_detailed(urban_areas, title=None, output_path=None, figsize=(15, 15)):
    """
    Create a more detailed visualization of urban areas with OSM context.
    
    Args:
        urban_areas: GeoDataFrame with urban area polygons
        title: Title for the map (optional)
        output_path: Path to save the output image (optional)
        figsize: Size of the figure as tuple (width, height)
    """
    if urban_areas.empty:
        print("No data to visualize")
        return
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get the bounds of the urban areas to determine map extent
    bounds = urban_areas.total_bounds  # (minx, miny, maxx, maxy)
    
    # Plot urban areas
    urban_areas.plot(ax=ax, facecolor='#FF9999', alpha=0.5, edgecolor='red', linewidth=1, label='Residential Areas')
    
    # Add basemap with more detailed tiles
    try:
        ctx.add_basemap(
            ax, 
            source=ctx.providers.CartoDB.Positron,
            zoom=14  # Higher zoom level for more detail
        )
    except Exception as e:
        print(f"Could not add detailed basemap: {e}")
        # Try fallback to regular basemap
        try:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)
        except Exception as e2:
            print(f"Could not add fallback basemap: {e2}")
    
    # Add scale bar
    ax.add_artist(ScaleBar(1))
    
    # Remove axis
    ax.set_axis_off()
    
    # Set title and legend
    if title:
        plt.title(title, fontsize=15)
    else:
        plt.title('Detailed Urban Areas Map', fontsize=15)
    plt.legend(loc='upper right')
    
    # Save figure if output path is provided
    if output_path:
        # Fix: Ensure output_path has a directory component
        directory = os.path.dirname(output_path)
        if directory:  # If there is a directory component
            os.makedirs(directory, exist_ok=True)
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Detailed map saved to {output_path}")
    
    # Show figure
    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download OpenStreetMap residential landuse data and visualize it')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--location', type=str, help='Location name to download data for')
    group.add_argument('--bbox', type=float, nargs=4, metavar=('NORTH', 'SOUTH', 'EAST', 'WEST'),
                      help='Bounding box coordinates (north, south, east, west)')
    
    parser.add_argument('--no-viz', action='store_true', help='Skip visualization')
    parser.add_argument('--detailed', action='store_true', help='Create detailed visualization')
    
    args = parser.parse_args()
    
    # Set default location for demonstration if running directly
    if len(os.sys.argv) <= 1:
        print("No arguments provided, using Manchester, UK as default location")
        location = "Manchester, UK"
        bbox = None
    else:
        location = args.location
        bbox = args.bbox
    
    # Download OSM data
    urban_areas = download_osm_data(location=location, bbox=bbox)
    
    # Output summary information
    if not urban_areas.empty:
        print("\nSummary of downloaded data:")
        print(f"Number of polygons: {len(urban_areas)}")
        print(f"Total area: {urban_areas.geometry.area.sum() / 1000000:.2f} km²")
        print(f"Coordinate reference system: {urban_areas.crs}")
        print(f"Columns available: {', '.join(urban_areas.columns)}")
    else:
        print("No data was downloaded.")
    
    # Create visualizations unless --no-viz flag is used
    if not args.no_viz and not urban_areas.empty:
        title = f"Urban Areas in {location}" if location else f"Urban Areas in Selected Bounding Box"
        
        if not args.detailed:
            # Create standard visualization
            output_path = f"urban_areas{'_' + location.replace(' ', '_').replace(',', '').lower() if location else '_bbox'}.png"
            visualize_urban_areas(urban_areas, title=title, output_path=output_path)
        else:
            # Create detailed visualization
            output_path = f"urban_areas_detailed{'_' + location.replace(' ', '_').replace(',', '').lower() if location else '_bbox'}.png"
            visualize_urban_areas_detailed(urban_areas, title=title, output_path=output_path)