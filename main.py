#!/usr/bin/env python3
import duckdb
import matplotlib.pyplot as plt
import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape, box, Point
import json
import numpy as np

def create_grid(bounds, cell_size=500):
    """
    Create a grid of points within the given bounds.
    
    Parameters:
    -----------
    bounds : tuple
        Bounds in the format (xmin, ymin, xmax, ymax)
    cell_size : float
        Size of each grid cell in meters
        
    Returns:
    --------
    GeoDataFrame with grid cells
    """
    xmin, ymin, xmax, ymax = bounds
    
    # Create grid points
    x_coords = np.arange(xmin, xmax, cell_size)
    y_coords = np.arange(ymin, ymax, cell_size)
    
    # Create a grid of points
    grid_points = []
    for x in x_coords:
        for y in y_coords:
            grid_points.append(Point(x, y))
    
    # Create grid cells (buffers around points)
    grid_cells = [point.buffer(cell_size/2) for point in grid_points]
    
    # Create GeoDataFrame with grid cells
    grid_gdf = gpd.GeoDataFrame({'geometry': grid_cells}, crs="EPSG:4326")
    
    return grid_gdf

def calculate_urban_areas(bbox=None, output_dir="output", use_local_data=False, local_data_path=None):
    """
    Calculate urban areas based on building density using DuckDB and GeoPandas.
    Optimized to use Parquet files for intermediate data storage.
    
    Parameters:
    -----------
    bbox : tuple, optional
        Bounding box in the format (xmin, ymin, xmax, ymax)
        Default is Liverpool, UK area (-2.39, 53.23, -2.07, 53.31)
    output_dir : str, optional
        Directory to store output parquet files
    use_local_data : bool, optional
        If True, use local data instead of fetching from Azure
    local_data_path : str, optional
        Path to local parquet files if use_local_data is True
        
    Returns:
    --------
    GeoDataFrame with urban area multipolygons
    """
    if bbox is None:
        bbox = (-2.39, 53.23, -2.07, 53.31)  # Default to Liverpool, UK
    
    xmin, ymin, xmax, ymax = bbox
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Connect to the DuckDB database with SSL configuration
    conn = duckdb.connect('building_data.duckdb')
    
    # Configure DuckDB for proper SSL handling
    conn.execute("INSTALL httpfs;")
    conn.execute("LOAD httpfs;")
        
    # Load spatial extension
    conn.execute("LOAD spatial;")
    
    # Create a table with buildings from Overture Maps
    print("Loading buildings data...")
    buildings_created = False
    
    try:
        if use_local_data and local_data_path:
            print(f"Using local data from: {local_data_path}")
            conn.execute(f"""
            CREATE OR REPLACE TABLE buildings AS (
              SELECT
                id,
                names.primary as primary_name,
                height,
                geometry
              FROM read_parquet('{local_data_path}')
              WHERE names.primary IS NOT NULL
              AND bbox.xmin BETWEEN {xmin} AND {xmax}
              AND bbox.ymin BETWEEN {ymin} AND {ymax}
            )
            """)
            buildings_created = True
        else:
            print("Fetching data from Azure Blob Storage...")
            try:
                # Try with az:// protocol first
                conn.execute(f"""
                CREATE OR REPLACE TABLE buildings AS (
                  SELECT
                    id,
                    names.primary as primary_name,
                    height,
                    geometry
                  FROM read_parquet('az://overturemapswestus2.blob.core.windows.net/release/2025-04-23.0/theme=buildings/type=building/*', 
                                  filename=true, 
                                  hive_partitioning=1)
                  WHERE names.primary IS NOT NULL
                  AND bbox.xmin BETWEEN {xmin} AND {xmax}
                  AND bbox.ymin BETWEEN {ymin} AND {ymax}
                )
                """)
                buildings_created = True
            except Exception as e:
                print(f"Error with az:// protocol: {str(e)}")
                # Try with https:// protocol instead
                try:
                    print("Trying with https:// protocol...")
                    conn.execute(f"""
                    CREATE OR REPLACE TABLE buildings AS (
                      SELECT
                        id,
                        names.primary as primary_name,
                        height,
                        geometry
                      FROM read_parquet('https://overturemapswestus2.blob.core.windows.net/release/2025-04-23.0/theme=buildings/type=building/*', 
                                      filename=true, 
                                      hive_partitioning=1)
                      WHERE names.primary IS NOT NULL
                      AND bbox.xmin BETWEEN {xmin} AND {xmax}
                      AND bbox.ymin BETWEEN {ymin} AND {ymax}
                    )
                    """)
                    buildings_created = True
                except Exception as e2:
                    print(f"Error with https:// protocol: {str(e2)}")
                    
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        if not buildings_created:
            print("Could not create buildings table. Using sample data instead.")
            # Create a simple sample table for testing
            conn.execute("""
            CREATE OR REPLACE TABLE buildings AS (
              SELECT 
                CAST(row_number() over() AS VARCHAR) as id, 
                'Building ' || CAST(row_number() over() AS VARCHAR) as primary_name,
                CAST(random() * 50 AS FLOAT) as height,
                ST_GeomFromText('POINT(' || (-2.2 + random()/10) || ' ' || (53.25 + random()/10) || ')') as geometry
              FROM range(1000)
            )
            """)
    
    # Save the filtered buildings to parquet for reuse
    print("Saving filtered buildings to parquet...")
    conn.execute(f"COPY buildings TO '{output_dir}/buildings.parquet' (FORMAT PARQUET)")
    
    # Get the building geometries as GeoJSON
    print("Extracting building geometries...")
    buildings_df = conn.execute("""
    SELECT 
        id, 
        ST_AsGeoJSON(geometry) as geom_json,
        ST_X(geometry) as x,
        ST_Y(geometry) as y
    FROM buildings
    """).fetchdf()
    
    # Calculate the bounds for grid creation
    min_x = buildings_df['x'].min()
    max_x = buildings_df['x'].max()
    min_y = buildings_df['y'].min()
    max_y = buildings_df['y'].max()
    
    print("Creating building density grid...")
    
    # Create a grid using GeoPandas
    grid = create_grid((min_x, min_y, max_x, max_y), cell_size=0.005)  # ~500m at this latitude
    
    # Convert buildings to GeoDataFrame
    buildings_gdf = gpd.GeoDataFrame({
        'id': buildings_df['id'],
        'geometry': [shape(json.loads(geom)) for geom in buildings_df['geom_json']]
    }, crs="EPSG:4326")
    
    # Count buildings in each grid cell
    print("Counting buildings in grid cells...")
    building_counts = []
    for idx, cell in grid.iterrows():
        count = buildings_gdf.intersects(cell.geometry).sum()
        if count > 5:  # Only consider cells with more than 5 buildings as urban
            building_counts.append({
                'geometry': cell.geometry,
                'building_count': count
            })
    
    # Create a GeoDataFrame with building density
    building_density_gdf = gpd.GeoDataFrame(building_counts, crs="EPSG:4326")
    
    # Save the building density grid to parquet
    print("Saving building density grid to parquet...")
    building_density_gdf.to_parquet(f"{output_dir}/building_density.parquet")
    
    # Dissolve adjacent cells to create urban area polygons
    print("Creating urban areas by dissolving adjacent cells...")
    
    # Fix 1: Use union_all() instead of deprecated unary_union
    dissolved_geom = building_density_gdf.geometry.union_all()
    
    urban_areas_gdf = gpd.GeoDataFrame({
        'geometry': [dissolved_geom],
        'total_buildings': [building_density_gdf['building_count'].sum()],
        'avg_density': [building_density_gdf['building_count'].mean()]
    }, crs="EPSG:4326")
    
    # Create detailed urban areas by density class
    print("Creating detailed urban areas by density class...")
    building_density_gdf['density_class'] = pd.cut(
        building_density_gdf['building_count'],
        bins=[5, 20, 50, float('inf')],
        labels=['low_density', 'medium_density', 'high_density']
    )
    
    # Dissolve by density class
    urban_areas_detailed_gdf = building_density_gdf.dissolve(by='density_class', aggfunc={
        'building_count': 'sum'
    }).reset_index()
    urban_areas_detailed_gdf = urban_areas_detailed_gdf.rename(columns={'building_count': 'total_buildings'})
    
    # Fix 2: Project to a local UTM zone for accurate area calculations
    # Find appropriate UTM zone for the area (Liverpool is approximately in UTM zone 30N)
    urban_areas_projected = urban_areas_detailed_gdf.to_crs("EPSG:32630")  # UTM zone 30N
    
    # Calculate density using the projected coordinates for accurate area
    urban_areas_projected['area_m2'] = urban_areas_projected.geometry.area
    urban_areas_detailed_gdf['avg_density'] = urban_areas_projected['total_buildings'] / urban_areas_projected['area_m2'] * 1000000  # Buildings per sq km
    
    # Save the results to parquet files
    print("Saving final results to parquet...")
    urban_areas_gdf.to_parquet(f"{output_dir}/urban_areas.parquet")
    urban_areas_detailed_gdf.to_parquet(f"{output_dir}/urban_areas_detailed.parquet")
    
    # Close the connection
    conn.close()
    
    return urban_areas_gdf, urban_areas_detailed_gdf

def main():
    # Try to use local data first if it exists
    local_buildings_path = 'output/buildings.parquet'
    use_local_data = os.path.exists(local_buildings_path)
    
    # Calculate urban areas for the default bounding box
    urban_areas, urban_areas_detailed = calculate_urban_areas(
        use_local_data=use_local_data,
        local_data_path=local_buildings_path if use_local_data else None
    )
    
    # Save the results to GeoJSON format
    urban_areas.to_file("urban_areas.geojson", driver="GeoJSON")
    urban_areas_detailed.to_file("urban_areas_detailed.geojson", driver="GeoJSON")
    
    # Plot the overall urban areas
    fig, ax = plt.subplots(figsize=(12, 10))
    urban_areas.plot(ax=ax, column='avg_density', legend=True, cmap='Reds', 
                     legend_kwds={'label': "Avg. Building Density (bldgs/km²)"})
    ax.set_title("Urban Areas Derived from Building Density")
    plt.savefig("urban_areas.png", dpi=300)
    
    # Plot the detailed urban areas classification
    fig2, ax2 = plt.subplots(figsize=(12, 10))
    urban_areas_detailed.plot(ax=ax2, column='density_class', legend=True, 
                              categorical=True, cmap='viridis')
    ax2.set_title("Urban Areas Classification by Building Density")
    plt.savefig("urban_areas_detailed.png", dpi=300)
    
    print(f"Found {len(urban_areas)} urban areas")
    print(f"Total buildings: {urban_areas['total_buildings'].sum()}")
    print("Results saved as:")
    print("- urban_areas.geojson and urban_areas.png (overall urban boundaries)")
    print("- urban_areas_detailed.geojson and urban_areas_detailed.png (density classifications)")
    print("- All intermediate data saved as parquet files in the 'output' directory")
    
    # Additional analysis using DuckDB with parquet files
    print("\nPerforming additional analysis using the saved parquet files...")
    conn = duckdb.connect('building_data.duckdb')
    
    # Configure DuckDB for proper SSL handling for this connection too
    conn.execute("INSTALL httpfs;")
    conn.execute("LOAD httpfs;")
    
    # Example of reusing the parquet files for further analysis
    try:
        result = conn.execute("""
        SELECT 
            density_class,
            sum(total_buildings) as buildings,
            ROUND(AVG(avg_density), 1) as avg_density_per_sqkm
        FROM read_parquet('output/urban_areas_detailed.parquet')
        GROUP BY density_class
        ORDER BY avg_density_per_sqkm DESC
        """).fetchdf()
        
        print("\nUrban density classification statistics:")
        print(result)
    except Exception as e:
        print(f"Error performing additional analysis: {str(e)}")
        
    conn.close()

if __name__ == "__main__":
    main()