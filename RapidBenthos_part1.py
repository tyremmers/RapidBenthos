#import libraries
import os
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from functools import reduce
from tqdm import tqdm
import rasterio
from rasterio.mask import mask
from PIL import Image, ImageFile
Image.MAX_IMAGE_PIXELS = None
from RB_fcn_part1 import Filter_segments, camera_point_from_segment_centerPoint,hexagrid
from osgeo import gdal
import numpy as np
import cv2
from datetime import datetime
def timestamp():
    return datetime.now().strftime('%Y-%m-%d')

from samgeo import SamGeo, tms_to_geotiff, get_basemaps#


ImageFile.LOAD_TRUNCATED_IMAGES = True

os.environ['CUDA_VISIBLE_DEVICES'] = '0'
import Metashape
#agisoft_LICENSE=5053@agisoft-lmgr.aims.gov.au

ts = timestamp()

#%% set orthomosaic path

ortho =  #input path to orthomosaic
out_folder = #input path to output folder
plot_id = #input plot|site name

print(plot_id)

#%% Polygonise Raster X2


#1 finest moving window


#Set kwargs and SAM parameters
device = 'cuda:0'

sam_kwargs = {
    'points_per_side': 128,
    'points_per_batch': 16,
    'pred_iou_thresh': 0.88,
    'stability_score_thresh': 0.94,
    'stability_score_offset': 1.0,
    'box_nms_thresh':  0.35,
    'crop_n_layers': 0,
    'crop_nms_thresh': 0.9,
    #'crop_overlap_ratio': 1500 / 1500,
    'crop_n_points_downscale_factor': 1,
    'min_mask_region_area': 1600,}

sam = SamGeo(
    model_type="vit_h",
    checkpoint="sam_vit_h_4b8939.pth",
    device = device,
    sam_kwargs=sam_kwargs,
)



mask_1 = os.path.join(out_folder, '{}.tif'.format(plot_id  + ts))
sam.generate(ortho, mask_1, batch=True, foreground=False, mask_multiplier=255, erosion_kernel=(3, 3), sample_size=(5400, 5400), bound=100)

#2 larger moving window

#Set kwargs and SAM parameters
device = 'cuda:0'

sam_kwargs = {
    'points_per_side': 200,
    'points_per_batch': 8,
    'pred_iou_thresh': 0.88,
    'stability_score_thresh': 0.94,
    'stability_score_offset': 1.0,
    'box_nms_thresh':  0.35,
    'crop_n_layers': 0,
    'crop_nms_thresh': 0.9,
    #'crop_overlap_ratio': 1500 / 1500,
    'crop_n_points_downscale_factor': 1,
    'min_mask_region_area': 1600,}

sam = SamGeo(
    model_type="vit_h",
    checkpoint="sam_vit_h_4b8939.pth",
    device = device,
    sam_kwargs=sam_kwargs,
)


mask_2 = os.path.join(out_folder, '{}.tif'.format(plot_id +  ts))
sam.generate(ortho, mask_2, batch=True, foreground=False, mask_multiplier=255, erosion_kernel=(3, 3), sample_size=(8400, 8400), bound=200)



#%%Gdal raster calculator to merge all raster
dir_path = #input path to lib/python3.8/site-packages/osgeo_utils/
file_name = "gdal_calc.py"
gdal_calc_path = os.path.join(dir_path, file_name)


Combined_seg_tif =  os.path.join(out_folder, '{}.tif'.format(plot_id +  ts))
calc_expr = '"A * B"'
typeof = '"Float32"'

# Generate string of process.
gdal_calc_str = 'python {0} -A {1} -B {2}  --outfile={3} --calc={4} --type={5} --hideNoData'
gdal_calc_process = gdal_calc_str.format(os.path.join(dir_path, file_name), mask_1, mask_2,
                                         Combined_seg_tif, calc_expr, typeof)
# Call process
os.system(gdal_calc_process)
print('gdal_calc_process done')


#%% vectorise tif
Combined_seg_shp = os.path.join(out_folder, '{}.gpkg'.format(plot_id +  ts))
sam.tiff_to_gpkg(Combined_seg_tif, Combined_seg_shp, simplify_tolerance=None)

#%%Filter segments function
SEG_shp_df_path=os.path.join(out_folder, '{}.shp'.format(plot_id + ts))
PTS_shp_df_path=os.path.join(out_folder, '{}.shp'.format(plot_id +  ts))
PTS_csv_df_path=os.path.join(out_folder, '{}.csv'.format(plot_id +  ts))



segments_df_filtered, segments_pts_df = Filter_segments(Combined_seg_shp, SEG_shp_df_path, PTS_shp_df_path, PTS_csv_df_path)


#%% Create hexagrid, merge with segments and export shp + csv
full_grid_path = os.path.join(out_folder, '{}.shp'.format(plot_id + ts))
clip_grid_path = os.path.join(out_folder, '{}.shp'.format(plot_id + ts))
hexagrid_RB_union_path_seg_shp = os.path.join(out_folder, '{}.shp'.format(plot_id + ts))
hexagrid_RB_union_pst_shp = os.path.join(out_folder, '{}.shp'.format(plot_id + ts))
hexagrid_RB_union_path_pts_csv = os.path.join(out_folder, '{}.csv'.format(plot_id + ts))


hexagird_union_shp, hexagrid_union_pts = hexagrid(ortho, 0.05, full_grid_path, SEG_shp_df_path, clip_grid_path, hexagrid_RB_union_path_seg_shp, hexagrid_RB_union_pst_shp, hexagrid_RB_union_path_pts_csv)


#%% Point from centroid no filter function
MetashapeProject_path = #input path to Metashape project
Chunk_number= #input chunk number
PhotoPath = #input path to underlying photos repository
OutputPath = os.path.join(out_folder, plot_id + '_{}.csv')

camera_uv = camera_point_from_segment_centerPoint(MetashapeProject_path, Chunk_number, PhotoPath, OutputPath, hexagrid_RB_union_path_pts_csv)

