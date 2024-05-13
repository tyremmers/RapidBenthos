
#from RB_fcn_part3 import Select_class_ReefCloud_pts

import numpy as np
import pandas as pd
import geopandas
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import Polygon
import tqdm
from functools import reduce
import multiprocess as mp
import matplotlib.gridspec as gridspec
from RB_fcn_part3 import Select_class_ReefCloud_pts, format_percent_cover,stack_catplot, rgb_to_hex, ColonyLevel_Segments


#%% Select_class_ReefCloud_pts
RB_centroid_csv = #input path to RapidBenthos part 1 polygon center point csv (hexagrid_RB_union_path_pts_csv)
RC_csv = #input path to ReefCloud annotation csv output
label_file = #input path to label file
polygon_file = #input path to RapidBenthos part 1 polygon shapefile (hexagrid_RB_union_path_seg_shp)
label_poygon_seg = #input path to RapidBenthos labeled polygon output shapefile
label_poygon_csv = #input path to RapidBenthos labeled polygon output CSV

gdf = Select_class_ReefCloud_pts(RB_centroid_csv, RC_csv, label_file, polygon_file, label_poygon_seg, label_poygon_csv)

#%% percent_cover_barGraph

RB_shp_df_classes_result_group_morpho, palette_RB = format_percent_cover(label_poygon_seg, label_file)
PercentCover = #input path to RapidBenthos percent cover csv
RB_shp_df_classes_result_group_morpho.to_csv(PercentCover)
RB_shp_df_classes_result_group_morpho['Morphology'] = RB_shp_df_classes_result_group_morpho['Morphology'].replace({'Columnar (CAAB 11 290915)': 'Columnar'})
Class_order = ['Acropora Corymbose (ACO)', 'Acropora', 'Branching_non_acropora', "Massive", 'Foliose', 'Encrusting', 'Columnar', 'Fire_Coral','Fungiidae', 'Soft_Coral', 'Sponge', 'Algae', 'Abiotic_substrate','Mobile_Biota', 'Markers', 'Unidentifiable' ]

fi2 = plt.figure(figsize=(6,3), dpi=600)
out_fig = # input path to community composition stacked bar graph
stack_catplot(x='Morphology', y='Percent_cover', cat='treatment', stack='label_set', data=RB_shp_df_classes_result_group_morpho, palette = palette_RB, out_fig = out_fig, order =Class_order)

#%% Colony level segments

ColonyLevelSegments_shp = #input path to RapidBenthos colony level segments shapefile

gdf_colonyLevel = ColonyLevel_Segments(label_poygon_seg, ColonyLevelSegments_shp)
