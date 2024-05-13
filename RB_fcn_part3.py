
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
#%% Select_class_ReefCloud_pts

def Select_class_ReefCloud_pts(RB_centroid_csv, RC_csv, label_file, polygon_file, label_poygon_seg, label_poygon_csv):

    #access and read files
    RB_centroid_csv = pd.read_csv(RB_centroid_csv)
    #sam_centroi = RB_centroid_csv[['segment_un', 'center_point_x', 'center_point_y']]

    df = pd.read_csv(RC_csv)

    label = pd.read_csv(label_file)
    label['label_set'] = label['label_set'].replace({'AC _Bran_OTH': 'AC_Bran_OTH'})

    #drop water and unid classes
    df_UNid = df.apply(lambda row: row[df['pred_desc'].isin(['Un_ID','WA'])]).index
    df.drop(df_UNid, inplace = True)
    df['pred_desc']= df['pred_desc'].replace({'G_STY ': 'G_STY'})

    df_all_sub_04 = df[df['pred_score']>=0.4] #filter 0.40 pred_score


    # append Moprhology classes
    df_classNum = df_all_sub_04.merge(label[['label_set', 'Morphology','ID_number','Morpho_number']], how='left', left_on='pred_desc', right_on='label_set')

    #for Hexagrid drop morpho soft_coral, branching_non_ac, Acropora (except ACS and ACX)
    Hexagrid_filter = df_classNum.loc[((df_classNum.SAM_centroid.str.contains('Hexa')) &
                                       ((df_classNum['Morphology'].isin(['Soft_Coral','Branching_non_acropora']))|
                                        (df_classNum['pred_desc'].isin(['ACO_MIL','ACO_TEN', 'ACO_OTH', 'ACD', 'AC_Bran_OTH', 'ACT']))))]
    df_classNum.drop(Hexagrid_filter.index, inplace = True) #, 'Foliose', 'Branching_non_acropora',

    # pre_filter morphology
    pred = 'Morphology'
    df = df_classNum
    #create series with unique values per mask

    most_occuring_morpho = df.groupby(['SAM_centroid'])[pred].agg(pd.Series.mode)


    most_occuring_morpho_df = most_occuring_morpho.to_frame().reset_index()

    df_drop_morpho = pd.merge(df, most_occuring_morpho_df,
                              how='inner',
                              on='SAM_centroid')

    df_t = df_drop_morpho.loc[df_drop_morpho.apply(lambda row: row.Morphology_x in row.Morphology_y, axis=1)]
    df_t = df_t.drop(columns = 'Morphology_y')
    df_t = df_t.rename(columns ={'Morphology_x':'Morphology'})
    # get label value with weighted formula
    pred = 'pred_desc' #'Morphology'
    df = df_t

    #create series with unique values per mask
    unique = df.groupby(['SAM_centroid']).size()
    most_occuring_class = df.groupby(['SAM_centroid'])[pred].agg(pd.Series.mode)
    count = df.groupby(['SAM_centroid']).size()
    max_pred_row = df.loc[df.groupby(["SAM_centroid"])["pred_score"].idxmax()]

    #get most occuring class and mean of prediction score
    occur = df.groupby(['SAM_centroid', pred]).size()
    mean = df.groupby(['SAM_centroid', pred])['pred_score'].mean()
    max_pred = df.groupby(['SAM_centroid', pred])['pred_score'].max()
    fctn = (occur * mean)

    #create dataframe with most occuring classes
    df_result = pd.DataFrame( columns= ['n_occur' , 'pred_score_mean', 'Weighted_score', 'max_pred'])
    df_result['n_occur'] = occur
    df_result['pred_score_mean'] = mean
    df_result['Weighted_score'] = fctn
    df_result["max_pred"] = max_pred
    Weighted_class_result = df_result[df_result['Weighted_score'] == df_result.groupby(level=0)['Weighted_score'].transform(max)]
    Weighted_class_result=Weighted_class_result.reset_index(level=[pred])

    #merge unique value and weighted classes
    df_merged =Weighted_class_result.merge(max_pred_row, on='SAM_centroid') # left_index=True, right_index=True)
    df_merged.reset_index(inplace=True)

    # Export file
    RB_SHP = geopandas.read_file(polygon_file)
    df_merged = df_merged.merge(RB_SHP, left_on='SAM_centroid', right_on = "segment_un")
    gdf = GeoDataFrame(df_merged, geometry='geometry')

    gdf.to_file(label_poygon_seg)
    gdf.to_csv(label_poygon_csv)

    return gdf

def format_percent_cover(label_polygon_seg, label_file):
    label = pd.read_csv(label_file)
    label['label_set'] = label['label_set'].replace({'AC _Bran_OTH': 'AC_Bran_OTH'})
    RB_shp_df = geopandas.read_file(label_poygon_seg)
    RB_shp_df['pred_desc_']= RB_shp_df.pred_desc_.replace({'G_STY ': 'G_STY'})
    RB_shp_df_classes = RB_shp_df.groupby(['pred_desc_'])['area'].sum()
    RB_shp_df_classes_result = pd.DataFrame(columns=['SAM_10X4_Total_area', 'SAM_10X4_Percent_cover'])
    RB_shp_df_classes_result['SAM_10X4_Total_area']= RB_shp_df_classes
    RB_shp_df_classes_result['SAM_10X4_Percent_cover']= RB_shp_df_classes/(RB_shp_df_classes.sum())*100

    RB_shp_df_classes_result_group = RB_shp_df_classes_result.copy()
    RB_shp_df_classes_result_group["class"] = RB_shp_df_classes_result_group.index
    RB_shp_df_classes_result_group = RB_shp_df_classes_result_group.reset_index()
    RB_shp_df_classes_result_group =RB_shp_df_classes_result_group.drop(['pred_desc_'], axis= 1)
    RB_shp_df_classes_result_group = RB_shp_df_classes_result_group.rename(columns={"SAM_10X4_Total_area": "Total_area", "SAM_10X4_Percent_cover": "Percent_cover"})
    label = label.dropna()
    RB_shp_df_classes_result_group_morpho = label.merge(RB_shp_df_classes_result_group, left_on='label_set', right_on='class', how='outer')
    RB_shp_df_classes_result_group_morpho = RB_shp_df_classes_result_group_morpho.dropna(subset = ['class'])
    RB_shp_df_classes_result_group_morpho['treatment']='RB'
    #create color list dictionary

    hex=[]

    for i in RB_shp_df_classes_result_group_morpho['RGB']:
        chunks = i.split(',')
        hex_1 =rgb_to_hex(int(chunks[0]), int(chunks[1]), int(chunks[2]))
        hex_1 = "#"+hex_1
        hex.append(hex_1)

    RB_shp_df_classes_result_group_morpho['hex']=hex

    color = zip(RB_shp_df_classes_result_group_morpho['label_set'].unique(), RB_shp_df_classes_result_group_morpho['hex'].unique())
    color_list = list(color)
    color_list = dict(color_list)

    label = label.drop(label[label.Morphology.isin([0])].index)
    label = label.dropna()
    hex=[]

    for i in label['RGB']:
        chunks = i.split(',')
        hex_1 =rgb_to_hex(int(chunks[0]), int(chunks[1]), int(chunks[2]))
        hex_1 = "#"+hex_1
        hex.append(hex_1)

    label['hex']=hex

    #plot multiple dataset per graph



    enmax_palette = label['hex']
    color_codes_wanted = label['label_set']
    cdict = dict(zip(color_codes_wanted, [mcolors.to_rgba(c) for c in enmax_palette]))

    palette_RB = mcolors.get_named_colors_mapping().update(cdict)

    return RB_shp_df_classes_result_group_morpho, palette_RB


def stack_catplot(x, y, cat, stack, data, palette, out_fig, order):


    f = plt.figure( figsize=(6,4), dpi=600)
    gs = gridspec.GridSpec(2,1, height_ratios=[1, 3])
    ax = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    # pivot the data based on categories and stacks
    df = data.pivot_table(values=y, index=[cat, x], columns=stack,
                          dropna=False, aggfunc='sum').fillna(0)
    df['Morph'] = df.index.get_level_values('Morphology')
    df['Morph'] = pd.Categorical(df['Morph'], order)
    df = df.sort_values('Morph')
    morph_class = df.Morph
    df=df.drop(columns = 'Morph')
    total_list = df.sum(axis=1).round(2)
    total_list = list(map(str, total_list))
    total_list =  list(map(lambda x: '<0.01' if x == '0.0' else x, total_list))
    #total_list = list(map(lambda x: x.replace('0.0', '<0.01'), total_list))


    ncat = data[cat].nunique()
    nx = data[x].nunique()
    nstack = data[stack].nunique()
    range_x = np.arange(nx)
    width = 0.8 / ncat # width of each bar

    for i, c in enumerate(data[cat].unique()):
        # iterate over categories, i.e., Conditions
        # calculate the location of each bar
        loc_x = (0.5 + i - ncat / 2) * width + range_x
        bottom = 0
        for j, s in enumerate(data[stack].unique()):
            #print(j, s)
            # iterate over stacks, i.e., Hosts
            # obtain the height of each stack of a bar
            height = df.loc[c][s].values
            # plot the bar, you can customize the color yourself
            ax.bar(x=loc_x, height=height, bottom=bottom, width=width, zorder=10, color=s)
            ax2.bar(x=loc_x, height=height, bottom=bottom, width=width, zorder=10, color=s)
            # change the bottom attribute to achieve a stacked barplot
            bottom += height


    ax.set_ylim(20, 100)  # outliers only
    ax2.set_ylim(0, 20)  # most of the data
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.set_visible(False)#ax.xaxis.tick_top()
    #ax.tick_params(labeltop='off')  # don't put tick labels at the top
    #ax2.xaxis.tick_bottom()
    d = .015
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    # make xlabel
    ax2.set_xticks(range_x)
    ax2.xaxis.tick_bottom()
    ax2.set_xticklabels(sorted(data[x].unique()), rotation=90)
    #ax.set_ylabel(y) # add y axis label
    ax2.set_xticklabels((morph_class), rotation=90)
    ax2.set_ylabel('     Percent cover')

    #labels = ct.tot.tolist()
    ax.bar_label(ax.containers[-1], labels=total_list, padding=3, fontsize=7 )
    ax2.bar_label(ax.containers[-1], labels=total_list, padding=3,fontsize=7 )

    f.savefig(out_fig, bbox_inches='tight')
    return f

def rgb_to_hex(r, g, b):
    return ('{:02X}' * 3).format(r, g, b)

def ColonyLevel_Segments(label_poygon_seg, ColonyLevelSegments_shp):
    #apply buffers based on morphology to merge
    gdf = geopandas.read_file(label_poygon_seg)
    gdf_NH = gdf.loc[~((gdf.SAM_centro.str.contains('Hexa')) & (~gdf['pred_desc_'].isin(['ACS', 'ACX', 'CB_OTH_C']))) #'G_POR'
                     | ((gdf['Morphology'].isin(['Acropora']) & (gdf['area'] <= 0.0006)))] #filterout Hexagrid
    gdf_buff = gdf_NH.loc[(gdf_NH['area'] <= 0.0025) | (~gdf_NH['Morphology'].isin(['Acropora', 'Branching_non_acropora', 'Soft_Coral']))]
    gdf_buff['geometry'] = gdf_buff.geometry.buffer(0.0004)
    gdf_Non_buff = gdf_NH[~gdf_NH.index.isin(gdf_buff.index)]
    gdf_buff_non_buff = geopandas.GeoDataFrame(pd.concat([gdf_buff, gdf_Non_buff]))
    gdf_buff_non_buff['SAM_centro_ID']= gdf_buff_non_buff.SAM_centro

    #merge based on class
    unique_ID = gdf_buff_non_buff['pred_desc_'].unique()
    colony_gdf = []
    for i in unique_ID:
        class_shp = gdf_buff_non_buff[gdf_buff_non_buff['pred_desc_'] ==i]
        class_shp_uu = class_shp.unary_union
        cuu = geopandas.GeoSeries(class_shp_uu)
        df2 = cuu.explode(ignore_index=True)
        df3= GeoDataFrame(geopandas.GeoSeries(df2))
        df3['Colony_ID']= i
        colony_gdf.append(df3)

    colony_gdf_c = GeoDataFrame(pd.concat(colony_gdf, ignore_index=True) )
    colony_gdf_c.rename(columns={ colony_gdf_c.columns[0]: "geometry" }, inplace = True)
    colony_gdf_c = colony_gdf_c.set_geometry("geometry")
    colony_gdf_c['area']=colony_gdf_c.area

    #fill holes
    colony_gdf_c['rings'] = colony_gdf_c['geometry'].apply(lambda geom: [i for i in geom.interiors])
    colony_gdf_c['to_fill'] = colony_gdf_c['rings'].apply(lambda rings: [Polygon(ring) for ring in rings if Polygon(ring).area < 0.001])
    colony_gdf_c_toFill= colony_gdf_c[(colony_gdf_c['to_fill'].apply(len) != 0)]


    def process_row(chunk):     # Function to process each row
        from functools import reduce
        return chunk.apply(lambda row: reduce(lambda geom1, geom2: geom1.union(geom2), [row['geometry']] + row['to_fill']), axis=1)

    num_processes = mp.cpu_count()      # Number of parallel processes
    chunk_size = len(colony_gdf_c_toFill) // num_processes
    chunks = [colony_gdf_c_toFill[i:i+chunk_size] for i in range(0, len(colony_gdf_c_toFill), chunk_size)]      # Split the DataFrame into chunks
    pool = mp.Pool(processes=num_processes)     # Create a multiprocessing Pool
    results = pool.map(process_row, chunks)     # Apply the function to each chunk in parallel

    # Close the pool to release resources
    pool.close()
    pool.join()

    # Concatenate the results
    final_result = pd.concat(results)

    colony_gdf_c_toFill['newgeom']= final_result

    colony_gdf_c_toFill['new_geometry'] = colony_gdf_c_toFill.apply(lambda row: reduce(lambda geom1, geom2: geom1.union(geom2), [row['geometry']] + row['to_fill']), axis=1)

    colony_gdf_c_NotToFill= colony_gdf_c[(colony_gdf_c['to_fill'].apply(len) == 0)]
    colony_gdf_c_NotToFill['newgeom']=colony_gdf_c_NotToFill['geometry']
    colony_gdf_c_merged = pd.concat([colony_gdf_c_NotToFill,colony_gdf_c_toFill], sort = True)
    colony_gdf_c_merged.drop(['rings', 'to_fill', 'geometry'], axis=1, inplace=True)
    colony_gdf_c_merged=colony_gdf_c_merged.rename(columns={'newgeom':'geometry'})
    colony_gdf_c_merged['area']=colony_gdf_c_merged.area
    colony_gdf_c_merged.drop(['new_geometry'], axis=1, inplace=True)

    colony_gdf_c_merged.to_file(ColonyLevelSegments_shp)
    return colony_gdf_c_merged
