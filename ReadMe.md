The official repository for the paper: "RapidBenthos – Automated segmentation and multi-view classification of coral reef communities from photogrammetric reconstruction".

If this repository contributes to your research, please consider citing the publication below.

```
Remmers, T., Boutros, N., Wyatt, M., Gordon, S., Toor, T., Roelfsema, C., Fabricius, K., Grech, A., Lechene, L., Ferrari, R. 2024. RapidBenthos – Automated segmentation and multi-view classification of coral reef communities from photogrammetric reconstruction.  
```

### Bibtex
```
insert bibtex here

```

The full paper can be accessed at: \[[Paper]()]. Manuscript in revision!

## Table of Contents
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Acknowledgements](#acknowledgements)

<a name="installation"></a>

## Installation
1. Create and activate a virtual/conda environment, then inside the root directory of the project run:

    ```conda env create --RapidBenthos -f environment.yml```
    
2. Download the segment-geospatial github repository and save it inside the root directory of the project:

    ```https://github.com/opengeos/segment-geospatial.git```

3. Download the dense_inference github repository and save it inside the root directory of the project:

    ```https://github.com/open-AIMS/dense_inference.git```

<a name="quick-start"></a>
## Quick Start
1. To segment the orthomosaic and extract the points to classify in the undelying field photo, fill all the requiered fields in the RapidBenthos_part1.py:

    ```ortho =  #input path to orthomosaic```

    ```out_folder = #input path to output folder```

    ```plot_id = #input plot|site name```

    ```dir_path = #input path to lib/python3.8/site-packages/osgeo_utils/```

    ```MetashapeProject_path = #input path to Metashape project```

    ```Chunk_number= #input chunk number```

    ```PhotoPath = #input path to underlying photos repository```

   Once all field are populated run RapidBenthos_part1.py.


3. To classify points from a csv (e.g. from the RapidBenthos_part1 outputs)

From the root directory of the project run:

    python scripts/inference_from_csv.py -d 'path/to/points.csv' -i 'path/to/image_directory' -e 'JPG'

For example, if the csv is `data/images/SAM_points_test.csv` and images are JPGs and are in `data/images/test_ims`

    python scripts/inference_from_csv.py -d 'data/images/SAM_points_test.csv' -i data/images/test_ims/ -e 'JPG'

would classify the csv points and save the results in `results/SAM_points/`

3. To classify the segments, Export communitly composition shapefile, extract percent cover, and export collony level segments fill all the requiered fields in the RapidBenthos_part3.py:

    ```RB_centroid_csv = #input path to RapidBenthos part 1 polygon center point csv (hexagrid_RB_union_path_pts_csv)```

    ```RC_csv = #input path to ReefCloud annotation csv output```

    ```label_file = #input path to label file```
 
    ```polygon_file = #input path to RapidBenthos part 1 polygon shapefile (hexagrid_RB_union_path_seg_shp)```

    ```label_poygon_seg = #input path to RapidBenthos labeled polygon output shapefile```

    ```label_poygon_csv = #input path to RapidBenthos labeled polygon output CSV```

    ```PercentCover = #input path to RapidBenthos percent cover csv```

    ```out_fig = #input path to community composition stacked bar graph```

    ```ColonyLevelSegments_shp = #input path to RapidBenthos colony level segments shapefile```

    Once all field are populated run RapidBenthos_part3.py.

<a name="acknowledgements"></a>
## Acknowledgements

This research was supported by the Reef Restoration and Adaptation Program funded by the partnership between the Australian Government’s Reef Trust and the Great Barrier Reef Foundation, the College of Science and Engineering at James Cook University and the Remote Sensing Research Centre of the School of the Environment at the University of Queensland. We acknowledge the Wulgurukaba and Bindal people as the Traditional Owners of the land this research was conducted on. More broadly, we acknowledge the Traditional Owners of the Great Barrier Reef area and their continuing connection to their land and sea country.
