import os
import argparse
import glob
import random
import numpy as np
import pandas as pd
import multiprocessing as mp
from math import sqrt, pow
from collections import namedtuple
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='popcon',
        description=('Generates a csv file containing'
                     ' population context features derived'
                     ' from feature-values and metadata csv files.'
                     )
    )
    parser.add_argument('source_dir', help='path to source directory')
    parser.add_argument('target_dir', help='path to destination directory')
    parser.add_argument('-m','--main_object', default="Cells", type=str,
                        help='main object which labels the feature files')
    parser.add_argument('-c','--centroid_object', default="Nuclei", type=str,
                        help='object from which centroids should be taken')
    parser.add_argument('-y','--image_size_y', default=2160, type=int,
                        help='x-dimension of image (pixels)')
    parser.add_argument('-x','--image_size_x', default=2560, type=int,
                        help='x-dimension of image (pixels)')
    parser.add_argument('-d','--downsample', default=2, type=int,
                        help='factor by which to downsample distances ' +
                        'to save memory/time in density calculation')

    return(parser.parse_args())


def get_local_density(
        features, image_shape, downsample_factor,
        global_centroid_name_y, global_centroid_name_x,
        object_type, radius=100):

    resized_image_shape = (int(image_shape[0] / downsample_factor),
                           int(image_shape[1] / downsample_factor))

    dot_image = np.zeros(shape=resized_image_shape, dtype=np.float64)
    feature_name = "LocalDensity_" + object_type + "_" + str(radius)

    # create a dot image of the features
    for row in features.itertuples():
        dot_image[
            int(getattr(row, global_centroid_name_y) / downsample_factor),
            int(getattr(row, global_centroid_name_x) / downsample_factor)] = 1

    # blur the dot image using a gaussian filter
    gaussian_image = gaussian_filter(
        dot_image,
        sigma=int(radius / downsample_factor))

    # measure the blurred_image
    d = []
    for row in features.itertuples():
        d.append({
            feature_name: gaussian_image[
                int(getattr(row, global_centroid_name_y) / downsample_factor),
                int(getattr(row, global_centroid_name_x) / downsample_factor)],
            'mapobject_id': row.mapobject_id})

    return(pd.DataFrame(d))


def get_local_crowding(
        features, image_shape,
        global_centroid_name_y, global_centroid_name_x,
        object_type, deterministic=False):

    if (deterministic):
        random.seed(123)

    maximal_distance = sqrt(
        pow(float(image_shape[0]),2) + pow(float(image_shape[1]),2))
    feature_name = "LocalCrowding_" + object_type

    lcc = []
    for row in features.itertuples():
        random_centroids_y = np.random.randint(
            0,image_shape[0],size=len(features)-1)
        random_centroids_x = np.random.randint(
            0,image_shape[1],size=len(features)-1)

        # find distances from current cell to random cells
        random_distances = cdist(
            np.column_stack(
                [getattr(row, global_centroid_name_y),
                 getattr(row, global_centroid_name_x)]
            ),
            np.column_stack(
                [random_centroids_y,
                 random_centroids_x]
            ),
            metric='euclidean')

        # find distances from current cell to real cells
        real_distances = cdist(
            np.column_stack(
                [getattr(row, global_centroid_name_y),
                 getattr(row, global_centroid_name_x)]
            ),
            np.column_stack(
                [getattr(features.drop(row.Index), global_centroid_name_y),
                 getattr(features.drop(row.Index), global_centroid_name_x)]
            ),
            metric='euclidean')

        # invert distances to get "crowding"
        inverse_random_distances = np.divide(
            np.full_like(random_distances, fill_value=maximal_distance),
            random_distances)
        inverse_real_distances = np.divide(
            np.full_like(real_distances, fill_value=maximal_distance),
            real_distances)

        lcc.append({
            feature_name:
                np.sum(inverse_real_distances) - np.sum(inverse_random_distances),
            'mapobject_id': row.mapobject_id})

    return(pd.DataFrame(lcc))


def convert_local_to_global_centroids(
        local_features, local_metadata, image_shape,
        centroid_name_y, centroid_name_x):

    global_centroids = namedtuple(
        'global_centroids', ['df','centroid_name_y','centroid_name_x']
    )

    global_centroids.centroid_name_y = centroid_name_y.replace(
        'Morphology_Local','Morphology_Global')
    global_centroids.centroid_name_x = centroid_name_x.replace(
        'Morphology_Local','Morphology_Global')

    local = local_features.merge(local_metadata, on='mapobject_id')
    d = []
    for row in local.itertuples():
        d.append({
            'mapobject_id': row.mapobject_id,
            global_centroids.centroid_name_y: getattr(row, centroid_name_y) +
                row.well_pos_y * image_shape[0],
            global_centroids.centroid_name_x: getattr(row, centroid_name_x) +
                row.well_pos_x * image_shape[1]
        })

    global_centroids.df = pd.DataFrame(d)
    return(global_centroids)


def calculate_popcon_features(
        target_dir,
        features_filename, metadata_filename,
        centroid_name_y, centroid_name_x,
        image_size_y=2560, image_size_x=2160,
        downsample_factor=1,
        object_name="Cells"):

    image_shape = (image_size_y, image_size_x)
    metadata = pd.read_csv(
        metadata_filename,
        usecols=['mapobject_id', 'well_pos_y', 'well_pos_x'])
    features = pd.read_csv(
        features_filename,
        usecols=['mapobject_id', centroid_name_y, centroid_name_x])

    well_shape = (
        int(image_size_y * (1 + metadata.well_pos_y.max())),
        int(image_size_x * (1 + metadata.well_pos_x.max()))
    )

    global_coordinates = convert_local_to_global_centroids(
        features, metadata, image_shape,
        centroid_name_y, centroid_name_x
    )

    local_cell_density_100 = get_local_density(
        global_coordinates.df,
        well_shape, downsample_factor,
        global_coordinates.centroid_name_y,
        global_coordinates.centroid_name_x,
        object_name,radius=100)
    local_cell_density_200 = get_local_density(
        global_coordinates.df,
        well_shape, downsample_factor,
        global_coordinates.centroid_name_y,
        global_coordinates.centroid_name_x,
        object_name,radius=200)
    local_cell_density_300 = get_local_density(
        global_coordinates.df,
        well_shape, downsample_factor,
        global_coordinates.centroid_name_y,
        global_coordinates.centroid_name_x,
        object_name,radius=300)
    local_cell_crowding = get_local_crowding(
        global_coordinates.df, well_shape,
        global_coordinates.centroid_name_y,
        global_coordinates.centroid_name_x,
        object_name,deterministic=True,reps=4)

    popcon_features = local_cell_density_100.merge(
        local_cell_density_200,on='mapobject_id')
    popcon_features = popcon_features.merge(
        local_cell_density_300,on='mapobject_id')
    popcon_features = popcon_features.merge(
        local_cell_crowding, on='mapobject_id')
    popcon_features = popcon_features.merge(
        global_coordinates.df, on='mapobject_id')

    f = os.path.basename(features_filename).replace('feature-values','popcon')
    output_file = os.path.join(os.path.expanduser(target_dir),f)

    popcon_features.to_csv(output_file,index=False,index_label=False)

    return


def calculate_popcon_features_star(args):
    return calculate_popcon_features(*args)


def main(args):

    features_files = [os.path.basename(full_path) for full_path in glob.glob(args.source_dir + '/*' + args.main_object + '_feature-values.csv')]
    metadata_files = [f.replace('feature-values','metadata') for f in features_files]

    features_paths = [os.path.join(args.source_dir,f) for f in features_files]
    metadata_paths = [os.path.join(args.source_dir,f) for f in metadata_files]

    if (args.centroid_object == args.main_object):
        local_centroid_name_y = 'Morphology_Local_Centroid_y'
        local_centroid_name_x = 'Morphology_Local_Centroid_x'
    else:
        local_centroid_name_y = args.centroid_object + '_Morphology_Local_Centroid_y'
        local_centroid_name_x = args.centroid_object + '_Morphology_Local_Centroid_x'

    zipped_args = zip(
        [args.target_dir for f in features_paths],
        features_paths,
        metadata_paths,
        [local_centroid_name_y for f in features_paths],
        [local_centroid_name_x for f in features_paths],
        [args.image_size_y for f in features_paths],
        [args.image_size_x for f in features_paths],
        [args.downsample for f in features_paths],
        [args.centroid_object for f in features_paths]
    )

    # use a multi-processing pool to get the work done
    pool = mp.Pool()
    pool.map(
        calculate_popcon_features_star,
        zipped_args
    )
    pool.close()
    pool.join()

    return

if __name__ == "__main__":
    args = parse_arguments()
    mp.freeze_support()
    main(args)
