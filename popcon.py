import os
import argparse
import glob
import time
import random
import logging
import numpy as np
import pandas as pd
import multiprocessing as mp
import euclidean_distance
from scipy import integrate
from collections import namedtuple
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from MultiProcessingLog import MultiProcessingLog

logger = logging.getLogger()


def get_local_density(
        features, image_shape, downsample_factor,
        global_centroid_name_y, global_centroid_name_x,
        object_type, radius=100):
    """
    Returns a DataFrame containing the mapobject_id and the object
    density with column label "LocalDensity" + object_type + radius
    """

    resized_image_shape = (int(image_shape[0] / downsample_factor),
                           int(image_shape[1] / downsample_factor))
    logger.debug('resized_image_shape: {}'.format(resized_image_shape))

    dot_image = np.zeros(shape=resized_image_shape, dtype=np.float64)
    feature_name = "LocalDensity_" + object_type + "_" + str(radius)

    for row in features.itertuples():
        logger.debug('Generating dot_image for object {}'.format(row.Index))
        dot_image[
            int(getattr(row, global_centroid_name_y) / downsample_factor),
            int(getattr(row, global_centroid_name_x) / downsample_factor)] = 1

    logger.debug('Blurring dot_image with radius {}'.format(radius))
    gaussian_image = gaussian_filter(
        dot_image,
        sigma=int(radius / downsample_factor))

    d = []
    for row in features.itertuples():
        logger.debug('Measuring blurred image for object {} at radius {}'.format(row.Index, radius))
        d.append({
            feature_name: gaussian_image[
                int(getattr(row, global_centroid_name_y) / downsample_factor),
                int(getattr(row, global_centroid_name_x) / downsample_factor)],
            'mapobject_id': row.mapobject_id})
    return(pd.DataFrame(d))


def mean_euclidean_distance(x_i, y_i, x_max, y_max):
    """
    Returns the mean euclidean distance from a point (x_i,y_i) to
    all other points in the rectangle with vertices (0,0),(x_max,y_max)
    """

    d = integrate.dblquad(lambda x, y: euclidean_distance.metric(x,y,x_i,y_i),
                          0, y_max,
                          0, x_max)

    area = x_max * y_max
    return d[0] / area


def get_crowding(
        features, image_shape,
        global_centroid_name_y, global_centroid_name_x,
        object_type):
    """
    Returns a DataFrame containing the mapobject_id and the object
    crowding with column label "Crowding" + object_type
    """

    maximal_distance = np.sqrt(float(image_shape[0])**2 +
                               float(image_shape[1])**2)

    crowding_list = []
    for row in features.itertuples():
        logger.debug('Calculating crowding for cell {}'.format(row.Index))
        # find distances from current cell to all other cells
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

        # find expected mean distance to all other points in the domain
        expected_mean_distance = mean_euclidean_distance(
            getattr(row, global_centroid_name_x),
            getattr(row, global_centroid_name_y),
            int(image_shape[1]),
            int(image_shape[0]))

        # invert real distances to get crowding
        inverse_real_distances = np.divide(
            np.ones_like(real_distances),
            real_distances)

        # subtract the expected crowding and normalise by maximal distance
        crowding = maximal_distance * \
            (np.mean(inverse_real_distances) -
                (1.0 / expected_mean_distance))

        crowding_list.append({
            "Crowding_" + object_type: crowding,
            'mapobject_id': row.mapobject_id})

    return(pd.DataFrame(crowding_list))


# Note that this function is not currently called by main because it
# is not deterministic and generates large random errors due to the
# numerical 1/distance step.
#
# TODO: implement replicates and averaging for local crowding
def get_local_cell_crowding(
        features, image_shape,
        global_centroid_name_y, global_centroid_name_x,
        object_type, deterministic=False, reps=1):
    """
    Returns a DataFrame containing the mapobject_id and the object
    local crowding with column label "LCC" + object_type.
    """
    if (deterministic):
        random.seed(123)

    maximal_distance = np.sqrt(float(image_shape[0])**2 +
                               float(image_shape[1])**2)
    feature_name = "LCC_" + object_type

    lcc = []
    for row in features.itertuples():
        random_centroids_y = np.random.randint(
            0, image_shape[0], size=len(features) - 1)
        random_centroids_x = np.random.randint(
            0, image_shape[1], size=len(features) - 1)

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
            feature_name: (np.sum(inverse_real_distances) -
                           np.sum(inverse_random_distances)),
            'mapobject_id': row.mapobject_id})

    return(pd.DataFrame(lcc))


def convert_local_to_global_centroids(
        local_features, local_metadata, image_shape,
        centroid_name_y, centroid_name_x):

    logger.debug('Converting to global coordinates')
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


def calculate_edge(features, metadata,
        centroid_name_y, centroid_name_x,
        image_size_y=2160, image_size_x=2560,
        object_name="Nuclei",
        edge_expansion = 200,
        edge_site_range = 1):

    # TODO: Implement a downsample_factor
    # TODO: Add a global option => don't calculate per site, but for the whole
    # well at the same time. Would save computation, but be expensive in RAM
    # Edge detection: Run massive dilation based on centroids of objects
    # Using centroids, because loading masks is slow and cell segmentation is
    # not always available

    # Parameters:
    # How far out is an edge searched for each site? Defaults to 1 => looks at the 3x3 grid around the site

    edge_measurements = pd.DataFrame(columns = ['mapobject_id', 'DistanceToEdge', 'isEdge'])

    # Function gets the parameters for a whole well, calculates site by site to avoid huge memory usage
    metadata = metadata.assign(well_pos_combined = zip(metadata['well_pos_y'], metadata['well_pos_x']))
    existing_sites = list(metadata['well_pos_combined'].unique())

    for site in existing_sites:
        logger.info('Edge detection for site {}'.format(site))
        surrounding_sites = []
        min_y = site[0]
        min_x = site[1]
        max_y = site[0]
        max_x = site[1]
        for y in range(-edge_site_range,edge_site_range +1):
            for x in range(-edge_site_range,edge_site_range +1):
                potential_site = (site[0] + y, site[1] + x)
                if potential_site in existing_sites:
                    surrounding_sites.append(potential_site)
                    # Check if the current site is a new min or max in x or y direction
                    if potential_site[0] < min_y:
                        min_y = potential_site[0]
                    if potential_site[1] < min_x:
                        min_x = potential_site[1]
                    if potential_site[0] > max_y:
                        max_y = potential_site[0]
                    if potential_site[1] > max_x:
                        max_x = potential_site[1]

        # Create a binary numpy array
        surrounding_size = ((max_y - min_y + 1)*image_size_y, (max_x - min_x + 1)* image_size_x)
        local_surrounding = np.zeros(surrounding_size, dtype=bool)
        # Go through each image and set the centroids to True
        # Need to calculate the position in the local_surrounding image depending on where an image is relative to the others
        for sub_site in surrounding_sites:
            x_shift = (sub_site[1] - min_x) * image_size_x
            y_shift = (sub_site[0] - min_y) * image_size_y
            subsite_centroids = features.loc[metadata['well_pos_combined'] == sub_site]
            for row_index in range(subsite_centroids.shape[0]):
                y_pos = int(subsite_centroids[centroid_name_y].iloc[row_index] + y_shift)
                x_pos = int(subsite_centroids[centroid_name_x].iloc[row_index] + x_shift)
                # print(subsite_centroids[centroid_name_y].iloc[row_index], subsite_centroids[centroid_name_x].iloc[row_index])
                # print(y_pos, x_pos)
                local_surrounding[y_pos, x_pos] = True

        # Do a dilation of the points to fill in holes between cells
        from scipy.ndimage.morphology import binary_dilation, binary_closing, distance_transform_edt
        local_surrounding = binary_dilation(local_surrounding, iterations=edge_expansion)
        local_surrounding = binary_closing(local_surrounding, structure=np.ones((15,15)))

        # Distance transform currently treats edge as 0 => funny biases in dense regions
        local_surrounding = distance_transform_edt(local_surrounding)

        site_centroids = features.loc[metadata['well_pos_combined'] == site]
        site_centroids = site_centroids.assign(DistanceToEdge=0)
        site_centroids = site_centroids.assign(isEdge=0)
        x_shift = (site[1] - min_x) * image_size_x
        y_shift = (site[0] - min_y) * image_size_y
        for row_index in range(site_centroids.shape[0]):
            y_pos = int(site_centroids[centroid_name_y].iloc[row_index] + y_shift)
            x_pos = int(site_centroids[centroid_name_x].iloc[row_index] + x_shift)

            # Write values back to the data frame. For TissueMaps: Write to database here?
            col_index = site_centroids.columns.get_loc('DistanceToEdge')
            site_centroids.iloc[row_index, col_index] = max(local_surrounding[y_pos, x_pos] - edge_expansion, 0)

            col_index2 = site_centroids.columns.get_loc('isEdge')
            site_centroids.iloc[row_index, col_index2] = int((local_surrounding[y_pos, x_pos] - edge_expansion) <= 0)

        edge_measurements = edge_measurements.append(site_centroids[['mapobject_id', 'DistanceToEdge', 'isEdge']], ignore_index=False, verify_integrity=False, sort=None)

    return edge_measurements


def calculate_popcon_features(
        target_dir,
        features_filename, metadata_filename,
        centroid_name_y, centroid_name_x,
        image_size_y=2160, image_size_x=2560,
        downsample_factor=1,
        object_name="Cells",
        calculate_density=True,radii=[100,200,300],
        find_edge=True, edge_expansion = 200, edge_site_range = 1):

    logger.info('Calculating popcon features')
    logger.info('features_filename: %s', features_filename)
    logger.info('metadata_filename: %s', metadata_filename)
    logger.info('target_dir: %s', target_dir)
    logger.info('downsample_factor: %d', downsample_factor)
    logger.info('object_name: %s', object_name)
    logger.info('radii: {}'.format(radii))
    logger.info('centroid_name_y: %s',centroid_name_y)
    logger.info('centroid_name_x: %s',centroid_name_x)


    image_shape = (image_size_y, image_size_x)
    metadata = pd.read_csv(
        metadata_filename,
        usecols=['mapobject_id', 'well_pos_y', 'well_pos_x'])
    features = pd.read_csv(
        features_filename,
        usecols=['mapobject_id', centroid_name_y, centroid_name_x])

    if calculate_density:
        logger.debug('Extracting well_shape from metadata_file:')
        well_shape = (
            int(image_size_y * (1 + metadata.well_pos_y.max())),
            int(image_size_x * (1 + metadata.well_pos_x.max()))
        )
        logger.info('Well shape is {}'.format(well_shape))

        global_coordinates = convert_local_to_global_centroids(
            features, metadata, image_shape,
            centroid_name_y, centroid_name_x
        )

        crowding = get_crowding(
            global_coordinates.df, well_shape,
            global_coordinates.centroid_name_y,
            global_coordinates.centroid_name_x,
            object_name)

        logger.info('Finished crowding calculation for: {}'.format(features_filename))

        popcon_features = global_coordinates.df.merge(
            crowding,on='mapobject_id')

        #lcc = get_local_cell_crowding(
        #    global_coordinates.df, well_shape,
        #    global_coordinates.centroid_name_y,
        #    global_coordinates.centroid_name_x,
        #    object_name)

        #popcon_features = popcon_features.merge(
        #    lcc,on='mapobject_id')

        for radius in radii:
            logger.debug(
                'At Radius {}, calculating object density for: {}'.format(radius, features_filename)
            )
            local_density = get_local_density(
                global_coordinates.df,
                well_shape, downsample_factor,
                global_coordinates.centroid_name_y,
                global_coordinates.centroid_name_x,
                object_name,radius=int(radius))

            logger.info('At Radius {}, finished density calculation for: {}'.format(radius,features_filename))

            popcon_features = popcon_features.merge(
                local_density,on='mapobject_id')

        f = os.path.basename(features_filename).replace('feature-values','popcon')
        output_file = os.path.join(os.path.expanduser(target_dir),f)

        popcon_features.to_csv(output_file,index=False,index_label=False)

    if find_edge:
        logger.info('edge_expansion: %d', edge_expansion)
        logger.info('edge_site_range: %d', edge_site_range)

        edge_measurements = calculate_edge(features, metadata, centroid_name_y,
                                           centroid_name_x, image_size_y,
                                           image_size_x, object_name,
                                           edge_expansion, edge_site_range)
        f = os.path.basename(features_filename).replace('feature-values','edge')
        output_file = os.path.join(os.path.expanduser(target_dir),f)

        edge_measurements.to_csv(output_file,index=False,index_label=False)

    return


def calculate_popcon_features_star(args):
    return calculate_popcon_features(*args)


def main(args):

    features_files = [
        os.path.basename(full_path) for full_path in
        glob.glob(
            args.source_dir +
            '/*' +
            args.main_object +
            '_feature-values.csv')]
    metadata_files = [f.replace('feature-values','metadata')
                      for f in features_files]

    features_paths = [os.path.join(args.source_dir,f) for f in features_files]
    metadata_paths = [os.path.join(args.source_dir,f) for f in metadata_files]

    if (args.centroid_object == args.main_object):
        local_centroid_name_y = 'Morphology_Local_Centroid_y'
        local_centroid_name_x = 'Morphology_Local_Centroid_x'
    else:
        local_centroid_name_y = args.centroid_object + '_Morphology_Local_Centroid_y'
        local_centroid_name_x = args.centroid_object + '_Morphology_Local_Centroid_x'

    if not os.path.exists(args.target_dir):
        logger.info('Creating popcon results folder {}'.format(args.target_dir))
        os.makedirs(args.target_dir)

    # use a multi-processing pool to get the work done
    # this requires zipping the arguments and passing to a
    # "star" function

    zipped_args = zip(
        [args.target_dir for f in features_paths],
        features_paths,
        metadata_paths,
        [local_centroid_name_y for f in features_paths],
        [local_centroid_name_x for f in features_paths],
        [args.image_size_y for f in features_paths],
        [args.image_size_x for f in features_paths],
        [args.downsample for f in features_paths],
        [args.centroid_object for f in features_paths],
        [args.calculate_density for f in features_paths],
        [args.radii for f in features_paths],
        [args.find_edge for f in features_paths],
        [args.edge_expansion for f in features_paths],
        [args.edge_site_range for f in features_paths]
    )

    # find out how many cores are available
    if args.parallel:
        try:
            n_cpus = min(len(zipped_args), mp.cpu_count())
            logger.info('Starting a parallel pool with {} processes'.format(n_cpus))
        except NotImplementedError:
            n_cpus = 1
            logger.warning('Error in cpu_count, not executing parallel processes')
    else:
        n_cpus = 1
        logger.info('Not running in parallel mode')

    start_time = time.time()

    # start the pool of workers
    pool = mp.Pool(processes=n_cpus)
    pool.map(
        calculate_popcon_features_star,
        zipped_args
    )
    pool.close()
    pool.join()

    end_time = time.time()

    logger.info('Processing took {:0.2f} s'.format(end_time - start_time))

    return


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='popcon',
        description=('Generates a csv file containing'
                     ' population context features derived'
                     ' from feature-values and metadata csv files.'
                     )
    )
    parser.add_argument('--verbose', '-v', action='count')
    parser.add_argument('--parallel', action='store_true', default=False)
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
    parser.add_argument('--calculate_density', action='store_true',
                        help='Whether density measurements should be made')
    parser.add_argument('-r','--radii', nargs='*', default=[],
                        help='List of radii (pixels) for density calculations')
    parser.add_argument('-d','--downsample', default=2, type=int,
                        help='factor by which to downsample distances ' +
                        'to save memory/time in density calculation')
    parser.add_argument('--find_edge', action='store_true',
                        help='Whether edge detection should be performed')
    parser.add_argument('-exp','--edge_expansion', default=200, type=int,
                        help='Amount of expansion of centroids for edge' +
                        'detection (in pixels)')
    parser.add_argument('-er','--edge_site_range', default=1, type=int,
                        help='What radius of surrounding sites is used for' +
                        'edge detection. 1 equals 3x3 sites')

    return(parser.parse_args())


def setup_logger(args):
    global logger

    logfile = os.path.join(
        os.getcwd(),
        'popcon-' + time.strftime('%Y%m%d-%H%M%S') + '.log')
    mp_log = MultiProcessingLog(logfile, 'w', 0, 0)
    formatter = logging.Formatter(
        '%(asctime)s [%(thread)d] %(funcName)s %(levelname)s: %(message)s')
    mp_log.setFormatter(formatter)
    logger.addHandler(mp_log)

    logger.setLevel(logging.INFO)
    if args.verbose > 0:
        logger.setLevel(logging.DEBUG)
    print 'Logging to {} at level {}'.format(
        logfile,logging.getLevelName(logger.getEffectiveLevel()))
    return


if __name__ == "__main__":
    args = parse_arguments()
    setup_logger(args)
    mp.freeze_support()
    main(args)
