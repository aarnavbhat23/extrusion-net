import numpy as np
import pandas as pd
import ast
from glob import glob
import os
import tifffile as tiff
from tqdm import tqdm
from .interactive_analysis_utils import zoomed_image

def interpolate(x_to_interp:np.array,y_to_interp:np.array, im_crop:np.array,radius:int,df_roi:pd.DataFrame):
    # Define the new x values
    new_x_values = np.arange(0, im_crop.shape[0], 1)

    interpolation_values = df_roi.frame.unique()
    # Perform the interpolation
    interpolated_values_unbleached_x = np.interp(new_x_values, interpolation_values,x_to_interp )
    interpolated_values_unbleached_y = np.interp(new_x_values, interpolation_values,y_to_interp )


    height, width = im_crop[0,...].shape
    mean_list_interp = []

    for (x,y,i) in zip(interpolated_values_unbleached_x,interpolated_values_unbleached_y,new_x_values):    
        y_indices, x_indices = np.ogrid[:height, :width]

        # Create a binary mask where the pixels inside the circle are True
        mask_unbleached = (x_indices - x)**2 + (y_indices - y)**2 <= radius**2

        # Use the mask to index into your image and extract the pixel values
        pixels_in_unbleached = im_crop[i,...][mask_unbleached]

        # Actualize the coordinate of the circles to be able to plot them in the next frame
        # compute the mean intensity of the pixels inside the circles
        mean_list_interp.append(np.mean(pixels_in_unbleached))
    return np.array(mean_list_interp)

def format_log(file_path:str) -> tuple:
    """
    Function to extract the relevant information from the log file

    Args:
    file_path (str): path to the log file

    Returns:
    tuple: a tuple containing the following information:
    - cell_coord: the coordinates of the cells analyzed (list of tuple)
    - bbox: the size of the zoomed image (int)
    - radius_unbleached: the radius of the unbleached circular ROI  (int)
    - radius_bleached: the radius of the bleached circular ROI (int)
    - im_path: the path to the image (str)
    """


    # open the log file
    with open(file_path, 'r') as log_file:
        content = log_file.read()

    # create a list of strings that can be used to extract information from the log file
    log = content.split('\n')
    # get the coordinates of the cells

    cell_coord = [x for x in log if 'Coordinates of cells analyzed' in x]
    cell_coord = cell_coord[0].split('[')[1].split(']')[0]
    # Convert the string to a tuple
    cell_coord = ast.literal_eval(cell_coord)

    # get the size of the zoomed image
    bbox  = [x for x in log if 'Size of bbox zoom' in x]
    bbox = int(bbox[0].split(' ')[-1])
    
    # get the frame actualization
    frame_actualization = [x for x in log if 'Frame actualization' in x]
    frame_actualization = int(frame_actualization[0].split(' ')[-1])

    # get the radius of the unbleached spot
    radius_unbleached = [x for x in log if 'Radius unbleach spot' in x]
    radius_unbleached = int(radius_unbleached[0].split(' ')[-1])
    
    # get the radius of the bleached spot
    radius_bleached = [x for x in log if 'Radius bleach spot' in x]
    radius_bleached = int(radius_bleached[0].split(' ')[-1])

    # get the image path
    im_path = [x for x in log if 'Loading image' in x]
    im_path = im_path[0].split(' ')[-1]
    
    return cell_coord, bbox, radius_unbleached, radius_bleached, im_path


def format_coord_df(file_path:str) -> tuple:
    """
    Function to format the ROI dataframe and interpolate the coordinates of the unbleached and bleached spots

    Args:
    file_path: str: path to the ROI dataframe

    Returns:
    df_roi: pd.DataFrame: dataframe containing the ROI information
    x_un: np.ndarray: x coordinate of the unbleached spot (interpolate), np.array shape: (n_frames, n_nuclei)
    y_un: np.ndarray: y coordinate of the unbleached spot (interpolated), np.array shape: (n_frames, n_nuclei)
    x_bleached: np.ndarray: x coordinate of the bleached spot (interpolated), np.array shape: (n_frames, n_nuclei)
    y_bleached: np.ndarray: y coordinate of the bleached spot (interpolated), np.array shape: (n_frames, n_nuclei)

    """
    df_roi = pd.read_csv(file_path)
    # get the coordinate of the cells in x and y for all frames

    df_roi['unbleached'] = df_roi['unbleached'].apply(ast.literal_eval)
    df_roi['bleached'] = df_roi['bleached'].apply(ast.literal_eval)

    shape = (len(df_roi.frame.unique()), len(df_roi.nucleus.unique()))
    x_un,y_un,x_bleached,y_bleached = list(map(lambda _: np.zeros(shape), range(4)))

    for nucleus in df_roi.nucleus.unique():
        df_r = df_roi[df_roi['nucleus'] == nucleus]
        x_un[:,nucleus] = [x[0] for x in df_r['unbleached']]
        y_un[:,nucleus] = [x[1] for x in df_r['unbleached']]
        x_bleached[:,nucleus] = [x[0] for x in df_r['bleached']]
        y_bleached[:,nucleus] = [x[1] for x in df_r['bleached']]

    return df_roi, x_un, y_un, x_bleached, y_bleached


def parse_files(file_path:str) -> tuple:
    """
    Parse files in the given directory and return the paths of the latest log file, ROI file, raw file,
    background intensity file, and phluorophore bleaching intensity file.

    Args:
        file_path (str): The directory path where the files are located.

    Returns:
        tuple: A tuple containing the paths of the latest log file, ROI file, raw file,
        background intensity file, and phluorophore bleaching intensity file.
    """
    # get the last log file
    list_log = glob(f'{file_path}/*.log')
    # sort the list by date
    list_log.sort(key=lambda x: os.path.getmtime(x))
    list_log = list_log[-1]

    # get the last ROI file
    list_ROI = glob(f'{file_path}/*_ROI_*.csv')
    list_ROI.sort(key=lambda x: os.path.getmtime(x))
    list_ROI = list_ROI[-1]

    # get the last raw file
    list_raw = glob(f'{file_path}/*_raw_*.csv')
    list_raw.sort(key=lambda x: os.path.getmtime(x))
    list_raw = list_raw[-1]

    # get the last background intensity file
    list_back = glob(f'{file_path}/*_mask_back_*.npy')
    list_back.sort(key=lambda x: os.path.getmtime(x))
    list_back = list_back[-1]

    # get the last phluorophore bleaching intensity file 
    list_bleach = glob(f'{file_path}/*_intensity_bleach_*.npy')
    list_bleach.sort(key=lambda x: os.path.getmtime(x))
    list_bleach = list_bleach[-1]

    return list_log, list_ROI, list_raw, list_back,list_bleach

def format_background(file_path:str) -> np.ndarray:
    """
    This function takes a file path as input and opens the last background file specified by the file path.
    It then calculates the average background intensity over all the circles for each frame.
    
    Args:
        file_path (str): The file path of the background file.
        
    Returns:
        background_int (np.ndarray): A 1D array containing the mean intensity for each frame.
    """
    # open the last background file
    background = np.load(file_path)
    # average the background intensity over all the circle for each frame
    background_int = np.median(background,axis=1)
    # this outputs a 1D array with the mean intensity for each frame 
    return background_int


def perform_analysis(file_path:str) -> pd.DataFrame:
    """
    Perform analysis on the data from the given file path.

    This function reads various data files, performs interpolation and mean_list_interp_bleached_flat calculations for each nucleus, and returns the results in a DataFrame.

    Args:
    file_path (str): The path to the directory containing the data files.

    Returns:

    pd.DataFrame: A DataFrame containing the mean_list_interp_bleached_flat values, interpolated values for bleached and unbleached regions, background intensity, unfrap cell values, nucleus number, and frame number for each nucleus and frame.
    """
    # parse_files
    list_log, list_ROI, list_raw, list_back,list_bleach = parse_files(file_path)

    # format the various df

    cell_coord, bbox, radius_unbleached, radius_bleached, im_path = format_log(list_log)
    background_int = format_background(list_back)
    df_roi, x_un, y_un, x_bleached, y_bleached = format_coord_df(list_ROI)

    ## fluorophore bleaching
    unfrap_cell = np.load(list_bleach)
    ## average the background intensity over all the cells for each frame to get I_bleaching(t)
    unfrap_cell = np.mean(unfrap_cell, axis=0)

    # interpolate the data
    im = tiff.imread(im_path)
    shape = (np.shape(im)[0], len(df_roi.nucleus.unique()))
    mean_list_interp_unbleached,mean_list_interp_bleached = list(map(lambda _: np.zeros(shape), range(2)))
    
    for nucleus in tqdm(df_roi.nucleus.unique()):
        if len(df_roi.nucleus.unique()) == 1:
            #only one nuclei anlysez thus cell_coord is not a list but a tuple
            c = [int(x) for x in cell_coord]
        else:
            c = [int(x) for x in cell_coord[nucleus]]

        _ = None
        im_crop = zoomed_image(_,im, c, bbox)
        mean_list_interp_unbleached[:,nucleus] = interpolate(x_un[:,nucleus],y_un[:,nucleus],im_crop,radius_unbleached,df_roi)
        mean_list_interp_bleached[:,nucleus] = interpolate(x_bleached[:,nucleus],y_bleached[:,nucleus],im_crop,radius_bleached,df_roi)
        
    mean_list_interp_bleached_flat = mean_list_interp_bleached.T.ravel()
    mean_list_interp_unbleached_flat = mean_list_interp_unbleached.T.ravel()

    # Create a DataFrame for easier readability
    df = pd.DataFrame({
        'mean_list_interp_bleached': mean_list_interp_bleached_flat,
        'mean_list_interp_unbleached': mean_list_interp_unbleached_flat,
        'background_int': np.repeat(background_int, len(df_roi.nucleus.unique())),
        'unfrap_cell': np.repeat(unfrap_cell, len(df_roi.nucleus.unique())),
        'nucleus': np.repeat(range(len(df_roi.nucleus.unique())),np.shape(im)[0]),
        'frame': np.tile(range(np.shape(im)[0]), len(df_roi.nucleus.unique()))
    })
    return  df

def concat_runs(pattern:list, condition:str,path:str = '../runs/',nuclei_to_exclude:pd.DataFrame='../scratchpad/data/nuclei_wrong.csv')->pd.DataFrame:
    """
    Concatenate all the runs from the path and return the dataframe

    Args:

    pattern: list
        List of strings to match the files
    
    condition: str
        Experimental condition to match the nuclei to exclude

    path: str
        Path to the runs

    nuclei_to_exclude: pd.DataFrame
        DataFrame containing the nuclei to exclude from the analysis

    
    """
    df_interp_d = []
    list_files = os.listdir(path)
    nuclei_to_exclude = pd.read_csv(nuclei_to_exclude)
    nuclei_condition = nuclei_to_exclude[condition].dropna(how='all').values.tolist()

    for l,i in enumerate(list_files):
        counter = 0
        for pat in pattern:
            if pat in i:
                counter +=1 
        if counter >=2:
            path_df = path+i
            try:
                print(f'Processing the file {i}')
                df_interp = perform_analysis(path_df)
                df_interp['nucleus'] = df_interp['nucleus'] + (10*l)
                df_interp['replicate'] = i
                df_interp = df_interp[~df_interp['nucleus'].isin(nuclei_condition)]
                df_interp_d.append(df_interp)
            except IndexError as e:
                print(f'Error, the file {i} could not be opened, might not contain data or is not analized yet.')
                continue

    df_interp_d = pd.concat(df_interp_d)
    return df_interp_d