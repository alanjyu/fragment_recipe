"""
Description: Output continental fragment (terrane) size.
Created by: Alan J Yu
Environment: Paraview 5.11+; Python 3.9+

To run this script in Paraview, 
open the Python Shell via "View" -> "Python Shell",
then click "Run Script" in the shell menu.
"""

import os
import numpy as np
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


################################################# SETTINGS #################################################
ROOT_PATH = f'--' # change this to the ASPECT output folder
OUTPUT_PATH = f'--' # change this to your desired output folder
OUTPUT_NAME = r'results.csv'

OVERWRITE_SAVED_PROGRESS = True # set this to false if you want checkpoints
START_TIME, END_TIME = 25, 120
############################################################################################################


def process_timestep(i, j, model, step):
    """
    Load the ASPECT .pvtu file and check for CONTINENTAL_BREAKUP.
    """
    global CONTINENTAL_BREAKUP  # declare CONTINENTAL_BREAKUP as global

    formatted_timesteps = np.char.zfill(step.astype(str), 5)
    path = f'{ROOT_PATH}/{model}/solution/solution-{formatted_timesteps}.pvtu'

    # append the above path if data file exists
    if os.path.isfile(path):
        print(f'{model}: processing timestep {step}')
        src = XMLPartitionedUnstructuredGridReader(
            FileName=path, registrationName=f'{model}_{step}')

        if not CONTINENTAL_BREAKUP:
            # see if mantle has risen to the surface
            # query 1: astheno > 0.001
            # query 2: T is min (273K)
            SetActiveSource(src)

            QuerySelect(QueryString='(astheno > 0.001) & (T == min(T))',
                        FieldType='POINT', InsideOut=0)

            ast = ExtractSelection(registrationName=f'{model}_{step}_ast')
            ast_fetch = paraview.servermanager.Fetch(ast)
            number_of_points = ast_fetch.GetNumberOfPoints()

            if number_of_points > 0:
                CONTINENTAL_BREAKUP = True
                breakup_times[i] = step
                print(f'{model}: continent breaks apart at timestep {step}!')
                get_fragment_width(i, j, model, step, src)

            # delete ast
            Delete(ast)
            del ast

        # delete source after processing
        Delete(src)
        del src
    else:
        print(f'{model}: missing timestep {step}.')


def get_fragment_width(i, j, model, step, src):
    """
    Calculate the difference of x_min and x_max of the fragment selection.
    """
    # select the surficial points of the terrane
    SetActiveSource(src)  

    # Query 1: upper_2 >= 0.5
    # Query 2: T is min (273K)
    # Query 3: horizontal velocity > 0
    QuerySelect(
        QueryString='(upper_2 >= 0.5) & (T == min(T)) & (velocity[:,0] > 0)', FieldType='POINT', InsideOut=0)

    # extract the selected points
    terrane = ExtractSelection(registrationName=f'{model}_{step}_surf')

    SetActiveSource(terrane)

    # append location attributes to the selected points
    locAttr = AppendLocationAttributes(
        registrationName=f'{model}_{step}_attr', Input=terrane)

    SetActiveSource(locAttr)

    # calculate the terrane range and size based on the min and max of pos_x
    try:
        terrane_range = locAttr.PointData.GetArray(
            'PointLocations').GetRange()
        results_array[i][j] = abs(terrane_range[0]-terrane_range[1])
        fragment_sizes[i] = results_array[i][j]

        print(f'Fragment width: {fragment_sizes[i]}')

        export_as_csv(output_data)
    except AttributeError:
        print(f'{model}: cannot get range at timestep {step}')

    # delete extracted data after processing
    Delete(terrane)
    del terrane

    Delete(locAttr)
    del locAttr


def export_as_csv(output_data):
    np.savetxt(f'{OUTPUT_PATH}/{OUTPUT_NAME}',
               output_data,
               delimiter=',',
               fmt='%s')
    print(f'output: {OUTPUT_PATH}/{OUTPUT_NAME}')


CONTINENTAL_BREAKUP = False # initial state of continental breakup. Do not change.

timesteps = np.arange(START_TIME, END_TIME + 1, 1) # list of timesteps (int array)
savepath = f'{OUTPUT_PATH}/{OUTPUT_NAME}'

print(f'root path: {ROOT_PATH}')

# create the output folder if it doesn't exist
os.makedirs(OUTPUT_PATH, exist_ok=True)

# get a list of all models in the root folder
detected_models = []  # list of model names

for item in os.listdir(f'{ROOT_PATH}/'):
    # filter out unrelated files in the root folder
    if os.path.isdir(f'{ROOT_PATH}/{item}'):
        detected_models.append(item)
    
print(detected_models)


if os.path.isfile(savepath):
    savedata = np.genfromtxt(savepath, delimiter=',', skip_header=1, dtype=str)

    if not OVERWRITE_SAVED_PROGRESS:
        processed_models = savedata[:, 0]
        processed_widths = savedata[:, 1].astype(float)
        processed_times = savedata[:, 2].astype(float)

        unprocessed_models = np.setdiff1d(detected_models, processed_models)

        print("Processed Models:", processed_models)
        print("Unprocessed Models:", unprocessed_models)

output_data = []
output_data.append(['model', 'fragment width', 'breakup time'])  # append column headers for the newly created output data


# get the fragment width from each model
if not OVERWRITE_SAVED_PROGRESS:
    results_array = np.zeros((len(unprocessed_models), len(timesteps)))
    fragment_sizes = np.zeros(len(unprocessed_models))
    breakup_times = np.zeros(len(unprocessed_models))

    for i, model in enumerate(processed_models):
        output_data.append([model, processed_widths[i], processed_times[i]])

    for i, model in enumerate(unprocessed_models):
        CONTINENTAL_BREAKUP = False
        for j, step in enumerate(timesteps):
            process_timestep(i, j, model, step)

        output_data.append([model, fragment_sizes[i], breakup_times[i]])
elif OVERWRITE_SAVED_PROGRESS:
    results_array = np.zeros((len(detected_models), len(timesteps)))
    fragment_sizes = np.zeros(len(detected_models))
    breakup_times = np.zeros(len(detected_models))

    for i, model in enumerate(detected_models):
        CONTINENTAL_BREAKUP = False
        for j, step in enumerate(timesteps):
            process_timestep(i, j, model, step)

        output_data.append([model, fragment_sizes[i], breakup_times[i]])


# export data as csv
export_as_csv(output_data)
