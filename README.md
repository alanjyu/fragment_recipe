# Supplementary material for "A recipe for exotic continental fragment formation: automated data analysis of rift models"

This is a repository of the ASPECT input parameter files and Python scripts used for the initial model setup and model analysis. 

The [overview.xlsx](https://github.com/alanjyu/fragment_recipe/blob/main/overview.xlsx) provides a list of tested models and their indices, as well as calculations for the initial geotherm and material constant under various strength parameters ($f$).

The Python scripts used to calculate the initial geotherm and yield strength of the continent are stored in [py/](https://github.com/alanjyu/fragment_recipe/tree/main/py). The ParaView Python script used to detect continental breakup and analyze the widths of continental fragment is stored in [pvpy/](https://github.com/alanjyu/fragment_recipe/tree/main/pvpy). The ASPECT input parameter file for the reference model, as well as its movie render and log files, is stored in [prm/](https://github.com/alanjyu/fragment_recipe/tree/main/prm).


## Tested environment

- ParaView 5.11
- Python 3.9


## Automating data analysis using ParaView and Python

ParaView is equipped with powerful visualization and data filtering tools. One of its most valuable features is the ability to trace your analytical workflows into Python code, which can then be automated using Python's scripting capabilities.

To use this feature in ParaView:

1. Navigate to **Tools** → **Start Trace**. You can skip rendering components if camera adjustments are unnecessary.
2. Perform your desired analytical procedures on a single output file.
3. Once done, go to **Tools** → **Stop Trace**. This will generate Python code capturing your workflow.
4. Implement a For loop to automate the process for all your output files. An example of this implementation is shown in [pv_analyze.py](https://github.com/alanjyu/fragment_recipe/blob/main/pvpy/pv_analyze.py).


To run this script in ParaView:

1. Open the Python Shell via **View** → **Python Shell**.
2. Click **Run Script** in the shell menu and select the Python script.


For further scripting guidance, refer to the official [ParaView Python documentation](https://www.paraview.org/paraview-docs/nightly/python/).
