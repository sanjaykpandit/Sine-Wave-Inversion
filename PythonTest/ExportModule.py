import os
from shutil import copyfile

copy_from = r"C:\Users\Sanjay Pandit\Documents\Study\PhD\Resistivity Instrument\Cpp\CosFit\x64\Release\CosFit.pyd"
copy_to = [r'C:\Users\Sanjay Pandit\Documents\Study\PhD\Resistivity Instrument\PythonApp\ImpedanceMeter\model\inversion\CosFit.pyd',
           r'C:\Users\Sanjay Pandit\Documents\Study\PhD\Resistivity Instrument\PythonApp\ImpedanceMeter\model\inversion\CosFit.pyd']
for c in copy_to:
    copyfile(copy_from, c)