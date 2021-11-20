# PhaseFinder

## Installation and Setup:

### S4 for Windows:
- Install Miniconda or Anaconda [here](https://docs.conda.io/projects/conda/en/4.6.1/user-guide/install/windows.html)
- Open the Anaconda Prompt and run `where conda` to find the path to conda (ex. `C:\Users\USERNAME\Miniconda3\Scripts\conda.exe`)
- Add to PATH
    - Edit environmental variables
    - Put as first in the list to prevent potential conflicts
- Open up Command Prompt and run the following:
    - `conda init cmd.exe`
    - `conda config --add channels conda-forge`
    - `conda config --set ssl_verify no`
    - `conda upgrade conda`
    - if you already have python installed, run `conda create -n s4 python`. otherwise, run `conda create -n s4 python=3.7`
    - `conda activate s4`
    - `conda install -c paulgoulain s4`
- Test the installation by running:
    - `python`
    - `import S4 as S4` (capitalization matters)
    - if nothing happens it was installed right
    - Ctrl + Z to exit
- Install necessary python packages by running `pip install numpy matplotlib scipy tqdm`

### S4 MacOS* (Note: the python version of S4 doesn't work on M1 Macs yet):
- Install Miniconda or Anaconda [here](https://docs.conda.io/projects/conda/en/4.6.1/user-guide/install/macos.html)
- Open up Terminal and run the following:
    - `conda config --add channels conda-forge`
    - `conda upgrade conda`
    - if you already have python installed, run `conda create -n s4 python`. otherwise, run `conda create -n s4 python=3.7`
    - `conda activate s4`
    - `conda install -c paulgoulain s4`
- Test the installation by running:
    - `python`
    - `import S4 as S4` (capitalization matters)
    - if nothing happens it was installed right
    - ^Z to exit
- Install necessary python packages by running `pip3 install numpy matplotlib scipy tqdm`

<!-- ### MacOS Building From Scratch (for M1 Macbooks)
- Install the Command Line Tools: `xcode-select --install`
- Install Homebrew (package manager): `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`
- `git clone https://github.com/phoebe-p/S4`
- `cd S4`
- `brew install openblas`
- `brew install fftw`
- `brew install suite-sparse`
- `brew install wget`
- `make boost`
- `pip3 install numpy`
- edit `# CHOLMOD_INC = -I/usr/local/include` -> `CHOLMOD_INC = -I/opt/homebrew/opt/suite-sparse/include` if needed
- edit "S4.h" and move Simulation_SaveSolution and Simulation_LoadSolution out of Internal functions section
- edit `CFLAGS = -Wall -O3 -m64 -mcpu=lightning -mtune=native -msse3 -msse2 -msse -fPIC`
- `make -f Makefile.osx s4_pyext` -->

## Usage
- start by verifying that you are in the s4 environment (run `conda activate s4`)
- run the program by typing `python main.py` (Windows) or `python3 main.py` (MacOS)
- input requested information
- if any information is incorrect, terminate the program and start over (Ctrl + C or ^C)
- all data files will be saved to the specified data directory
    - data directory will be in the phase-shift directory
- selected parameters will be saved in `R_selected_parameters.txt` or ``T_selected_parameters.txt`` in the data directory (R for reflection, T for transmission)
    - for reference, raw data will be in a separate folder inside the directory too (to reduce clutter)
- simulation time for set intervals will be saved in `time_logs.txt` in the logs directory (for reference)
- at the end of the program, the graphs of selected parameter values will be generated using matplotlib
    - each set of parameters will produce a separate graph
    - use Ctrl + W or Cmd + W to close current graph
    - next graph will automatically pop up once previous graph is closed

## Folders and Files [not updated (yet)]:
- [archived](https://github.com/aaronzr1/phase-shift/tree/master/archived): sample simulation files
    - [lua_simu14](https://github.com/aaronzr1/phase-shift/tree/master/archived/lua_simu14): sample lua simulation files
    - [python_simu14](https://github.com/aaronzr1/phase-shift/tree/master/archived/python_simu14): sample python simulation files
- [logs](https://github.com/aaronzr1/phase-shift/tree/master/logs): miscellaneous log files
    - [time_logs.txt](https://github.com/aaronzr1/phase-shift/tree/master/logs/time_logs.txt): records time taken for simulation (only most recent info is saved)
    - [compare.py](https://github.com/aaronzr1/phase-shift/tree/master/logs/compare.py): simple python script that compares 2 txt files
- [data](https://github.com/aaronzr1/phase-shift/tree/master/data): contains data files generated from simulation code (directory name is customizable)

## Misc Implementation Details:
- E<sub>reflected</sub> = E<sub>Z<sub>0</sub></sub> - E<sub>incident</sub>
- Phase Shift: atan2(ratioI, ratioR)
    - ratioI = (reference/reflected).imag (imaginary component)
    - ratioR = (reference/reflected).real (real component)
- frequency = 1000/wavelength
- np.linspace was used instead of something like np.arange to prevent some float precision errors

## TODO:
- add input details
- progress bar
- add alternate method: run from command line and use sys.argv
- more detailed installation guide
