# ceq -- Library for Equilibrium Chemistry Calculations

- Author: Nick Gibbons (n.gibbons(at)uq.edu.au)

- References:

    "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"\
    NASA Reference Publication 1311, October 1995\
    Sanford Gordon and Bonnie J. McBride

    "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"\
    NASA/TP - 2002-211556, September 2002\
    Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

- Build Requirements

    + python3
    + python3-numpy
    + gcc

- Build Instructions

    To build and install in '$HOME/ceq' type:\
    $ cd source\
    $ make all\
    $ make install


    Installation in a custom directory is accomplished by the INSTALL_DIR variable:\
    $ make install INSTALL_DIR='/path/to/install/dir'


    Once installed, add the location to your PYTHONPATH in your .bashrc file:\
    export PYTHONPATH=${PYTHONPATH}:'/path/to/install/dir'


- Use Instructions

    The code is accessed through a python script that handles reading of data, memory management, and initialisation. All of this data is then passed to a set of c routines that do most of the actual work. See examples in the "tests" directory.

```python
from numpy import array
import pyeq

species = ['CO2', 'CO', 'O2']
ceq = pyeq.EqCalculator(species)

X0 = array([1.0, 0.0, 0.0])
X1 = ceq.pt(p=10e3, T=2500.0, Xs0=X0)
print(X1)

>>> array([ 0.65897178,  0.22735215,  0.11367607])
```

- Licensing

    This program is MIT licensed, which allows you do to almost anything with it except pretend that you wrote it or mess with the license itself. See the mit.txt file for details.

