# equilibrium-c -- Build Instructions for Windows

- Build Requirements

    eqc has the same requirements on windows as linux, you will need Python3 and Numpy, as well as the GNU C Compiler (GCC). Please make sure to install either the 32 bit GCC, or the 64 bit GCC, depending on what kind of system you have.

- Build Steps:

    1. Click the green "Code" button on the github page and select "Download ZIP"
    
    2. Unzip the resulting file in a directory of your choosing, for example, 'C:Program Files/eqc'
    
    3. Open powershell in this location and navigate to the 'source' directory, then run .\compile.ps1
    
    4. Check that the file libeqc.dll has been created correctly.

- Installation Steps:

    Once the code has been built, you will need to add the source directory to your Path and Pythonpath. This can be accomplished as follows:

    1.  Open the Start menu and search for "Environment Variables."
    
    2.  Click on "Edit the system environment variables."
    
    3.  In the System Properties window, click on the "Environment Variables" button.
    
    4.  Under "System variables," scroll down and find the "Path" variable.
    
    5.  Select the "Path" variable and click on the "Edit" button.
    
    6.  Click on the "New" button and enter the path to the directory that you want to add.
    
    7.  Click on the "OK" button to close the Edit Environment Variable window.
    
    8.  Under "System variables," scroll down and find the "PYTHONPATH" variable. If the "PYTHONPATH" variable does not exist, you will need to create it by clicking on the "New" button and entering a name of
    YTHONPATH" and a value of the directory that you want to add.
    
    9.  Select the "PYTHONPATH" variable and click on the "Edit" button.
    
    10. Click on the "New" button and enter the path to the directory that you want to add.
    
    11. Click on the "OK" button to close the Edit Environment Variable window.
    
    12. Click on the "OK" button to close the Environment Variables window.
    
    13. Click on the "OK" button to close the System Properties window.

