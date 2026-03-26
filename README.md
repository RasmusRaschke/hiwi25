# Welcome!
If you wonder where to find specific things, this readme should make things clear. I will go shortly over each of the folders.
## [Code](code/)
In this place, you find all the code I have written so far. 
- The [archive](code/archive/) contains old code not needed in the forseeable future.
- The folder [cpp](code/cpp/) contains the up-to-date code. You can find [source code](code/cpp/src/) as well as many [examples](code/cpp/simulations/) in there. The source code is structured as follows: The subdirectory [magsphere](code/cpp/src/magsphere/) contains the main code of this project, written in C++. You can download the source code and compile it yourself with ```g++ -o magsphere.out -I/usr/include/eigen3 magsphere.cpp``` (you may need to adapt the Eigen path). Then, modify the [input file](code/cpp/src/magsphere/input) to your liking and execute the simulation by ```./magsphere.out < input```.
To analyse the data, you find the (well-documented) python package [simutils](code/cpp/src/simutils/). To install it, simply go to the [parent directory](code/cpp/src) and execute ```pip install .```. Since the package uses numpy-docstrings, you can call help() on any function to obtain extensive information.
- The [animation directory](code/cpp/src/animation/) is still work in progress. It may contain code to animate the output of the main program in [manim](https://www.manim.community/) one day.

## [Inkscape](inkscape/)
Here you can find the images I created with inkscape. Those are mostly overview images of the system's physics.

## [Script](script/)
This folder contains the complete theoretical background and all my derivations for the main system.

## [Topology](topology/)
This contains a quick detour I did about the topology of the [sine surface](https://www.mathcurve.com/surfaces.gb/sinus/sinus.shtml) and the [Seifert-van Kampen Theorem](https://en.wikipedia.org/wiki/Seifert%E2%80%93Van_Kampen_theorem).