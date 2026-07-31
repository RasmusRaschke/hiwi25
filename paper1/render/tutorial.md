# Redering with Manim

To render a video or a freeze frame on the HPC with manim, follow these steps:

1. Load ```anaconda3``` to install the ```manim``` package with ```pip``` in a suitable environment. Do not forget to replace the name of the environment in the job file (mine is called ```visual```).
2. Use ```manim.cfg``` to modify general options such as background color and default resolution.
3. Copy some ```magsphere.out```-output which you want to animate in this folder and add the filenames to the dictionary in ```main.py```. 
4. Edit ```main.py``` as indicated in the documentation.
5. Modify the ```manim``` line in ```job.sh```: If you want a test animation, use ```manim -pql main.py Anim```; for the high quality final animation, use ```manim -pqh main.py Anim```; for a freeze frame, use ```manim -s -r 3840,2160 --format=png main.py Anim```
6. Run the job with ```job.sh```.
7. The output will be saved in the ```media``` folder with the subdictionary indicating the kind of output and the resolution.
8. Ignore that the job technically fails, this is just because ```manim``` calles ```xdg-open``` at the end.