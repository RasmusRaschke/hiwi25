#!/bin/bash

g++ -o sim.out -I/usr/include/eigen3 sim.cpp

./sim.out < input_rolling_down
python3 visualizer.py --input input_rolling_down
python3 other_solution.py --input input_rolling_down

./sim.out < input_osc 
python3 visualizer.py --input input_osc
python3 other_solution.py --input input_osc