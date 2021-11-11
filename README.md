# VLSI_optimization
Optimization for Very Large Scale Integration using Constraint Programming, SAT and SMT.

Project developed as part of the course of Combinatorial Decision Making and Optimization of the University of Bologna.

The instances are not provided, but you can create your own by following the structure explained in the following paragraph.

## Problem Description
Given a fixed-width plate and a list of rectangular objects, decide how to place them on the plate so that the height of the final plate is minimized. 
Initially, each circuit must be placed in a fixed orientation with respect to the others, therefore it cannot be rotated.

The input files are txt files consinsting of multiple lines of integers, as shown below:
- width of the plate
- number of objects
- x, y (horizontal and vertical dimensions for object 0)
- x, y (horizontal and vertical dimensions for object 1)
...

## Usage
In each folder there's a README.md file with a description of the steps required.
