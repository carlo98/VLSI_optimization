# Constraint Programming VLSI
The following paragraphs show how the given problem was encoded in MiniZinc and the results obtained through these encodings.

## Problem Description, main constraints and objective function
The code is available in "model_0.mzn", the constraints used are those that describe the problem; given a 2Darray, called "dims", that represent the width and the length of each circuit and another 2Darray, called "sol", that represent the position of left-bottom corner of each rectangle, and "w" the width of the enclosing rectangle, we have:

1) The variables of "sol" are in the range 0-max(max\_width, max\_height), where max\_width = min(w, sum(c in CIRCUITS)(dims[c, 1])) and max\_height=sum(c in CIRCUITS)(dims[c, 2]);

2) The objective function to be minimized is the maximum of the array "heights", with heights[i] = sol[i, 2]+dims[i, 2];

3) For each pair of rectangles (c1, c2) the following constraints need to hold, in order to avoid intersections: sol[c1, 1]+dims[c1, 1]<=sol[c2, 1] \/ sol[c2, 1]+dims[c2, 1]<=sol[c1, 1] \/ sol[c1, 2]+dims[c1, 2]<=sol[c2, 2] \/ sol[c2, 2]+dims[c2, 2]<=sol[c1, 2];

4) Each rectangle c needs to stay inside the given width, so sol[c, 1]+dims[c, 1]<=w.

## Best Model
The code is available in "model_2.mzn", the paragraphs below describe the encodings of all the suggestions available on the project description.

### Horizontal and Vertical sum
We know that the for each circuit c sol[c, 1] should be between 0 and the maximum width showed above minus the length of the circuit, dims[c, 1]; the same can be said for the vertical coordinate, sol[c, 2], which should be between 0 and the maximum height, described in the above paragraph, minus its height, so to recap:

constraint forall(c in CIRCUITS)(sol\_tmp[c, 1] <= max_width - dims[per[c], 1]);
constraint forall(c in CIRCUITS)(sol\_tmp[c, 2] <= max_height - dims[per[c], 2]);

### Global constraints
One can rewrite the constraints used to avoid intersection using the global constraint "diffn", moreover two implicit constraints could be added using cumulative both on the horizontal and vertical axis, as shown below:

constraint diffn(sol\_hor, sol\_ver, dim\_hor, dim\_ver);  % Avoiding intersections
constraint cumulative(sol\_ver, dim\_ver, dim\_hor, max\_width);  % Cumulative on y
constraint cumulative(sol\_hor, dim\_hor, dim\_ver, obj);  % Cumulative on x

### Symmetry breaking constraint
This problem has various symmetries, each rectangle disposition could be flipped on the x or y axis and the position of each rectangle with equal shape in one or both axis could be changed with that of similar rectangles; in order to restrict the search space the following symmetry constraints have been used:

1) Lexicographic constraints on relative order between rectangles, in order to avoid flips on x or y axis;
2) Lexicographic constraints on both axis for rectangles with equal shapes;
3) Lexicographic constraints on a single axis for rectangles with only one of the two dimensions equal.

### Search
The best results have been obtained by sorting the circuits, based on their areas, from biggest to smallest and then searching with restarts, input order and random domains.

### Rotations
The model is provided in "model_rotation.mzn", the following changes have been made:
1) 
