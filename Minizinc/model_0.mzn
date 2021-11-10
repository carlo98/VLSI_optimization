%% Inputs
int: w;  % Width of plate
int: n;  % Number of circuits
set of int: CIRCUITS = 1..n; 
array[CIRCUITS, 1..2] of int: dims;  % Dimensions xy for each circuit

%% Variables
int: max_width = min(w, sum(c in CIRCUITS)(dims[c, 1]));
int: max_height = sum(c in CIRCUITS)(dims[c, 2]);
array[CIRCUITS, 1..2] of var 0..max(max_width, max_height): sol;  % xy positions of left-bottom vertex for each circuit
array[CIRCUITS] of var int: heights;  % Height for each circuit

%% Constraints
constraint forall(c in CIRCUITS)(heights[c] == sol[c, 2]+dims[c, 2]);  % Consistence for heights
constraint forall(c1, c2 in CIRCUITS where c1!=c2)(avoid_intersection(sol, dims, c1, c2));  % Avoiding intersections
constraint forall(c in CIRCUITS)(sol[c, 1]+dims[c, 1]<=w);  % Forcing width

var int: obj = max(heights);

%% Search
ann: search_ann;
search_ann = int_search(sol, dom_w_deg, indomain_random);

%% Solve & output
solve
  :: search_ann
  minimize obj;
  
output [
  "obj: ", show(obj), "\n",
  "Solution: ", 
  show([sol[c, 1] | c in CIRCUITS]),
  show([sol[c, 2] | c in CIRCUITS])
];


%% Utilities
predicate avoid_intersection(array[CIRCUITS, 1..2] of var int: sol, array[CIRCUITS, 1..2] of int: dims, int: c1, int: c2) = 
    sol[c1, 1]+dims[c1, 1]<=sol[c2, 1] \/ sol[c2, 1]+dims[c2, 1]<=sol[c1, 1] \/ sol[c1, 2]+dims[c1, 2]<=sol[c2, 2] \/ sol[c2, 2]+dims[c2, 2]<=sol[c1, 2]
