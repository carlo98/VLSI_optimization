from minizinc import Instance, Model, Solver
import matplotlib.pyplot as plt
import random
import matplotlib.patches as patches
import numpy as np
import datetime


def plot_result(result):
    dims = result["dims_v"]
    w = result["w_v"]
    ys = np.array(result["sol"])[:,1]
    xs = np.array(result["sol"])[:,0]

    lims = (0, max(w, max([ys[i]+dims[i][1] for i in range(len(ys))])))

    fig1 = plt.figure(figsize=(10, 10))
    ax1 = fig1.add_subplot(111, aspect='equal')
    for i in range(len(xs)):
        ax1.add_patch(patches.Rectangle((xs[i], ys[i]), dims[i][0], dims[i][1], edgecolor='black', facecolor=random.choice(['r', 'g', 'y'])))

    plt.xticks(range(w+1))
    plt.yticks(range(max([ys[i]+dims[i][1] for i in range(len(ys))])+1))
    plt.show()


def read_instance(instance_id):
    filepath = "../instances/ins-" + str(instance_id) + ".txt"
    with open(filepath, "r") as f_in:
        f = f_in.readlines()
        for i in range(len(f)):
            if not f[i][-1].isnumeric():
                f[i] = f[i][:-1]
                
        W = int(f[0])
        n = int(f[1])
        dims_line = f[2].split(" ")    
        dims = [[int(dims_line[0]), int(dims_line[1])]]
        for i in range(1, int(f[1])-1):
            dims_line = f[2 + i].split(" ")
            dims.append([int(dims_line[0]), int(dims_line[1])])
        dims_line = f[-1].split(" ")
        dims.append([int(dims_line[0]), int(dims_line[1])])
        
    return dims, W, n


def solve_instance(instance_id):
    # Load model from file
    model = Model("./model_rotation.mzn")
    print("Model Rotation")
    
    # Find the MiniZinc solver configuration for Gecode
    solver = Solver.lookup("gecode")

    # Assign data
    model.add_file("./instances/ins-" + str(instance_id) + ".dzn")

    # Create an Instance of model for the solver
    instance = Instance(solver, model)
    
    timeout = datetime.timedelta(seconds=300)
    return instance.solve(timeout=timeout)
    

def solve_all(max_instance):
    results = []
    times = []
    for i in range(1, max_instance+1):
        print("Solving: ", i)
        res = solve_instance(i)
        print("\tObj: ", res["objective"])
        print("\tTime: ", res.statistics["solveTime"].total_seconds())
        print("\tFailures: ", res.statistics["failures"])
        print("Solved: ", i)

        dims, W, n = read_instance(i)

        f = open("out_rotation/out-" + str(i) + ".txt", "w")
        f.write(str(W) + " " + str(res["objective"]) + "\n")
        f.write(str(n) + "\n")
        for i in range(n):
            if i < n-1:
                f.write(str(dims[i][0]) + " " + str(dims[i][1]) + " " + str(res["sol"][i][0]) + " " + str(res["sol"][i][1]) + "\n")
            else:
                f.write(str(dims[i][0]) + " " + str(dims[i][1]) + " " + str(res["sol"][i][0]) + " " + str(res["sol"][i][1]))
        f.close()


if __name__ == "__main__":
    solve_all(40)

