import matplotlib as plt
import math
import sys



def f(x):
    return (1.0 / math.tan(x))


if __name__ == "__main__":
    Xi = [1.0, 1.9, 2.8, 3.7, 4.6]
    Yi = [2.4142, 1.0818, 0.50953, 0.11836, -0.24008]
    Xt = 2.66666667
    filename = "data/d2.txt"
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) != 1:
        raise ValueError("Wrong number of args")
    
    with open(filename, "w") as f:
        f.write(" ".join(str(x) for x in Xi) + "\n")
        f.write(" ".join(str(y) for y in Yi) + "\n")
        f.write(str(Xt) + "\n")


