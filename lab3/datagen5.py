import matplotlib as plt
import math
import sys


if __name__ == "__main__":
    X0 = -2
    X1 = 2
    h = 0.5
    filename = "data/d5.txt"
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) != 1:
        raise ValueError("Wrong number of args")
    
    with open(filename, "w") as f:
        f.write(str(X0) + "\n")
        f.write(str(X1) + "\n")
        f.write(str(h) + "\n")
