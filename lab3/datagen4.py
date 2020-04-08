import math
import sys


if __name__ == "__main__":
    Xi = [0.0, 0.1, 0.2, 0.3, 0.4]
    Yi = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]
    Xt = 0.2
    n = 2
    filename = "data/d4.txt"
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) != 1:
        raise ValueError("Wrong number of args")
    
    with open(filename, "w") as f:
        f.write(" ".join(str(x) for x in Xi) + "\n")
        f.write(" ".join(str(y) for y in Yi) + "\n")
        f.write(str(Xt) + "\n")
        f.write(str(n) + "\n")