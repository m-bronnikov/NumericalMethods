import matplotlib as plt
import math
import sys



def f(x):
    return (1.0 / math.tan(x))


if __name__ == "__main__":
    Xi = [math.pi/8.0, 5.0*math.pi/16.0, 3*math.pi/8.0, math.pi/2.0]
    Xt = math.pi / 3.0
    Yi = list(map(f, Xi))
    filename = "data/d1.txt"
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) != 1:
        raise ValueError("Wrong number of args")
    
    with open(filename, "w") as f:
        f.write(" ".join(str(x) for x in Xi) + "\n")
        f.write(" ".join(str(y) for y in Yi) + "\n")
        f.write(str(Xt) + "\n")


