import numpy as np

transition_radii = np.array([1233.0, 3482.0, 5702.0, 6371.0])
layer_densities = np.array([13077.0, 12772.0, 12142.0, 9941.0, 5589.0, 4443.0])
# layer_densities = np.array([12772., 9941., 4443., 1060.])


def Do():
    file = open("/home/os18o068/Documents/PHD/Projects/Planets/Data/PREM_Data.txt")
    data = [[], []]

    # Read in density data from file
    for line in file:
        dat = line.split(";")
        dat[0] = dat[0].replace(",", ".")
        dat[1] = dat[1].replace(",", ".")
        data[0].append(float(dat[0]))
        data[1].append(float(dat[1]))

        # Extract layer transitions

    return np.asarray(data)
