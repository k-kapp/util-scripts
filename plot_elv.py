import argparse
import numpy as np
import matplotlib.pyplot as plt
import png

def read_png(filename):
    """
    Read a PNG file and return as an array. Note that this assumes the image
    is grayscale and that the data is contained in the first channel
    
    Parameters
    ---------------------
    filename : str
        PNG filename (full path)
    """

    w, h, data, info = png.Reader(filename).read()
    return np.array(list(data))[::info["planes"]]

def read_elv(filename, notherparams=0):
    """
    Read a file containing elevation or canopy height information.
    This file can be either in PNG format, ELV, or TXT format.
    Extension must be either .png or .txt

    Parameters
    ---------------------
    filename : str
        PNG, ELV, TXT filename
    notherparams: int
        If filename is a TXT or ELV file, then this indicates how many additional parameters from the file
        should be read. If notherparams = 1, then the step size is read, and if notherparams = 2, then
        the latitude is also read, in addition to step size
    """

    if filename.lower().endswith(".png"):
        return read_png(filename)
    if notherparams < 0 or notherparams > 2:
        raise ValueError("Invalid notherparams argument given to read_elv: must be between 0 and 2, inclusive")
    with open(filename, "r") as infile:
        arr = []
        for lnum, line in enumerate(infile):
            line = line.strip()
            line = line.split(" ")
            line = [el for el in line if len(el) > 0]
            if len(line) > 0:
                if lnum == 0:
                    width, height = int(line[0]), int(line[1])
                    if notherparams == 2:
                        try:
                            step, lat = float(line[2]), float(line[3])
                        except IndexError:
                            raise ValueError("file {} does not contain step and latitude parameters".format(filename))
                    elif notherparams == 1:
                        try:
                            step = float(line[2])
                        except IndexError:
                            raise ValueError("file {} does not contain a step parameter".format(filename))
                    elif notherparams == 0:
                        if len(line) > 2:
                            step = float(line[2])
                            print("NOTE: step parameter not requested, but present in file. step size: {}".format(step))
                else:
                    try:
                        line = [float(el) for el in line]
                    except ValueError:
                        for el in line:
                            try:
                                float(el)
                            except ValueError:
                                print("Can't convert '{}' to float".format(el))
                                raise
                        raise
                    #if lnum < height / 2:
                    #    line = [el for el in line]
                    arr.append(line)
        arr = np.array(arr).astype(np.float32)
        if arr.shape != (height, width):
            arr = arr.reshape((height, width))
    if notherparams == 2:
        return arr, step, lat
    elif notherparams == 1:
        return arr, step
    else:
        return arr

def write_elv(filename, arr, step=None, lat=None):
    """
    Write an array containing elevation or canopy height data to an ELV file.

    Parameters
    --------------------
    filename : str
        Filepath of file to write to.
    arr : np.array
        Numpy array containing data to be written
    step : float
        Size of each cell containing data (optional)
    lat : float
        Latitude parameter (optional)
    """

    h, w = arr.shape
    with open(filename, "w+") as outfile:
        if step is None and lat is None:
            wh_str = "{} {}\n".format(w, h)
        elif step is not None and lat is None:
            wh_str = "{} {} {}\n".format(w, h, step)
        elif step is not None and lat is not None:
            wh_str = "{} {} {} {}\n".format(w, h, step, lat)
        else:
            raise ValueError("lat cannot be given without step, due to positionality")
        outfile.write(wh_str)
        for rownum in range(h):
            rowvals = arr[rownum,:]
            vals_str = ["{} " for _ in range(w)]
            vals_str[-1] = vals_str[-1][:-1]
            vals_str += "\n"
            vals_str = "".join(vals_str)
            vals_str = vals_str.format(*list(rowvals))
            outfile.write(vals_str)

def transpose_elv(filename, outfilename):
    arr = read_elv(filename)
    arrt = arr.transpose()
    write_elv(outfilename, arrt)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("elv_filenames", type=str, nargs="+")
    arg_parser.add_argument("--titles", type=str, nargs="+")
    arg_parser.add_argument("--scales", type=str, nargs="+")

    a = arg_parser.parse_args()

    nfiles = len(a.elv_filenames)

    if a.scales is None:
        a.scales = []
    while len(a.scales) < nfiles * 2:
        a.scales.append("None")

    noname_count = 1
    if a.titles:
        titles = []
        for t in a.titles:
            titles.append(t)
        while len(titles) < nfiles:
            titles.append("Untitled" + str(noname_count))
            noname_count += 1

    scales = a.scales[:]

    for i in range(nfiles * 2):
        try:
            num = float(scales[i])
            scales[i] = num
        except ValueError:
            scales[i] = None

    arrs = [_ for _ in range(nfiles)]

    for i in range(nfiles):
        elv_filename = a.elv_filenames[i]
        arrs[i] = read_elv(elv_filename)

    fig, axes = plt.subplots(1, nfiles)

    if nfiles == 1:
        vmin, vmax = scales[0], scales[1]
        axes.imshow(arrs[0], vmin=vmin, vmax=vmax)
        if a.titles:
            axes.title.set_text(titles[0])
    else:
        for i in range(nfiles):
            vmin = scales[i * 2]
            vmax = scales[i * 2 + 1]
            axes[i].imshow(arrs[i], vmin=vmin, vmax=vmax)
            if a.titles:
            	axes[i].title.set_text(titles[i])
    plt.show()
