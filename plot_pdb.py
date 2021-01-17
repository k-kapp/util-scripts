import matplotlib.pyplot as plt
from read_pdb import pdb_content
import argparse
import png
import numpy as np
import copy

def get_chm_data(chm_filename, plot_rect):
    """
    Get CHM data from file

    Parameters
    ----------------
    chm_filename : str
        CHM filename
    plot_rect : list
        A list of floats describing rectangle in which to include trees, in
        the format [xmin, ymin, xmax, ymax]
    """
    if chm_filename is not None:
        if chm_filename.endswith(".png"):
            w, h, chm_data, info = png.Reader(chm_filename).read()
            chm_data = np.array(list(chm_data))
            chm_data = chm_data[:,::info["planes"]]
        else:
            with open(chm_filename, "r") as infile:
                lines = []
                for lnum, line in enumerate(infile):
                    line = line.strip()
                    if len(line) > 0:
                        line = line.split(" ")
                        if lnum == 0:
                            w, h = int(line[0]), int(line[1])
                        else:
                            line = [float(el) for el in line]
                            lines.append(line)
                chm_data = np.array(lines)
        if plot_rect is not None:
            chm_data = chm_data[plot_rect[1] : plot_rect[3], plot_rect[0] : plot_rect[2]]

    else:
        chm_data = None
    return chm_data


def process_titles(titlesarg, npdbs):
    """
    Ensure that each plot has a title (by default it assigns 'Untitled<int>')

    Parameters
    ----------------------
    titlesarg : list
        A list of titles passed by user as argument to program
    """
    notitle_count = 1
    titles = None
    if titlesarg:
        titles = copy.copy(titlesarg)
        while len(titlesarg) < npdbs:
            titles.append("Untitled" + str(notitle_count))
    return titles

def process_colors(color_list, nspecies):
    """
    Create a valid color list, i.e. so that each species has a corresponding color

    Parameters
    -----------------
    color_list : list
        Existing list of colors, passed by user, if any
    nspecies : int
        Number of species
    """
    color_list = copy.deepcopy(color_list)
    if color_list is not None and nspecies != len(color_list):
        color_list = [color_list[i % len(color_list)] for i in range(nspecies)]
    elif color_list is None:
        color_list = [(np.random.random(), np.random.random(), np.random.random()) for _ in range(nspecies)]
    return color_list


class PdbData:
    """
    Class that imports, processes and stores relevant info on canopy trees

    Parameters
    ------------------
    pdb_filenames : list
        List of PDB filenames
    plot_rect : list
        A list of floats describing a rectangle in the format [xmin, ymin, xmax, ymax]
    """
    def __init__(self, pdb_filenames, plot_rect):
        pdbs = self._get_pdb_data(pdb_filenames)
        self._process_data(pdbs, plot_rect)


    def _get_pdb_data(self, pdb_filenames):
        """
        Import pdb data and compute some basic stats, like the number of pdbs, list of unique species, etc.

        Parameters
        ------------------
        pdb_filenames : Refer to class docstring
        """
        pdbs = [pdb_content(pdbfname, normalize_distances=True) for pdbfname in pdb_filenames]
        self.npdbs = len(pdbs)
        self.allspecies = set()
        allspecies_list = [pdb.get_species_ids() for pdb in pdbs]
        for allsp in allspecies_list:
            self.allspecies = self.allspecies.union(allsp)
        self.allspecies = list(self.allspecies)
        self.nspecies = len(self.allspecies)
        return pdbs


    def _process_data(self, pdbs, plot_rect):
        """
        Create plot-friendly data from PDB file contents

        Parameters
        ---------------
        pdbs : list
            List of pdb_content objects, which summarise pdb file contents
        plot_rect : refer to this class' docstring for info on this param
        """
        nplants = 0
        maxr = 0
        maxh = 0
        self.radii = []
        self.heights = []
        self.plnt_xs, self.plnt_ys = [ [] for _ in range(self.npdbs)], [ [] for _ in range(self.npdbs)]
        self.plnt_colors = [ [] for _ in range(self.npdbs)]
        if plot_rect is None:
            plot_rect = [-float("inf"), -float("inf"), float("inf"), float("inf")]
        for pdbidx, pdb in enumerate(pdbs):
            for plnt in pdb.all_plants():
                if plnt.x < plot_rect[0] or plnt.x > plot_rect[2] or plnt.y < plot_rect[1] or plnt.y > plot_rect[3]:
                    continue
                self.plnt_xs[pdbidx].append(plnt.x)
                self.plnt_ys[pdbidx].append(plnt.y)
                self.plnt_colors[pdbidx].append(plnt.specie)
                if plnt.r > maxr:
                    maxr = plnt.r
                if plnt.h > maxh:
                    maxh = plnt.h
                self.radii.append(plnt.r)
                self.heights.append(plnt.h)
                nplants += 1
        print("max radius: {}".format(maxr))
        print("max height: {}".format(maxh))
        print("Number of plants: {}".format(nplants))
        
def plot_stats(pdb_data):
    """
    Plot canopy tree radiuses or heights, as requested

    Parameters
    -----------------
    pdb_data : PdbData
        PdbData instance containing relevant canopy tree info
    """
    plt.hist(pdb_data.radii, rwidth=0.9)
    plt.title("Plant radiuses")
    plt.show()
    plt.hist(pdb_data.heights, rwidth=0.9)
    plt.title("Plant heights")
    plt.show()

def plot(chm_data, pdb_data, titles, color_list, plot_rect, marker_size):
    """
    Plot desired data

    Parameters
    -------------------
    chm_data : np.array
        Canopy height data as a numpy array
    pdb_data : PdbData
        PdbData class instance containing relevant data for canopy trees
    titles : list
        list of titles
    color_list : list
        Valid list of colors (each species has a color)
    plot_rect : list
        Refer to PdbData class docstring for info on this param
    marker_size : int
        Size of dot that will be used to represent each tree
    """
    marker_size = marker_size if marker_size is not None else 1
    color_list = process_colors(color_list, pdb_data.nspecies)
    specie_colors = {spec_id: color for spec_id, color in zip(pdb_data.allspecies, color_list)}
    npdbs = pdb_data.npdbs

    fig, ax = plt.subplots(1, npdbs)

    for idx in range(npdbs):
        colarr = np.array(pdb_data.plnt_colors[idx])
        unique, counts = np.unique(colarr, return_counts=True)
        for v, c in zip(unique, counts):
            print("Species {}, count: {}".format(v, c))
        plnt_plot_colors = [specie_colors[spec_id] for spec_id in pdb_data.plnt_colors[idx]]

        if chm_data is not None:
            extent = None
            if plot_rect is not None:
                extent = [plot_rect[0], plot_rect[2], plot_rect[3], plot_rect[1]]
            if npdbs == 1:
                ax.imshow(chm_data, extent=extent)
            else:
                ax[idx].imshow(chm_data, extent=extent)
        if npdbs == 1:
            ax.scatter(pdb_data.plnt_xs[idx], pdb_data.plnt_ys[idx], s=marker_size, c=plnt_plot_colors)
        else:
            ax[idx].scatter(pdb_data.plnt_xs[idx], pdb_data.plnt_ys[idx], s=marker_size, c=plnt_plot_colors)
        if chm_data is None:
            if npdbs == 1:
                ax.set_ylim(ax.get_ylim()[::-1])
            else:
                ax[idx].set_ylim(ax[idx].get_ylim()[::-1])
        if titles:
            if npdbs == 1:
                ax.title.set_text(titles[0])
            else:
                for i in range(npdbs):
                    ax[i].title.set_text(titles[i])
    plt.show()


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("pdb_filenames", nargs="+")
    arg_parser.add_argument("--chm_filename")
    arg_parser.add_argument("--marker_size", type=int)
    arg_parser.add_argument("--color_list", nargs="+", type=str)
    arg_parser.add_argument("--plot_stats", action="store_true")
    arg_parser.add_argument("--titles", type=str, nargs="+")
    arg_parser.add_argument("--plot_rect", type=int, nargs=4, help="Rectangle in which plants are allowed: X1 Y1 X2 Y2")
    a = arg_parser.parse_args()


    chm_data = get_chm_data(a.chm_filename, a.plot_rect)
    pdb_data = PdbData(pdb_filenames=a.pdb_filenames, 
                        plot_rect=a.plot_rect)
    titles = process_titles(a.titles, pdb_data.npdbs)

    if a.plot_stats:
        plot_stats(pdb_data)

    print(len(pdb_data.allspecies))
    print(pdb_data.allspecies)


    plot(chm_data, pdb_data, titles, a.color_list, a.plot_rect, a.marker_size)
