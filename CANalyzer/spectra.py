"""
This is meant to be called as a library from a Jupyter notebook.
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import subprocess
import linecache
from itertools import chain

import CANalyzer.utilities

class Spectra():
    def __init__(self):
        self.figheight = 30
        self.figwidth = 40
        self.axiswidth = 10
        self.majorxtick_fontsize = 75
        self.majorytick_fontsize = 0
        self.tickpad = 15
        self.legend_font = 80
        self.legend_ncol = 1
        self.legend_colspacing = 1
        self.x_decimal = 0
        self.y_decimal = 0
        self.xbins = None
        self.ybins = None
        self.display_xlim = None
        self.display_ylim = None
        self.linewidth = 10
        self.colors = list(mcolors.CSS4_COLORS.keys())
        self.linestyles = ["-"]
        self.labels = []
        self.legend_loc = "upper left"
        self.xlabel = "Energy (eV)"
        self.ylabel = None
        self.axislabel_font = 90
        self.alpha = 0.5

        self.redshift = 0
        self.spectra_yscale = 1
        self.sticks_yscale = 1
        self.npoints = 5000
        self.spectra = {}
        self.exp_spectra = None
        self.thresh = 1E-6

    def gaussian(self, energy, broadening, osc_str, xspace):
        ones = np.ones(xspace.shape)
        return (osc_str / (np.sqrt(2 * np.pi * broadening ** 2))) * np.exp(
            -(xspace - ones * energy) ** 2 / (2 * broadening ** 2))


    def lorentzian(self, energy, broadening, osc_str, xspace):
        ones = np.ones(xspace.shape)
        d = np.pi * (xspace - ones * energy) ** 2 + broadening ** 2
        lorentzian = broadening / d
        return osc_str * lorentzian


    def make_spectrum(self, title, xstart, xend, roots, os, broadening, lineshape=None):
        line_list = []
        xspace = np.linspace(xstart, xend, self.npoints)
        os_list = os
        for i in range(os_list.shape[0]):
            root = roots[i]
            oscstr = os[i]
            if abs(oscstr) < self.thresh:
                continue
            if root >= xstart and root <= xend:
                if lineshape:
                    line_list.append(lineshape(root, broadening, oscstr, xspace))
                else:
                    line_list.append(self.lorentzian(root, broadening, oscstr, xspace))
        spectrum = np.zeros(self.npoints)
        for i in line_list:
            spectrum += i
        self.spectra[title] = spectrum
        self.labels.append(title)


    def read_exp_spectrum(self, csvfile, xlabel, ylabel):
        exp_csv = pd.read_csv(csvfile)
        exp_x = exp_csv[xlabel]
        exp_y = exp_csv[ylabel]
        self.exp_spectra = (exp_x, exp_y)


    def plot(self, spectra_names, xstart, xend, plotname=None, show=False, vline=None):
        plt.rc('axes', linewidth=self.axiswidth)
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"

        fig, ax = plt.subplots()
        fig.set_figwidth(self.figwidth)
        fig.set_figheight(self.figheight)
        fig.canvas.draw()

        xspace = np.linspace(xstart, xend, self.npoints)

        ncolors = len(self.colors)
        nlinestyles = len(self.linestyles)
        for i in range(len(spectra_names)):
            s = spectra_names[i]
            y = self.spectra[s]
            ax.plot(xspace - self.redshift, y * self.spectra_yscale, label=s, alpha=self.alpha, linewidth=self.linewidth,
                    color=self.colors[i%ncolors], linestyle=self.linestyles[i%nlinestyles])

        if self.exp_spectra:
            ax.plot(self.exp_spectra[0], self.exp_spectra[1], label="Exp", alpha=self.alpha, linewidth=self.linewidth,
                    color="black", linestyle='--')

        if vline and self.display_ylim:
            try:
                try:
                    for x in vline:
                        ax.vlines(x, self.display_ylim[0]*self.spectra_yscale, self.display_ylim[1]*self.spectra_yscale, colors="gray", linestyle='--', linewidth=self.linewidth)
                except:
                    ax.vlines(vline, self.display_ylim[0]*self.spectra_yscale, self.display_ylim[1]*self.spectra_yscale, colors="gray", linestyle='--', linewidth=self.linewidth)
            except:
                print("Invalid argument for vline. Must be float, int, or list of floats ort ints.")
                pass

        ax.set_xlabel(self.xlabel, fontsize=self.axislabel_font, fontweight='bold')
        if self.ylabel:
            ax.set_ylabel(self.ylabel, fontsize=self.axislabel_font, fontweight='bold')

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(self.majorxtick_fontsize)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(self.majorytick_fontsize)
            tick.label1.set_fontweight('bold')
        ax.tick_params(axis='both', which='major', pad=self.tickpad)
        ax.legend(loc=self.legend_loc, prop={'weight': 'bold', 'size': self.legend_font}, frameon=False,
                  ncol=self.legend_ncol, columnspacing=self.legend_colspacing)
        #ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(f'{x:,.2f}'))
        #ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(f'{x:,.2f}'))
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(f'%.{self.y_decimal}f'))
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(f'%.{self.x_decimal}f'))

        if self.xbins:
            plt.locator_params(axis='x', nbins=self.xbins)
        if self.ybins:
            plt.locator_params(axis='y', nbins=self.ybins)

        if self.display_xlim:
            plt.xlim(self.display_xlim)
        if self.display_ylim:
            plt.ylim(self.display_ylim)

        if plotname:
            plt.savefig(plotname, dpi=200, format="png", bbox_inches='tight')

        if show:
            plt.show()

        return fig, ax


