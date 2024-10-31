from CANalyzer.tddft_spectra import TDDFT_spectra

plot = TDDFT_spectra("cepn-star.log", "cepn-star.fchk", None)

plot.purge()

plot.occupied_spaces = {"Occ 1" : (344, 347), "Occ 2" : (342, 343), "Occ 3" : (341, 341)}
plot.virtual_spaces = {r"t$_1$ (4f)" : (348, 350), r"a$_1$ (4f)" : (351, 351), r"t$_2$ (4f)" : (352, 354)}

plot.decompose_byexcat()

cat_names = []
energy = plot.energy
for icat in range(plot.ncat):
    os = plot.decomp_byexcat[:, icat]
    fromcat = plot.categories[icat][0]
    tocat = plot.categories[icat][1]
    if fromcat != "Other Occ" and tocat != "Other Virt":
        label = f"{fromcat} -> {tocat}"
        cat_names.append(label)
        plot.make_spectrum(label, 1, 7, energy, os*4, 0.5)

plot.make_spectrum("Full Spectrum", 1, 7, energy, plot.oscstr, 0.5)
cat_names.append("Full Spectrum")

plot.spectra_yscale = 2
plot.legend_ncol = 3
plot.legend_font = 50
plot.display_xlim = (2, 5.2)
plot.display_ylim = (0, 4)
plot.xbins = 5
plot.x_decimal = 2
plot.alpha=1
plot.colors = ["red", "maroon", "salmon", "green", "lime", "mediumspringgreen", "blue", "teal", "royalblue", "darkorange", "gold", "tan", "magenta", "purple", "violet"]
plot.plot(cat_names, 1, 7, plotname="cepn-star-byexcat.png")

