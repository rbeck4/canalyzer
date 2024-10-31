from CANalyzer.ci_spectra import CI_spectra

plot = CI_spectra("cq-ci.out", "cq-ci.bin", None)

plot.spaces = {r"p$_{1/2}$" : (1, 2), r"p$_{3/2}$" : (3, 6), r"d$_{3/2}$" : (7, 12), r"d$_{5/2}$" : (13, 16)}
space_names = list(plot.spaces.keys())

plot.decompose_byspaces()

energy = plot.energy[0, :]*27.2114
for space in range(len(space_names)):
    os = plot.decomp_byspaces[0, :, space]
    plot.make_spectrum(space_names[space], 450, 495, energy, os, 0.2)

plot.make_spectrum("Full Spectrum", 450, 495, energy, plot.oscstr[0], 0.2)

#space_names += ["Full Spectrum"]

plot.spectra_yscale = 1

plot.plot(space_names, 450, 490, plotname="test.png")

