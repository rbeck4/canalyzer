#Pray, pray to St. Isidore; may he have mercy on your poor soul for trying to
#run something that I've written.
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy
import string
import sys
sys.path.append("../")
from CANalyzer.ci_spectra import CI_spectra

font = { 
    'size' : '20', 
    'weight' : 'bold'
    }
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 4.0
matplotlib.rcParams['hatch.linewidth'] = 2.0
matplotlib.rcParams['lines.linewidth'] = 2.0
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
har2ev = 27.2114

enMin = 0
enMax = 5
res = .001
nPoints = int((enMax-enMin)/res)

plots = []
subPlots = []
energy = []
subEnergy = []

fileList = [\
            "outputFiles/Na_ci1.out", \
           ]

subList  = [[\
            "outputFiles/Na_ci2.out", \
            ],[ \
            ]]


fchkList = [\
            "outputFiles/Na_ci1.bin", \
           ]

subFchk = [[\
            "outputFiles/Na_ci2.bin", \
            ],[ \
            ]]


names    = [\
            "NaTest", \
           ]

shifts   = [\
            0, \
            0, \
           ]

for i, fl in enumerate(fileList):
  print("WORKING, ", fl, " [", i+1, "/", len(fileList), "]")
  plots.append(CI_spectra(fl, fchkList[i], None))
  
  if len(subList[i]) > 0:
    subs = []
    for j, sbfl in enumerate(subList[i]):
      print("SUB WORKING, ", sbfl, " [", j+1, "/", len(subList[i]), "]")
      subs.append(CI_spectra(sbfl, subFchk[i][j], None))
    subPlots.append(subs)
  else:
    subPlots.append(None)

for i in range(len(fileList)):
  plots[i].decompose_byorbital()
  energy.append((plots[i].energy[0,:])*har2ev)
  subs = []
  for j in range(len(subList[i])):
    subPlots[i][j].decompose_byorbital()
    subs.append((subPlots[i][j].energy[0,:])*har2ev)
  subEnergy.append(subs)
  print("Plots: ", names[i], " max: ", max(energy[-1]))

space_names  = list(plots[0].spaces)
for i in range(len(fileList)):
  for space in range(len(space_names)):
    os = plots[i].decomp_byorbital[0, :, space]
    totEn = energy[i]
    for j in range(len(subList[i])):
      tmpEn = np.asarray(subEnergy[i][j])
      idxEn = np.where(tmpEn > totEn[-1])
      totEn = np.concatenate([totEn, tmpEn[idxEn]])
      tmpOs = subPlots[i][j].decomp_byorbital[0, :, space]
      os = np.concatenate([os, tmpOs[idxEn]])
    plots[i].make_spectrum(space_names[space], shifts[i]-5, max(totEn), totEn, os, 0.14, npoints=nPoints)

space_names += ["Full Spectrum"]

fullEn = [[0] for x in range(len(fileList))]
fullOs = [[0] for x in range(len(fileList))]
for i in range(len(fileList)):
  fullOs[i] = plots[i].oscstr[0]
  fullEn[i] = energy[i]
  for j in range(len(subList[i])):
      tmpOs = subPlots[i][j].oscstr[0]
      tmpEn = np.asarray(subEnergy[i][j])
      idxEn = np.where(tmpEn>fullEn[i][-1])
      fullEn[i] = np.concatenate([fullEn[i],tmpEn[idxEn]])
      fullOs[i] = np.concatenate([fullOs[i],tmpOs[idxEn]])
  plots[i].make_spectrum("Full Spectrum", shifts[i]-5, max(fullEn[i]), fullEn[i], fullOs[i], 0.14, npoints=nPoints)
  #plots[i].spectra_yscale = 1

#PLOTTING:
fig, ax = plt.subplots(len(fileList), sharex='col', figsize=[15,15], gridspec_kw={'hspace': 0},squeeze=False)

for e,p in enumerate(plots):
  for i, space in enumerate(space_names):
    y = p.spectra[space]
    x = np.linspace(shifts[e]-5, max(fullEn[e]), len(y))
    if "Full" in space:
      ax[e,0].plot(x-shifts[e], y, label=str(space), color='k', linewidth=2)
    else:
      ax[e,0].plot(x-shifts[e], y, label=str(space), color=colors[i], linewidth=4)
  
  ax[e,0].set_yticks([])
  ax[e,0].set_xlim([enMin,enMax])
  ax[e,0].axhline(y=0, xmin=enMin, xmax=enMax, color='k')
  #handles, labels = ax[e,0].get_legend_handles_labels()
  ax[e,0].legend(loc='upper center', handlelength=0.5, frameon=False, 
      prop=dict(size=15), ncol=6)

  ax[e,0].text(0.01, 0.9, names[e], ha='left', va='top', transform=ax[e,0].transAxes)
  ax[e,0].text(0.99, 0.01, "+%.2f eV" %(shifts[e]), ha='right', va='bottom', transform=ax[e,0].transAxes)
  
  #Final Energy
  ax[e,0].axvline(x=fullEn[e][-1]-shifts[e], ymin=-1, ymax=1, color='yellow', linestyle='-.')
  

fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.ylabel("Intensity (a.u.)", fontweight='bold')
plt.xlabel("Energy (eV)", fontweight='bold')
fig.tight_layout()
plt.savefig("CQ_byOrbital.png")
  
