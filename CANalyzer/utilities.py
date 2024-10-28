def onlytype(l, type):
    """
    Filters list, keeping only characters that can be converted to said type, and saves it in that type.
    :param l: list
    :param type: type
    :return: list
    """
    new_l = []
    for x in l:
        try:
            new_l.append(type(x))
        except:
            pass
    return new_l


def write_fchk(startstr, nri, matrix, matsize, oldfchk, fchk):
    import numpy as np
    import os
    import subprocess

    startline = int(str(subprocess.check_output(f"grep -n '{startstr}' {oldfchk}", shell=True)).split(" ")[
                        0].split("'")[1].split(":")[0])

    os.system(f"sed -n '1,{startline}p' {oldfchk} >> {fchk}")

    totalsize = matsize * nri
    if nri == 1:
        flat = matrix.flatten(order='F')
    if nri == 2:
        flatreal = np.real(matrix).flatten(order='F')
        flatimag = np.imag(matrix).flatten(order='F')
        flat = np.zeros(totalsize)
        flat[::2] = flatreal
        flat[1::2] = flatimag
    rows = int(np.ceil(totalsize / 5))
    remainder = totalsize % 5
    for i in range(rows):
        if i == rows - 1 and remainder != 0:
            line = flat[-remainder:]
        else:
            line = flat[5*i:5*i+5]
        strline = ["  {:.8E}".format(x) if x >= 0 else " {:.8E}".format(x) for x in line]
        finalstring = "".join(strline)
        with open(fchk, "a", newline='') as file:
            file.write(f"{finalstring}\n")

    startline2 = startline + rows + 1
    os.system(f"tail -n +{startline2} {oldfchk} >> {fchk}")


def eig(A):
    import numpy as np
    import numpy.linalg as npl
    evalue, evector = npl.eig(A)
    idx = evalue.argsort()[::-1]   
    evalue = evalue[idx]
    evector = evector[:,idx]
    return evalue, evector


def normalize_matrix(A):
    import numpy as np

    nrows, ncols = A.shape
    for j in range(ncols):
        vec = A[:, j] / np.linalg.norm(A[:, j])
        for i in range(nrows):
            A[i, j] = vec[i]
    return A


def gaussian(energy, broadening, osc_str, xspace):
    import numpy as np
    ones = np.ones(xspace.shape)
    return (osc_str / (np.sqrt(2*np.pi*broadening**2))) * np.exp(-(xspace-ones*energy)**2/(2*broadening**2))


def lorentzian(energy, broadening, osc_str, xspace):
    import numpy as np
    ones = np.ones(xspace.shape)
    d = np.pi * (xspace - ones * energy) ** 2 + broadening ** 2
    lorentzian = broadening / d
    return osc_str * lorentzian


"""
def plot_spectra(x, y, colors, labels, alphas, roots=None, os=None, redshift=0, yscale=1, yscale_sticks=1):
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    fig, ax = plt.subplots()
    fig.set_figheight(30)
    fig.set_figwidth(40)
    fig.canvas.draw()

    plt.rc('axes', linewidth=10)
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"

    for i in range(len(x)):
        color = colors[i]
        l = labels[i]
        a = alphas[i]
        ax.plot(x[i] - redshift, y[i] * yscale, c=color, label=l, alpha=a, linewidth=10)
        if roots and os and roots[i]:
            for n in range(len(roots[i])):
                root = roots[i][n] - redshift
                height = os[i][n]*yscale_sticks
                ax.plot((root, root), (0, height), c=color, linewidth=10)

    ax.set_xlabel('Energy (eV)', fontsize=90, fontweight='bold')
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(75)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(0)
        tick.label1.set_fontweight('bold')
    ax.tick_params(axis='both', which='major', pad=15)
    ax.legend(fontsize=30, loc='upper left', prop={'weight': 'bold', 'size': 80}, frameon=False, ncol=1,
              columnspacing=1)
    ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.2f}'))
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.2f}'))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(xstart, xend)
    plt.show()
"""




c = 137.035999074


OAM = {
    "S": 0,
    "P": 1,
    "D": 2,
    "F": 3,
    "G": 4,
    "H": 5,
    "I": 6
}

OAM2 = {
    0: "S",
    1: "P",
    2: "D",
    3: "F",
    4: "G",
    5: "H",
    6: "I"
}

ptable = {
1:"H",
2:"He",
3:"Li",
4:"Be",
5:"B",
6:"C",
7:"N",
8:"O",
9:"F",
10:"Ne",
11:"Na",
12:"Mg",
13:"Al",
14:"Si",
15:"P",
16:"S",
17:"Cl",
18:"Ar",
19:"K",
20:"Ca",
21:"Sc",
22:"Ti",
23:"V",
24:"Cr",
25:"Mn",
26:"Fe",
27:"Co",
28:"Ni",
29:"Cu",
30:"Zn",
31:"Ga",
32:"Ge",
33:"As",
34:"Se",
35:"Br",
36:"Kr",
37:"Rb",
38:"Sr",
39:"Y",
40:"Zr",
41:"Nb",
42:"Mo",
43:"Tc",
44:"Ru",
45:"Rh",
46:"Pd",
47:"Ag",
48:"Cd",
49:"In",
50:"Sn",
51:"Sb",
52:"Te",
53:"I",
54:"Xe",
55:"Cs",
56:"Ba",
57:"La",
58:"Ce",
59:"Pr",
60:"Nd",
61:"Pm",
62:"Sm",
63:"Eu",
64:"Gd",
65:"Tb",
66:"Dy",
67:"Ho",
68:"Er",
69:"Tm",
70:"Yb",
71:"Lu",
72:"Hf",
73:"Ta",
74:"W",
75:"Re",
76:"Os",
77:"Ir",
78:"Pt",
79:"Au",
80:"Hg",
81:"Tl",
82:"Pb",
83:"Bi",
84:"Po",
85:"At",
86:"Rn",
87:"Fr",
88:"Ra",
89:"Ac",
90:"Th",
91:"Pa",
92:"U",
93:"Np",
94:"Pu",
95:"Am",
96:"Cm",
97:"Bk",
98:"Cf",
99:"Es",
100:"Fm",
101:"Md",
102:"No",
103:"Lr",
104:"Rf",
105:"Db",
106:"Sg",
107:"Bh",
108:"Hs",
109:"Mt",
110:"Ds",
111:"Rg",
112:"Cn",
113:"Nh",
114:"Fl",
115:"Mc",
116:"Lv",
117:"Ts",
118:"Og",
}


