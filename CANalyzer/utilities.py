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


def eig(A, ifHermitian=False):
    import scipy.linalg as spl
    if ifHermitian:
        evalue, evector = spl.eigh(A)
    else:
        evalue, evector = spl.eig(A)
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


