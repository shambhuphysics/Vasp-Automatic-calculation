# The MIT License (MIT)

import numpy as np
from copy import deepcopy

# First part: Convert XDATCAR to XYZ format

with open('XDATCAR', 'r') as xdatcar, open('XDATCAR.xyz', 'w') as xyz, open('XDATCAR_fract.xyz', 'w') as xyz_fract:
    system = xdatcar.readline()
    scale = float(xdatcar.readline().strip())
    print(scale)

    # get lattice vectors
    a1 = np.array([float(s) * scale for s in xdatcar.readline().strip().split()])
    a2 = np.array([float(s) * scale for s in xdatcar.readline().strip().split()])
    a3 = np.array([float(s) * scale for s in xdatcar.readline().strip().split()])
    print(a1)
    print(a2)
    print(a3)

    # Save scaled lattice vectors
    with open('lattice.vectors', 'w') as lat_rec:
        lat_rec.write(f"{a1[0]} {a1[1]} {a1[2]}\n")
        lat_rec.write(f"{a2[0]} {a2[1]} {a2[2]}\n")
        lat_rec.write(f"{a3[0]} {a3[1]} {a3[2]}")

    # Read xdatcar
    element_names = xdatcar.readline().strip().split()
    element_dict = {}
    element_numbers = xdatcar.readline().strip().split()
    N = 0
    for t, num in enumerate(element_numbers):
        element_dict[element_names[t]] = int(num)
        N += int(num)
    print(element_dict)

    while True:
        line = xdatcar.readline()
        if not line:
            break
        xyz.write(f"{N}\ncomment\n")
        xyz_fract.write(f"{N}\ncomment\n")
        for el in element_names:
            for _ in range(element_dict[el]):
                p = xdatcar.readline().strip().split()
                coords = np.array([float(s) for s in p])
                cartesian_coords = coords[0] * a1 + coords[1] * a2 + coords[2] * a3
                xyz.write(f"{el} {cartesian_coords[0]} {cartesian_coords[1]} {cartesian_coords[2]}\n")
                xyz_fract.write(f"{el} {coords[0]} {coords[1]} {coords[2]}\n")

# Second part: MSD calculation

def MSD(xyz_file, L):
    a = L
    l = [np.sqrt(np.dot(v, v)) for v in a]  # basis vector lengths

    with open(xyz_file, 'r') as file, open("msd.out", 'w') as recorder, open("unwrapped.xyz", 'w') as coord_rec:
        origin_list = []  # Stores the origin as [element,[coords]]
        prev_list = []  # Stores the wrapped previous step
        unwrapped_list = []  # Stores the instantenous unwrapped

        msd = []  # Stores atom-wise MSD  Stores msd as [msd]
        msd_dict = {}  # Stores element-wise MSD
        msd_lattice = []
        msd_dict_lattice = {}

        element_list = []  # element list
        element_dict = {}  # number of elements stored

        N = int(file.readline())

        msd = [np.float64('0.0') for _ in range(N)]
        msd_lattice = [[0.0, 0.0, 0.0] for _ in range(N)]

        file.readline()
        step = 0

        while True:
            step += 1
            # Get and store the origin coordinates in origin_dict at first step
            if step == 1:
                for _ in range(N):
                    t = file.readline().strip().split()
                    element = t[0]
                    if element not in element_list:
                        element_list.append(element)
                    element_dict[element] = element_dict.get(element, 0) + 1
                    coords = np.array([float(s) for s in t[1:]])
                    origin_list.append([element, coords])
                # Copy the first set of coordinates as prev_dict and unwrapped
                unwrapped_list = deepcopy(origin_list)
                prev_list = deepcopy(origin_list)
                recorder.write("step " + " ".join(element_list) + "\n")

            # Read wrapped coordinates into wrapped_dict
            content = file.readline()
            if not content:
                print("\n---End of file---\n")
                break
            N = int(content)
            file.readline()
            wrapped_list = []  # Erase the previous set of coordinates
            for _ in range(N):
                t = file.readline().strip().split()
                element = t[0]
                coords = np.array([float(s) for s in t[1:]])
                wrapped_list.append([element, coords])

            coord_rec.write(f"{N}\ncomment\n")

            # Unwrap coodinates and get MSD
            for atom in range(N):
                msd[atom] = 0.0
                coord_rec.write(wrapped_list[atom][0])

                # decompose wrapped atom coordinates to onto lattice vectors:
                w1, w2, w3 = wrapped_list[atom][1]
                # decompose prev atom coordinates to onto lattice vectors:
                p1, p2, p3 = prev_list[atom][1]

                # get distance between periodic images and use the smallest one
                u1 = w1 - p1 - np.sign(w1 - p1) if abs(w1 - p1) > 0.5 else w1 - p1
                u2 = w2 - p2 - np.sign(w2 - p2) if abs(w2 - p2) > 0.5 else w2 - p2
                u3 = w3 - p3 - np.sign(w3 - p3) if abs(w3 - p3) > 0.5 else w3 - p3

                # add unwrapped displacements to unwrapped coords
                unwrapped_list[atom][1] += np.array([u1, u2, u3])

                uw = sum(unwrapped_list[atom][1][i] * a[i] for i in range(3))
                ol = sum(origin_list[atom][1][i] * a[i] for i in range(3))

                msd[atom] = np.linalg.norm(uw - ol) ** 2
                msd_lattice[atom] = [(uw[i] - ol[i]) ** 2 for i in range(3)]

                coord_rec.write(" " + np.array2string(uw, separator=' ', suppress_small=True)[1:-1] + "\n")

            prev_list = deepcopy(wrapped_list)
            # record msd
            recorder.write(f"{step} ")

            for el in element_list:
                msd_dict[el] = 0.0
                msd_dict_lattice[el] = [0., 0., 0.]

            for atom, (el, _) in enumerate(wrapped_list):
                msd_dict[el] += msd[atom] / element_dict[el]
                for i in range(3):
                    msd_dict_lattice[el][i] += msd_lattice[atom][i] / element_dict[el]

            for el in element_list:
                recorder.write(f"{msd_dict[el]} {' '.join(map(str, msd_dict_lattice[el]))} ")

            recorder.write("\n")
            if step % 10 == 0:
                print(step)

# Run the MSD calculator with XDATCAR_fract.xyz and the lattice vector defined above
lattice = [a1, a2, a3]
MSD("XDATCAR_fract.xyz", lattice)
