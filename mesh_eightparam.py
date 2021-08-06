import numpy as np
import trimesh

mesh = trimesh.load_mesh('eightparam.off')
print(mesh.is_watertight)

mesh.vertices += [10, 10, 10]

points = mesh.bounding_box.sample_volume(count=10)
print("bounding box vertices")
print(mesh.bounding_box.vertices)
print("bounding box bounds")
print(mesh.bounding_box.bounds)
print("bounding box extents")
print(mesh.bounding_box.extents)
print("signed distances -ve means outside")
print(trimesh.proximity.signed_distance(mesh, points))
print(points)
x = 64
y = 64
z = 64
gridSize = 64
alpha = 32
# earlier it was x y 32 and z 64
boxCordinates = np.zeros((x, y, z, 3))
ix = 9.757102
iy = 9.8957
iz = 9.4900
incrementX = 0.5 / gridSize
incrementY = 0.21 / alpha
incrementZ = 1 / gridSize
# print(increment)
for i in range(gridSize):
    for j in range(gridSize):
        for k in range(gridSize):
            boxCordinates[i][j][k][0] += (ix + i * incrementX)
            boxCordinates[i][j][k][1] += (iy + j * incrementY)
            boxCordinates[i][j][k][2] += (iz + k * incrementZ)
# print(boxCordinates)
points = np.reshape(boxCordinates, (x*y*z, 3))
# print(pts)

inside = trimesh.proximity.signed_distance(mesh, points)

f = open("eighty80alpha.txt", "a")
total = x*y*z
for i in range(total):
    if inside[i] < 0:
        f.write(str(0))
    else:
        f.write(str(1))
    # f.write(str(round(inside[i])))
    f.write("\n")
'''
for i in range(gridSize):
    for j in range(gridSize):
        for k in range(gridSize * 2):
            f.write(str(boxCordinates[i][j][k][0] ))
            f.write("\n")
            f.write(str(boxCordinates[i][j][k][1] ))
            f.write("\n")
            f.write(str(boxCordinates[i][j][k][2] ))
            f.write("\n")
f.close()
'''
