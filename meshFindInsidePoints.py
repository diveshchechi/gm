import numpy as np
import trimesh

mesh = trimesh.load_mesh('mushroom.off')
print(mesh.is_watertight)

mesh.vertices += [10, 10, 10]
mesh.export('mushTranslate.off')

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
x = 100
y = 100
z = 100
boxCordinates = np.zeros((x, y, z, 3))
ix = 9.200
iy = 9.1140
iz = 9.084106
gridSize = 100

incrementX = 1.8 / gridSize
incrementY = 1.50 / gridSize
incrementZ = 2.0 / gridSize
# print(increment)
for i in range(gridSize):
    for j in range(gridSize):
        for k in range(gridSize):
            boxCordinates[i][j][k][0] += (ix + i * incrementX)
            boxCordinates[i][j][k][1] += (iy + j * incrementY)
            boxCordinates[i][j][k][2] += (iz + k * incrementZ)
# print(boxCordinates)
points = np.reshape(boxCordinates, (x * y * z, 3))
# print(pts)

inside = trimesh.proximity.signed_distance(mesh, points)

f = open("mushroom100.txt", "a")
total = x * y * z
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
        for k in range(gridSize):
            f.write(str(boxCordinates[i][j][k][0]))
            f.write("\n")
            f.write(str(boxCordinates[i][j][k][1]))
            f.write("\n")
            f.write(str(boxCordinates[i][j][k][2]))
            f.write("\n")
f.close()
'''
