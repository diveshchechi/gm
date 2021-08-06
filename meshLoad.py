import numpy as np
import trimesh
mesh = trimesh.load_mesh('eightparam.off')
print(mesh.is_watertight)

mesh.vertices += [10 ,10, 10]

points = mesh.bounding_box.sample_volume(count=10)
print("bounding box vertices")
print(mesh.bounding_box.vertices)
print("bounding box bounds")
print(mesh.bounding_box.bounds)
print("bounding box extents")
print(mesh.bounding_box.extents)
print("signed distances -ve means outside")
print(trimesh.proximity.signed_distance(mesh, points))
print(np.shape(points))


# create a PointCloud object out of each (n,3) list of points
cloud_original = trimesh.points.PointCloud(points)
cloud_close    = trimesh.points.PointCloud(points)

# create a unique color for each point
cloud_colors = np.array([trimesh.visual.random_color() for i in points])

# set the colors on the random point and its nearest point to be the same
cloud_original.vertices_color = cloud_colors
cloud_close.vertices_color    = cloud_colors

# create a scene containing the mesh and two sets of points

scene = trimesh.Scene([mesh,
                       cloud_original,
                       cloud_close])

# show the scene wusing
scene.show()

