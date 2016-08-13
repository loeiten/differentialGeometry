from diff_geom import Arrow3D, transformation
# sympy imports
from sympy import lambdify
# Numerical imports
import numpy as np
from random import uniform
# Plot imports
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def rnd_num(from_nr = 0.1, to_nr = 2.0):
    """Random number generator"""
    return uniform(from_nr, to_nr)

# Create a figure
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")


# CARTESIAN SYSTEM
#==============================================================================
# Basis vectors
# -------------
# First we would like to plot the basis vectors in the cartesian system
# The basis vectors can be written like
# x_vec = [x1, x2, x3]
# y_vec = [y1, y2, y3]
# z_vec = [z1, z2, z3]
x_vec = np.array([1,0,0])
y_vec = np.array([0,1,0])
z_vec = np.array([0,0,1])
# We will place our basis in the origin, so
# x_start = y_start = z_start = [0, 0, 0]
origin = np.array([0,0,0])
# We will now write the basis vectors as
# x_start_in_origin = [(x_origin, x1), (y_origin, x2), (z_origin, x3)]
# y_start_in_origin = [(x_origin, y1), (y_origin, y2), (z_origin, y3)]
# z_start_in_origin = [(x_origin, z1), (y_origin, z2), (z_origin, z3)]
x_start_in_origin = [(origin[i], x_vec[i]) for i in range(len(x_vec))]
y_start_in_origin = [(origin[i], y_vec[i]) for i in range(len(y_vec))]
z_start_in_origin = [(origin[i], z_vec[i]) for i in range(len(z_vec))]
# Draw the basis vectors
x_ar= Arrow3D(x_start_in_origin[0], x_start_in_origin[1], x_start_in_origin[2],\
              ax, color="b")
y_ar= Arrow3D(y_start_in_origin[0], y_start_in_origin[1], y_start_in_origin[2],\
              ax, color="b")
z_ar= Arrow3D(z_start_in_origin[0], z_start_in_origin[1], z_start_in_origin[2],\
              ax, color="b")
ax.add_artist(x_ar)
ax.add_artist(y_ar)
ax.add_artist(z_ar)
# -------------


# A vector with a starting point different from the origin
# -------------
# Position vector which starts at the origin
R = np.array([rnd_num() for ind in range(3)])
# The vector A
A = np.array([rnd_num() for ind in range(3)])
# We want A to start from the endpoint of R (not in the origin). If we want
# to describe the same vector as starting from the origin, it would
# simply be A_from_the_origin = R + A
# If we say that R = [r1,r2,r3] and A = [a1,a2,a3], we can write the
# vector A_start_in_R as A_start_in_R = [(r1,a1),(r2,a2),(r3,a3)]
A_start_in_R = [(R[i], A[i]) for i in range(len(R))]

print('The A vector in cartesian coordinate system:')
print(A_start_in_R)

# Draw the vector arrow
A_ar = Arrow3D(A_start_in_R[0], A_start_in_R[1], A_start_in_R[2],\
               ax, color="g")
ax.add_artist(A_ar)
# -------------
#==============================================================================



# THE SPHERICAL MANIFOLD
#==============================================================================
# The spherical manifold can be obtianed from the spherical coordinates
# The cartesian coordinates written in a spherical coordinate system is
# stored in the class diff_geom.py. We will obtain the coordinates by
# creating an instance of the transformation class
cartesian_written_in_spherical =\
    transformation(from_coord = 'spherical', to_coord = 'cartesian')
# The coordinates are written symbolically with sympy.
cartesian_coordinates = cartesian_written_in_spherical.to_coords
spherical_coordinates = cartesian_written_in_spherical.from_coords
# As we would like to plot the sphere, it's a good idea to transform the
# symbols to something numerically as sympy's own plotting routine is
# somewhat limited. To do this, we first need to obtain the symbols
r     = spherical_coordinates[0]
theta = spherical_coordinates[1]
phi   = spherical_coordinates[2]
# The radius of the manifold is given by the length of R
radius = np.linalg.norm(R)
# First, we substitute the symbol r with the radius. Then, we make a
# numerical function of the rest of the cartesian coordinates written in
# terms of spherical coordinates
x = lambdify((theta,phi), cartesian_coordinates[0].subs(r, radius), 'numpy')
y = lambdify((theta,phi), cartesian_coordinates[1].subs(r, radius), 'numpy')
z = lambdify((theta,phi), cartesian_coordinates[2].subs(r, radius), 'numpy')

# To plot the mainfold, we need a mesh
# We let u correspond to theta and v correspond to phi
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
U, V = np.meshgrid(u,v)

# We put the meshgrid into the cartesian coordinate system written with
# spherical coordinates
x = x(U,V)
y = y(U,V)
z = z(U,V)

ax.plot_surface(x, y, z, rstride=4, cstride=4,\
                color='y', alpha=0.2, linewidth=0.1)
#==============================================================================



# SPHERICAL COORDINATE SYSTEM
#==============================================================================
# Obtaining the local spherical coordinate system
#------------------------------------------------
# We would now like to plot the local spherical coordinate system atop
# of the vector R. This can be done by transforming the basis vectors in
# the spherical coordinate system to the local cartesian coordinate
# system. The basis vectors in the spherical basis is written as
r_vec     = np.array([1,0,0])
theta_vec = np.array([0,1,0])
phi_vec   = np.array([0,0,1])
# These basis vectors needs to be right dotted with the (contravariant)
# transformation matrix in order to have them described in the local
# cartesian coordinate system.
# Let the instance calculate the transformation matrix
cartesian_written_in_spherical.get_transformation_matrix(do_simplify=True)
# As this is to be evaluated in the endpoint of R, we should rewrite
# {r, theta, phi} to {x, y, z} in order to insert the values of R
# directly
# We firstly save the transformation matrix
transformation_matrix_in_spherical_coords =\
    cartesian_written_in_spherical.trans_matrix
# Subsequently, as sympy treats the transformation matrix symbolically,
# we need to have the x, y and z as symbols. This could of course be
# done by importing sympy and manually declare the symbols, but due to
# laziness, we just obtain these from the transformation class
spherical_written_in_cartesian =\
    transformation(from_coord = 'cartesian', to_coord = 'spherical')
# So, the (old) spherical coordinates in
# transformation_matrix_in_spherical_coords needs to be replaced by the
# (new) cartesian coordinates. We can do that in the following way
old_coordinates = cartesian_written_in_spherical.from_coords
new_coordinates = spherical_written_in_cartesian.to_coords
old_to_new = zip(old_coordinates, new_coordinates)
transformation_matrix_in_cartesian_coords =\
    transformation_matrix_in_spherical_coords.subs(old_to_new)
# The transformation matrix needs to be evaluated in a point. Hence we
# need to put in some numbers. To do this, we first need to write the
# symbolical transformation matrix over to a numerical form. We can do
# this in the same way as we did for the cartesian case
# First we get the coordinates in order to a lambdify
cartesian_coordinates = spherical_written_in_cartesian.from_coords
x = cartesian_coordinates[0]
y = cartesian_coordinates[1]
z = cartesian_coordinates[2]
transformation_matrix = lambdify((x,y,z),\
                        transformation_matrix_in_cartesian_coords, 'numpy')
# Evaluate the matrix at the start of the vector
transformation_matrix = transformation_matrix(R[0], R[1], R[2])
# We are now ready to transform the spherical basis coordinates to
# cartesian coordinates. The spherical coordinate basis vectors reads
# The (contravariant) transformation is obtained by dotting the transformation
# matrix with the basis vectors.
r_vec_in_cartesian     = np.array(np.dot(transformation_matrix,r_vec))[0]
theta_vec_in_cartesian = np.array(np.dot(transformation_matrix,theta_vec))[0]
phi_vec_in_cartesian   = np.array(np.dot(transformation_matrix,phi_vec))[0]
# Note that we recasted the matrix back to an array by taking np.array and
# chosing the zeroth element.
# Note that this transformation gives the vectors in the local cartesian
# basis, positioned at the end of the R vector. However, we need to
# correlate this local origin to the origin of our first cartesian
# system. Therfore, in order to draw the arrows of the basis vectors of
# the local spherical coordinate system, we need to translate the coordinate
# system from [0,0,0] out to the endpoint of R
r_vec_in_cartesian_translated     = r_vec_in_cartesian + R
theta_vec_in_cartesian_translated = theta_vec_in_cartesian + R
phi_vec_in_cartesian_translated   = phi_vec_in_cartesian + R
# We make the vectors on the form
# x = [(x1_start, x1_end), (x2_start, ...), ...]
r_basis_start_in_R     = [(R[i], r_vec_in_cartesian_translated[i])\
                          for i in range(len(R))]
theta_basis_start_in_R = [(R[i], theta_vec_in_cartesian_translated[i])\
                          for i in range(len(R))]
phi_basis_start_in_R   = [(R[i], phi_vec_in_cartesian_translated[i])\
                          for i in range(len(R))]
# Draw the basis vectors
r_ar     = Arrow3D(\
    r_basis_start_in_R[0],\
    r_basis_start_in_R[1],\
    r_basis_start_in_R[2],\
    ax, color="r")
theta_ar = Arrow3D(\
    theta_basis_start_in_R[0],\
    theta_basis_start_in_R[1],\
    theta_basis_start_in_R[2],\
    ax, color="r")
phi_ar   = Arrow3D(\
    phi_basis_start_in_R[0],\
    phi_basis_start_in_R[1],\
    phi_basis_start_in_R[2],\
    ax, color="r")
ax.add_artist(r_ar)
ax.add_artist(theta_ar)
ax.add_artist(phi_ar)
#------------------------------------------------


# Transforming finding the vector A in spherical coordinates
#------------------------------------------------
# At last we would like to find the vector A in spherical coordinates.
# The vector A starts in the local cartesian coordinate system, so the
# vector in spherical coordinates simply becomes
A_in_spherical = np.array(np.dot(transformation_matrix,A))[0]
# Since this starts in the origin of the spherical coordinate system, we have
A_spherical_start_in_R = [(0.0, A_in_spherical[i]) for i in range(len(R))]

print('\n\nThe A vector in spherical coordinate system:')
print(A_spherical_start_in_R)
#------------------------------------------------
#==============================================================================

plt.show()
