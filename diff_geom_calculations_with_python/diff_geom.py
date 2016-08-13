from __future__ import division
# Analytic imports
from sympy import symbols, MatrixSymbol, Matrix, zeros,\
                  sin, cos, sqrt, acos, atan,\
                  diff, Derivative,\
                  simplify
# Plotting imports
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
# Numerical imports
import numpy as np

class Arrow3D(FancyArrowPatch):
    """Draws an arrow"""
    def __init__(self, xs, ys, zs, ax, color = 'k',\
                 mutation_scale=20, lw=1, arrowstyle="-|>", *args, **kwargs):
        """Takes xs, ys and zs as input, which all are on the form
        (start, end)."""
        # PosA and PosB are initialized to (0.0), and changed in draw
        FancyArrowPatch.__init__(self, (0,0), (0,0),\
                color = color, mutation_scale=20, lw=1, arrowstyle="-|>",\
                *args, **kwargs)
        self._verts3d = xs, ys, zs
        # Draw the starting point
        ax.scatter(xs[0],ys[0],zs[0],color=color,s=35)
        ax.scatter(xs[1],ys[1],zs[1],color=color,s=0)

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        # PosA and PosB changed
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class transformation:
    """Transformation from one coordinate system to the other"""
    def __init__(self, from_coord = None, to_coord = None):
        if from_coord == None or to_coord == None:
            raise TypeError ("'from_coord' and 'to_coord' must be set")

        # Decide coordinate transformations
        if from_coord == 'cartesian':
            x, y, z = symbols('x y z', real=True)
            if to_coord == 'spherical':
                r1 = symbols('r1', positive = True)
                theta1, phi1 = symbols('theta1 phi1', real=True)
                # Spherical written as cartesian
                # ===============================
                r1 = sqrt(x**2 + y**2 + z**2)
                theta1 = acos(z/sqrt(x**2 + y**2 + z**2))
                phi1 = atan(y/x)
                # ===============================
                self.from_coords = [x, y, z]
                self.to_coords   = [r1, theta1, phi1]
            else:
                self.errors(from_coord = from_coord, to_coord=to_coord)
        elif from_coord == 'spherical':
            r = symbols('r', positive = True)
            theta, phi  = symbols('theta phi', real=True)
            if to_coord == 'cartesian':
                x1, y1, z1 = symbols('x1 y1 z1', real=True)
                # Cartesian  written as spherical
                # ===============================
                x1=r*sin(theta)*cos(phi)
                y1=r*sin(theta)*sin(phi)
                z1=r*cos(theta)
                # ===============================
                self.from_coords = [r, theta, phi]
                self.to_coords   = [x1, y1, z1]
            else:
                self.errors(from_coord = from_coord, to_coord=to_coord)
        else:
            self.errors(from_coord = from_coord)

        # Set additional parameters
        self.trans_matrix = None
        self.jacobian = None

#{{{errors
    def errors(self, from_coord=None, to_coord=None):
        """Giving the errors"""
        from_coord_set = False
        to_coord_set = False
        if from_coord != None:
            if type(from_coord) != str:
                raise TypeError ('from_coord must be given as a string')
            from_coord_set=True
        if to_coord != None:
            if type(to_coord) != str:
                raise TypeError ('to_coord must be given as a string')
            to_coord_set = True

        if to_coord_set == True and from_coord_set == True:
            raise ValueError ("from_coord '" + from_coord + "' to " +\
                              "to_coord '" + to_coord + "' not implimented")
        elif to_coord_set == True and from_coord_set == False:
            raise ValueError ("to_coord '" + to_coord + "' not implimented")
        elif to_coord_set == False and from_coord_set == True:
            raise ValueError ("from_coord '" + from_coord + "' not implimented")
#}}}

    def get_transformation_matrix(self, do_simplify=False):
        """Returns the transformation matrix for contravariant
        transformations"""
        self.trans_matrix = Matrix.zeros(len(self.from_coords), len(self.to_coords))
        for from_ind, from_coord in enumerate(self.from_coords):
            for to_ind, to_coord in enumerate(self.to_coords):
                self.trans_matrix[to_ind, from_ind] =\
                    diff(to_coord, from_coord)
                if do_simplify:
                    simplify(self.trans_matrix[to_ind, from_ind])

    def get_jacobian(self, do_simplify=False):
        """Returns the jacobian of the transformation matrix for
        covariant transformations"""
        if self.trans_matrix == None:
            # Get the transformation matrix
            self.get_transformation_matrix(simplify)
        if do_simplify == False:
            self.jacobian = self.trans_matrix.det()
        else:
            self.jacobian = simplify(self.trans_matrix.det())
