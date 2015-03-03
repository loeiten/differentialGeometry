from diff_geom import transformation
from sympy import init_printing, pprint
from sympy import simplify
init_printing()

# Spherical to cartesian
print('\n'*4)
print('Spherical to cartesian:')
s_to_c = transformation(from_coord = 'spherical', to_coord = 'cartesian')
s_to_c.get_jacobian(do_simplify=True) 
print('Transformation matrix:')
pprint(s_to_c.trans_matrix)
print('\nJacobian:')
pprint(s_to_c.jacobian)

# Cartesian to shperical
print('\n'*4)
print('Cartesian to spherical:')
c_to_s = transformation(from_coord = 'cartesian', to_coord = 'spherical')
c_to_s.get_jacobian(do_simplify=True) 
print('Transformation matrix:')
pprint(c_to_s.trans_matrix)
print('\nJacobian:')
pprint(c_to_s.jacobian)

# Rewritten in spherical
print('\n'*4)
print('Rewritten in spherical:')
# Get the coordinates
old = c_to_s.from_coords
new = s_to_c.to_coords
old_to_new = zip(old,new)
print('We replace in the following way:')
pprint(old_to_new)
# Rewrite the matrix
rew_m = simplify(c_to_s.trans_matrix.subs(old_to_new))
# Rewrite the jacobian
rew_j = simplify(c_to_s.jacobian.subs(old_to_new))
print('Transformation matrix:')
pprint(rew_m)
print('\nJacobian:')
pprint(rew_j)
