# # отладочный файл

import sopromat as sp

frame = sp.frame([[0, 0], [2, 0]])
frame.add_a_holder(0, 0, type='hinge')
frame.add_a_holder(2, 0, 'joint')

frame.add_a_force(1, 0, 10, vec = [0,1], name='force', type='lin', len=1)




frame.basic_draw(size='1500x500')


