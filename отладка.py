# # отладочный файл

import sopromat as sp
import numpy as np

# frame = sp.twodframe([[0, 0], [0, .5], [.5, .5], [.5, 0], [1.5, 0], [1.5, -.5], [.5, -.5]], raster=50)
# frame.add_a_holder(0, 0, 'joint', [0,1])
# frame.add_a_holder(.5, -.5, 'hinge', [0, 1])
# frame.add_a_force(1.5, -.5, 7, vec=[1, 1], name='moment')
# frame.add_a_force(.5, .5, -2, [1, 0], 'force', 'lin', stopx=.5, stopy=0)

# frame=sp.twodframe([[0, 0], [0, .5], [.5, 1], [1, 1]], raster=50 )
# frame.add_a_holder(0, 0, 'hinge', [0, 1])
# frame.add_a_holder(1, 1, 'joint', [1, 0])
# frame.add_a_force(0, 0, 2, [1,0], 'force', 'lin', stopx=0, stopy=0.5)
# frame.add_a_force(0, .5, -7, vec=[1,1], name='moment')



# frame.set_size('1200x1000')
# frame.start_drawing(mesh=[3, 1])
# frame.basic_draw(pos = [0,0], text='Начальная система')
# frame.draw_forces(pos=[1, 0], type='продольные', text='Продольные усилия', scale=1)
# frame.draw_moments(pos = [2, 0], text='Эпюра моментов', scale=2)
# frame.print_it()


# Вариант 112. Момент 1 10 кн, q1 4 кн, a=0.8, схемы 4, 15, 25, 37
# P1 = 22. M2 = 13. a = 1.2. горизонтальные перемещения схема 15

frame=sp.twodframe([[1, 1], [0, 1], [0, 0], [3, 0], [3, 1]], raster=100)
frame.add_a_holder(1, 1, 'hinge', [0, 1])
frame.add_a_holder(0, 0, 'joint', [0, 1])
frame.add_a_holder(3, 0, 'hinge', [0, 1])
frame.add_a_force(3, 1, -21, [1, 0], 'force', 'point')
frame.add_a_force(3, 0, -15, [1, 0], 'force', 'lin', stopx=3, stopy=1)



frame.set_size('1200x1000')
frame.start_drawing(mesh=[2, 2])
frame.basic_draw(pos = [0,0], text='Начальная система')
frame.system_draw(pos=[1, 0], text='Эквивалентная система')
frame.draw_forces(pos = [0,1], type='продольные', text='Эпюра продольных сил', scale=0.2)
frame.draw_moments(pos = [1, 1], text='Эпюра моментов', scale=0.5)
frame.print_it()
################################################################################################################
# frame=sp.threedframe(coords=[[0,0,0], [1,0,0], [1,0,1], [1,1,1]])
# frame.add_a_holder(0, 0, 0, type='seal', vector=[0,1,0])

# frame.set_size('1200x1000')
# frame.start_drawing(mesh=[1, 1])
# frame.base_draw(pos=[0,0])
# frame.print_it()
