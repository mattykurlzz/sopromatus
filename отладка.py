# отладочный файл

import sopromat as sp

frame = sp.frame(((0,0), (0, 100), (200,100), (200,0)))
#frame.draw_raw()
frame.add_a_force(0, 100, 10, [1,0], 'force')
frame.add_a_force(100, 100, 20, 'clockwise', name='moment')
frame.add_a_force(200, 100, 10, [1,1], name='force')

frame.basic_draw()

