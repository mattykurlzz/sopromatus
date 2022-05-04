## тут я буду писать подключаемый файл

# для начала надо придумать как хранить плоскую балку

# балка хранится в виде точек. Кажая точка состоит из списка двух(для пространственной потом допишу) координат
from math import sqrt
import tkinter as tk
import numpy as np

def isgood(lst): 
    if type(lst) is (tuple or list): 
        for element in lst: 
            if type(element) is (tuple or list) and len(element) == 2: continue
            else: return False
        return True
    else: return False

def flatten(coord):
    coord = np.array(coord)
    return coord.flatten() 

def dims(coord): 
    return max(coord[:, 0]) - min(coord[:, 0]), max(coord[:, 1]) - min(coord[:, 1]) 



class frame: 



    coord = []
    dots_num = 0 #объявляю явно, чтобы не забыть, что у меня тут



    def __init__(self, coord):
        if isgood(coord): 
            self.coord = np.array(coord)
            self.dots_num = len(coord)
            return None
        else: 
            print('Bad coordinates')
            return -1



    def draw_raw(self, size = '900x900', mesh=[3,3], place=[1,1]): #рисование "сырой" рамы. рудимент, который оставлю для референса
        # mesh - переменная, хранящая разбиение экрана на области. По-умолчанию бьёт на 9 одинаковых квадратов, нумерация с единицы
        # place - место, куда в сетке вставляется рисунок. нумерация с нуля.
        wd, hg = list(map(int, size.split('x')))
        spacing = wd//200
        vspacing = hg//100
        frame_wd, frame_hg = dims(self.coord) #это размеры самой рамки. сейчас не используются, нужны будут потом для масштабирования картинки. Пока на доверии.

        printpos_x = (wd)//mesh[0]*place[0]+spacing
        printpos_y = (hg)//mesh[1]*(place[1]+1)-vspacing

# инициализация окна рисования тут нужна только если надо нарисовать сырую раму.
        window = tk.Tk()
        window.geometry(size)
        canvas = tk.Canvas(window, width=wd, height=hg, bg = 'white')

        for i in range(self.dots_num-1): 
            canvas.create_line(self.coord[i][0] + printpos_x, printpos_y-self.coord[i][1], self.coord[i+1][0]+printpos_x, printpos_y-self.coord[i+1][1], fill='black', width=900//300) #если этот кусок подогнать под 
            # рисование одной и той же рамы где угодно - можно использовать для шаблона, сверху которого рисуем эпюры

        canvas.pack()
        window.mainloop()



    def framedraw(self, canvas, wd, hg, mesh, place, spacing=10, vspacing=10): 
        printpos_x = (wd)//mesh[0]*place[0]+spacing
        printpos_y = (hg)//mesh[1]*(place[1]+1)-vspacing

        for i in range(self.dots_num-1): 
            canvas.create_line(self.coord[i][0] + printpos_x, printpos_y-self.coord[i][1], self.coord[i+1][0]+printpos_x, printpos_y-self.coord[i+1][1], fill='black', width=5) #если этот кусок подогнать под 
        # рисование одной и той же рамы где угодно - можно использовать для шаблона, сверху которого рисуем эпюры

    def factordraw(self, canvas, wd, hg, mesh, place, factor, spacing, vspacing): 
        frame_wd, frame_hg = dims(self.coord)
        truex = factor.xpos+(wd)//mesh[0]*place[0]+spacing
        truey = (hg)//mesh[1]*(place[1]+1)-vspacing - factor.ypos

        if factor.name == 'force' and factor.type=='point':
            canvas.create_line(truex, truey,
            truex+frame_wd//5*factor.vec[0], truey-frame_wd//5*factor.vec[1], 
            fill='red', width=4, arrow='last')

        if factor.name== 'moment': 
            canvas.create_line(truex, truey, 
            truex, truey-frame_wd//10,
            fill = 'blue', width=4)
            canvas.create_arc(truex-frame_wd//10, truey-frame_wd//10, 
            truex+frame_wd//10, truey+frame_wd//10, 
            start=90, extent=45*(-1 if factor.vec=='clockwise' else 1), 
            width=4, outline='blue', style=tk.ARC)
            canvas.create_line(truex+frame_wd//10/sqrt(2), truey-frame_wd//10/sqrt(2), 
            truex+frame_wd//10/sqrt(2)+5, truey-frame_wd//10/sqrt(2)+5,
            fill = 'blue', width=4, arrow='last')
            #переделать)0)0)0




    def basic_draw(self, size = '900x300', mesh=[3,1]): 
        wd, hg = list(map(int, size.split('x')))
        spacing = wd//20
        vspacing = hg//10
        frame_wd, frame_hg = dims(self.coord) #это размеры самой рамки. сейчас не используются, нужны будут потом для масштабирования картинки. Пока на доверии.

        window = tk.Tk()
        window.geometry(size)
        canvas = tk.Canvas(window, width=wd, height=hg, bg = 'white')

        self.framedraw(canvas, wd, hg, mesh, [0,0], spacing, vspacing)
        self.framedraw(canvas, wd, hg, mesh, [1,0], spacing, vspacing)
        self.framedraw(canvas, wd, hg, mesh, [2,0], spacing, vspacing) #отрисовывает три рамки потому что я так сказал

        self.factordraw(canvas, wd, hg, mesh, [0,0], self.forces[0], spacing, vspacing)
        self.factordraw(canvas, wd, hg, mesh, [0,0], self.forces[1], spacing, vspacing)
        self.factordraw(canvas, wd, hg, mesh, [0,0], self.forces[2], spacing, vspacing)

        canvas.pack()
        window.mainloop()


    
    forces = [] # тут хранятся все силы и моменты. Каждый элемент - переменная класса factor. 
    
    def add_a_force(self, x, y, mod, vec, name='force', type='point', len='0'): #метод добавляет приложенные к раме силы и моменты. x, y - координаты. mod - модуль. type - 'po
        # 'point' - точечно приложенная, 'lin' - равномерно распределённая. Если ввиодм не 'point' - надо указать участок len. name - тип: 'moment' или 'force'
        self.forces.append(factor(x, y, mod, vec, name, type, len))



class factor: 
    def __init__(self, x, y, mod, vec, name='force', type='point', len='0'): #vec - любой вектор из двух переменных. пока что я честно не буду их отсеивать, 
        #но для правильной работы по модулю они должнны быть от 0 до 1
        self.xpos = x
        self.ypos = y 
        self.mod = mod
        self.type = type
        self.name = name
        self.vec = vec # значения up, down, left, right
        if type == 'const': self.len = len
