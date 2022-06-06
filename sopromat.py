## тут я буду писать подключаемый файл

# для начала надо придумать как хранить плоскую балку

# балка хранится в виде точек. Кажая точка состоит из списка двух(для пространственной потом допишу) координат
from math import sqrt
import tkinter as tk
import numpy as np

def isgood(mas, type=2): #проверяет список точек на соответствие шаблонному виду
    return True


def flatten(coord): #из обычного списка делает np.array плоского вида. Возможно позже перейду полностью на них
    coord = np.array(coord)
    return coord.flatten() 

def dims(coord): #возвращает список из высоты графика по x и y. решение не самое элегантное, но пока что пойдёт с пивом
    return max(coord[:, 0]) - min(coord[:, 0]) + 0.001, max(coord[:, 1]) - min(coord[:, 1])  + 0.001

def check_zero_frac (num): 
    return num if num != 0 else 1


def linforce_checker(coord, start, finish, type): 
    if type == 1: 
        if coord <= start and finish >= start: return finish-start
        elif coord > start and coord <= finish: return coord-start
        elif coord > finish and finish >= start: return finish-start
        elif coord >= start and start >= finish: return start-finish
        elif coord < start and coord >= finish: return start-coord
        elif coord < finish and finish<=start: return start-finish
    if type == 2: 
        if coord <= start and finish >= start: return ((finish-start)/2 + start-coord)
        elif coord > start and coord <= finish: return -(coord-start)/2
        elif coord > finish and finish >= start: return -((finish-start)/2 + coord-finish)
        elif coord >= start and start >= finish: return -((start-finish)/2 + coord-start)
        elif coord < start and coord >= finish: return (start-coord)/2
        elif coord < finish and finish<=start: return (start-finish)/2 + finish-coord

def getsign(force, zero_is_ok=True): 
    if force=='unknown': return 1
    else: return (force)/abs(check_zero_frac(force)) if zero_is_ok else check_zero_frac(force)/abs(check_zero_frac(force)) 


class twodframe: #самый важный класс. Возможно стоит силы и стойки сделать подклассами(?)



    coord = [] #набор точек, которые мы задали
    dots_num = 0 #объявляю явно, чтобы не забыть, что у меня тут
    rastered_frame = [] #набор координат точек, по которым строится эпюра.
#   self.draw_scale
    forces = [] # тут хранятся все силы и моменты. Каждый элемент - переменная класса factor. 
    holders = [] #тут список всех стоек - переменных типа holder
    size = 0
    raster_num = 0

    window = 0
    canvas = 0

    solved_trigger = 0
    mesh = [0, 0]
    vspacing=0

    def __init__(self, coord, raster=100):
        if isgood(mas=coord): 
            self.coord = np.array(coord)
            self.dots_num = len(coord)
            self.make_a_raster(raster)
            self.raster_num = raster
            self.set_size()
        else: 
            print('Bad coordinates')

    
    def make_a_raster(self, numopoints=10): # разбивает рамку на набор точек - нужно для растрового рисования эпюр
        for i in range(len(self.coord)-1): 
            #print(self.coord[i], self.coord[i+1])
            xbeat=self.coord[i+1][0]-self.coord[i][0]
            ybeat=self.coord[i+1][1]-self.coord[i][1]
            for k in range(numopoints): self.rastered_frame.append([self.coord[i][0]+xbeat*(k/numopoints), self.coord[i][1]+ybeat*(k/numopoints)])
        #for k in range(numopoints): self.rastered_frame.append([self.coord[-1][0]+xbeat*(k/numopoints), self.coord[-1][1]+ybeat*(k/numopoints)])
        print(self.rastered_frame)
    
    def set_size(self, size = '900x300'): 
        self.size = size

    def getstep(self, dot0, dot1): 
        return ((dot1[0]-dot0[0])**2 + (dot1[1]-dot0[1])**2)**(1/2)
########################################################################################################################
## инициализация и всё связанное ## 

    def add_a_force(self, x, y, mod, vec=[0,0], name='force', type='point', stopx=0, stopy=0): #метод добавляет приложенные к раме силы и моменты. x, y - координаты. mod - модуль. type - 'po
        # 'point' - точечно приложенная, 'lin' - равномерно распределённая. Если ввиодм не 'point' - надо указать участок len. name - тип: 'moment' или 'force'
        self.forces.append(factor(x=x, y=y, mod=mod, vec=vec, name=name, type=type, stopx=stopx, stopy=stopy))

    def add_a_holder(self, x, y, type, vector=[0, 1], force_x='unknown', force_y='unknown', moment_cw='unknown'):#типы - hinge(2 реакции), joint(1 реакция), seal(3 реакции)
        if type=='joint': 
            if vector in ([0,1], [0, -1]): force_x, moment_cw = 0, 0
            if vector in ([1,0], [-1, 0]): force_y, moment_cw = 0, 0
        if type=='hinge': moment_cw=0
        self.add_a_force(x, y, force_x, [1,0], 'force') if force_x != 0 else 0
        self.add_a_force(x, y, force_y, [0,1], 'force') if force_y != 0 else 0
        self.add_a_force(x, y, moment_cw, 'clockwise', 'moment') if moment_cw != 0 else 0

        self.holders.append(holder(x, y, type=type, vec=vector))
    


##############################################################  графические методы  #####################################

    def start_drawing(self, mesh = [3, 3]):
        wd, hg = list(map(int, self.size.split('x')))
        self.window = tk.Tk()
        self.window.geometry(self.size)
        self.canvas = tk.Canvas(self.window, width=wd, height=hg, bg = 'white')
        self.mesh = mesh

    def basic_draw(self, pos=[1,0], text=''): #шаблон рисования простой картинки
        wd, hg = list(map(int, self.size.split('x')))
        spacing = wd//20
        self.vspacing = -hg//3
        self.getscale(self.size)
        printpos_x = (wd)//self.mesh[0]*pos[0]+spacing
        printpos_y = (hg)//self.mesh[1]*(pos[1]+1)+self.vspacing

        

        self.framedraw(self.canvas, wd, hg, pos, spacing) # три метода переезжают: теперь они будут отрисовываться мануально, для шаблонов можно создать метод
        # self.framedraw(canvas, wd, hg, mesh, [1,0], spacing, vspacing)
        # self.framedraw(canvas, wd, hg, mesh, [2,0], spacing, vspacing) #отрисовывает три рамки потому что я так сказал

        for element in self.holders: 
            self.holderdraw(self.canvas, wd, hg, pos, element, spacing)
            # self.holderdraw(canvas, wd, hg, mesh, [1,0], element, spacing, vspacing)
            # self.holderdraw(canvas, wd, hg, mesh, [2,0], element, spacing, vspacing)

        for element in self.forces: 
            self.factordraw(self.canvas, wd, hg, pos, element, spacing) if element.mod != 0 else 0
            # self.factordraw(canvas, wd, hg, mesh, [1,0], element, spacing, vspacing) if element.name=='force' else 0
            # self.factordraw(canvas, wd, hg, mesh, [2,0], element, spacing, vspacing)

        #for unkn_force in self.ret_unknown(): self.determine_a_force(unkn_force)

        [print(element.mod, element.xpos, element.ypos, element.type, element.name ) for element in self.forces] # это для отладки, потом удалить

        self.canvas.create_text(printpos_x+dims(self.coord)[0]/2*self.draw_scale, printpos_y+dims(self.coord)[0]*self.draw_scale*0.2, anchor='center', font='Arial', text=text, fill='black')

    def system_draw(self, pos, text=''): 
        wd, hg = list(map(int, self.size.split('x')))
        spacing = wd//20
        printpos_x = (wd)//self.mesh[0]*pos[0]+spacing
        printpos_y = (hg)//self.mesh[1]*(pos[1]+1)+self.vspacing

        self.framedraw(self.canvas, wd, hg, pos, spacing)
        for element in self.forces: 
            self.factordraw(self.canvas, wd, hg, pos, element, spacing, '1') if element.mod != 0 else 0

        self.canvas.create_text(printpos_x+dims(self.coord)[0]/2*self.draw_scale, printpos_y+dims(self.coord)[0]*self.draw_scale*0.2, anchor='center', font='Arial', text=text, fill='black')


    def draw_moments (self, pos=[2,0], text='', scale=1):
        self.solve() if not self.solved_trigger else None
        wd, hg = list(map(int, self.size.split('x')))
        spacing = wd//20
        self.momentdraw(self.canvas, wd, hg, pos, spacing, text, scale)

    def draw_forces (self, pos=[1, 0], type='вертикальные', text='', scale=1): 
        self.solve() if not self.solved_trigger else None
        wd, hg = list(map(int, self.size.split('x')))
        spacing = wd//20
        self.forcesdraw(self.canvas, wd, hg, pos, spacing, type, text, scale)


    def getscale(self, size): #метод выясняет, как надо растянуть рамку, чтобы она нормально выглядела в картинке нашего размера и записывает в self.draw_scale
        wd, hg = list(map(int, size.split('x')))
        frame_wd, frame_hg = dims(self.coord)
        if wd/hg < frame_wd/frame_hg: #ширина оказывается важнее
            self.draw_scale = wd/(4*frame_wd)
        else: 
            self.draw_scale = hg/(4*frame_hg)

        #print('this is getscale. the scale is {}'.format(self.draw_scale)) #

    def framedraw(self, canvas, wd, hg, place, spacing=10): 
        printpos_x = (wd)//self.mesh[0]*place[0]+spacing
        printpos_y = (hg)//self.mesh[1]*(place[1]+1)+self.vspacing

        for i in range(self.dots_num-1): 
            canvas.create_line(self.coord[i][0]*self.draw_scale + printpos_x, printpos_y-self.coord[i][1]*self.draw_scale, self.coord[i+1][0]*self.draw_scale+printpos_x, printpos_y-self.coord[i+1][1]*self.draw_scale, fill='black', width=3) #если этот кусок подогнать под 
        # рисование одной и той же рамы где угодно - можно использовать для шаблона, сверху которого рисуем эпюры

    def factordraw(self, canvas, wd, hg, place, factor, spacing, unk_sign=''): #рисует в place ячейке один конкретный фактор
        frame_wd, frame_hg = dims(self.coord)
        frame_wd, frame_hg = frame_wd*self.draw_scale, frame_hg*self.draw_scale
        truex = factor.xpos*self.draw_scale+(wd)//self.mesh[0]*place[0]+spacing
        truey = (hg)//self.mesh[1]*(place[1]+1)+self.vspacing-factor.ypos*self.draw_scale


        if factor.name == 'force' and factor.type=='point':
            canvas.create_line(truex, truey,
            truex+frame_wd//5*factor.vec[0]*getsign(factor.mod), truey-frame_wd//5*factor.vec[1]*getsign(factor.mod), 
            fill='red', width=4, arrow='last')

            canvas.create_text((truex-(factor.vec[0]==0)*frame_wd//15)+(factor.vec[1]==0)*frame_wd//15*getsign(factor.mod), truey-frame_wd//12*getsign(factor.mod), anchor='center', font='Arial', text=str(round(abs(factor.mod), 1) if factor.mod!='unknown' else unk_sign), fill='red')

        if factor.name== 'moment': 
            canvas.create_line(truex, truey, 
            truex, truey-frame_wd//10,
            fill = 'purple', width=4)
            canvas.create_arc(truex-frame_wd//10, truey-frame_wd//10, 
            truex+frame_wd//10, truey+frame_wd//10, 
            start=90, extent=-45*getsign(factor.mod), 
            width=4, outline='purple', style=tk.ARC)
            canvas.create_line(truex+frame_wd//10/sqrt(2)*getsign(factor.mod), truey-frame_wd//10/sqrt(2), 
            truex+(frame_wd//10/sqrt(2)+5)*getsign(factor.mod), truey-frame_wd//10/sqrt(2)+5,
            fill = 'purple', width=4, arrow='last')
            canvas.create_text(truex-frame_wd//10*getsign(factor.mod), truey-frame_wd//8, anchor='center', font='Arial', text=str(round(abs(factor.mod), 2) if factor.mod!='unknown' else ''), fill='purple')

        
        if factor.name == 'force' and factor.type=='lin':
            canvas.create_line(truex, truey, truex-factor.vec[0]*frame_wd//20*getsign(factor.mod), truey+factor.vec[1]*frame_wd//20*getsign(factor.mod), fill='green', width=3) #
            canvas.create_line(truex-factor.vec[0]*frame_wd//20*getsign(factor.mod)-factor.vec[1]*factor.len*self.draw_scale, truey+factor.vec[1]*frame_wd//20*getsign(factor.mod)+factor.vec[0]*factor.len*self.draw_scale,
            truex-factor.vec[1]*factor.len*self.draw_scale, truey+factor.vec[0]*factor.len*self.draw_scale,
            fill='green', width=3)

            canvas.create_line(truex-factor.vec[0]*frame_wd//20*getsign(factor.mod), truey+factor.vec[1]*frame_wd//20*getsign(factor.mod), 
            truex-factor.vec[0]*frame_wd//20*getsign(factor.mod)-factor.vec[1]*factor.len*self.draw_scale, truey+factor.vec[1]*frame_wd//20*getsign(factor.mod)+factor.vec[0]*factor.len*self.draw_scale,
            fill='green', width=3) 

            canvas.create_text(truex-factor.vec[0]*frame_wd//10*getsign(factor.mod)-factor.vec[1]*factor.len*self.draw_scale//2, truey+factor.vec[1]*frame_wd//10*getsign(factor.mod)+factor.vec[0]*factor.len*self.draw_scale//2
            , anchor='center', font='Arial', text=str(round(abs(factor.mod), 2) if factor.mod!='unknown' else ''), fill='green')


            for i in range(int(frame_wd//80)): 
                canvas.create_line(truex-factor.vec[0]*frame_wd//20*getsign(factor.mod)-factor.vec[1]*(factor.len*i/int(frame_wd//80))*self.draw_scale, truey+factor.vec[1]*frame_wd//20*getsign(factor.mod)+factor.vec[0]*(factor.len*i/int(frame_wd//80))*self.draw_scale,
                truex-factor.vec[1]*(factor.len*i/int(frame_wd//80))*self.draw_scale, truey+factor.vec[0]*(factor.len*i/int(frame_wd//80))*self.draw_scale,
                fill='green', width=3, arrow='last', arrowshape="6 5 2")

    def holderdraw(self, canvas, wd, hg, place, hold, spacing): 
        truex = hold.xpos*self.draw_scale+(wd)//self.mesh[0]*place[0]+spacing
        truey = (hg)//self.mesh[1]*(place[1]+1)+self.vspacing - hold.ypos*self.draw_scale

        if hold.type=='hinge': 
            canvas.create_polygon(truex, truey, truex+wd//120, truey+wd//120, truex-wd//120, truey+wd//120, outline='black', width=3, fill='white')

            canvas.create_oval(truex-wd//280, truey-wd//280, 
            truex+wd//280, truey+wd//280, 
            width=2, outline='black', fill='white')

            canvas.create_line(truex-wd//90, truey+wd//120, truex+wd//90, truey+wd//120, width=3, fill='black')
            canvas.create_line(truex-wd//90, truey+wd//120, truex-wd//90-wd//180, truey+wd//120+wd//180, width=3, fill='black')
            canvas.create_line(truex, truey+wd//120, truex-wd//180, truey+wd//120+wd//180, width=3, fill='black')
            canvas.create_line(truex+wd//90, truey+wd//120, truex+wd//90-wd//180, truey+wd//120+wd//180, width=3, fill='black')
        
        if hold.type=='joint': 
            canvas.create_line(truex-wd//90, truey+wd//120*hold.vec[1], truex+wd//90, truey+wd//120*hold.vec[1], width=3, fill='black')
            canvas.create_line(truex-wd//90, truey+wd//120*hold.vec[1], truex-wd//90-wd//180, truey+wd//120*hold.vec[1]+wd//180*hold.vec[1], width=3, fill='black')
            canvas.create_line(truex, truey+wd//120*hold.vec[1], truex-wd//180, truey+wd//120*hold.vec[1]+wd//180*hold.vec[1], width=3, fill='black')
            canvas.create_line(truex+wd//90, truey+wd//120*hold.vec[1], truex+wd//90-wd//180, truey+wd//120*hold.vec[1]+wd//180*hold.vec[1], width=3, fill='black')

            canvas.create_oval(truex-wd//280, truey-wd//280, 
            truex+wd//280, truey+wd//280, 
            width=2, outline='black', fill='white')

            canvas.create_line(truex, truey, truex-wd//120*hold.vec[0], truey+wd//120*hold.vec[1], width=3, fill='black')

            canvas.create_oval(truex-wd//120*hold.vec[0]-wd//280, truey+wd//120*hold.vec[1]-wd//280, 
            truex-wd//120*hold.vec[0]+wd//280, truey+wd//120*hold.vec[1]+wd//280, 
            width=2, outline='black', fill='white')

        if hold.type=='seal': 
            canvas.create_line(truex+wd//90*hold.vec[1], truey+wd//90*hold.vec[0], truex-wd//90*hold.vec[1], truey-wd//90*hold.vec[0], width=3, fill='black')

            canvas.create_line(truex+wd//90*hold.vec[1], truey+wd//90*hold.vec[0], truex+wd//90*hold.vec[1]-wd//180*(1 if hold.vec[0] != -1 else -1), truey+wd//90*hold.vec[0]+wd//180*(1 if hold.vec[1] != -1 else -1), width=3, fill='black')
            canvas.create_line(truex, truey, truex-wd//180*(1 if hold.vec[0] != -1 else -1), truey+wd//180*(1 if hold.vec[1] != -1 else -1), width=3, fill='black')
            canvas.create_line(truex-wd//90*hold.vec[1], truey-wd//90*hold.vec[0], truex-wd//90*hold.vec[1]-wd//180*(1 if hold.vec[0] != -1 else -1), truey-wd//90*hold.vec[0]+wd//180*(1 if hold.vec[1] != -1 else -1), width=3, fill='black')

    def momentdraw(self, canvas, wd, hg, place, spacing, text, scale): #рисует эпюру моментов 
        #print('momentdraw is active')
        #print('моменты будут создаваться силами величиной {}'.format([element.mod for element in self.forces]))
        #print('сейчас я выведу все точки маршрута: {}\nА теперь я буду выводить моменты'.format(self.rastered_frame))
        printpos_x = (wd)//self.mesh[0]*place[0]+spacing
        printpos_y = (hg)//self.mesh[1]*(place[1]+1)+self.vspacing

        prev_dot = 0

        counts = []

        ctr = 0

        prevdrawdot = self.rastered_frame[0][0]*self.draw_scale + printpos_x, printpos_y-self.rastered_frame[0][1]*self.draw_scale
        for dot in self.rastered_frame: 
            if prev_dot == 0: 
                prev_dot = dot
                continue
            for factor in self.forces: 
                if factor not in counts and((factor.xpos >= prev_dot[0] and factor.xpos <= dot[0]) or (factor.xpos <= prev_dot[0] and factor.xpos >= dot[0]) ) and ( (factor.ypos >= prev_dot[1] and factor.ypos <= dot[1])or (factor.ypos <= prev_dot[1] and factor.ypos >= dot[1])): 
                    counts.append(factor)

            summom = sum([self.moment_count(factor, dot[0], dot[1]) for factor in counts])

            canvas.create_line(prevdrawdot[0], prevdrawdot[1], 
            dot[0]*self.draw_scale-self.draw_scale**-0.01*summom*scale*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale**-.01*summom*scale*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))), fill='blue', width = wd//500)
            

            prevdrawdot = [dot[0]*self.draw_scale-self.draw_scale**-.01*summom*scale*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale**-.01*summom*scale*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)))]

            if (ctr+2)%(round(self.raster_num//5)) == 0: canvas.create_line(dot[0]*self.draw_scale + printpos_x, printpos_y - dot[1]*self.draw_scale, prevdrawdot[0], prevdrawdot[1], fill='blue', width = wd//500)
            
            if (ctr)%(self.raster_num) == 0 and round(summom, 1)!=0: 
                canvas.create_text(dot[0]*self.draw_scale+printpos_x + wd//50, -dot[1]*self.draw_scale+printpos_y-wd//90, anchor='center', font='Arial', text=str(round(summom, 2)), fill='black')


            prev_dot = dot

            ctr+=1

            print('сумма моментов - {}, конкретные моменты - {}, создаваемые силами: {}'.format(summom, [self.moment_count(factor, dot[0], dot[1]) for factor in counts], [element.mod for element in counts]))
        self.basic_draw(place)

        canvas.create_text(printpos_x+dims(self.coord)[0]/2*self.draw_scale, printpos_y+dims(self.coord)[0]*self.draw_scale*0.2, anchor='center', font='Arial', text=text, fill='black')

    def forcesdraw(self, canvas, wd, hg, place, spacing, type, text, scale): # рисует эпюру сил - можно выбрать "продольные", "поперечные", "вертикальные" и "горизонтальные"
        printpos_x = (wd)//self.mesh[0]*place[0]+spacing
        printpos_y = (hg)//self.mesh[1]*(place[1]+1)+self.vspacing

        counts = []
        prev_dot = 0
        ctr = 0

        prevdrawdot = self.rastered_frame[0][0]*self.draw_scale + printpos_x, printpos_y-self.rastered_frame[0][1]*self.draw_scale
        for dot in self.rastered_frame: 
            if prev_dot == 0: 
                prev_dot = dot
                continue
            for factor in self.forces: 
                if factor not in counts and( (factor.xpos >= prev_dot[0] and factor.xpos <= dot[0]) or (factor.xpos <= prev_dot[0] and factor.xpos >= dot[0]) ) and ( (factor.ypos >= prev_dot[1] and factor.ypos <= dot[1])or (factor.ypos <= prev_dot[1] and factor.ypos >= dot[1]) ): 
                    counts.append(factor)
                    
            if type == 'вертикальные': 
                sumforce = sum([self.getforce(fac, dot[0], dot[1], vec2=[0,1]) for fac in counts])
            if type == 'горизонтальные': 
                sumforce = sum([self.getforce(fac, dot[0], dot[1], vec2=[1,0]) for fac in counts])
            if type == 'поперечные': 
                sumforce = sum([(dot[0]-prev_dot[0]*self.getforce(fac, dot[0], dot[1], vec2=[0,1])/(dot[0]-prev_dot[0]) + dot[1]-prev_dot[1]*self.getforce(fac, dot[0], dot[1], vec2=[1,0])/(dot[1]-prev_dot[1] if (dot[1]-prev_dot[1]) != 0 else 1)) for fac in counts])
            if type == 'продольные': 
                sumforce = sum([(dot[1]-prev_dot[1])*self.getforce(fac, dot[0], dot[1], vec2=[0,1])/check_zero_frac(dot[1]-prev_dot[1]) + (dot[0]-prev_dot[0])*self.getforce(fac, dot[0], dot[1], vec2=[1,0])/check_zero_frac(dot[0]-prev_dot[0]) for fac in counts])
            

            canvas.create_line(prevdrawdot[0], prevdrawdot[1], 
            dot[0]*self.draw_scale+self.draw_scale**0.1*sumforce*scale*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale+self.draw_scale**0.1*sumforce*scale*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))), fill='orange', width = wd//500)
            
            
            prevdrawdot = [dot[0]*self.draw_scale+self.draw_scale**0.1*sumforce*scale*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale+self.draw_scale**0.1*sumforce*scale*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt(check_zero_frac((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)))]

            if (ctr+2)%(round(self.raster_num/5)) == 0: canvas.create_line(dot[0]*self.draw_scale + printpos_x, printpos_y - dot[1]*self.draw_scale, prevdrawdot[0], prevdrawdot[1], fill='orange', width = wd//500)
            
            if (ctr)%(self.raster_num) == 0 and round(sumforce, 1)!=0: 
                canvas.create_text(dot[0]*self.draw_scale+printpos_x + wd//50, -dot[1]*self.draw_scale+printpos_y-wd//120, anchor='center', font='Arial', text=str(round(sumforce, 2)), fill='black')
            

            #print('\n', prevdrawdot, dot)
            #print('сумма сил - {}, создаваемая силами: {}'.format(sumforce, [self.getforce(fac, dot[0], dot[1], vec2=[0,1]) for fac in counts]))
            prev_dot = dot

            ctr+=1
        self.basic_draw(place)
        canvas.create_text(printpos_x+dims(self.coord)[0]/2*self.draw_scale, printpos_y+dims(self.coord)[0]*self.draw_scale*0.2, anchor='center', font='Arial', text=text, fill='black')

    def print_it(self):
        self.canvas.pack()
        self.window.mainloop()

##############################  methods working with forces  ########################################

    def getforce(self, factor, x, y, vec2=[0,1]): # в данной точке выясняем, какую силу данный фактор. vec2 - это единичный(для правильной работы - по длине) вектор, который определяет направление. 
        # то есть для получения вертикальной составляющей нужно указать vec2=[0,1], для горизонтальной - [1,0]
        if factor.name == 'force' and factor.type == 'point': 
            return factor.mod*factor.vec[0]*vec2[0] + factor.mod*factor.vec[1]*vec2[1]
        elif factor.name == 'force' and factor.type == 'lin': 
            Y = factor.mod*factor.vec[1]*linforce_checker(y, factor.ypos, factor.stopy, 1)*vec2[1]
            X = factor.mod*factor.vec[0]*linforce_checker(x, factor.xpos, factor.stopx, 1)*vec2[0]
            return X+Y
        else: return 0


    def ret_unknown(self):  # возвращает список сил с значением модуля "unknown"
        ret = []
        for element in self.forces: 
            if element.mod=='unknown': ret.append(element)
        return ret

    def ret_known(self): #позвращает список сил с известным значением модуля
        ret = []
        for element in self.forces: 
            if element.mod!='unknown': ret.append(element)

        print(ret)
        return ret

############################  moments  ####################################################

    
    def moment_count(self, force, x, y, is_unknown=False):# находит значение момента от передаваемого фактора в передаваемой точке
        force.mod = 1 if is_unknown else force.mod 
        if force.name=='moment': return force.mod
        elif force.type=='point': 
            #print('сила расположена в [{}, {}], в точке [{}, {}] создает момент {}'.format(force.xpos, force.ypos, x, y, force.mod*(x-force.xpos)*force.vec[1]-force.mod*(y-force.ypos)*force.vec[0]))
            #print('moment {} from point force'.format(force.mod*(x-force.xpos)*force.vec[1]-force.mod*(y-force.ypos)*force.vec[0]))
            return force.mod*(x-force.xpos)*force.vec[1]-force.mod*(y-force.ypos)*force.vec[0]
        else: 
            mx = -force.mod*force.vec[1]*linforce_checker(x, force.xpos, force.stopx, 1)*linforce_checker(x, force.xpos, force.stopx, 2)
            my = -force.mod*force.vec[0]*linforce_checker(y, force.ypos, force.stopy, 1)*linforce_checker(y, force.ypos, force.stopy, 2)
            print(mx, 'its momentcounters msg, {} is mod, {} is vec and {} is len making it {}, where mx is {} and my is {}'.format(force.mod, force.vec, linforce_checker(y, force.ypos, force.stopy, 2), mx+my, mx, my))
            
            return mx+my

##############################  solving  ######################################################################

    def solve(self): #эта система занимается тем, что решает систему и записывает значения факторов в них, поэтому ей не надо ничего писать.\ В V0.0.2 переехала и теперь её надо руками вызывать пользователю, проверок нет ни на что
        if not self.iscomplex():
            #в соучае если система статически определима, можно посчитать через систему из 3 уравнений, где надо стоставить уравнение моментов в начале рамы, 
            dot = self.coord[0] #это записываю точку, в которой всё считаю                # уравнение вертикальных сил и уравнение горизонтальных сил
            known_forces = self.ret_known()
            unknown_forces = self.ret_unknown() #в переменной unknown forces храню все силы с неизвестными модулями. Они на следующем этапе буут перезаписаны как 
 # как единичные, а к концу функции в них будут записаны найденные значения
            mat = np.array([[self.moment_count(element, dot[0], dot[1], is_unknown=True) for element in unknown_forces]])   #это я создал матрицу, которую будем решать. Сейчас там по одному коэффициенту неизвестного(одна строчка)
            vec = np.array([-1*sum([self.moment_count(element, dot[0], dot[1]) for element in known_forces])]) #это вектор. Там сейчас сумма моментов от известных сил на минус 1

            mat = np.append(mat, [[element.mod*element.vec[1] for element in unknown_forces]], axis=0) #ДОБАВЛЯЮ вторую строчку к матрице - вертикальные силы. находятся как неизвестная сила(равна одному) на y-составляющую вектора направления
            vec = np.append(vec, -1*sum(element.mod*element.vec[1]*(abs(element.len) if element.type=='lin' else 1)*(element.name=='force') for element in known_forces)) # добавляю сумму вертикальных сил на минус один
            
            mat = np.append(mat, [[element.mod*element.vec[0] for element in unknown_forces]], axis=0) #то же самое для горизонтальных
            vec = np.append(vec, -1*sum(element.mod*element.vec[0]*(abs(element.len) if element.type=='lin' else 1)*(element.name=='force') for element in known_forces))

            print('[{}] - this is matrix'.format(mat))
            print('[{}] - this is vector'.format(vec))
            answer = np.linalg.solve(mat, vec)
            print('answer is {}'.format)
            for i in range(len(unknown_forces)): unknown_forces[i].mod = answer[i]
        else: ################################################## the complex frame option
            print('it\'s hard but ill manage')

            known_forces = self.ret_known()
            unknown_forces = self.ret_unknown()
            pos = [self.holders[0].xpos, self.holders[0].ypos]
            
            int = np.array([[0.0 for _ in unknown_forces] for _ in unknown_forces])
            vec = np.array([0.0 for _ in unknown_forces])
            for i in range(len(unknown_forces)):
                prevdot=[0,0]
                for dot in self.rastered_frame: 
                    selfmoment = self.moment_count(unknown_forces[i], dot[0], dot[1], True)
                    known_forces_moment = sum([self.moment_count(element, dot[0], dot[1]) for element in known_forces])
                    for j in range(len(unknown_forces)): 
                        int[i][j] += (selfmoment*self.moment_count(unknown_forces[j], dot[0], dot[1], True)+self.moment_count(unknown_forces[i], prevdot[0], prevdot[1], True)*self.moment_count(unknown_forces[j], dot[0], dot[1], True))*self.getstep(dot, prevdot)*1/2
                    vec[i] += -(selfmoment*known_forces_moment+self.moment_count(unknown_forces[i], prevdot[0], prevdot[1], True)*sum([self.moment_count(element, prevdot[0], prevdot[1]) for element in known_forces]))*self.getstep(dot, prevdot)*1/2
                    prevdot = dot
                    #print(int, selfmoment)
                    print(self.getstep(dot, prevdot))
            print(int, '\n', vec)

            
            intcopy = int
            for i in np.arange(len(int)): 
                if sum(intcopy[i]) == 0:
                    print('ok')
                    int = np.delete(int, i, axis=0)
                    int = np.delete(int, i, axis=1)
                    vec = np.delete(vec, i)
                    unknown_forces[i].mod = 0
                    unknown_forces = np.delete(unknown_forces, i)
                    i += -1
            print(vec)
            print(int)
            answer=np.linalg.solve(int, vec)

            print('\n', answer)

            for i in range(len(unknown_forces)): unknown_forces[i].mod = answer[i]

            



        self.solved_trigger=1



    def iscomplex(self, num=False): # выясняет, является ли данная система статически неопределимой. Если неопределима - вернёт True
        unknown_num = 0
        for element in self.holders: 
            if element.type == 'joint': unknown_num+=1
            elif element.type == 'hinge': unknown_num+=2
            else: unknown_num+=3
            print('Это  iscomplex. Я встретил стойку типа {} и записал, что связей уже {}'.format(element.type, unknown_num))
        return unknown_num > 3 if not num else unknown_num-3

class threedframe: 
    coord = [] #точки рамки
    rater_num=0 #точность разбиения
    dots_num=0 # len of the coord var
    size = '' # screen size
    wd, hg= 0,0 # screen dims
    window = 0 # contains the window
    canvas=0 #canvas in the window
    mesh=0 # сетка разбиения картинки
    hspacing, wspacing=0,0 #spacing between drawings
    zdrawdegr=np.pi/6
    rastered_frame=[]
    forces=[]
    holders=[]
    
    def __init__(self, coords, raster=10):

        

        if isgood(coords, 3): self.coord = np.array(coords)
        else: return False
        self.raster_num = raster
        self.dots_num = len(coords)
        self.make_a_raster()
        self.set_size()
        self.getscale(self.size)


    def getscale(self, size): #метод выясняет, как надо растянуть рамку, чтобы она нормально выглядела в картинке нашего размера и записывает в self.draw_scale
        wd, hg = list(map(int, size.split('x')))
        self.frame_wd, self.frame_hg = dims(np.delete(self.coord, 2, axis=1))
        if wd/hg < self.frame_wd/self.frame_hg: #ширина оказывается важнее
            self.draw_scale = wd/(4*self.frame_wd)
        else: 
            self.draw_scale = hg/(4*self.frame_hg)

    def make_a_raster(self): 
        for i in range(len(self.coord)-1):  
            xbeat=self.coord[i+1][0]-self.coord[i][0]
            ybeat=self.coord[i+1][1]-self.coord[i][1]
            zbeat=self.coord[i+1][1]-self.coord[i][1]
            for k in range(self.raster_num): self.rastered_frame.append([self.coord[i][0]+xbeat*(k/self.raster_num), self.coord[i][1]+ybeat*(k/self.raster_num), self.coord[i][1]+zbeat*(k/self.raster_num)])
    def set_size(self, size = '900x300'): 
        self.size = size
        self.getscale(size)

    def start_drawing(self, mesh = [1, 1]):
        self.wd, self.hg = list(map(int, self.size.split('x')))
        self.window = tk.Tk()
        self.window.geometry(self.size)
        self.canvas = tk.Canvas(self.window, width=self.wd, height=self.hg, bg = 'white')
        self.mesh = mesh
        self.hspacing, self.wspacing= self.hg//(10*mesh[1]), self.wd//(10*mesh[0]) 

    def framedraw(self, pos, text=''): 
        printpos_x = (self.wd)//self.mesh[0]*pos[0]+self.wspacing #просто координаты левого верхнего угла бокса, в котором рисуем
        printpos_y = (self.hg)//self.mesh[1]*(pos[1]+0.5)

        for i in range(self.dots_num-1): 
            self.canvas.create_line(self.coord[i][0]*self.draw_scale + self.coord[i][2]*self.draw_scale*np.cos(self.zdrawdegr) + printpos_x, printpos_y-self.coord[i][1]*self.draw_scale-self.coord[i][2]*self.draw_scale*np.sin(self.zdrawdegr), 
            self.coord[i+1][0]*self.draw_scale+self.coord[i+1][2]*self.draw_scale*np.cos(self.zdrawdegr)+printpos_x, printpos_y-self.coord[i+1][1]*self.draw_scale-self.coord[i+1][2]*self.draw_scale*np.sin(self.zdrawdegr), fill='black', width=3) #если этот кусок подогнать под 
            print()# рисование одной и той же рамы где угодно - можно использовать для шаблона, сверху которого рисуем эпюры
            
    def add_a_force(self, x, y, z, mod, vec=[0,0], name='force', type='point', stopx=0, stopy=0): #метод добавляет приложенные к раме силы и моменты. x, y - координаты. mod - модуль. type - 'po
        # 'point' - точечно приложенная, 'lin' - равномерно распределённая. Если ввиодм не 'point' - надо указать участок len. name - тип: 'moment' или 'force'
        self.forces.append(factor(x=x, y=y, z=z,  mod=mod, vec=vec, name=name, type=type, stopx=stopx, stopy=stopy))

    def add_a_holder(self, x, y, z, type, vector=[0, 1], force_x='unknown', force_y='unknown', force_z='unknown', moment_cw='unknown'):#типы - hinge(2 реакции), joint(1 реакция), seal(3 реакции)
        if type=='joint': 
            if vector in ([0,1,0], [0, -1,0]): force_x, moment_cw, force_z = 0, 0, 0
            if vector in ([1,0,0], [-1, 0,0]): force_y, moment_cw, force_z = 0, 0, 0
            if vector in ([0,0,1], [0,0,-1]): force_y, moment_cw, force_x = 0, 0, 0
        if type=='hinge': moment_cw=0
        self.add_a_force(x, y, z, force_x, [1,0], 'force') if force_x != 0 else 0
        self.add_a_force(x, y, z, force_y, [0,1], 'force') if force_y != 0 else 0
        self.add_a_force(x, y, z, moment_cw, 'clockwise', 'moment') if moment_cw != 0 else 0

        self.holders.append(holder(x, y, z, type=type, vec=vector))

    def base_draw(self, pos, text=''): 
        self.framedraw(pos, text)

        for holder in self.holders: 
            self.holderdraw(pos, holder)

    def holderdraw(self, pos, holder): 
        printpos_x = (self.wd)//self.mesh[0]*pos[0]+self.wspacing #просто координаты левого верхнего угла бокса, в котором рисуем
        printpos_y = (self.hg)//self.mesh[1]*(pos[1]+0.5)

        self.canvas.create_polygon(printpos_x+(holder.xpos-np.cos(self.zdrawdegr)*self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos-np.sin(self.zdrawdegr)*self.frame_hg/20)*self.draw_scale+20, 
        printpos_x+(holder.xpos+np.cos(self.zdrawdegr)*self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos+np.sin(self.zdrawdegr)*self.frame_hg/20)*self.draw_scale+20,
        printpos_x+(holder.xpos+np.cos(self.zdrawdegr)*self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos+np.sin(self.zdrawdegr)*self.frame_hg/20)*self.draw_scale-20,
        printpos_x+(holder.xpos-np.cos(self.zdrawdegr)*self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos-np.sin(self.zdrawdegr)*self.frame_hg/20)*self.draw_scale-20, 
        outline='black', width=3, fill='')

        print(printpos_x+(holder.xpos-self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos-self.frame_hg/20)*self.draw_scale, 
        printpos_x+(holder.xpos+self.frame_wd/20)*self.draw_scale, printpos_y-(holder.ypos+self.frame_hg/20)*self.draw_scale)

    def print_it(self):
        self.canvas.pack()
        self.window.mainloop()

class factor: 
    def __init__(self, x, y, z=0, mod=0, vec=[0, 0], name='force', type='point', stopx = 0, stopy = 0): #vec - любой вектор из двух переменных. пока что я честно не буду их отсеивать, 
        #но для правильной работы по модулю они должнны быть от 0 до 1
        #если сила распределённая, указывается точка начала, направление действия и точка окончания - пока что только равномерная сила (V 0.0.2). проверок нет.
        self.xpos = x
        self.ypos = y 
        self.mod = mod
        self.type = type
        self.name = name
        self.vec = vec # значения up, down, left, right. Если указываем момент, принимает clockwise или counterclockwise
        self.stopx = stopx
        self.stopy = stopy
        self.len = getsign(x-stopx, False)*getsign(y-stopy, False)*sqrt((x-stopx)**2 + (y-stopy)**2) if type=='lin' else 0
        
    


class holder: 
    def __init__(self, x, y, z=0, type='hinge', vec=[0,1]):
        self.xpos=x
        self.ypos=y
        self.xpos=z
        self.type=type
        self.vec=vec