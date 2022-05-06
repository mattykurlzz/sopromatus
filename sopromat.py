## тут я буду писать подключаемый файл

# для начала надо придумать как хранить плоскую балку

# балка хранится в виде точек. Кажая точка состоит из списка двух(для пространственной потом допишу) координат
from math import sqrt
import tkinter as tk
import numpy as np

def isgood(lst): #проверяет список точек на соответствие шаблонному виду
    if type(lst) is tuple or list: 
        for element in lst: 
            if type(element) is tuple or list and len(element) == 2: continue
            else: return False
        return True
    else: return False

def flatten(coord): #из обычного списка делает np.array плоского вида. Возможно позже перейду полностью на них
    coord = np.array(coord)
    return coord.flatten() 

def dims(coord): #возвращает список из высоты графика по x и y. решение не самое элегантное, но пока что пойдёт с пивом
    return max(coord[:, 0]) - min(coord[:, 0]) + 0.001, max(coord[:, 1]) - min(coord[:, 1])  + 0.001



class frame: #самый важный класс. Возможно стоит силы и стойки сделать подклассами(?)



    coord = [] #набор точек, которые мы задали
    dots_num = 0 #объявляю явно, чтобы не забыть, что у меня тут
    rastered_frame = [] #набор координат точек, по которым строится эпюра.
#   self.draw_scale
    forces = [] # тут хранятся все силы и моменты. Каждый элемент - переменная класса factor. 
    holders = [] #тут список всех стоек - переменных типа holder


    def __init__(self, coord):
        if isgood(coord): 
            self.coord = np.array(coord)
            self.dots_num = len(coord)
            self.make_a_raster(100) #стоит дать выбор челикам, но )))
        else: 
            print('Bad coordinates')

    def getscale(self, size): #метод выясняет, как надо растянуть рамку, чтобы она нормально выглядела в картинке нашего размера и записывает в self.draw_scale
        wd, hg = list(map(int, size.split('x')))
        frame_wd, frame_hg = dims(self.coord)
        if wd/hg < frame_wd/frame_hg: #ширина оказывается важнее
            self.draw_scale = wd/(4*frame_wd)
        else: 
            self.draw_scale = hg/(4*frame_hg)

        #print('this is getscale. the scale is {}'.format(self.draw_scale)) #

    def framedraw(self, canvas, wd, hg, mesh, place, spacing=10, vspacing=90): 
        printpos_x = (wd)//mesh[0]*place[0]+spacing
        printpos_y = (hg)//mesh[1]*(place[1]+1)-vspacing

        for i in range(self.dots_num-1): 
            canvas.create_line(self.coord[i][0]*self.draw_scale + printpos_x, printpos_y-self.coord[i][1]*self.draw_scale, self.coord[i+1][0]*self.draw_scale+printpos_x, printpos_y-self.coord[i+1][1]*self.draw_scale, fill='black', width=3) #если этот кусок подогнать под 
        # рисование одной и той же рамы где угодно - можно использовать для шаблона, сверху которого рисуем эпюры

    def factordraw(self, canvas, wd, hg, mesh, place, factor, spacing, vspacing): #рисует в place ячейке один конкретный фактор
        frame_wd, frame_hg = dims(self.coord)
        frame_wd, frame_hg = frame_wd*self.draw_scale, frame_hg*self.draw_scale
        truex = factor.xpos*self.draw_scale+(wd)//mesh[0]*place[0]+spacing
        truey = (hg)//mesh[1]*(place[1]+1)-vspacing - factor.ypos*self.draw_scale

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
            canvas.create_line(truex+frame_wd//10/sqrt(2)*(1 if factor.vec=='clockwise' else -1), truey-frame_wd//10/sqrt(2), 
            truex+(frame_wd//10/sqrt(2)+5)*(1 if factor.vec=='clockwise' else -1), truey-frame_wd//10/sqrt(2)+5,
            fill = 'blue', width=4, arrow='last')
        
        if factor.name == 'force' and factor.type=='lin':
            canvas.create_line(truex, truey, truex-factor.vec[0]*frame_wd//20, truey+factor.vec[1]*frame_wd//20, fill='green', width=3)
            canvas.create_line(truex-factor.vec[0]*frame_wd//20+factor.vec[1]*factor.len*self.draw_scale, truey+factor.vec[1]*frame_wd//20-factor.vec[0]*factor.len*self.draw_scale,
            truex+factor.vec[1]*factor.len*self.draw_scale, truey-factor.vec[0]*factor.len*self.draw_scale,
            fill='green', width=3)

            canvas.create_line(truex-factor.vec[0]*frame_wd//20, truey+factor.vec[1]*frame_wd//20, 
            truex-factor.vec[0]*frame_wd//20+factor.vec[1]*factor.len*self.draw_scale, truey+factor.vec[1]*frame_wd//20-factor.vec[0]*factor.len*self.draw_scale,
            fill='green', width=3) 

            for i in range(int(frame_wd//80)): 
                canvas.create_line(truex-factor.vec[0]*frame_wd//20+factor.vec[1]*(factor.len*i/int(frame_wd//80))*self.draw_scale, truey+factor.vec[1]*frame_wd//20-factor.vec[0]*(factor.len*i/int(frame_wd//80))*self.draw_scale,
                truex+factor.vec[1]*(factor.len*i/int(frame_wd//80))*self.draw_scale, truey-factor.vec[0]*(factor.len*i/int(frame_wd//80))*self.draw_scale,
                fill='green', width=3, arrow='last', arrowshape="6 5 2")

    def getforce(self, factor, x, y, vec2=[0,1]): # в данной точке выясняем, какую силу данный фактор. vec2 - это единичный(для правильной работы - по длине) вектор, который определяет направление. 
        # то есть для получения вертикальной составляющей нужно указать vec2=[0,1], для горизонтальной - [1,0]
        if factor.name == 'force' and factor.type == 'point': 
            return factor.mod*factor.vec[0]*vec2[0] + factor.mod*factor.vec[1]*vec2[1]
        elif factor.name == 'force' and factor.type == 'lin': 
            Y = factor.mod*factor.vec[1]*abs(factor.len-(factor.xpos+factor.len-x if x<=factor.xpos+factor.len and x>=factor.xpos else 0))*vec2[1]
            X = factor.mod*factor.vec[0]*abs(factor.len-(factor.ypos-factor.len-y if y>factor.xpos and factor.xpos+factor.len else 0))*vec2[0]
            return X+Y
        else: return 0


    def holderdraw(self, canvas, wd, hg, mesh, place, hold, spacing, vspacing): 
        truex = hold.xpos*self.draw_scale+(wd)//mesh[0]*place[0]+spacing
        truey = (hg)//mesh[1]*(place[1]+1)-vspacing - hold.ypos*self.draw_scale

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


    def momentdraw(self, canvas, wd, hg, mesh, place, spacing, vspacing): #рисует эпюру моментов 
        print('momentdraw is active')
        print('моменты будут создаваться силами величиной {}'.format([element.mod for element in self.forces]))
        print('сейчас я выведу все точки маршрута: {}\nА теперь я буду выводить моменты'.format(self.rastered_frame))
        printpos_x = (wd)//mesh[0]*place[0]+spacing
        printpos_y = (hg)//mesh[1]*(place[1]+1)-vspacing

        prev_dot = 0

        counts = []

        ctr = 0

        prevdrawdot = self.rastered_frame[0][0]*self.draw_scale + printpos_x, printpos_y-self.rastered_frame[0][1]*self.draw_scale
        for dot in self.rastered_frame: 
            if prev_dot == 0: 
                prev_dot = dot
                continue
            for factor in self.forces: 
                if factor not in counts and( (factor.xpos >= prev_dot[0] and factor.xpos <= dot[0]) or (factor.xpos <= prev_dot[0] and factor.xpos >= dot[0]) ) and ( (factor.ypos >= prev_dot[1] and factor.ypos <= dot[1])or (factor.ypos <= prev_dot[1] and factor.ypos >= dot[1]) ): 
                    counts.append(factor)

            summom = sum([self.moment_count(factor, dot[0], dot[1]) for factor in counts])

            canvas.create_line(prevdrawdot[0], prevdrawdot[1], 
            dot[0]*self.draw_scale-self.draw_scale//25*summom*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale//25*summom*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)), fill='blue', width = wd//500)
            

            prevdrawdot = [dot[0]*self.draw_scale-self.draw_scale//25*summom*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale//25*summom*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))]

            if ctr%5 == 0: canvas.create_line(dot[0]*self.draw_scale + printpos_x, printpos_y - dot[1]*self.draw_scale, prevdrawdot[0], prevdrawdot[1], fill='blue', width = wd//500)
            
            
            prev_dot = dot

            ctr+=1

            print('сумма моментов - {}, конкретные моменты - {}, создаваемые силами: {}'.format(summom, [self.moment_count(factor, dot[0], dot[1]) for factor in counts], [element.mod for element in counts]))

    def forcesdraw(self, canvas, wd, hg, mesh, place, spacing, vspacing, type='вертикальные'): # рисует эпюру сил - можно выбрать "продольные", "поперечные", "вертикальные" и "горизонтальные"
        printpos_x = (wd)//mesh[0]*place[0]+spacing
        printpos_y = (hg)//mesh[1]*(place[1]+1)-vspacing

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
                sumforce = sum([dot[1]-prev_dot[1]*self.getforce(fac, dot[0], dot[1], vec2=[0,1])/(dot[1]-prev_dot[1]) + dot[0]-prev_dot[0]*self.getforce(fac, dot[0], dot[1], vec2=[1,0])/(dot[0]-prev_dot[0]) for fac in counts])
            
            canvas.create_line(prevdrawdot[0], prevdrawdot[1], 
            dot[0]*self.draw_scale-self.draw_scale//25*sumforce*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale//25*sumforce*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)), fill='orange', width = wd//500)
            

            prevdrawdot = [dot[0]*self.draw_scale-self.draw_scale//25*sumforce*np.sin((np.pi/2)*(dot[1]-prev_dot[1])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2)) + printpos_x, 
            printpos_y - dot[1]*self.draw_scale-self.draw_scale//25*sumforce*np.sin((np.pi/2)*(dot[0]-prev_dot[0])/sqrt((dot[1]-prev_dot[1])**2 + (dot[0]-prev_dot[0])**2))]

            if ctr%5 == 0: canvas.create_line(dot[0]*self.draw_scale + printpos_x, printpos_y - dot[1]*self.draw_scale, prevdrawdot[0], prevdrawdot[1], fill='orange', width = wd//500)
            
            print('сумма сил - {}, создаваемая силами: {}'.format(sumforce, [self.getforce(fac, dot[0], dot[1], vec2=[0,1]) for fac in counts]))
            prev_dot = dot

            ctr+=1
            


    def basic_draw(self, size = '900x300', mesh=[3,1]): #шаблон рисования простой трёхэлементной картинки
        wd, hg = list(map(int, size.split('x')))
        spacing = wd//20
        vspacing = hg//5
        self.getscale(size)

        window = tk.Tk()
        window.geometry(size)
        canvas = tk.Canvas(window, width=wd, height=hg, bg = 'white')

        self.framedraw(canvas, wd, hg, mesh, [0,0], spacing, vspacing)
        self.framedraw(canvas, wd, hg, mesh, [1,0], spacing, vspacing)
        self.framedraw(canvas, wd, hg, mesh, [2,0], spacing, vspacing) #отрисовывает три рамки потому что я так сказал

        for element in self.holders: 
            self.holderdraw(canvas, wd, hg, mesh, [0,0], element, spacing, vspacing)
            self.holderdraw(canvas, wd, hg, mesh, [1,0], element, spacing, vspacing)
            self.holderdraw(canvas, wd, hg, mesh, [2,0], element, spacing, vspacing)

        for element in self.forces: 
            self.factordraw(canvas, wd, hg, mesh, [0,0], element, spacing, vspacing)
            self.factordraw(canvas, wd, hg, mesh, [1,0], element, spacing, vspacing) if element.name=='force' else 0
            self.factordraw(canvas, wd, hg, mesh, [2,0], element, spacing, vspacing)

        #for unkn_force in self.ret_unknown(): self.determine_a_force(unkn_force)
        self.solve()

        [print(element.mod, element.xpos, element.ypos ) for element in self.forces]

        self.momentdraw(canvas, wd, hg, mesh, [2,0], spacing, vspacing)
        self.forcesdraw(canvas, wd, hg, mesh, [1,0], spacing, vspacing, type='вертикальные')

        canvas.pack()
        window.mainloop()

    def ret_unknown(self):  # возвращает список сил с значением модуля "unknown"
        ret = []
        for element in self.forces: 
            if element.mod=='unknown': ret.append(element)
        return ret

    def ret_known(self): #позвращает список сил с известным значением модуля
        ret = []
        for element in self.forces: 
            if element.mod!='unknown': ret.append(element)
        return ret
    

    def add_a_force(self, x, y, mod, vec, name='force', type='point', len='0'): #метод добавляет приложенные к раме силы и моменты. x, y - координаты. mod - модуль. type - 'po
        # 'point' - точечно приложенная, 'lin' - равномерно распределённая. Если ввиодм не 'point' - надо указать участок len. name - тип: 'moment' или 'force'
        self.forces.append(factor(x, y, mod, vec, name, type, len))

    def add_a_holder(self, x, y, type, vector=[0, 1], force_x='unknown', force_y='unknown', moment_cw='unknown'):#типы - hinge(2 реакции), joint(1 реакция), seal(3 реакции)
        if type=='joint': 
            if vector==[0,1]: force_x, moment_cw = 0, 0
            if vector==[1,0]: force_y, moment_cw = 0, 0
        if type=='hinge': moment_cw=0
        self.add_a_force(x, y, force_x, [1,0], 'force') if force_x != 0 else 0
        self.add_a_force(x, y, force_y, [0,1], 'force') if force_y != 0 else 0
        self.add_a_force(x, y, moment_cw, 'clockwise', 'moment') if moment_cw != 0 else 0

        self.holders.append(holder(x, y, type, vector))
        
    def moment_count(self, force, x, y, is_unknown=False):# находит значение момента от передаваемого фактора в передаваемой точке
        force.mod = 1 if is_unknown else force.mod 
        if force.name=='moment': return force.mod if force.vec=='clockwise' else -1*force.mod
        elif force.type=='point': 
            #print('сила расположена в [{}, {}], в точке [{}, {}] создает момент {}'.format(force.xpos, force.ypos, x, y, force.mod*(x-force.xpos)*force.vec[1]-force.mod*(y-force.ypos)*force.vec[0]))
            return force.mod*(x-force.xpos)*force.vec[1]-force.mod*(y-force.ypos)*force.vec[0]
        else: 
            mx = force.mod*force.vec[1]*abs(force.len-(force.xpos+force.len-x if x<=force.xpos+force.len and x>=force.xpos else 0))*(x-force.xpos+(force.len-(force.xpos+force.len-x if x<=force.xpos+force.len and x>=force.xpos else 0))/2)
            my = force.mod*force.vec[0]*abs(force.len+(force.ypos-force.len-y if y>=force.ypos and y<=force.ypos+force.len else 0))*(force.ypos-y-(force.len-(force.ypos-force.len-y if y>=force.ypos and y<=force.ypos+force.len else 0))/2)
            return mx+my

    def make_a_raster(self, numopoints=10): # разбивает рамку на набор точек - нужно для растрового рисования эпюр
        for i in range(len(self.coord)-1): 
            print(self.coord[i], self.coord[i+1])
            xbeat=self.coord[i+1][0]-self.coord[i][0]
            ybeat=self.coord[i+1][1]-self.coord[i][1]
            for k in range(numopoints+1): self.rastered_frame.append([self.coord[i][0]+xbeat*(k/numopoints), self.coord[i][1]+ybeat*(k/numopoints)])
        #for k in range(numopoints): self.rastered_frame.append([self.coord[-1][0]+xbeat*(k/numopoints), self.coord[-1][1]+ybeat*(k/numopoints)])
    

    def integrate_moments_mult(self, fac1, fac2, EI=1, backwards = False): #возвращает интеграл произведения двух моментов - заменяет правила Симпсона и Верещагина для статически неопределимых задач
        rast_frame = self.rastered_frame if not backwards else self.rastered_frame[::-1]
        sum=0
        prevdot = rast_frame[0]
        counts = []
        for dot in rast_frame: 
            for factor in [fac1, fac2]:
                if len(counts) == 2: continue
                elif factor.xpos>=prevdot[0] and factor.xpos<=dot[0] and factor.ypos>=prevdot[1] and factor.ypos<=dot[1]: counts.append(factor)
            sum += self.moment_count(fac1, dot[0], dot[1])*self.moment_count(fac2, dot[0], dot[1])*sqrt((dot[0]-prevdot[0])**2+(dot[1]-prevdot[1])**2)
            prevdot = dot
        print('this is Integrator. your forces are {} and {}. Their integral is {}'.format(fac1.mod, fac2.mod, sum))
        return sum/EI

    def determine_a_force(self, force_to_det): #находит силу в случае статически неопределимой задачи. Возможно стоит переписать под систему уравнений.
        force_to_det.mod = 1
        moments_multed=[]
        for factor in self.forces: 
            if factor.mod=='unknown' or factor == force_to_det: continue 
            else: 
                moments_multed.append(self.integrate_moments_mult(force_to_det, factor))
        force_to_det.mod = sum(moments_multed)/self.integrate_moments_mult(force_to_det, force_to_det) if self.integrate_moments_mult(force_to_det, force_to_det) != 0 else 0
            
        return force_to_det.mod

    def solve(self): #эта система занимается тем, что решает систему и записывает значения факторов в них, поэтому ей не надо ничего писать.
        if not self.iscomplex():
            #в соучае если система статически определима, можно посчитать через систему из 3 уравнений, где надо стоставить уравнение моментов в начале рамы, 
            dot = self.coord[0] #это записываю точку, в которой всё считаю                # уравнение вертикальных сил и уравнение горизонтальных сил
            known_forces = self.ret_known()
            unknown_forces = self.ret_unknown() #в переменной unknown forces храню все силы с неизвестными модулями. Они на следующем этапе буут перезаписаны как 
 # как единичные, а к концу функции в них будут записаны найденные значения
            mat = np.array([[self.moment_count(element, dot[0], dot[1], is_unknown=True) for element in unknown_forces]])   #это я создал матрицу, которую будем решать. Сейчас там по одному коэффициенту неизвестного(одна строчка)
            vec = np.array([-1*sum([self.moment_count(element, dot[0], dot[1]) for element in known_forces])]) #это вектор. Там сейчас сумма моментов от известных сил на минус 1

            mat = np.append(mat, [[element.mod*element.vec[1] for element in unknown_forces]], axis=0) #ДОБАВЛЯЮ вторую строчку к матрице - вертикальные силы. находятся как неизвестная сила(равна одному) на y-составляющую вектора направления
            vec = np.append(vec, -1*sum(element.mod*element.vec[1] for element in known_forces)) # добавляю сумму вертикальных сил на минус один
            
            mat = np.append(mat, [[element.mod*element.vec[0] for element in unknown_forces]], axis=0) #то же самое для горизонтальных
            vec = np.append(vec, -1*sum(element.mod*element.vec[0] for element in known_forces))

            print(mat)
            print(vec)
            answer = np.linalg.solve(mat, vec)
            for i in range(len(unknown_forces)): unknown_forces[i].mod = answer[i]
        else: 
            print('это статически неопределимая задача. Я пока не умею такие решать! Iscomplex вернула {}'.format(self.iscomplex()))


    def iscomplex(self): # выясняет, является ли данная система статически неопределимой. Если неопределима - вернёт True
        unknown_num = 0
        for element in self.holders: 
            if element.type == 'joint': unknown_num+=1
            elif element.type == 'hinge': unknown_num+=2
            else: unknown_num+=3
            print('Это  iscomplex. Я встретил стойку типа {} и записал, что связей уже {}'.format(element.type, unknown_num))
        return unknown_num > 3


class factor: 
    def __init__(self, x, y, mod, vec, name='force', type='point', len='0'): #vec - любой вектор из двух переменных. пока что я честно не буду их отсеивать, 
        #но для правильной работы по модулю они должнны быть от 0 до 1
        #если сила распределённая, указывается точка начала, направление действия и точка окончания - пока что только равномерная сила (V 0.0.2). проверок нет.
        self.xpos = x
        self.ypos = y 
        self.mod = mod
        self.type = type
        self.name = name
        self.vec = vec # значения up, down, left, right. Если указываем момент, принимает clockwise или counterclockwise
        if type == 'lin': self.len = len
    


class holder: 
    def __init__(self, x, y, type, vec=[0,1]):
        self.xpos=x
        self.ypos=y
        self.type=type
        self.vec=vec