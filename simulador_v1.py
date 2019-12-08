#	Developed by José Juan Galera Rodríguez. josgalrod on Github.

#	"Simulador de sistemas mecánicos" is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with "Simulador de sistemas mecánicos".  If not, see <https://www.gnu.org/licenses/>

from tkinter import *
import numpy as np
import math as math
import matplotlib.pyplot as plt
import pygame as pygame
from tkinter import messagebox

raiz = Tk()
raiz.title("Simulador")

miFrame = Frame(raiz)
miFrame.pack(fill="both", expand=True)
miFrame.config(width=800, height=600)

m_eq = DoubleVar()
k_eq = DoubleVar()
c_eq = DoubleVar()
x0 = DoubleVar()
v0 = DoubleVar()
tipoFuerza = ""
fuerza = False
cteF = DoubleVar()
wF = DoubleVar()
pendF = DoubleVar()
aF = DoubleVar()
bF = DoubleVar()

# ---------FUNCIONES---------#


def asignaValores(masa, k, c, x, v):

    global m_eq
    global k_eq
    global c_eq
    global x0
    global v0
    global guardado

    m_eq = float(masa)
    k_eq = float(k)
    c_eq = float(c)
    x0 = float(x)
    v0 = float(v)


def borraDatos():

	global fuerza
	global tipoFuerza
	global cteF
	global pendF
	global wF

	fuerza = False
	tipoFuerza = ""
	cteF = 0
	pendF = 0
	wF = 0
	print("Fuerzas borradas")


def eligeFuerza():
	fuerzas = Toplevel(raiz)
	fuerzas.title("Tipo de fuerza")
	frameFuerzas = Frame(fuerzas)
	frameFuerzas.pack()

	selecFuerza = Label(frameFuerzas, text="Seleccione la forma de la fuerza aplicada:")
	selecFuerza.grid(row=0, column=0, sticky="w")
	fuerzaCte = Button(frameFuerzas, text="Escalón", width=6, command=lambda:ponFuerza("escalon"))
	fuerzaCte.grid(row=1, column=0, sticky="w")
	fuerzaRampa = Button(frameFuerzas, text="Rampa", width=6, command=lambda:ponFuerza("rampa"))
	fuerzaRampa.grid(row=2, column=0, sticky="w")
	fuerzaCoseno = Button(frameFuerzas, text="Senoidal", width=6, command=lambda:ponFuerza("senoidal"))
	fuerzaCoseno.grid(row=3, column=0, sticky="w")


def ponFuerza(tipo):

	global tipoFuerza
	global fuerza
	global cteF
	global wF
	global pendF

	fuerzasEscribe = Toplevel(raiz)
	fuerzasEscribe.title("Tipo de fuerza")
	frameFuerzasE = Frame(fuerzasEscribe)
	frameFuerzasE.pack()

	if tipo == "escalon":

		constante = StringVar()
		tipoLabel = Label(frameFuerzasE, text="F(t)=C")
		tipoLabel.grid(row=0, column=0, sticky="w")
		cteLabel = Label(frameFuerzasE, text="Valor C:")
		cteLabel.grid(row=1, column=0, sticky="w")
		fuerzaCte = Entry(frameFuerzasE, textvariable=cteF)
		fuerzaCte.grid(row=1, column=1, sticky="w")
		cteAsignar = Button(frameFuerzasE, text="OK", command=lambda:guardaDato(fuerzaCte.get(), "constante"))
		cteAsignar.grid(row=1, column=2)
		tipoFuerza = "escalon"

	if tipo == "rampa":

		pendiente = StringVar()
		tipoLabel = Label(frameFuerzasE, text="F(t)=B*t")
		tipoLabel.grid(row=0, column=0, sticky="w")
		pendLabel = Label(frameFuerzasE, text="Valor B:")
		pendLabel.grid(row=1, column=0, sticky="w")
		fuerzaPend = Entry(frameFuerzasE, textvariable=pendF)
		fuerzaPend.grid(row=1, column=1, sticky="w")
		pendAsignar = Button(frameFuerzasE, text="OK", command=lambda:guardaDato(fuerzaPend.get(), "pendiente"))
		pendAsignar.grid(row=1, column=2)
		tipoFuerza = "rampa"

	if tipo == "senoidal":

		tipoLabel = Label(frameFuerzasE, text="F(t)=A*cos(w*t)+B*sen(w*t)")
		tipoLabel.grid(row=0, column=0, sticky="w")
		cteALabel = Label(frameFuerzasE, text="Valor A:")
		cteALabel.grid(row=1, column=0, sticky="w")
		fuerzaCteA = Entry(frameFuerzasE, textvariable=aF)
		fuerzaCteA.grid(row=1, column=1, sticky="w")
		cteAAsignar = Button(frameFuerzasE, text="OK", command=lambda: guardaDato(fuerzaCteA.get(), "ctecos"))
		cteAAsignar.grid(row=1, column=2)
		cteBLabel = Label(frameFuerzasE, text="Valor B:")
		cteBLabel.grid(row=2, column=0, sticky="w")
		fuerzaCteB = Entry(frameFuerzasE, textvariable=bF)
		fuerzaCteB.grid(row=2, column=1, sticky="w")
		cteBAsignar = Button(frameFuerzasE, text="OK", command=lambda: guardaDato(fuerzaCteB.get(), "ctesen"))
		cteBAsignar.grid(row=2, column=2)
		omLabel = Label(frameFuerzasE, text="Valor w:")
		omLabel.grid(row=3, column=0, sticky="w")
		fuerzaOm = Entry(frameFuerzasE, textvariable=wF)
		fuerzaOm.grid(row=3, column=1, sticky="w")
		omAsignar = Button(frameFuerzasE, text="OK", command=lambda: guardaDato(fuerzaOm.get(), "omega"))
		omAsignar.grid(row=3, column=2)
		tipoFuerza = "senoidal"


def guardaDato(dato, tipo):

	global cteF
	global wF
	global pendF
	global fuerza
	global aF
	global bF

	if tipo == "constante":
		cteF = float(dato)

	elif tipo == "omega":
		wF = float(dato)

	elif tipo == "pendiente":
		pendF = float(dato)

	elif tipo == "ctecos":
		aF = float(dato)

	elif tipo == "ctesen":
		bF = float(dato)

	fuerza = True


def muestraDatos():

	global m_eq
	global k_eq
	global c_eq
	global x0
	global v0
	global wF

	wn = DoubleVar()
	xi = DoubleVar()
	wd = DoubleVar()
	tau = DoubleVar()
	fdc = DoubleVar()

	wn = math.sqrt(k_eq/m_eq)
	xi = c_eq/(2*m_eq*wn)


	if xi == 0:
		wd = wn

	elif xi < 1:
		wd = wn*math.sqrt(1-pow(xi,2))

	elif xi >= 1:
		wd = 0

	datos = Toplevel(raiz)
	datos.title("Datos")
	frameDatos = Frame(datos)
	frameDatos.pack()

	if fuerza:
		if xi == 0:
			wd = wn
		if xi <= 1:
			wd = wn*math.sqrt(1-pow(xi,2))
		if tipoFuerza == "senoidal":
			tau = wF/wn
			try:
				fdc = math.pow(math.sqrt((1-tau**2)**2+(2*xi*tau**2)),-1)
			except ValueError:
				fdc = StringVar()
				fdc = "Infinito"

	valorWn = Label(frameDatos, text="wn:")
	valorWn.grid(row=0, column=0)
	valorWd = Label(frameDatos, text="wd:")
	valorWd.grid(row=1, column=0)
	valorXi = Label(frameDatos, text="xi:")
	valorXi.grid(row=2, column=0)
	if tipoFuerza == "senoidal":
		valorTau = Label(frameDatos, text="tau:")
		valorTau.grid(row=3, column=0)
		valorFDC = Label(frameDatos, text="FDC:")
		valorFDC.grid(row=4, column=0)

	entryWn = Entry(frameDatos, textvariable=wn)
	entryWn.grid(row=0, column=1)
	entryWn.insert(0, round(wn, 5))
	entryWd = Entry(frameDatos, textvariable=wd)
	entryWd.grid(row=1, column=1)
	entryWd.insert(0, round(wd, 5))
	entryXi = Entry(frameDatos, textvariable=xi)
	entryXi.grid(row=2, column=1)
	entryXi.insert(0, round(xi, 5))
	if tipoFuerza == "senoidal":
		entryTau = Entry(frameDatos, textvariable=tau)
		entryTau.grid(row=3, column=1)
		entryTau.insert(0, round(tau, 5))
		entryFDC = Entry(frameDatos, textvariable=fdc)
		entryFDC.grid(row=4, column=1)
		entryFDC.insert(0, round(fdc, 5))


def eligePlot():
	plots = Toplevel(raiz)
	plots.title("Plots")
	framePlots = Frame(plots)
	framePlots.pack()

	selecPlot = Label(framePlots, text="Seleccione el plot que desea ver:")
	selecPlot.grid(row=0, column=0, sticky="w")
	plotPos = Button(framePlots, text="x-t", width=6, command=lambda: plotPosicion())
	plotPos.grid(row=1, column=0, sticky="w")
	plotTauFDC = Button(framePlots, text="tau-fdc", width=6, command=lambda: plotFDC())
	plotTauFDC.grid(row=2, column=0, sticky="w")


def plotFDC():
	global m_eq
	global k_eq
	global c_eq
	global x0
	global v0
	global fuerza
	global tipoFuerza
	global cteF
	global pendF
	global wF
	global inicF

	wn = math.sqrt(k_eq/m_eq)
	xi = c_eq/(2*m_eq*wn)

	try:
		def fdc(tau):
			return math.pow(math.sqrt((1-tau**2)**2+(2*xi*tau**2)), -1)

		tau_output = np.arange(0, 5, 0.01)
		plt.xlabel('tau') 
		plt.plot(tau_output, [fdc(i) for i in tau_output], label="FDC")
		plt.legend()
		plt.show()

	except ValueError:
		messagebox.showwarning("Atención", "Para xi=0, el FDC tiende a infinito con tau=1")


def plotPosicion():
	global m_eq
	global k_eq
	global c_eq
	global x0
	global v0
	global fuerza
	global tipoFuerza
	global cteF
	global pendF
	global wF

	wn = math.sqrt(k_eq/m_eq)
	xi = c_eq/(2*m_eq*wn)

	if fuerza is False:

		if xi == 0:

			def x(t):
				return x0*math.cos(wn*t)+v0/wn*math.sin(wn*t)

			def v(t):
				return -x0*wn*math.sin(wn*t)+v0*math.cos(wn*t)

		elif xi < 1:
			wd = wn*math.sqrt(1-pow(xi, 2))

			def x(t):
				return pow(math.e, -xi*wn*t)*((x0*math.cos(wd*t)+(v0+xi*wn*x0)/wd*math.sin(wd*t)))
			
			def v(t):
				return (-xi*wn)*pow(math.e, -xi*wn*t)*((x0*math.cos(wd*t)+(v0+xi*wn*x0)/wd*math.sin(wd*t)))+pow(math.e, -xi*wn*t)*((-wd*x0*math.sin(wd*t)+(v0+xi*wn*x0)*math.cos(wd*t)))

		elif xi == 1:

			def x(t):
				return pow(math.e, -wn*t)*(x0+(v0+wn*x0)*t)

			def v(t):
				return (-wn)*pow(math.e, -wn*t)*(x0+(v0+wn*x0)*t)+pow(math.e, -wn*t)*(v0+wn*x0)

		elif xi > 1:
			r1 = wn*(-xi+math.sqrt(math.pow(xi, 2)-1))
			r2 = wn*(-xi-math.sqrt(math.pow(xi, 2)-1))

			def x(t):
				return (x0/2*(1+xi/(math.sqrt(math.pow(xi, 2)-1)))+v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r1*t)+(x0/2*(1-xi/(math.sqrt(math.pow(xi, 2)-1)))-v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r2*t)
			
			def v(t):
				return r1*(x0/2*(1+xi/(math.sqrt(math.pow(xi, 2)-1)))+v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r1*t)+r2*(x0/2*(1-xi/(math.sqrt(math.pow(xi, 2)-1)))-v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r2*t)

	else:

		if tipoFuerza == "escalon":

			if xi == 0:

				def x(t):
					return cteF/k_eq*(1-math.cos(wn*t))

				def v(t):
					return cteF/k_eq*wn*math.sin(wn*t)

			else:

				wd = wn*math.sqrt(1-pow(xi, 2))

				def x(t):
					return math.pow(math.e, -xi*wn*t)*((x0-cteF/k_eq)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-xi*wn*(cteF/k_eq)/wd)*math.sin(wd*t))+cteF/k_eq
				
				def v(t):
					return (-xi*wn)*math.pow(math.e, -xi*wn*t)*((x0-cteF/k_eq)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-xi*wn*(cteF/k_eq)/wd)*math.sin(wd*t))+math.pow(math.e, -xi*wn*t)*(-wd*(x0-cteF/k_eq)*math.sin(wd*t)+((v0+xi*wn*x0)-xi*wn*(cteF/k_eq))*math.cos(wd*t))

		elif tipoFuerza == "rampa":

			if xi == 0:

				def x(t):
					return pendF/(k_eq*wn)*(wn*t-math.sin(wn*t))
				
				def v(t):
					return pendF/k_eq*(1-math.cos(wn*t))

			else:

				wd = wn*math.sqrt(1-pow(xi, 2))
				R = pendF/k_eq
				S = -c_eq*float(pendF)/(k_eq**2)

				def x(t):
					return math.pow(math.e, -xi*wn*t)*((x0-S)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-(R+xi*wn*S)/wd)*math.sin(wd*t))+R*t+S
				
				def v(t):
					return -xi*wn*math.pow(math.e, -xi*wn*t)*((x0-S)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-(R+xi*wn*S)/wd)*math.sin(wd*t))+math.pow(math.e, -xi*wn*t)*(-wd*(x0-S)*math.sin(wd*t)+((v0+xi*wn*x0)-(R+xi*wn*S))*math.cos(wd*t))+R

		elif tipoFuerza == "senoidal":

			tau = wF/wn
			xc_est = aF/k_eq
			xs_est = bF/k_eq
			fi = math.atan(xc_est/xs_est)

			if xi == 0:

				if tau <= 0.9 or tau >= 1.1:

					def x(t):
						return (x0-xc_est/(1-tau**2))*math.cos(wn*t)+(v0/wn-tau*xs_est/(1-tau**2))*math.sin(wn*t)+(xc_est/(1-tau**2))*math.cos(wF*t)+(xs_est/(1-tau**2))*math.sin(wF*t)
					
					def v(t):
						return -wn*(x0-xc_est/(1-tau**2))*math.sin(wn*t)+wn*(v0/wn-tau*xs_est/(1-tau**2))*math.cos(wn*t)-wF*(xc_est/(1-tau**2))*math.sin(wF*t)+wF*(xs_est/(1-tau**2))*math.cos(wF*t)

				elif tau == 1:
					
					def x(t):
						return (x0-bF*t/(2*m_eq*wn))*math.cos(wn*t)+(v0/wn+bF/(2*k_eq)+(aF*t)/(2*m_eq*wn))*math.sin(wn*t)
					
					def v(t):
						return (bF/(2*m_eq*wn))*math.cos(wn*t)-(x0-bF*t/(2*m_eq*wn))*wn*math.sin(wn*t)+wn*(v0/wn+bF/(2*k_eq)+(aF*t)/(2*m_eq*wn))*math.cos(wn*t)+aF/(2*m_eq*wn)*math.sin(wn*t)

				elif tau > 0.9 and tau < 1:

					def x(t):
						return 2*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.sin((wF-wn)*t/2)
					
					def v(t):
						return -(wF*wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.sin(((wF+wn)*t+math.pi)/2)*math.sin((wF-wn)*t/2)+(wF-wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.cos((wF-wn)*t/2)

				elif tau > 1 and tau < 1.1:
					
					def x(t):
						return 2*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.sin((wF-wn)*t/2)					
					def v(t):
						return -(wF*wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.sin(((wF+wn)*t+math.pi)/2)*math.sin((wF-wn)*t/2)+(wF-wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.cos((wF-wn)*t/2)

			elif xi < 1:

				wd = wn*math.sqrt(1-pow(xi, 2))
				fi = math.atan(-2*xi*tau/(1-tau**2))
				xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
				vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
				
				def x(t):
					return math.pow(math.e, -xi*wn*t)*(((x0-xp0)*math.cos(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)/wd*math.sin(wd*t)))+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)
				
				def v(t):
					return -xi*wn*math.pow(math.e, -xi*wn*t)*((-wd*(x0-xp0)*math.sin(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)*math.cos(wd*t)))+math.pow(math.e, -xi*wn*t)*((-wd*(x0-xp0)*math.sin(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)*math.cos(wd*t)))+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)

			elif xi == 1:

				fi = math.atan(-2*xi*tau/(1-tau**2))
				xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
				vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
				A = x0-xp0
				B = (x0-xp0)*xi*wn+v0-vp0

				def x(t):
					return math.pow(math.e,-xi*wn*t)*(A+B*t)+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)

				def v(t):
					return -xi*wn*math.pow(math.e,-xi*wn*t)*(A+B*t)+math.pow(math.e,-xi*wn*t)*B+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)

			elif xi > 1:

				fi = math.atan(-2*xi*tau/(1-tau**2))
				xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
				vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
				A = 1/2*((x0-xp0)-(v0-vp0)/(wn*xi*math.sqrt(xi**2-1)))
				B = 1/2*((x0-xp0)+(v0-vp0)/(wn*xi*math.sqrt(xi**2-1)))

				def x(t):
					return math.pow(math.e, -xi*wn*t)*(A*math.pow(math.e,-wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)

				def v(t):
					return -xi*wn*math.pow(math.e, -xi*wn*t)*(A*math.pow(math.e,-wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+math.pow(math.e, -xi*wn*t)*wn*math.sqrt(xi**2-1)*(-A*math.pow(math.e,-wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)

	t_output = np.arange(0, 50, 0.01)
	plt.xlabel('tiempo') 
	plt.plot(t_output, [x(i) for i in t_output], label="Posicion")
	plt.plot(t_output, [v(i) for i in t_output], label="Velocidad")
	plt.legend()
	plt.show()


def ejecutaSimulador():

	global m_eq
	global k_eq
	global c_eq
	global x0
	global v0
	global fuerza
	global tipoFuerza
	global cteF
	global pendF
	global wF

	dimension_ventana = (1200, 800)
	coord_eq_masa = 400
	coord_eq_resorte = 0
	coord_eq_amortiguador = 0
	coord_barra = (500, 0)
	coord_suelo = (0, 350)
	ancho_resorte = 400
	ancho_amortiguador = 400
	background = (255, 255, 255)
	fps = 60

	dt = 1/(2*float(fps))
	wn = math.sqrt(k_eq/m_eq)
	xi = c_eq/(2*m_eq*wn)
	t = 0

	ventana = pygame.display.set_mode(dimension_ventana)
	pygame.display.set_caption("Simulación")
	ventana.fill(background)
	reloj = pygame.time.Clock()

	suelo = pygame.image.load('suelo.png')
	barra = pygame.image.load('barra.png').convert_alpha()
	masa_0 = pygame.image.load('masa.png')
	resorte_0 = pygame.image.load('resorte.png').convert_alpha()
	amortiguador_0 = pygame.image.load('amortiguador.png').convert_alpha()

	x = float(x0)
	v = float(v0)
	while True:
		for evento in pygame.event.get():
			if evento.type == pygame.QUIT:
				pygame.quit()
				break

		masa = pygame.transform.scale2x(masa_0)
		resorte = pygame.transform.scale(resorte_0, (int(abs(ancho_resorte+x)), 90))
		amortiguador = pygame.transform.scale(amortiguador_0, (int(abs(ancho_amortiguador+x)), 90))

		ventana.blit(resorte, (coord_eq_resorte, 250))
		ventana.blit(masa, (coord_eq_masa+x, 150))
		ventana.blit(barra, coord_barra)
		ventana.blit(suelo, coord_suelo)
		ventana.blit(amortiguador, (coord_eq_amortiguador, 150))

		pygame.display.flip()
		ventana.fill(background)
		ventana.fill(background)

		if fuerza is False:
			if xi == 0:
				t += dt
				x = x0*math.cos(wn*t)+v0/wn*math.sin(wn*t)
				v = -x0*wn*math.sin(wn*t)+v0*math.cos(wn*t)

			elif xi < 1:
				t += dt
				wd = wn*math.sqrt(1-pow(xi, 2))
				x = pow(math.e, -xi*wn*t)*((x0*math.cos(wd*t)+(v0+xi*wn*x0)/wd*math.sin(wd*t)))
				v = (-xi*wn)*pow(math.e, -xi*wn*t)*((x0*math.cos(wd*t)+(v0+xi*wn*x0)/wd*math.sin(wd*t)))+pow(math.e, -xi*wn*t)*((-wd*x0*math.sin(wd*t)+(v0+xi*wn*x0)*math.cos(wd*t)))

			elif xi == 1:
				t += dt
				x = pow(math.e, -wn*t)*(x0+(v0+wn*x0)*t)
				v = (-wn)*pow(math.e, -wn*t)*(x0+(v0+wn*x0)*t)+pow(math.e, -wn*t)*(v0+wn*x0)

			elif xi > 1:
				t += dt
				r1 = wn*(-xi+math.sqrt(math.pow(xi, 2)-1))
				r2 = wn*(-xi-math.sqrt(math.pow(xi, 2)-1))

				x = (x0/2*(1+xi/(math.sqrt(math.pow(xi, 2)-1)))+v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r1*t)+(x0/2*(1-xi/(math.sqrt(math.pow(xi, 2)-1)))-v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r2*t)
				v = r1*(x0/2*(1+xi/(math.sqrt(math.pow(xi, 2)-1)))+v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r1*t)+r2*(x0/2*(1-xi/(math.sqrt(math.pow(xi, 2)-1)))-v0/(2*wn*math.sqrt(math.pow(xi, 2)-1)))*math.pow(math.e, r2*t)

		else:
			if tipoFuerza == "escalon":
				if xi == 0:
					t += dt
					x = cteF/k_eq*(1-math.cos(wn*t))
					v = cteF/k_eq*wn*math.sin(wn*t)
				else:
					t += dt
					wd = wn*math.sqrt(1-pow(xi, 2))
					x = math.pow(math.e, -xi*wn*t)*((x0-cteF/k_eq)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-xi*wn*(cteF/k_eq)/wd)*math.sin(wd*t))+cteF/k_eq
					v = (-xi*wn)*math.pow(math.e, -xi*wn*t)*((x0-cteF/k_eq)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-xi*wn*(cteF/k_eq)/wd)*math.sin(wd*t))+math.pow(math.e, -xi*wn*t)*(-wd*(x0-cteF/k_eq)*math.sin(wd*t)+((v0+xi*wn*x0)-xi*wn*(cteF/k_eq))*math.cos(wd*t))

			elif tipoFuerza == "rampa":
				if xi == 0:
					t += dt
					x = pendF/(k_eq*wn)*(wn*t-math.sin(wn*t))
					v = pendF/k_eq*(1-math.cos(wn*t))
				else:
					t += dt
					wd = wn*math.sqrt(1-pow(xi, 2))
					R = pendF/k_eq
					S = -c_eq*pendF/(k_eq**2)
					x = math.pow(math.e, -xi*wn*t)*((x0-S)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-(R+xi*wn*S)/wd)*math.sin(wd*t))+R*t+S
					v = -xi*wn*math.pow(math.e, -xi*wn*t)*((x0-S)*math.cos(wd*t)+((v0+xi*wn*x0)/wd-(R+xi*wn*S)/wd)*math.sin(wd*t))+math.pow(math.e, -xi*wn*t)*(-wd*(x0-S)*math.sin(wd*t)+((v0+xi*wn*x0)-(R+xi*wn*S))*math.cos(wd*t))+R

			elif tipoFuerza == "senoidal":

				tau = wF/wn
				xc_est = aF/k_eq
				xs_est = bF/k_eq
				fi = math.atan(xc_est/xs_est)

				if xi == 0:
					if tau <= 0.9 or tau >= 1.1:
						t += dt
						x = (x0-xc_est/(1-tau**2))*math.cos(wn*t)+(v0/wn-tau*xs_est/(1-tau**2))*math.sin(wn*t)+(xc_est/(1-tau**2))*math.cos(wF*t)+(xs_est/(1-tau**2))*math.sin(wF*t)
						v = -wn*(x0-xc_est/(1-tau**2))*math.sin(wn*t)+wn*(v0/wn-tau*xs_est/(1-tau**2))*math.cos(wn*t)-wF*(xc_est/(1-tau**2))*math.sin(wF*t)+wF*(xs_est/(1-tau**2))*math.cos(wF*t)

					elif tau == 1:
						t += dt
						x = (x0-bF*t/(2*m_eq*wn))*math.cos(wn*t)+(v0/wn+bF/(2*k_eq)+(aF*t)/(2*m_eq*wn))*math.sin(wn*t)
						v = (bF/(2*m_eq*wn))*math.cos(wn*t)-(x0-bF*t/(2*m_eq*wn))*wn*math.sin(wn*t)+wn*(v0/wn+bF/(2*k_eq)+(aF*t)/(2*m_eq*wn))*math.cos(wn*t)+aF/(2*m_eq*wn)*math.sin(wn*t)

					elif  tau > 0.9 and tau < 1:
						t += dt
						x = 2*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.sin((wF-wn)*t/2)					
						v = -(wF*wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.sin(((wF+wn)*t+math.pi)/2)*math.sin((wF-wn)*t/2)+(wF-wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.cos((wF-wn)*t/2)

					elif  tau < 1.1 and tau > 1:
						t += dt
						x = 2*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.sin((wF-wn)*t/2)					
						v = -(wF*wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.sin(((wF+wn)*t+math.pi)/2)*math.sin((wF-wn)*t/2)+(wF-wn)*math.sqrt(xc_est**2+xs_est**2)/(1-tau**2)*math.cos(((wF+wn)*t+2*fi)/2)*math.cos((wF-wn)*t/2)

				elif xi < 1:

					wd = wn*math.sqrt(1-pow(xi, 2))
					fi = math.atan(-2*xi*tau/(1-tau**2))
					xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
					vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
					
					t += dt
					x = math.pow(math.e, -xi*wn*t)*(((x0-xp0)*math.cos(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)/wd*math.sin(wd*t)))+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)
					v = -xi*wn*math.pow(math.e, -xi*wn*t)*((-wd*(x0-xp0)*math.sin(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)*math.cos(wd*t)))+math.pow(math.e, -xi*wn*t)*((-wd*(x0-xp0)*math.sin(wd*t)+(v0+xi*wn*x0-vp0-xi*wn*xp0)*math.cos(wd*t)))+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)

				elif xi == 1:

					fi = math.atan(-2*xi*tau/(1-tau**2))
					xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
					vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
					A = x0-xp0
					B = (x0-xp0)*xi*wn+v0-vp0

					t += dt
					x = math.pow(math.e, -xi*wn*t)*(A+B*t)+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)
					v = -xi*wn*math.pow(math.e, -xi*wn*t)*(A+B*t)+math.pow(math.e, -xi*wn*t)*B+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)


				elif xi > 1:

					fi = math.atan(-2*xi*tau/(1-tau**2))
					xp0 = math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(fi)
					vp0 = math.sqrt(xc_est**2+xs_est**2)*wF/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)
					A = 1/2*((x0-xp0)-(v0-vp0)/(wn*xi*math.sqrt(xi**2-1)))
					B = 1/2*((x0-xp0)+(v0-vp0)/(wn*xi*math.sqrt(xi**2-1)))

					t += dt
					x = math.pow(math.e, -xi*wn*t)*(A*math.pow(math.e, -wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.sin(wF*t+fi)
					v = -xi*wn*math.pow(math.e, -xi*wn*t)*(A*math.pow(math.e, -wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+math.pow(math.e, -xi*wn*t)*wn*math.sqrt(xi**2-1)*(-A*math.pow(math.e,-wn*math.sqrt(xi**2-1))+B*math.pow(math.e,+wn*math.sqrt(xi**2-1)))+wF*math.sqrt(xc_est**2+xs_est**2)/math.sqrt((1-tau**2)**2+(2*xi*tau)**2)*math.cos(wF*t+fi)


# --------DISEÑO GRÁFICO----------#

barraMenu=Menu(raiz)
raiz.config(menu=barraMenu, width=300, height=300)

archivoMenu=Menu(barraMenu, tearoff=0)
archivoMenu.add_command(label="Salir", command=lambda:avisoSalir())

archivoAyuda=Menu(barraMenu, tearoff=0)
archivoAyuda.add_command(label="Licencia", command=lambda:avisoLicencia())
archivoAyuda.add_command(label="Acerca de", command=lambda:infoAdicional())

barraMenu.add_cascade(label="Archivo", menu=archivoMenu)
barraMenu.add_cascade(label="Ayuda", menu=archivoAyuda)

def infoAdicional():
	messagebox.showinfo("Simulador", "Desarrollado por José Juan Galera Rodríguez")

def avisoLicencia():
	messagebox.showinfo("Licencia", "Producto bajo licencia GNU GPL v3, lea el archivo license.txt")

def avisoSalir():
	valor = messagebox.askokcancel("Salir", "¿Desea salir de la aplicación?")
	if valor == True:
		raiz.destroy()


cuadroMasa = Entry(miFrame, textvariable=m_eq)
cuadroMasa.grid(row=0, column=1, sticky="w", padx=7, pady=5)
cuadroMasa.config(justify="center")

cuadroK = Entry(miFrame, textvariable=k_eq)
cuadroK.grid(row=1, column=1, sticky="w", padx=7, pady=5)
cuadroK.config(justify="center")

cuadroC = Entry(miFrame, textvariable=c_eq)
cuadroC.grid(row=2, column=1, sticky="w", padx=7, pady=5)
cuadroC.config(justify="center")

cuadroFuerza = Button(miFrame, text="Elegir fuerza", command=lambda: eligeFuerza())
cuadroFuerza.grid(row=3, column=1)

cuadroX = Entry(miFrame, textvariable=x0)
cuadroX.grid(row=4, column=1, sticky="w", padx=7, pady=5)
cuadroX.config(justify="center")

cuadroV = Entry(miFrame, textvariable=v0)
cuadroV.grid(row=5, column=1, sticky="w", padx=7, pady=5)
cuadroV.config(justify="center")

# Labels
masaLabel = Label(miFrame, text="Masa eq:")
masaLabel.grid(row=0, column=0, sticky="e", padx=10, pady=5)

kLabel = Label(miFrame, text="K eq:")
kLabel.grid(row=1, column=0, padx=7, pady=5)

cLabel = Label(miFrame, text="C eq:")
cLabel.grid(row=2, column=0, padx=7, pady=5)

fuerzaLabel = Label(miFrame, text="F(t):")
fuerzaLabel.grid(row=3, column=0, padx=7, pady=5)

xLabel = Label(miFrame, text="x0:")
xLabel.grid(row=4, column=0, padx=7, pady=5)

vLabel = Label(miFrame, text="v0:")
vLabel.grid(row=5, column=0, padx=7, pady=5)

botonEjecutar = Button(miFrame, text="SIMULAR", width=8, command=lambda: ejecutaSimulador())
botonEjecutar.grid(row=6, column=1)
botonEjecutar.config(bg="#a9a9a9")

botonEjecutar = Button(miFrame, text="PLOT", width=8, command=lambda: eligePlot())
botonEjecutar.grid(row=6, column=2)
botonEjecutar.config(bg="#a9a9a9")

botonGuardar = Button(miFrame, text="GUARDAR", width=8, command=lambda: asignaValores(cuadroMasa.get(), cuadroK.get(), cuadroC.get(), cuadroX.get(), cuadroV.get()))
botonGuardar.grid(row=2, column=2)

botonBorrarF = Button(miFrame, text="BORRAR F", width=8, command=lambda: borraDatos())
botonBorrarF.grid(row=3, column=2)

botonDatos = Button(miFrame, text="DATOS", width=8, command=lambda: muestraDatos())
botonDatos.grid(row=6, column=0)
botonDatos.config(bg="#a9a9a9")

raiz.mainloop()
