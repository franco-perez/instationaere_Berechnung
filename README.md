#!/usr/bin/python
# coding: utf8

#1) cambio: float durch Double ersetzen

#debugger: aufruf mit: " python -m pdb Inst_Schrank_final.py", next:n, continue:c, breakpoint:b, werte ausgeben: p
### simulacion hecha con la interpolacion de franco y datenschaltschrank igual que para el caso estacionario

### exigencias del programa

#1) die .txt Datei soll folgendes heißen "Kaelteleistung_messungen_782", die Zahl 782, kommt aus Daten_Schaltshcrank wird per angabe definiert.


## corregir : 
#1) Qkuehl_average , the average is not correctly calculated as Messung

#1) reparar el Luftwerte, graficar nudo del aire y no de la placa, revisar que sensores son los adecuados
#3) probar que tanto cambian las sol. con dif sol ec. diferenciales. die sollen nicht viel abweichen.

#slices indexing?
#https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.indexing.html
#x[1:2, 1:3]


####notas notes anmerkungen
### 1) os.path.abspath(__file__) => gib mir das Verzeichnis des aktuellen Dokumentes



########## https://github.com/scipy/scipy/blob/master/scipy/integrate/_ivp/ivp.py

############################## de la bibliografia 1 odeint vector
#### odeint uses an Adams integrator for non-stiff problems, and a
##backwards differentiation method for stiff problem.
## it is our problem stiff or non stiff.. difficult to determine. 

figurename= 'franco'
Stat_Berechnung = 0 #wahl 0 aus, 1 eingeschaltet.
Ventilator=0 #wert für initialisierung gleich 0 
init_time = 0 #Variabeln die erlaubt der Kontrolle  der Ausschaltung des Dachkuehlgeraets
Dachkuehlgeraet_Qkuehl=1226 #watt
dichte_p=2700 # (kg/m3) Platte ist aus Aluminium
dichte_g=7860 # (kg/m3) Gehaeuse ist aus Stahl.



import pdb

#importieren mathematischer Operatoren
import numpy as np
import string
import math
import re 
import csv
import time
import datetime
import os
import Luftwerte_Messungen
import shutil


from time import gmtime, strftime


#from os import*
from Luftwerte_Messungen import*
from numpy import *
from scipy.optimize import *
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

######################################################3
### benötigte Dokumente
### 1. Daten_Schaltschrank.txt
### 2. Sw.csv
### 3. mcp2.csv
### 4. Kaelteleistung_messungen.csv
######################################################3



##############################################################################################
##############################################################################################
##############################################################################################

dateiname = "Daten_Schaltschrank"

#Funktion, um Werte aus Textdatei auszulesen
def bekomme_wert (string):
	laenge = len(string)
	index=zeile.find(':')
	wert=string[index+1:laenge]
	return wert

#Werte aus Textdatei suchen "Daten_Schaltschrank"
datei = open(dateiname).readlines() #Daten_Schaltschrank öffnen und lesen.
laenge_datei=len(datei)
for i in range(0, laenge_datei):
	zeile=re.sub('[\s+]', '', datei[i])
	if zeile.find('Schrankhöhe[m]:') != -1:
		h=bekomme_wert(zeile)
	elif zeile.find('Schrankbreite[m]:') != -1:
		b=bekomme_wert(zeile)
	elif zeile.find('Schranktiefe[m]:') != -1:
		t=bekomme_wert(zeile)
	elif zeile.find('Wanddicke[m]:') != -1:
		s=bekomme_wert(zeile)
	elif zeile.find('Plattenhöhe[m]:') != -1:
		hp=bekomme_wert(zeile)
	elif zeile.find('Plattenbreite[m]:') != -1:
		bp=bekomme_wert(zeile)
	elif zeile.find('Plattentiefe[m]:') != -1:
		tp=bekomme_wert(zeile)
	elif zeile.find('xp[m]:') != -1:
		xp=bekomme_wert(zeile)
	elif zeile.find('yp[m]:') != -1:
		yp=bekomme_wert(zeile)
	elif zeile.find('zp[m]:') != -1:
		zp=bekomme_wert(zeile)
	elif zeile.find('cp[J/(kgK)]:') != -1:
		cp=bekomme_wert(zeile)
	elif zeile.find('Luftdichte[kg/m3]:') != -1:
		dichte=bekomme_wert(zeile)
	elif zeile.find('Cs[W*m^(-2)*K^(-4)]:') != -1:
		Cs=bekomme_wert(zeile) 
	elif zeile.find('EmissionsgradGehäuse:') != -1:
		eg=bekomme_wert(zeile)
	elif zeile.find('EmissionsgradPlatte:') != -1:
		ep=bekomme_wert(zeile)
	elif zeile.find('WärmeleitfähigkeitGehäuse:') != -1:
		lambda_g=bekomme_wert(zeile)
	elif zeile.find('WärmeleitfähigkeitPlatte:') != -1:
		lambda_p=bekomme_wert(zeile)
	elif zeile.find('Volumenstrom[m^3/h]:') != -1:
		v_dot=bekomme_wert(zeile)
	elif zeile.find('WärmestromKühlgerät[W]:') != -1:
		QKKK=bekomme_wert(zeile)
	elif zeile.find('Zone1:') != -1:
		QV1=bekomme_wert(zeile)
	elif zeile.find('Zone2:') != -1:
		QV2=bekomme_wert(zeile)
	elif zeile.find('Zone3:') != -1:
		QV3=bekomme_wert(zeile)
	elif zeile.find('Umgebungstemperatur[K]:') != -1:
		Tamb=bekomme_wert(zeile)
	elif zeile.find('Automatic_on_off:') != -1:
		Automatic_on_off=bekomme_wert(zeile)
	elif zeile.find('Dateilaenge:') != -1:
		NumOfData=bekomme_wert(zeile)



	
# string zu double konvertieren

h=double(h)
b=double(b)
t=double(t)
s=double(s)
hp=double(hp)
tp=double(tp)
bp=double(bp)
xp=double(xp)
yp=double(yp)
zp=double(zp)
cp=double(cp)
dichte=double(dichte)
Cs=double(Cs)
eg=double(eg)
ep=double(ep)
lambda_g=double(lambda_g)
lambda_p=double(lambda_p)
v_dot=double(v_dot)
QKKK=double(QKKK) # es wird benötigt nur beim Stationären Fall
Tamb=double(Tamb) #Tamb in kelvin
QV1=double(QV1)
QV2=double(QV2)
QV3=double(QV3)
Automatic_on_off=float(Automatic_on_off)
NumOfData=int(NumOfData)


## Ordner erstellen mit Name der Simulation.
#os.makedir()

#sim_date=datetime.datetime.today().strftime('%Y-%m-%d')
sim_date=datetime.datetime.today().strftime('%d_%b_%Y_%H:%M:%S')




#work_path= os.getcwd()
#gotodir = "cd /home/fabian/Dokumente/studienarbeit/"
#os.chdir(work_path)
#new_dir=os.getcwd()
#print new_dir
#os.path.dirname(file_path)


#os.makedirs(Results_carpet)
my_path=os.getcwd() ### Aktuelles Verzeichnis der Berechnungsdokumente.

## Ergebnisse in einem neuen Ordner speichern 
## for overwrite a folder
results_directory = r''+ str(my_path)+'/results_'+ str(sim_date) 


## wenn der Folder da ist, wird der Weg durch die if schleife, wenn nicht, wird der if schleife übersprungen und geht es weiter zum Ordner erstellen
if os.path.exists(results_directory):
    shutil.rmtree(results_directory)
os.makedirs(results_directory) 


workpath=os.path.abspath(__file__) # der aktuellen Pfad bis zum aktuellen Dokument.  '/home/fabian/Dokumente/TEST_studienarbeiten/test_old/Inst_Schrank_final.py'










if Stat_Berechnung == 1:
	
	global choice
	global ENDTIME
	global NDATEN	
	choice = 3			
	ENDTIME=20
	NDATEN= 101
	mcp = np.ones((35,2))


else:
	global choice
	choice = 1


if Stat_Berechnung == 0: 
	with open('mcp2_x3.csv') as csvfile:
		reader = csv.reader(csvfile)
		reader.next()
	
		mcp = np.zeros((35,2)) #weil die Tabelle in dem Dokument 37 zeilen hat. Python fängt von 0 an zu zahlen.

		i=0
		for row in reader:
			mcpr = np.genfromtxt(row)
			mcp[i] = mcpr
			i=i+1
        
		for i in range (0,35):
			for j in range(0,2):
				mcp[i][j]=float(mcp[i][j])





#print 'type mcp = ', type(mcp)

#### folgende Dokumente werden benötigt
#print mcp
#print np.shape(mcp)
#print type(mcp[1][1])
## Kaelteleistung wird vorgegeben, d.h. die enstprechende Kaelteleistung wird aus einem Excel-Datei eingelesen.
# see the dimension of the matrix
# wc -l nameofFile.csv
 
#NumOfData=8 


#with open('Kaelteleistung_messungen_3468.csv') as csvfile:
#with open('Kaelteleistung_messungen_782.csv') as csvfile:
#with open('Kaelteleistung_messungen_7103.csv') as csvfile:
#with open('Kaelteleistung_messungen_8.csv') as csvfile:
with open('Kaelteleistung_messungen_'+str(NumOfData)+'.csv') as csvfile:
    reader = csv.reader(csvfile)
    reader.next()
	# Title
	# 	: one index is lost by reader.next()

    aux1= NumOfData - 1
    messung = np.zeros((aux1,12)) #weil die Tabelle in dem Dokument 37 zeilen hat. Python fängt von 0 an zu zahlen.

    i=0
    for row in reader:
        messungr = np.genfromtxt(row)
	messung[i] = messungr
	i=i+1
    for i in range (0,aux1):
        for j in range(0,12):
			messung[i][j]=float(messung[i][j])

#print messung 
#print np.shape(messung)
#print type(messung[1][1])
####################################
#Qkuehl = messung[:,11] #gemessene Kaelteleistung in der Zeile 12, Python faengt mit 0 anzuzaehlen.
#print np.transpose(Qkuehl)
#print np.shape(Qkuehl)
#print type(Qkuehl[0])

#print 'length of Qkuehl=',len(Qkuehl)
#####print 'longitud',len(Qkuehl)
indexmatritze= int(NumOfData-2)
#print indexmatritze
#print "wert from Matrix[{:d}][{:d}]: {:.10f}\nlength of Qkuehl: {:d}".format(indexmatritze,2,messung[indexmatritze][2], len(Qkuehl)) 
#gemessene Kaeltelestung in der Zeile 12, Python faengt mit 0 anzuzaehlen.
## Zeitdiskretisierung 
#zeit=np.linspace(0,30.,1000)
#tzeit=messung[:,1]

#x[1:2, 1:3]

zeit_nichtfein = messung[:,2] #sie ist benutzt im bekomme_qkuehl
Qkuehl_nicht_interpoliert = messung[:,11]
Qkuehl_average_stationar = np.average (Qkuehl_nicht_interpoliert[651:6004])
Qkuehl_average_instat_proportional = np.average (Qkuehl_nicht_interpoliert[651:NumOfData-1]) # von scan 652 bis scan 781
print "Qkuehl_average_instat_proportional = ",Qkuehl_average_instat_proportional

kkk=400
###############################################################################
######## Zeit
###############################################################################




if choice == 1: # Die Zeitdiskretisierung passt sich an die Abtastrate der Messdaten
	end_time = messung[indexmatritze][2]
	zeit_eval =np.linspace(init_time, end_time, indexmatritze+1) 
	Qkuehl = messung[:,11]


elif choice == 2: #die Berechnung wird durchgeführt mit der Abtastrate der vorgegebenen Messung
	end_time = messung[indexmatritze][2]	
	Qkuehl = messung[:,11]
	zeit_eval = messung[:,2]

elif choice == 3: #stationaer
	end_time = ENDTIME
	zeit_eval = np.linspace(init_time,end_time,NDATEN)
#	Qkuehl = messung[:,11]



###################### bekomme_qkuehl

if Stat_Berechnung == 0: #stationäre Berechnung ausgeschaltet.

	def bekomme_qkuehl(zeit,Automatic_on_off,T23,Ventilator):
#	def bekomme_qkuehl(zeit,Automatic_on_off,T23,Ventilator,Qkuehl):        
		index = np.searchsorted(zeit_nichtfein,zeit)
		index_qkuehl = index-1
		
		if index_qkuehl < 0:
			index_qkuehl = 0
		
		
		
#		Qkuehl = np.interp(zeit, zeit_nichtfein, Qkuehl_nicht_interpoliert)
#		np.insert(Qkuehl, index,Qkuehl[index])
#		print 'holaaaaaaaaQkuehl',Qkuehl	
		
		if Automatic_on_off == 0:	
			# print 'Taktung aus'
			None
		else:

#			print 'Taktung an'

			if  Ventilator ==0: # if ventilator ist aus.
				if T23>310.65: # if T23 > 37,5C
					Qkuehl[index_qkuehl]=Dachkuehlgeraet_Qkuehl
					print "[index_qkuehl] =",[index_qkuehl]
					print "Qkuehl[index_qkuehl] = ", Qkuehl[index_qkuehl]
					print "Qkuehl =", Qkuehl
					Ventilator=1
				else:
					Qkuehl[index_qkuehl]=0
			else:
				if T23<305.65: #if T23< 32.5
					Qkuehl[index_qkuehl]=0
					Ventilator=0
				else:
					Qkuehl[index_qkuehl]=Dachkuehlgeraet_Qkuehl

		return Qkuehl[index_qkuehl],Ventilator

	

















# Stoffwertetabelle auslesen "sw.txt"
with open('sw.csv') as csvfile:
	reader = csv.reader(csvfile)
	reader.next()

	sw = np.zeros((36,5)) #weil die Tabelle in dem Dokument 37 zeilen hat. Python fängt von 0 an zu zahlen.

	i=0
	for row in reader:
		swr = np.genfromtxt(row)
		sw[i] = swr
		i=i+1
        
	for i in range (0,36):
		for j in range(0,5):
			sw[i][j]=float(sw[i][j])
#	print 'sw',sw
#	print 'sw0',sw[:,0]
#	print 'sw[:,1]',sw[:,1]

#Zonen des Schaltschrankes Unterteilen.



#Zonen Schaltschrank
hz=h/3

#Zonen Platte
hp2=hz
hp3=hz-zp
hp1=hp-hp2-hp3


#Schrank innen
hi = h-s*2 # innere Gesamthöhe des Schaltschranks
hzi = hz-s # innere Höhe der oberen und unteren Zone (innere Höhe der mittleren Zone bleibt hz)
bi = b-s*2 # innere Breite des Schaltschranks
ti = t-s*2 # innere Tiefe des Schaltschranks

tiv = t-xp-tp-s # Tiefe innen vor Platte
#print "tiefe innen vor Platte ", tiv
#print (t)
#print (xp)
#print (tp)
#print (s)
tih = xp-s # Tiefe innen hinter Platte
ypl = yp-s # Abstand Platte zu Gehäuse links
ypr = bi-ypl-bp # Abstand Platte zu Gehäuse rechts
zpu = zp-s # Abstand Platte zum Boden
zpo = hi-hp-zpu #Abstand Platte zur Decke



# Massenstrom
m_dot = v_dot/3600*dichte # bei Umgebungsdruck und Umbegungstemperatur


# Strecken von Plattenmitte bis Schrank vorne Bzw. Schrank hinten
tv = tiv+s+tp/2 # vorne
th = tih+s+tp/2 # hinten

# Flächen des Gehäuses (außen)
Aov = (tiv+s+tp/2)*b # Fläche oben auf deckel (auch Boden) von Mitte der Platte bis vorne zum Schaltschrank
Aoh = (tih+s+tp/2)*b # Fläche oben auf deckel (auch Boden) von Mitte der Platte bis hinten zum Schaltschrank
Av = hz*b # Fläche der jeweiligen Zonen vorne vom Schaltschrank
Ah = Av # Fläche der jeweiligen Zonen hinten vom Schaltschrank
Asv = hz*tv # Fläche der jeweiligen Zonen rechts und links vorne vom Schaltschrank
Ash = hz*th # Fläche der jeweiligen Zonen rechts und links hinten vom Schanltschrank



# Flächen des Gehäuses (innen)
# es sind 27 Flächen, d.h. 27 Temperaturknoten.
# die Flächen sind Innenflächen

A = np.zeros((1,28)) #siehe Seite 20 Bachelorarbeit, Deborah Bilgic
A[0][0] = tiv*bi #tiefe interne vorne * breite interne.
A[0][1] = hzi*tiv
A[0][2] = hz*tiv  
A[0][3] = hzi*tiv
A[0][4] = hzi*tiv
A[0][5] = hz*tiv
A[0][6] = hzi*tiv
A[0][7] = tiv*bi
A[0][8] = bi*hzi
A[0][9] = bi*hz
A[0][10] = bi*hzi

A[0][11] = tih*bi
A[0][12] = hzi*tih
A[0][13] = hz*tih
A[0][14] = hzi*tih
A[0][15] = hzi*tih
A[0][16] = hz*tih
A[0][17] = hzi*tih
A[0][18] = tih*bi
A[0][19] = bi*hzi
A[0][20] = bi*hz
A[0][21] = bi*hzi

# Flächen der Platte
A[0][22] = hp1*bp
A[0][23] = hp2*bp
A[0][24] = hp3*bp
A[0][25] = hp1*bp
A[0][26] = hp2*bp
A[0][27] = hp3*bp
#print "valor de A",A 			






### Massen der Gehäuse Aufteilung
masse = np.zeros((1,35)) #siehe Seite 20 Bachelorarbeit, Deborah Bilgic

## Masse der Gehaeuse vorne
masse[0][0] = (tv*b*s)*dichte_g #tiefe vorne*breite*wanddicke
masse[0][1] = (hzi*tv*s)*dichte_g
masse[0][2] = (hz*tv*s)*dichte_g  
masse[0][3] = (hzi*tv*s)*dichte_g
masse[0][4] = (hzi*tv*s)*dichte_g
masse[0][5] = (hz*tv*s)*dichte_g
masse[0][6] = (hzi*tv*s)*dichte_g
masse[0][7] = (tv*b*s)*dichte_g
masse[0][8] = (b*hzi*s)*dichte_g	
masse[0][9] = (b*hz*s)*dichte_g
masse[0][10] = (b*hzi*s)*dichte_g

## Masse der Gehaeuse hinten
masse[0][11] = (th*b*s)*dichte_g #tiefe interne hinten * breite *wanddicke
masse[0][12] = (hzi*th*s)*dichte_g
masse[0][13] = (hz*th*s)*dichte_g  
masse[0][14] = (hzi*th*s)*dichte_g
masse[0][15] = (hzi*th*s)*dichte_g
masse[0][16] = (hz*th*s)*dichte_g
masse[0][17] = (hzi*th*s)*dichte_g
masse[0][18] = (th*b*s)*dichte_g
masse[0][19] = (b*hzi*s)*dichte_g	
masse[0][20] = (b*hz*s)*dichte_g
masse[0][21] = (b*hzi*s)*dichte_g
mgehaeuse = masse[0][0:22].sum()
print 'mgehaeuse= ',mgehaeuse

# Massen der Platte
# masse platte vorne
masse[0][22] = (hp1*bp*tp)*dichte_p 
masse[0][23] = (hp2*bp*tp)*dichte_p
masse[0][24] = (hp3*bp*tp)*dichte_p
#masse platte hinten
masse[0][25] = (hp1*bp*tp/2)*dichte_p
masse[0][26] = (hp2*bp*tp/2)*dichte_p
masse[0][27] = (hp3*bp*tp/2)*dichte_p
mplatte = masse[0][22:28].sum()
print 'mplatte (kg) = ',mplatte


# Massen des Luftes
#masse Luft hinten
masse[0][28] = (Ash*bi)*dichte #knote 29. Tlein se asume que es un cubo de volumen.
masse[0][29] = (Ash*bi)*dichte
masse[0][30] = (Ash*bi)*dichte
masse[0][31] = (Ash*bi)*dichte
#Masse Luft vorne
masse[0][32] = (Asv*bi)*dichte
masse[0][33] = (Asv*bi)*dichte
masse[0][34] = (Asv*bi)*dichte
mluft=masse[0][28:35].sum()
print 'type masse', type(masse)
print 'mluft (Kg) =', mluft

print 'masse = ',masse.transpose()


koko = 400

# Interpolation der Stoffwerte franco
def ip(j,NrSe,Tb): # j = Zählvariable zum Einordnen in die Zeile, NrSe = Wahl der Spalte mit jeweiligen Stoffeigenschaften
	wert = np.interp(Tb,sw[:,0],sw[:,NrSe])
	return wert



# Wärmeübergangskoeffizient

# freie Konvektion bei einer senkrechten Wand
def alpha_frei(lc,Tw1,Tw2,Tw3): # lc = charakteristische Länge
#	print 'Tw1=',Tw1
#	print 'Tw2=',Tw2
#	print 'Tw3=',Tw3

	Tw = (Tw1+Tw2+Tw3)/3
	
	Tetha_b = (Tw+Tamb)/2-273.15	
	Tdiff = Tamb-Tw

	
	if Tdiff<0:
		Tdiff = Tdiff*(-1)
	
	j=0 # Zählvariable zum einordnen in das Stoffwerte-Array
	for i in range (0,35):
		if Tetha_b>sw[i][0]:
			j=i

	# Werte aus der Stoffwertetabelle entnehmen und interpolieren
	viskositaet = ip(j,3,Tetha_b)*10**(-7)		# Viskosität
	pr = ip(j,4,Tetha_b)				# Prandtl-Zahl
	lambda_konv = ip(j,2,Tetha_b)*10**(-3)		# Wärmeleitfähigkeit
	betha = ip(j,1,Tetha_b)*10**(-3)		# Wärmeausdehnungskoeffizient
	
	# Grashof-Zahl
	
	gr = 9.81*betha*Tdiff*lc**3/viskositaet**2

	# Raileigh-Zahl
	ra = gr*pr

	# Nu-Korrelation nach Churchill/Chu
	f1 = (1+(0.492/pr)**double(9.0/16.0))**(-double(16.0/9.0))
	nu = (0.825+0.387*(ra*f1)**double(1.0/6.0))**2


	# Wärmeübergangskoeffizient
	alpha = nu*lambda_konv/lc

	return alpha


# Nummerierung der Bereichsfälle
vorne=1 
hinten=2


# Mischkonvektion bei einer längsangeströmten ebenen Platte
def alpha_misch(bereich,lc,Tw1,Tw2,Tw3,Tl1,Tl2,Tl3): # Bereich: vorne/hinten; lc = charakteristische Länge
	
#	print'Tl1,Tl2,Tl3',Tl1,Tl2,Tl3
#	print'Tw1,Tw2,Tw3',Tw1,Tw2,Tw3
	Tl = (Tl1+Tl2+Tl3)/3
	Tw = (Tw1+Tw2+Tw3)/3
	
	Tetha_b = (Tw+Tl)/2-273.15
#	print 'Tw',Tw
#	print 'Tl',Tl
#	print 'Theta_b',Tetha_b
	Tdiff = Tl-Tw
	if Tdiff<0:
		Tdiff =Tdiff*(-1)		
	
	j=0 # Zählvariable zum einordnen in das Stoffwerte-Array
	for i in range (0,35):
		if Tetha_b>sw[i][0]:
			j=i

	# Geschwindigkeit
	if bereich == 1:
		w = v_dot/(tiv*bi)/3600
#		print 'bereich 1 alpha misch'	
	#print 'w = ' +str(w)

	if bereich == 2:
		w = v_dot/(tih*bi)/3600
#		print 'bereich 2 alpha_misch'

	# Werte aus der Stoffwertetabelle entnehmen und interpolieren
	viskositaet = ip(j,3,Tetha_b)*10**(-7)		# Viskosität
	pr = ip(j,4,Tetha_b)				# Prandtl-Zahl
	lambda_konv = ip(j,2,Tetha_b)*(10)**(-3)	# Wärmeleitfähigkeit
	betha = ip(j,1,Tetha_b)*10**(-3)		# Wärmeausdehnungskoeffizient


	# Erzwungene Konvektion

	# Reynolds-Zahl
	re = w*lc/viskositaet
#	print 'viskositaet=',viskositaet
#	print 'w=',w			
#	print 're=',re
	# Nu-Korrelation, erzwungen,laminar und turbulent
	nu_lam = 0.664*math.sqrt(re)*pr**double(1.0/3.0)
	nu_turb = 0.037*re**0.8*pr/(1+2.443*re**(-0.1)*(pr**double(2.0/3.0)-1))
	nu_erzw = math.sqrt(nu_lam**2+nu_turb**2)
#	print nu_lam
#	print nu_turb
#	print nu_erzw



	#freie Konvektion

	# Grashof-Zahl
	gr = 9.81*betha*Tdiff*lc**3/viskositaet**2
	
	# Raileigh-Zahl
	ra = gr*pr

	# Nu-Korrelation nach Churchill/Chu
	f1 = (1+(0.492/pr)**double(9.0/16.0))**(-double(16.0/9.0))
	nu_frei = (0.825+0.387*(ra*f1)**double(1.0/6.0))**2

	nu_verhaeltnis=nu_frei/nu_erzw
	gr_re_verhaeltnis=gr/(re**2)
	#print 'gr/re^2'+str(gr_re_verhaeltnis)

	# Mischkonvektion
	if bereich == 1:
		nu_misch = np.abs(nu_erzw**3+nu_frei**3)**double(1.0/3.0)
	if bereich == 2:
		nu_misch = np.abs(nu_erzw**3-nu_frei**3)**double(1.0/3.0)
	


	# Wärmeübergangskoeffizient
	if ((nu_verhaeltnis<0.8 and bereich==2) or bereich==1) and (gr_re_verhaeltnis<10 or gr_re_verhaeltnis>0.1):
		alpha = nu_misch*lambda_konv/lc
		#print 'Mischkonvektion'
	elif gr_re_verhaeltnis<1:
		alpha = nu_erzw*lambda_konv/lc
		#print 'erzwungene Konvektion'
	else:
		alpha = nu_frei*lambda_konv/lc
		#print 'freie Konvektion'
	return alpha


# Sichtfaktoren

# 1. Anwendungsfall: Anwendung bei senkrecht zueinander angrenzenden Platten mit einer gemeinsamen Seite
	# a = gemeinsame Seite/Länge; b = Breite A1; c = Breite A2
def sichtfaktor(a,b,c):


	B = b/a
	C = c/a

	if C == 0:
		sf1 = 1/(math.pi*B)*(B*math.atan(1/B)-math.sqrt(B**2)*math.atan(1/math.sqrt(B**2))\
		+(0.25)*(B**2*np.log((1+B**2)*B**2/(B**2*(1+B**2)))\
		-np.log((1+B**2)/((B**2+1)))))

	else:
		sf1 = 1/(math.pi*B)*(B*math.atan(1/B)+C*math.atan(1/C)-math.sqrt(B**2+C**2)*math.atan(1/math.sqrt(B**2+C**2))\
		+(0.25)*(B**2*np.log((1+B**2+C**2)*B**2/((B**2+C**2)*(1+B**2)))\
		+C**2*np.log((1+B**2+C**2)*C**2/((B**2+C**2)*(1+C**2)))-np.log((1+B**2+C**2)/((B**2+1)*(1+C**2)))))
	return sf1


# 2. Anwendungsfall: Anwendung bei senkrecht zueinander liegenden Platten, dessen gleich langen Seiten nicht direkt miteinander angrenzen
	# gs = gemeinsame Seite/Länge;
	# B_2 = Breite von gemeinsamer Seite bis Ende der Rechteckfläche A1;
	# B_1 = Breite von gemeinsamer Seite bis Anfang der Rechteckfläche A1;
	# C_2 = Breite von gemeinsamer Seite bis Ende der Rechteckfläche A2;
	# C_1 = Breite von gemeinsamer Seite bis Anfang der Rechteckfläche A2
def sichtfaktor2(gs,B_2,B_1,C_2,C_1):
	if B_1 == 0:
		sf2 = B_2*(sichtfaktor(gs,B_2,C_2)-sichtfaktor(gs,B_2,C_1))/B_2

	else:
		sf2 = (B_2*(sichtfaktor(gs,B_2,C_2)-sichtfaktor(gs,B_2,C_1))-B_1*(sichtfaktor(gs,B_1,C_2)\
		-sichtfaktor(gs,B_1,C_1)))/(B_2-B_1)
	return sf2


# 3. Anwendungsfall: Anwendung bei parallelen, gleich großen, gegenüberliegenden Rechteckflächen
	# a = Abstand der Flächen
	# b = Breite
	# c = Länge
def sichtfaktor3(a,b,c):

	B = b/a
	C = c/a

	sf3 = 1/math.pi*(1/(B*C)*np.log((1+B**2)*(1+C**2)/(1+B**2+C**2))-2/B*math.atan(C)-2/C*math.atan(B)\
	+2/C*math.sqrt(1+C**2)*math.atan(B/math.sqrt(1+C**2))+2/B*math.sqrt(1+B**2)*math.atan(C/math.sqrt(1+B**2)))
	return sf3

# 4. Anwendungsfall: Anwendung bei parallelen Rechteckflächen unterschiedlicher Größe und/oder nicht 					direkt gegenüberliegend; "geometrie" = parallel (p) oder rechtwinklig (rw)

# Nummerierung der Geometriefälle
p=1 # parallel
rw=2 # rechtwinklig




def sichtfaktor4 (geometrie,z,x1,x2,y1,y2,u1,u2,v1,v2):
	sf4 = 1/((x2-x1)*(y2-y1))*summe(geometrie,z,x1,x2,y1,y2,u1,u2,v1,v2)
	return sf4

def summe(geometrie,z,x1,x2,y1,y2,u1,u2,v1,v2):
	sum=0
	for l in range(1,3):
		for k in range(1,3):
			for j in range(1,3):
				for i in range (1,3):
					if i==1:
						x=x1
					if i==2:
						x=x2
					if j==1:
						y=y1
					if j==2:
						y=y2
					if k==1:
						v=v1
					if k==2:
						v=v2
					if l==1:
						u=u1
					if l==2:
						u=u2

					sum = sum + (-1)**(i+j+k+l)*G(geometrie,z,x,y,u,v)

	return sum


def G(geometrie,z,x,y,u,v):
	if geometrie == 1:
		G = 1/(2*math.pi)*((y-v)*((x-u)**2+z**2)**0.5*math.atan((y-v)/((x-u)**2+z**2)**0.5)\
		+(x-u)*((y-v)**2+z**2)**0.5*math.atan((x-u)/((y-v)**2+z**2)**0.5)\
		-z**2/2*np.log((x-u)**2+(y-v)**2+z**2))

	if geometrie == 2:
		if y-v != 0:
			if x+u != 0:
				K = (y-v)/(x**2+u**2)**0.5
				G = 1/(2*math.pi)*((y-v)*(x**2+u**2)**0.5*math.atan(K)\
				-0.25*((x**2+u**2)*np.log(1+K**2)-(y-v)**2*np.log(1+1/K**2)))
			if x+u == 0:
				G = 1/(2*math.pi)*(0.25*(y-v)**2*np.log(1))
		if y-v == 0:
			G = 1/(2*math.pi)*(-0.25*((x**2+u**2)*np.log(1)))


	return G







phi = np.zeros((28,28))


# 1. Anwendungsfall:

# Zone 1
phi[8][0] = sichtfaktor(bi,hzi,tiv)
phi[0][8] = phi[8][0]*A[0][8]/A[0][0]
phi[19][11] = sichtfaktor(bi,hzi,tih)
phi[11][19] = phi[19][11]*A[0][19]/A[0][11]
phi[0][1] = sichtfaktor(tiv,bi,hzi)
phi[0][4] = phi[0][1]
phi[1][0] = sichtfaktor(tiv,hzi,bi)
phi[4][0] = phi[1][0]
phi[11][12] = sichtfaktor(tih,bi,hzi)
phi[11][15] = phi[11][12]
phi[12][11] = sichtfaktor(tih,hzi,bi)
phi[15][11] = phi[12][11]

phi[8][1] = sichtfaktor(hzi,bi,tiv)
phi[1][8] = phi[8][1]*A[0][8]/A[0][1]
phi[8][4] = phi[8][1]
phi[4][8] = phi[1][8]
phi[19][12] = sichtfaktor(hzi,bi,tih)
phi[12][19] = phi[19][12]*A[0][19]/A[0][12]
phi[19][15] = phi[19][12]
phi[15][19] = phi[12][19]


# Zone 2

phi[9][2] = sichtfaktor(hz,bi,tiv)
phi[2][9] = phi[9][2]*A[0][9]/A[0][2]
phi[9][5] = phi[9][2]
phi[5][9] = phi[2][9]
phi[20][13] = sichtfaktor(hz,bi,tih)
phi[13][20] = phi[20][13]*A[0][20]/A[0][13]
phi[20][16] = phi[20][13]
phi[16][20] = phi[13][20]


# Zone 3

phi[10][3] = sichtfaktor(hzi,bi,tiv)
phi[3][10] = phi[10][3]*A[0][10]/A[0][3]
phi[6][10] = phi[3][10]
phi[10][6] = phi[10][3]
phi[21][14] = sichtfaktor(hzi,bi,tih)
phi[14][21] = phi[21][14]*A[0][21]/A[0][14]
phi[21][17] = phi[21][14]
phi[17][21] = phi[14][21]

phi[10][7] = sichtfaktor(bi,hzi,tiv)
phi[7][10] = phi[10][7]*A[0][10]/A[0][7]
phi[21][18] = sichtfaktor(bi,hzi,tih)
phi[18][21] = phi[21][18]*A[0][21]/A[0][18]
phi[7][3] = sichtfaktor(tiv,bi,hzi)
phi[7][6] = phi[7][3]
phi[3][7] = sichtfaktor(tiv,hzi,bi)
phi[6][7] = phi[3][7]
phi[18][14] = sichtfaktor(tih,bi,hzi)
phi[18][17] = phi[18][14]
phi[14][18] = sichtfaktor(tih,hzi,bi)
phi[17][18] = phi[14][18]


# 2. Anwendungsfall:

# Deckel und Gehäuse vorne

phi[9][0] = sichtfaktor2(bi,hzi+hz,hzi,tiv,0)
phi[0][9] = phi[9][0]*A[0][9]/A[0][0]
phi[10][0] = sichtfaktor2(bi,hi,hzi+hz,tiv,0)
phi[0][10] = phi[10][0]*A[0][10]/A[0][0]

# Deckel und Gehäuse hinten

phi[20][11] = sichtfaktor2(bi,hzi+hz,hzi,tih,0)
phi[11][20] = phi[20][11]*A[0][20]/A[0][11]
phi[21][11] = sichtfaktor2(bi,hi,hzi+hz,tih,0)
phi[11][21] = phi[21][11]*A[0][21]/A[0][11]

# Boden und Gehäuse vorne

phi[7][9] = phi[0][9]
phi[9][7] = phi[9][0]
phi[7][8] = phi[0][10]
phi[8][7] = phi[10][0]

phi[18][20] = phi[11][20]
phi[20][18] = phi[20][11]
phi[19][18] = phi[21][11]
phi[18][19] = phi[11][21]

# Deckel vorne und Gehäuse rechts und links
phi[2][0] = sichtfaktor2(tiv,hz+hzi,hzi,bi,0)
phi[0][2] = phi[2][0]*A[0][2]/A[0][0]
phi[0][5] = phi[0][2]
phi[5][0] = phi[2][0]
phi[3][0] = sichtfaktor2(tiv,hi,hz+hzi,bi,0)
phi[0][3] = phi[3][0]*A[0][3]/A[0][0]
phi[0][6] = phi[0][3]
phi[6][0] = phi[3][0]

# Deckel hinten und Gehäuse rechts und links
phi[13][11] = sichtfaktor2(tih,hz+hzi,hzi,bi,0)
phi[11][13] = phi[13][11]*A[0][13]/A[0][11]
phi[11][16] = phi[11][13]
phi[16][11] = phi[13][11]
phi[14][11] = sichtfaktor2(tih,hi,hz+hzi,bi,0)
phi[11][14] = phi[14][11]*A[0][14]/A[0][11]
phi[11][17] = phi[11][14]
phi[17][11] = phi[14][11]

# Boden vorne und Gehäuse rechts und links

phi[7][2] = phi[0][2]
phi[2][7] = phi[2][0]
phi[7][5] = phi[0][5]
phi[5][7] = phi[5][0]
phi[7][1] = phi[0][3]
phi[1][7] = phi[3][0]
phi[7][4] = phi[0][6]
phi[4][7] = phi[6][0]

# Boden hinten und Gehäuse rechts und links

phi[18][13] = phi[11][13]
phi[13][18] = phi[13][11]
phi[18][16] = phi[11][16]
phi[16][18] = phi[16][11]
phi[18][12] = phi[11][14]
phi[12][18] = phi[14][11]
phi[18][15] = phi[11][17]
phi[15][18] = phi[17][11]

 # Platte vorne und Gehäuse rechts und links

phi[23][5] = sichtfaktor2(hz,ypl+bp,ypl,tiv,0)
phi[5][23] = phi[23][5]*A[0][23]/A[0][5]

phi[23][2] = sichtfaktor2(hz,ypr+bp,ypr,tiv,0)
phi[2][23] = phi[23][2]*A[0][23]/A[0][2]


 # Platte hinten und Gehäuse rechts und links

phi[26][16] = sichtfaktor2(hz,ypl+bp,ypl,tih,0)
phi[16][26] = phi[26][16]*A[0][26]/A[0][16]

phi[26][13] = sichtfaktor2(hz,ypr+bp,ypr,tih,0)
phi[13][26] = phi[26][13]*A[0][26]/A[0][13]


# 3. Anwendungsfall:

# Parallele Zonen der Wand rechts und links vorne

phi[1][4] = sichtfaktor3(bi,hzi,tiv)
phi[4][1] = phi[1][4]
phi[2][5] = sichtfaktor3(bi,hz,tiv)
phi[5][2] = phi[2][5]
phi[6][3] = sichtfaktor3(bi,hzi,tiv)
phi[3][6] = phi[6][3]

# Deckel und Boden vorne
phi[0][7] = sichtfaktor3(hi,tiv,bi)
phi[7][0] = phi[0][7]


# Parallele Zonen der Wand rechts und links hinten

phi[12][15] = sichtfaktor3(bi,hzi,tih)
phi[15][12] = phi[12][15]
phi[13][16] = sichtfaktor3(bi,hz,tih)
phi[16][13] = phi[13][16]
phi[17][14] = sichtfaktor3(bi,hzi,tih)
phi[14][17] = phi[17][14]

# Deckel und Boden vorne
phi[11][18] = sichtfaktor3(hi,tih,bi)
phi[18][11] = phi[11][18]


# 4. Anwendungsfall(parallel):
	# jeweils Einführung eines KOS innerhalb der parallelen Ebenen,Ursprung liegt immer auf der x-Achse des alten KOS
		# 1. Sichtfaktoren zwischen vorderem Gehäuse (x,y-KOS) und Platte (u,v-KOS)
# Gehäuse 8
phi[8][22] = sichtfaktor4(p,tiv,hz*2,h-s,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[22][8] = phi[8][22]*A[0][8]/A[0][22]
phi[8][23] = sichtfaktor4(p,tiv,hz*2,h-s,s,b-s,hz,hz*2,yp,yp+bp)
phi[23][8] = phi[8][23]*A[0][8]/A[0][23]
phi[8][24] = sichtfaktor4(p,tiv,hz*2,h-s,s,b-s,zp,hz,yp,yp+bp)
phi[24][8] = phi[8][24]*A[0][8]/A[0][24]

# Gehäuse 9
phi[9][22] = sichtfaktor4(p,tiv,hz,hz*2,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[22][9] = phi[9][22]*A[0][9]/A[0][22]
phi[9][23] = sichtfaktor4(p,tiv,hz,hz*2,s,b-s,hz,hz*2,yp,yp+bp)
phi[23][9] = phi[9][23]*A[0][9]/A[0][23]
phi[9][24] = sichtfaktor4(p,tiv,hz,hz*2,s,b-s,zp,hz,yp,yp+bp)
phi[24][9] = phi[9][24]*A[0][9]/A[0][24]

# Gehäuse 10
phi[10][22] = sichtfaktor4(p,tiv,s,hz,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[22][10] = phi[10][22]*A[0][10]/A[0][22]
phi[10][23] = sichtfaktor4(p,tiv,s,hz,s,b-s,hz,hz*2,yp,yp+bp)
phi[23][10] = phi[10][23]*A[0][10]/A[0][23]
phi[10][24] = sichtfaktor4(p,tiv,s,hz,s,b-s,zp,hz,yp,yp+bp)
phi[24][10] = phi[10][24]*A[0][10]/A[0][24]


		# 1. Sichtfaktoren zwischen hinterem Gehäuse (u,v-KOS) und Platte (x,y-KOS)
# Gehäuse 19
phi[19][25] = sichtfaktor4(p,tih,hz*2,h-s,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[25][19] = phi[19][25]*A[0][19]/A[0][25]
phi[19][26] = sichtfaktor4(p,tih,hz*2,h-s,s,b-s,hz,hz*2,yp,yp+bp)
phi[26][19] = phi[19][26]*A[0][19]/A[0][26]
phi[19][27] = sichtfaktor4(p,tih,hz*2,h-s,s,b-s,zp,hz,yp,yp+bp)
phi[27][19] = phi[19][27]*A[0][19]/A[0][27]

# Gehäuse 20
phi[20][25] = sichtfaktor4(p,tih,hz,hz*2,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[25][20] = phi[20][25]*A[0][20]/A[0][25]
phi[20][26] = sichtfaktor4(p,tih,hz,hz*2,s,b-s,hz,hz*2,yp,yp+bp)
phi[26][20] = phi[20][26]*A[0][20]/A[0][26]
phi[20][27] = sichtfaktor4(p,tih,hz,hz*2,s,b-s,zp,hz,yp,yp+bp)
phi[27][20] = phi[20][27]*A[0][20]/A[0][27]

# Gehäuse 21
phi[21][25] = sichtfaktor4(p,tih,s,hz,s,b-s,hz*2,hz*2+hp1,yp,yp+bp)
phi[25][21] = phi[21][25]*A[0][21]/A[0][25]
phi[21][26] = sichtfaktor4(p,tih,s,hz,s,b-s,hz,hz*2,yp,yp+bp)
phi[26][21] = phi[21][26]*A[0][21]/A[0][26]
phi[21][27] = sichtfaktor4(p,tih,s,hz,s,b-s,zp,hz,yp,yp+bp)
phi[27][21] = phi[21][27]*A[0][21]/A[0][27]

# Sichtfaktoren zwischen Gehäuse rechts und links vorne

phi[1][5] = sichtfaktor4(p,bi,0,tiv,hzi+hz,hi,0,tiv,hzi,hzi+hz)
phi[5][1] = phi[1][5]*A[0][1]/A[0][5]

phi[1][6] = sichtfaktor4(p,bi,0,tiv,hzi+hz,hi,0,tiv,0,hzi)
phi[6][1] = phi[1][6]*A[0][1]/A[0][6]

phi[2][4] = sichtfaktor4(p,bi,0,tiv,hzi,hzi+hz,0,tiv,hzi+hz,hi)
phi[4][2] = phi[2][4]*A[0][2]/A[0][4]

phi[2][6] = sichtfaktor4(p,bi,0,tiv,hzi,hzi+hz,0,tiv,0,hzi)
phi[6][2] = phi[2][6]*A[0][2]/A[0][6]

phi[3][4] = sichtfaktor4(p,bi,0,tiv,0,hzi,0,tiv,hzi+hz,hi)
phi[4][3] = phi[3][4]*A[0][3]/A[0][4]

phi[3][5] = sichtfaktor4(p,bi,0,tiv,0,hzi,0,tiv,hzi,hzi+hz)
phi[5][3] = phi[3][5]*A[0][3]/A[0][5]


# Sichtfaktoren zwischen Gehäuse rechts und links hinten

phi[12][16] = sichtfaktor4(p,bi,0,tih,hzi+hz,hi,0,tih,hzi,hzi+hz)
phi[16][12] = phi[12][16]*A[0][12]/A[0][16]

phi[12][17] = sichtfaktor4(p,bi,0,tih,hzi+hz,hi,0,tih,0,hzi)
phi[17][12] = phi[12][17]*A[0][12]/A[0][17]

phi[13][15] = sichtfaktor4(p,bi,0,tih,hzi,hzi+hz,0,tih,hzi+hz,hi)
phi[15][13] = phi[13][15]*A[0][13]/A[0][15]

phi[13][17] = sichtfaktor4(p,bi,0,tih,hzi,hzi+hz,0,tih,0,hzi)
phi[17][13] = phi[13][17]*A[0][13]/A[0][17]

phi[14][15] = sichtfaktor4(p,bi,0,tih,0,hzi,0,tih,hzi+hz,hi)
phi[15][14] = phi[14][15]*A[0][14]/A[0][15]

phi[14][16] = sichtfaktor4(p,bi,0,tih,0,hzi,0,tih,hzi,hzi+hz)
phi[16][14] = phi[14][16]*A[0][14]/A[0][16]


# 4. Anwendungsfall(rechtwinklig):
	# jeweils Einführung eines KOS innerhalb der rechtwinkligen Ebenen,Ursprung liegt innen zwischen vorderem und linken 			gehäuse unten
		#z=0
		# 1. Sichtfaktoren zwischen vorderem Gehäuse und Gehäuse links und rechts
# Gehäuse 8

phi[8][5] = sichtfaktor4(rw,0,0,bi,hzi+hz,hi,0,tiv,hzi,hzi+hz)
phi[8][2] = phi[8][5]
phi[5][8] = phi[8][5]*A[0][8]/A[0][5]
phi[2][8] = phi[5][8]

phi[8][6] = sichtfaktor4(rw,0,0,bi,hzi+hz,hi,0,tiv,0,hzi)
phi[8][3] = phi[8][6]
phi[6][8] = phi[8][6]*A[0][8]/A[0][6]
phi[3][8] = phi[6][8]

# Gehäuse 9

phi[9][4] = sichtfaktor4(rw,0,0,bi,hzi,hzi+hz,0,tiv,hzi+hz,hi)
phi[9][1] = phi[9][4]
phi[4][9] = phi[9][4]*A[0][9]/A[0][4]
phi[1][9] = phi[4][9]

phi[9][6] = sichtfaktor4(rw,0,0,bi,hzi,hzi+hz,0,tiv,0,hzi)
phi[9][3] = phi[9][6]
phi[6][9] = phi[9][6]*A[0][9]/A[0][6]
phi[3][9] = phi[6][9]

# Gehäuse 10
phi[10][1] = phi[8][3]
phi[10][4] = phi[8][6]
phi[1][10] = phi[3][8]
phi[4][10] = phi[6][8]

phi[10][2] = phi[8][2]
phi[10][5] = phi[8][5]
phi[2][10] = phi[2][8]
phi[5][10] = phi[5][8]



		# 2. Sichtfaktoren zwischen hinterem Gehäuse und Gehäuse links und rechts
# Gehäuse 19

phi[19][16] = sichtfaktor4(rw,0,0,bi,hz+hzi,hi,0,tih,hzi,hzi+hz)
phi[19][13] = phi[19][16]
phi[16][19] = phi[19][16]*A[0][19]/A[0][16]
phi[13][19] = phi[16][19]

phi[19][17] = sichtfaktor4(rw,0,0,bi,hz+hzi,hi,0,tih,0,hzi)
phi[19][14] = phi[19][17]
phi[17][19] = phi[19][17]*A[0][19]/A[0][17]
phi[14][19] = phi[17][19]

# Gehäuse 20

phi[20][15] = sichtfaktor4(rw,0,0,bi,hzi,hzi+hz,0,tih,hzi+hz,hi)
phi[20][12] = phi[20][15]
phi[15][20] = phi[20][15]*A[0][20]/A[0][15]
phi[12][20] = phi[15][20]

phi[20][17] = sichtfaktor4(rw,0,0,bi,hzi,hzi+hz,0,tih,0,hzi)
phi[20][14] = phi[20][17]
phi[17][20] = phi[20][17]*A[0][20]/A[0][17]
phi[14][20] = phi[17][20]

# Gehäuse 21
phi[21][12] = phi[19][14]
phi[21][15] = phi[19][17]
phi[12][21] = phi[14][19]
phi[15][21] = phi[17][19]

phi[21][13] = phi[19][13]
phi[21][16] = phi[19][16]
phi[13][21] = phi[13][19]
phi[16][21] = phi[16][19]


		# 3. Sichtfaktoren zwischen platte vorne und Gehäuse links und rechts/Deckel/Boden
# Platte 1
phi[22][0] = sichtfaktor4(rw,0,zpo,hzi,ypl,bp+ypl,0,tiv,0,bi) # KOS an obere Kante links des Schranks
phi[0][22] = phi[22][0]*A[0][22]/A[0][0]

phi[22][7] = sichtfaktor4(rw,0,hzi+hz,hi-zpo,ypl,bp+ypl,0,tiv,0,bi) # KOS an untere Kante links des Schranks
phi[7][22] = phi[22][7]*A[0][22]/A[0][7]


phi[22][4] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tiv,hzi+hz,hi)
phi[4][22] = phi[22][4]*A[0][22]/A[0][4]

phi[22][1] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tiv,hzi+hz,hi)
phi[1][22] = phi[22][1]*A[0][22]/A[0][1]


phi[22][5] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tiv,hzi,hzi+hz)
phi[5][22] = phi[22][5]*A[0][22]/A[0][5]

phi[22][2] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tiv,hzi,hzi+hz)
phi[2][22] = phi[22][2]*A[0][22]/A[0][2]


phi[22][6] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tiv,0,hzi)
phi[6][22] = phi[22][6]*A[0][22]/A[0][6]

phi[22][3] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tiv,0,hzi)
phi[3][22] = phi[22][3]*A[0][22]/A[0][3]

# Platte 2
phi[23][0] = sichtfaktor4(rw,0,hzi,hz+hzi,ypl,bp+ypl,0,tiv,0,bi) # KOS an obere Kante links des Schranks
phi[0][23] = phi[23][0]*A[0][23]/A[0][0]

phi[23][7] = sichtfaktor4(rw,0,hzi,hz+hzi,ypl,bp+ypl,0,tiv,0,bi) # KOS an untere Kante links des Schranks
phi[7][23] = phi[23][7]*A[0][23]/A[0][7]


phi[23][4] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi,hzi+hz,0,tiv,hzi+hz,hi)
phi[4][23] = phi[23][4]*A[0][23]/A[0][4]

phi[23][1] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi,hzi+hz,0,tiv,hzi+hz,hi)
phi[1][23] = phi[23][1]*A[0][23]/A[0][1]


phi[23][6] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi,hzi+hz,0,tiv,0,hzi)
phi[6][23] = phi[23][6]*A[0][23]/A[0][6]

phi[23][3] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi,hzi+hz,0,tiv,0,hzi)
phi[3][23] = phi[23][3]*A[0][23]/A[0][3]


# Platte 3
phi[24][0] = sichtfaktor4(rw,0,hz+hzi,hi-zpu,ypl,bp+ypl,0,tiv,0,bi) # KOS an obere Kante links des Schranks
phi[0][24] = phi[24][0]*A[0][24]/A[0][0]

phi[24][7] = sichtfaktor4(rw,0,zpu,hzi,ypl,bp+ypl,0,tiv,0,bi) # KOS an untere Kante links des Schranks
phi[7][24] = phi[24][7]*A[0][24]/A[0][7]


phi[24][4] = sichtfaktor4(rw,0,ypl,bp+ypl,zpu,hzi,0,tiv,hzi+hz,hi)
phi[4][24] = phi[24][4]*A[0][24]/A[0][4]

phi[24][1] = sichtfaktor4(rw,0,ypr,bp+ypr,zpu,hzi,0,tiv,hzi+hz,hi)
phi[1][24] = phi[24][1]*A[0][24]/A[0][1]


phi[24][5] = sichtfaktor4(rw,0,ypl,bp+ypl,0,hzi,0,tiv,hzi,hzi+hz)
phi[5][24] = phi[24][5]*A[0][24]/A[0][5]

phi[24][2] = sichtfaktor4(rw,0,ypr,bp+ypr,0,hzi,0,tiv,hzi,hzi+hz)
phi[2][24] = phi[24][2]*A[0][24]/A[0][2]


phi[24][6] = sichtfaktor4(rw,0,ypl,bp+ypl,0,hzi,0,tiv,0,hzi)
phi[6][24] = phi[24][6]*A[0][24]/A[0][6]

phi[24][3] = sichtfaktor4(rw,0,ypr,bp+ypr,0,hzi,0,tiv,0,hzi)
phi[3][24] = phi[24][3]*A[0][24]/A[0][3]


# 3. Sichtfaktoren zwischen platte vorne und Gehäuse links und rechts/Deckel/Boden
# Platte 1
phi[25][11] = sichtfaktor4(rw,0,zpo,hzi,ypl,bp+ypl,0,tih,0,bi) # KOS an obere Kante links des Schranks
phi[11][25] = phi[25][11]*A[0][25]/A[0][11]

phi[25][18] = sichtfaktor4(rw,0,hzi+hz,hi-zpo,ypl,bp+ypl,0,tih,0,bi) # KOS an untere Kante links des Schranks
phi[18][25] = phi[25][18]*A[0][25]/A[0][18]


phi[25][15] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tih,hzi+hz,hi)
phi[15][25] = phi[25][15]*A[0][25]/A[0][15]

phi[25][12] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tih,hzi+hz,hi)
phi[12][25] = phi[25][12]*A[0][25]/A[0][12]


phi[25][16] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tih,hzi,hzi+hz)
phi[16][25] = phi[25][16]*A[0][25]/A[0][16]

phi[25][13] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tih,hzi,hzi+hz)
phi[13][25] = phi[25][13]*A[0][25]/A[0][13]


phi[25][17] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi+hz,hi-zpo,0,tih,0,hzi)
phi[17][25] = phi[25][17]*A[0][25]/A[0][17]

phi[25][14] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi+hz,hi-zpo,0,tih,0,hzi)
phi[14][25] = phi[25][14]*A[0][25]/A[0][14]

# Platte 2
phi[26][11] = sichtfaktor4(rw,0,hzi,hzi+hz,ypl,bp+ypl,0,tih,0,bi) # KOS an obere Kante links des Schranks
phi[11][26] = phi[26][11]*A[0][26]/A[0][11]

phi[26][18] = sichtfaktor4(rw,0,hzi,hzi+hz,ypl,bp+ypl,0,tih,0,bi) # KOS an untere Kante links des Schranks
phi[18][26] = phi[26][18]*A[0][26]/A[0][18]


phi[26][15] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi,hzi+hz,0,tih,hzi+hz,hi)
phi[15][26] = phi[26][15]*A[0][26]/A[0][15]

phi[26][12] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi,hzi+hz,0,tih,hzi+hz,hi)
phi[12][26] = phi[26][12]*A[0][26]/A[0][12]


phi[26][17] = sichtfaktor4(rw,0,ypl,bp+ypl,hzi,hzi+hz,0,tih,0,hzi)
phi[17][26] = phi[26][17]*A[0][26]/A[0][17]

phi[26][14] = sichtfaktor4(rw,0,ypr,bp+ypr,hzi,hzi+hz,0,tih,0,hzi)
phi[14][26] = phi[26][14]*A[0][26]/A[0][14]

# Platte 3
phi[27][11] = sichtfaktor4(rw,0,hzi+hz,hi-zpu,ypl,bp+ypl,0,tih,0,bi) # KOS an obere Kante links des Schranks
phi[11][27] = phi[27][11]*A[0][27]/A[0][11]


phi[27][18] = sichtfaktor4(rw,0,zpu,hzi,ypl,bp+ypl,0,tih,0,bi) # KOS an untere Kante links des Schranks
phi[18][27] = phi[27][18]*A[0][27]/A[0][18]


phi[27][15] = sichtfaktor4(rw,0,ypl,bp+ypl,zpu,hzi,0,tih,hzi+hz,hi)
phi[15][27] = phi[27][15]*A[0][27]/A[0][15]

phi[27][12] = sichtfaktor4(rw,0,ypr,bp+ypr,zpu,hzi,0,tih,hzi+hz,hi)
phi[12][27] = phi[27][12]*A[0][27]/A[0][12]


phi[27][16] = sichtfaktor4(rw,0,ypl,bp+ypl,0,hzi,0,tih,hzi,hzi+hz)
phi[16][27] = phi[27][16]*A[0][27]/A[0][16]

phi[27][13] = sichtfaktor4(rw,0,ypr,bp+ypr,0,hzi,0,tih,hzi,hzi+hz)
phi[13][27] = phi[27][13]*A[0][27]/A[0][13]


phi[27][17] = sichtfaktor4(rw,0,ypl,bp+ypl,0,hzi,0,tih,0,hzi)
phi[17][27] = phi[27][17]*A[0][27]/A[0][17]

phi[27][14] = sichtfaktor4(rw,0,ypr,bp+ypr,0,hzi,0,tih,0,hzi)
phi[14][27] = phi[27][14]*A[0][27]/A[0][14]


#Angleichen der Sichtfaktoren
n=0
for i in range (0,27):
	n=n+1
	for j in range (0,n):
		phi[i+1][j]=(phi[i+1][j]*A[0][i+1]+phi[j][i+1]*A[0][j])/(2*A[0][i+1])
		phi[j][i+1]=(phi[i+1][j]*A[0][i+1]+phi[j][i+1]*A[0][j])/(2*A[0][j])
# Reziprokregel
for j in range (0,28):
	for i in range (0,28):
		rp1 = phi[i][j]*A[0][i]
		rp2 = phi[j][i]*A[0][j]
		if rp1 <> rp2:
			None
#			phi[j][i]=rp1/A[0][j]
#			rp1 = phi[i][j]*A[0][i]
#			rp2 = phi[j][i]*A[0][j]
#			print ("rp1-rp2",rp1-rp2)
#			print '[' + str(i) + ']'+'[' + str(j) + ']'
#			print' '



#### Summenregel
phi_summen = np.zeros((28,1))
###print(phi_summen[0][0])


for j in range (0,28):
	for i in range (0,28):
		phi_summen[j][0] = phi_summen[j][0] + phi[j][i]
#print "phi_summen = ", phi_summen

for i in range(0,28):
	phi[i][i]=1-phi_summen[i]

phi_summen = np.zeros((28,1))
for j in range (0,28):
	for i in range (0,28):
		phi_summen[j][0] = phi_summen[j][0] + phi[j][i]

#print phi
#for i in range (0,28):
#	print 'Phi' + str(i) + ' = ' + str(phi_summen [i])
####print (phi_summen)





######### double 
tv=double(tv)
s =double(s)
t =double(t)
b =double(b)
hz=double(hz)
lambda_g=double(lambda_g)

t1 = time.time()
print ' Zeit für die Initialisierung',t1



def system(zeit,v):
	T = np.zeros((35)) 
	#Gehäuse
	T[0] = v[0]
	T[1] = v[1]
	T[2] = v[2]
	T[3] = v[3]
	T[4] = v[4]
	T[5] = v[5]
	T[6] = v[6]
	T[7] = v[7]
	T[8] = v[8]
	T[9] = v[9]
	T[10] = v[10]
	T[11] = v[11]
	T[12] = v[12]
	T[13] = v[13]
	T[14] = v[14]
	T[15] = v[15]
	T[16] = v[16]
	T[17] = v[17]
	T[18] = v[18]
	T[19] = v[19]
	T[20] = v[20]
	T[21] = v[21]

	#Platte
	T[22] = v[22]
	T[23] = v[23]
	T[24] = v[24]

	T[25] = v[25]
	T[26] = v[26]
	T[27] = v[27]
	
	# Luft 

	T[28]= v[28]  #Tlein 
	T[29] = v[29] #Tl0   
	T[30] = v[30] #Tl1   
	T[31] = v[31] #Tl2   
	T[32] = v[32] #Tl3   
	T[33] = v[33] #Tl4   
	T[34] = v[34] #Tl5   



	global alpha_misch_wi_v 
	global alpha_misch_wi_vr 
	global alpha_misch_wi_vl 
	global alpha_misch_pv 
	global alpha_misch_wi_hr 
	global alpha_misch_wi_hl 
	global alpha_misch_wi_h 
	global alpha_misch_ph 
	global alpha_frei_aussen_vr 
	global alpha_frei_aussen_vl 
	global alpha_frei_aussen_v 
	global alpha_frei_aussen_hr 
	global alpha_frei_aussen_hl 
	global alpha_frei_aussen_h 
	
	

	# Mischkonvektion
	alpha_misch_wi_v = alpha_misch(vorne,hi,T[8],T[9],T[10],T[32],T[33],T[34])
	alpha_misch_wi_vr = alpha_misch(vorne,hi,T[1],T[2],T[3],T[32],T[33],T[34])
	alpha_misch_wi_vl = alpha_misch(vorne,hi,T[4],T[5],T[6],T[32],T[33],T[34])
	alpha_misch_pv = alpha_misch(vorne,hp,T[22],T[23],T[24],T[32],T[33],T[34])
	alpha_misch_wi_hr = alpha_misch(hinten,hi,T[12],T[13],T[14],T[29],T[30],T[31])
	alpha_misch_wi_hl = alpha_misch(hinten,hi,T[15],T[16],T[17],T[29],T[30],T[31])
	alpha_misch_wi_h = alpha_misch(hinten,hi,T[19],T[20],T[21],T[29],T[30],T[31])
	alpha_misch_ph = alpha_misch(hinten,hp,T[22],T[23],T[24],T[29],T[30],T[31])

#	# freie Konvektion
	alpha_frei_aussen_vr = alpha_frei(h,T[1],T[2],T[3])
	alpha_frei_aussen_vl = alpha_frei(h,T[4],T[5],T[6])
	alpha_frei_aussen_v = alpha_frei(h,T[8],T[9],T[10])
	alpha_frei_aussen_hr = alpha_frei(h,T[12],T[13],T[14])
	alpha_frei_aussen_hl = alpha_frei(h,T[15],T[16],T[17])
	alpha_frei_aussen_h = alpha_frei(h,T[19],T[20],T[21])

	


######## Q strahlung
	global summ_Q_str
	global Q_str
	Q_str=np.zeros((28,1),dtype=double)
	kronecker_delta = np.eye(28, dtype=double)
	A_Matrix = np.zeros((28,28), dtype=double)
	C_Matrix = np.zeros((28,), dtype=double)
#	print 'phi',phi
#	print 'phi dimension',np.shape(phi)

	for i in range (0,28):
		for j in range (0, 28):
			#emissionsgrad gehäuse
			e = eg
			if i >=22:
				#emissionsgrad platte
				e = ep
			A_Matrix[i,j] = kronecker_delta[i][j]-(1-e)*phi[i][j]
		C_Matrix[i] = e*Cs*T[i]**4

	x= np.linalg.solve(A_Matrix,C_Matrix) # Helligkeit
	
	for j in range (0,28): 
		e = eg
		if j >=22:
			e = ep
		Q_str[j] = A[0][j]*e/(1-e)*(Cs*T[j]**4-x[j])
	summ_Q_str=Q_str[:].sum()	
#	print 'helligkeit',x
#	print 'helligkeit dim',np.shape(x)	
#	print 'Q_str',Q_str
#	print 'Q_str_dimension',np.shape(Q_str)


############### Qgehaeuse

	#Berechnung der Wärmeströme über Gehäuaseoberfläche
	global Q_gehaeuse
	global summe_Q_gehaeuse
	laenge = 22
	Q_gehaeuse = np.zeros((laenge,), dtype=double)
	Q_gehaeuse[0] = 0
	Q_gehaeuse[1] = alpha_frei_aussen_vr*Asv*(Tamb-T[1])+eg*Asv*Cs*(Tamb**4-T[1]**4)
	Q_gehaeuse[2] = alpha_frei_aussen_vr*Asv*(Tamb-T[2])+eg*Asv*Cs*(Tamb**4-T[2]**4)
	Q_gehaeuse[3] = alpha_frei_aussen_vr*Asv*(Tamb-T[3])+eg*Asv*Cs*(Tamb**4-T[3]**4)
	Q_gehaeuse[4] =	alpha_frei_aussen_vl*Asv*(Tamb-T[4])+eg*Asv*Cs*(Tamb**4-T[4]**4)
	Q_gehaeuse[5] = alpha_frei_aussen_vl*Asv*(Tamb-T[5])+eg*Asv*Cs*(Tamb**4-T[5]**4)
	Q_gehaeuse[6] = alpha_frei_aussen_vl*Asv*(Tamb-T[6])+eg*Asv*Cs*(Tamb**4-T[6]**4)
	Q_gehaeuse[7] = 0
	Q_gehaeuse[8] = alpha_frei_aussen_v*Av*(Tamb-T[8])+eg*Av*Cs*(Tamb**4-T[8]**4)
	Q_gehaeuse[9] = alpha_frei_aussen_v*Av*(Tamb-T[9])+eg*Av*Cs*(Tamb**4-T[9]**4)
	Q_gehaeuse[10] = alpha_frei_aussen_v*Av*(Tamb-T[10])+eg*Av*Cs*(Tamb**4-T[10]**4)
	Q_gehaeuse[11] = 0
	Q_gehaeuse[12] = alpha_frei_aussen_hr*Ash*(Tamb-T[12])+eg*Ash*Cs*(Tamb**4-T[12]**4)
	Q_gehaeuse[13] = alpha_frei_aussen_hr*Ash*(Tamb-T[13])+eg*Ash*Cs*(Tamb**4-T[13]**4)
	Q_gehaeuse[14] = alpha_frei_aussen_hr*Ash*(Tamb-T[14])+eg*Ash*Cs*(Tamb**4-T[14]**4)
	Q_gehaeuse[15] = alpha_frei_aussen_hl*Ash*(Tamb-T[15])+eg*Ash*Cs*(Tamb**4-T[15]**4)
	Q_gehaeuse[16] = alpha_frei_aussen_hl*Ash*(Tamb-T[16])+eg*Ash*Cs*(Tamb**4-T[16]**4)
	Q_gehaeuse[17] = alpha_frei_aussen_hl*Ash*(Tamb-T[17])+eg*Ash*Cs*(Tamb**4-T[17]**4)
	Q_gehaeuse[18] = 0			
	Q_gehaeuse[19] = alpha_frei_aussen_h*Ah*(Tamb-T[19])+eg*Ah*Cs*(Tamb**4-T[19]**4)
	Q_gehaeuse[20] = alpha_frei_aussen_h*Ah*(Tamb-T[20])+eg*Ah*Cs*(Tamb**4-T[20]**4)
	Q_gehaeuse[21] = alpha_frei_aussen_h*Ah*(Tamb-T[21])+eg*Ah*Cs*(Tamb**4-T[21]**4)

	summe_Q_gehaeuse = Q_gehaeuse[:].sum()
#	print 'Q_gehaeuse=',Q_gehaeuse
#	print 'Summe Q_gehauese', summe_Q_gehaeuse

	
	


	global dT
	dT = np.zeros((35)) 


	dT[0] = (-Q_str[0]\
		+lambda_g*b*s*(T[8]-T[0])/(0.5*(tv+hz))+lambda_g*b*s*(T[11]-T[0])/(0.5*t)\
		+lambda_g*tv*s*(T[4]-T[0])/(0.5*(b+hz))+lambda_g*tv*s*(T[1]-T[0])/(0.5*(b+hz)))/(mcp[0][0]*mcp[0][1])



	dT[1] = (alpha_frei_aussen_vr*Asv*(Tamb-T[1])\
		+alpha_misch_wi_vr*A[0][1]*(T[34]-T[1])\
		-Q_str[1]\
		+eg*Asv*Cs*(Tamb**4-T[1]**4)\
		+lambda_g*tv*s*(T[0]-T[1])/(0.5*(b+hz))+lambda_g*tv*s*(T[2]-T[1])/hz\
		+lambda_g*hz*s*(T[8]-T[1])/(0.5*(b+tv))+lambda_g*hz*s*(T[12]-T[1])/(0.5*t))/(mcp[1][0]*mcp[1][1])

	dT[2] = (alpha_frei_aussen_vr*Asv*(Tamb-T[2])\
		+alpha_misch_wi_vr*A[0][2]*(T[33]-T[2])\
		-Q_str[2]\
		+eg*Asv*Cs*(Tamb**4-T[2]**4)\
		+lambda_g*tv*s*(T[1]-T[2])/hz+lambda_g*tv*s*(T[3]-T[2])/hz\
		+lambda_g*hz*s*(T[9]-T[2])/(0.5*(b+tv))+lambda_g*hz*s*(T[13]-T[2])/(0.5*t))/(mcp[2][0]*mcp[2][1])

	dT[3] = (alpha_frei_aussen_vr*Asv*(Tamb-T[3])\
		+alpha_misch_wi_vr*A[0][3]*(T[32]-T[3])\
		-Q_str[3]\
		+eg*Asv*Cs*(Tamb**4-T[3]**4)\
		+lambda_g*tv*s*(T[2]-T[3])/hz+lambda_g*tv*s*(T[7]-T[3])/(0.5*(hz+b))\
		+lambda_g*hz*s*(T[10]-T[3])/(0.5*(b+tv))+lambda_g*hz*s*(T[14]-T[3])/(0.5*t))/(mcp[3][0]*mcp[3][1])



	dT[4] =	(alpha_frei_aussen_vl*Asv*(Tamb-T[4])\
		+alpha_misch_wi_vl*A[0][4]*(T[34]-T[4])\
		-Q_str[4]\
		+eg*Asv*Cs*(Tamb**4-T[4]**4)\
		+lambda_g*tv*s*(T[0]-T[4])/(0.5*(b+hz))+lambda_g*tv*s*(T[5]-T[4])/hz\
		+lambda_g*hz*s*(T[8]-T[4])/(0.5*(b+tv))+lambda_g*hz*s*(T[15]-T[4])/(0.5*t))/(mcp[4][0]*mcp[4][1])

	dT[5] = (alpha_frei_aussen_vl*Asv*(Tamb-T[5])\
		+alpha_misch_wi_vl*A[0][5]*(T[33]-T[5])\
		-Q_str[5]\
		+eg*Asv*Cs*(Tamb**4-T[5]**4)\
		+lambda_g*tv*s*(T[4]-T[5])/hz+lambda_g*tv*s*(T[6]-T[5])/hz\
		+lambda_g*hz*s*(T[9]-T[5])/(0.5*(b+tv))+lambda_g*hz*s*(T[16]-T[5])/(0.5*t))/(mcp[5][0]*mcp[5][1])

	dT[6] = (alpha_frei_aussen_vl*Asv*(Tamb-T[6])\
		+alpha_misch_wi_vl*A[0][6]*(T[32]-T[6])\
		-Q_str[6]\
		+eg*Asv*Cs*(Tamb**4-T[6]**4)\
		+lambda_g*tv*s*(T[5]-T[6])/hz+lambda_g*tv*s*(T[7]-T[6])/(0.5*(hz+b))\
		+lambda_g*hz*s*(T[10]-T[6])/(0.5*(b+tv))+lambda_g*hz*s*(T[17]-T[6])/(0.5*t))/(mcp[6][0]*mcp[6][1])

	dT[7] = (-Q_str[7]\
		+lambda_g*b*s*(T[10]-T[7])/(0.5*(tv+hz))+lambda_g*b*s*(T[18]-T[7])/(0.5*t)\
		+lambda_g*tv*s*(T[3]-T[7])/(0.5*(b+hz))+lambda_g*tv*s*(T[6]-T[7])/(0.5*(b+hz)))/(mcp[7][0]*mcp[7][1])

	dT[8] = (alpha_frei_aussen_v*Av*(Tamb-T[8])\
		+alpha_misch_wi_v*A[0][8]*(T[34]-T[8])\
		-Q_str[8]\
		+eg*Av*Cs*(Tamb**4-T[8]**4)\
		+lambda_g*b*s*(T[0]-T[8])/(0.5*(tv+hz))+lambda_g*b*s*(T[9]-T[8])/hz\
		+lambda_g*hz*s*(T[4]-T[8])/(0.5*(b+tv))+lambda_g*hz*s*(T[1]-T[8])/(0.5*(b+tv)))/(mcp[8][0]*mcp[8][1])

	dT[9] = (alpha_frei_aussen_v*Av*(Tamb-T[9])\
		+alpha_misch_wi_v*A[0][9]*(T[33]-T[9])\
		-Q_str[9]\
		+eg*Av*Cs*(Tamb**4-T[9]**4)\
		+lambda_g*b*s*(T[8]-T[9])/hz+lambda_g*b*s*(T[10]-T[9])/hz\
		+lambda_g*hz*s*(T[5]-T[9])/(0.5*(b+tv))+lambda_g*hz*s*(T[2]-T[9])/(0.5*(b+tv)))/(mcp[9][0]*mcp[9][1]) #Qstr_[9] wird als Qnetto,9 im Abschlussarbeit geschrieben, weil es als abkürzung für Q netto strahlung steht.

	dT[10] = (alpha_frei_aussen_v*Av*(Tamb-T[10])\
		+alpha_misch_wi_v*A[0][10]*(T[32]-T[10])\
		-Q_str[10]\
		+eg*Av*Cs*(Tamb**4-T[10]**4)\
		+lambda_g*b*s*(T[9]-T[10])/hz+lambda_g*b*s*(T[7]-T[10])/(0.5*(hz+tv))\
		+lambda_g*hz*s*(T[6]-T[10])/(0.5*(b+tv))+lambda_g*hz*s*(T[3]-T[10])/(0.5*(b+tv)))/(mcp[10][0]*mcp[10][1])

	dT[11] = (-Q_str[11]\
		+lambda_g*b*s*(T[19]-T[11])/(0.5*(th+hz))+lambda_g*b*s*(T[0]-T[11])/(0.5*t)\
		+lambda_g*th*s*(T[15]-T[11])/(0.5*(b+hz))+lambda_g*th*s*(T[12]-T[11])/(0.5*(b+hz)))/(mcp[11][0]*mcp[11][1])##2mal T15

	dT[12] = (alpha_frei_aussen_hr*Ash*(Tamb-T[12])\
		+alpha_misch_wi_hr*A[0][12]*(T[29]-T[12])\
		-Q_str[12]\
		+eg*Ash*Cs*(Tamb**4-T[12]**4)\
		+lambda_g*th*s*(T[11]-T[12])/(0.5*(b+hz))+lambda_g*th*s*(T[13]-T[12])/hz\
		+lambda_g*hz*s*(T[19]-T[12])/(0.5*(b+th))+lambda_g*hz*s*(T[1]-T[12])/(0.5*t))/(mcp[12][0]*mcp[12][1])

	dT[13] = (alpha_frei_aussen_hr*Ash*(Tamb-T[13])\
		+alpha_misch_wi_hr*A[0][13]*(T[30]-T[13])\
		-Q_str[13]\
		+eg*Ash*Cs*(Tamb**4-T[13]**4)\
		+lambda_g*th*s*(T[12]-T[13])/hz+lambda_g*th*s*(T[14]-T[13])/hz\
		+lambda_g*hz*s*(T[20]-T[13])/(0.5*(b+th))+lambda_g*hz*s*(T[2]-T[13])/(0.5*t))/(mcp[13][0]*mcp[13][1])

	dT[14] = (alpha_frei_aussen_hr*Ash*(Tamb-T[14])\
		+alpha_misch_wi_hr*A[0][14]*(T[31]-T[14])\
		-Q_str[14]\
		+eg*Ash*Cs*(Tamb**4-T[14]**4)\
		+lambda_g*th*s*(T[13]-T[14])/hz+lambda_g*th*s*(T[18]-T[14])/(0.5*(hz+b))\
		+lambda_g*hz*s*(T[21]-T[14])/(0.5*(b+th))+lambda_g*hz*s*(T[3]-T[14])/(0.5*t))/(mcp[14][0]*mcp[14][1])

	dT[15] = (alpha_frei_aussen_hl*Ash*(Tamb-T[15])\
		+alpha_misch_wi_hl*A[0][15]*(T[29]-T[15])\
		-Q_str[15]\
		+eg*Ash*Cs*(Tamb**4-T[15]**4)\
		+lambda_g*th*s*(T[16]-T[15])/hz+lambda_g*th*s*(T[11]-T[15])/(0.5*(hz+b))\
		+lambda_g*hz*s*(T[19]-T[15])/(0.5*(b+th))+lambda_g*hz*s*(T[4]-T[15])/(0.5*t))/(mcp[15][0]*mcp[15][1])

	dT[16] = (alpha_frei_aussen_hl*Ash*(Tamb-T[16])\
		+alpha_misch_wi_hl*A[0][16]*(T[30]-T[16])\
		-Q_str[16]\
		+eg*Ash*Cs*(Tamb**4-T[16]**4)\
		+lambda_g*th*s*(T[15]-T[16])/hz+lambda_g*th*s*(T[17]-T[16])/hz\
		+lambda_g*hz*s*(T[20]-T[16])/(0.5*(b+th))+lambda_g*hz*s*(T[5]-T[16])/(0.5*t))/(mcp[16][0]*mcp[16][1])

	dT[17] = (alpha_frei_aussen_hl*Ash*(Tamb-T[17])\
		+alpha_misch_wi_hl*A[0][17]*(T[31]-T[17])\
		-Q_str[17]\
		+eg*Ash*Cs*(Tamb**4-T[17]**4)\
		+lambda_g*th*s*(T[16]-T[17])/hz+lambda_g*th*s*(T[18]-T[17])/(0.5*(hz+b))\
		+lambda_g*hz*s*(T[21]-T[17])/(0.5*(b+th))+lambda_g*hz*s*(T[6]-T[17])/(0.5*t))/(mcp[17][0]*mcp[17][1])



	dT[18] = (-Q_str[18]\
		+lambda_g*th*s*(T[17]-T[18])/(0.5*(hz+b))+lambda_g*th*s*(T[14]-T[18])/(0.5*(hz+b))\
		+lambda_g*b*s*(T[21]-T[18])/(0.5*(hz+th))+lambda_g*b*s*(T[7]-T[18])/(0.5*t))/(mcp[18][0]*mcp[18][1])

	dT[19] = (alpha_frei_aussen_h*Ah*(Tamb-T[19])\
		+alpha_misch_wi_h*A[0][19]*(T[29]-T[19])\
		-Q_str[19]\
		+eg*Ah*Cs*(Tamb**4-T[19]**4)\
		+lambda_g*hz*s*(T[15]-T[19])/(0.5*(th+b))+lambda_g*hz*s*(T[12]-T[19])/(0.5*(th+b))\
		+lambda_g*b*s*(T[11]-T[19])/(0.5*(hz+th))+lambda_g*b*s*(T[20]-T[19])/hz)/(mcp[19][0]*mcp[19][1])

	dT[20] = (alpha_frei_aussen_h*Ah*(Tamb-T[20])\
		+alpha_misch_wi_h*A[0][20]*(T[30]-T[20])\
		-Q_str[20]\
		+eg*Ah*Cs*(Tamb**4-T[20]**4)\
		+lambda_g*hz*s*(T[16]-T[20])/(0.5*(th+b))+lambda_g*hz*s*(T[13]-T[20])/(0.5*(th+b))\
		+lambda_g*b*s*(T[19]-T[20])/hz+lambda_g*b*s*(T[21]-T[20])/hz)/(mcp[20][0]*mcp[20][1])

	dT[21] = (alpha_frei_aussen_h*Ah*(Tamb-T[21])\
		+alpha_misch_wi_h*A[0][21]*(T[31]-T[21])\
		-Q_str[21]\
		+eg*Ah*Cs*(Tamb**4-T[21]**4)\
		+lambda_g*hz*s*(T[17]-T[21])/(0.5*(th+b))+lambda_g*hz*s*(T[14]-T[21])/(0.5*(th+b))\
		+lambda_g*b*s*(T[20]-T[21])/hz+lambda_g*b*s*(T[18]-T[21])/(0.5*(hz+th)))/(mcp[21][0]*mcp[21][1])

#	Platte dt 22, 23, 24, 25, 26, 27

		
	dT[22] = (alpha_misch_pv*A[0][22]*(T[34]-T[22])\
		-Q_str[22]\
		+lambda_p*(bp*tp/2)*(T[23]-T[22])/(0.5*(hp1+hp2))+lambda_p*(bp*hp1)*(T[25]-T[22])/(0.5*tp)\
		+QV1/2)/(mcp[22][0]*mcp[22][1])
	     

	dT[23] = (alpha_misch_pv*A[0][23]*(T[33]-T[23])\
		-Q_str[23]\
		+lambda_p*(bp*tp/2)*(T[22]-T[23])/(0.5*(hp1+hp2))+lambda_p*(bp*tp/2)*(T[24]-T[23])/(0.5*(hp2+hp3))\
		+lambda_p*(bp*hp2)*(T[26]-T[23])/(0.5*tp)\
		+QV2/2)/(mcp[23][0]*mcp[23][1])
	     

	dT[24] = (alpha_misch_pv*A[0][24]*(T[32]-T[24])\
		-Q_str[24]\
		+lambda_p*(bp*tp/2)*(T[23]-T[24])/(0.5*(hp2+hp3))+lambda_p*(bp*hp3)*(T[27]-T[24])/(0.5*tp)\
		+QV3/2)/(mcp[24][0]*mcp[24][1])

	
	dT[25] = (alpha_misch_ph*A[0][25]*(T[29]-T[25])\
		-Q_str[25]\
		+lambda_p*(bp*tp/2)*(T[26]-T[25])/(0.5*(hp1+hp2))+lambda_p*(hp1*bp)*(T[22]-T[25])/(0.5*tp)\
		+QV1/2)/(mcp[25][0]*mcp[25][1])
	     

	dT[26] = (alpha_misch_ph*A[0][26]*(T[30]-T[26])\
		-Q_str[26]\
		+lambda_p*(bp*tp/2)*(T[25]-T[26])/(0.5*(hp1+hp2))+lambda_p*(bp*tp/2)*(T[27]-T[26])/(0.5*(hp2+hp3))\
		+lambda_p*(bp*hp2)*(T[23]-T[26])/(0.5*tp)\
		+QV2/2)/(mcp[26][0]*mcp[26][1])
	     

	dT[27] = (alpha_misch_ph*A[0][27]*(T[31]-T[27])\
		-Q_str[27]\
		+lambda_p*(bp*tp/2)*(T[26]-T[27])/(0.5*(hp2+hp3))+lambda_p*(bp*hp3)*(T[24]-T[27])/(0.5*tp)\
		+QV3/2)/(mcp[27][0]*mcp[27][1])


	# T28 = Tlein, T29 = Tl0, T30=Tl1, T31=Tl2,T32=Tl3,T33=Tl4, T34=Tl5	


	if Stat_Berechnung == 1:	
		dT[28] = (m_dot*cp*(T[34]-T[28])-QKKK)/(mcp[28][0]*mcp[28][1])
		
	if Stat_Berechnung ==0:
		global Ventilator
		Q, Ventilator = bekomme_qkuehl(zeit,Automatic_on_off,T[23],Ventilator)	
		dT[28] = (m_dot*cp*(T[34]-T[28])-Q)/(mcp[28][0]*mcp[28][1])

	dT[29] = (m_dot*cp*(T[28]-T[29])\
		+alpha_misch_wi_hr*A[0][12]*(T[12]-T[29])+alpha_misch_wi_hl*A[0][15]*(T[15]-T[29])\
		+alpha_misch_wi_h*A[0][19]*(T[19]-T[29])+alpha_misch_ph*A[0][22]*(T[22]-T[29]))/(mcp[29][0]*mcp[29][1])

	dT[30] = (m_dot*cp*(T[29]-T[30])\
		+alpha_misch_wi_hr*A[0][13]*(T[13]-T[30])+alpha_misch_wi_hl*A[0][16]*(T[16]-T[30])\
		+alpha_misch_wi_h*A[0][20]*(T[20]-T[30])+alpha_misch_ph*A[0][23]*(T[23]-T[30]))/(mcp[30][0]*mcp[30][1])
		
	dT[31] = (m_dot*cp*(T[30]-T[31])\
		+alpha_misch_wi_hr*A[0][14]*(T[14]-T[31])+alpha_misch_wi_hl*A[0][17]*(T[17]-T[31])\
		+alpha_misch_wi_h*A[0][21]*(T[21]-T[31])+alpha_misch_ph*A[0][24]*(T[24]-T[31]))/(mcp[31][0]*mcp[31][1])

	dT[32] = (m_dot*cp*(T[31]-T[32])\
		+alpha_misch_wi_vr*A[0][3]*(T[3]-T[32])+alpha_misch_wi_vl*A[0][6]*(T[6]-T[32])\
		+alpha_misch_wi_v*A[0][10]*(T[10]-T[32])+alpha_misch_pv*A[0][24]*(T[24]-T[32]))/(mcp[32][0]*mcp[32][1])

	dT[33] = (m_dot*cp*(T[32]-T[33])\
		+alpha_misch_wi_vr*A[0][2]*(T[2]-T[33])+alpha_misch_wi_vl*A[0][5]*(T[5]-T[33])\
		+alpha_misch_wi_v*A[0][9]*(T[9]-T[33])+alpha_misch_pv*A[0][23]*(T[23]-T[33]))/(mcp[33][0]*mcp[33][1])

	dT[34] = (m_dot*cp*(T[33]-T[34])\
		+alpha_misch_wi_vr*A[0][1]*(T[1]-T[34])+alpha_misch_wi_vl*A[0][4]*(T[4]-T[34])\
		+alpha_misch_wi_v*A[0][8]*(T[8]-T[34])+alpha_misch_pv*A[0][22]*(T[22]-T[34]))/(mcp[34][0]*mcp[34][1])

	global t2
	t2 = time.time()
	print 'die Berechnungszeit beträgt =',t2-t1,'sekunden'


	return dT



### Randbedingungen /boundary conditions(BC)

BC = np.ones ((35)) 
for i in range(35):

	BC[i]=Tamb #anfangsttemp ist die Umgebungstemperatur
	


#################### 2. Option T_anfang= T_stationäre Berechnung

#BC = empty (32) #vector de ecuaciones
#BCdim=size(BC)
#print 'BCdim=', BCdim
##
#for i in range (0,32):
#	BC[i] = zneu[i]

#BCshape=np.shape(BC),
#print 'BCshape=',BCshape
#print 'BC=',BC


#print 'zeit complete=', zeit
#print 'zeit[0:t_ende]=',zeit[0:t_ende] ## wenn zeit[35] = zeile 36, python faengt von 0 anzuzahlen

#sol=odeint(system,BC,zeit[0:t_ende])#,full_output=1) 

#sol = solve_ivp(system, [0,340], BC,t_eval= zeit[0:t_ende])
if Stat_Berechnung ==1:
	sol= solve_ivp(system, [init_time, end_time+1], BC, method='BDF',t_eval= zeit_eval)#, events= event(t,T))
elif Stat_Berechnung == 0:
	sol= solve_ivp(system, [init_time, end_time+1], BC, method='LSODA',t_eval= zeit_eval)
#sol= solve_ivp(system, [init_time, end_time+1], BC, method='RK45')

#t_eval= zeit_eval) # fuer stationaeren Fall

zeit_solver =sol.t
Temperature = sol.y #Temperature in Kelvin
Status=sol.status

#print 'Temperature dimension=',np.shape(Temperature)
#print 'Temperature solver =',Temperature

#print 'sol.t zeit solver',zeit_solver
#print 'sol.y',Temperature
	
#print 'sol.y before transpose', temperature
#print 'dim vor transpose', np.shape(Temperature)
Temperature=Temperature.transpose() ## Temperature in kelvin
Temp_ult = Temperature[-1:,:]
#Temp_mitte = Temperature[Aux_mitte/2,:]

#temp in celsius umwandeln
Aux_mitte = len(Temperature)
#print 'dimension Temeperature',np.shape(Temperature) ## deberia ser de ndaten filas por 35 columnas
Temperature_celsius=np.zeros((np.shape(Temperature))) ### Temperature in Celsius
#print Temperature_celsius

for i in range (0,Aux_mitte):
	for j in range (0,35):
		Temperature_celsius[i][j]=Temperature[i][j]-273.15
#print Temperature
#print Temperature_celsius
#print 'dimension temperature celsius',np.shape(Temperature_celsius)
#print 'Temperature_celsius[0][8]',Temperature_celsius[0][8]
Temp_ult_celsius = Temperature_celsius[-1:,:]



#Temp_ult_celsius=np.zeros((1,35))

#for i in range(0,35):
#	Temp_ult_celsius[0][i]=Temp_ult[0][i]-273.15

# print 'Temperature dimension=',np.shape(Temperature)
#print 'Temperature solver =',Temperature

#print 'dim sol.y after transpose',np.shape(temperature)
#print 'dim nach transpose', np.shape(Temperature)





#sol=odeint(system,BC,messung[:,1])
#print 'tcur',sol[1]['tcur']
#print zeit[0:t_ende]


###########################################################3
################### Ergebnisse drücken print results
######################################################3
#print 'zeit_eval',zeit_eval

print 'alpha_misch_wi_v = ' + str(alpha_misch_wi_v)
print 'alpha_misch_wi_vr = ' + str(alpha_misch_wi_vr)
print 'alpha_misch_wi_vl = ' + str(alpha_misch_wi_vl)
print 'alpha_misch_pv = ' + str(alpha_misch_pv)
print 'alpha_misch_wi_hr = ' + str(alpha_misch_wi_hr)
print 'alpha_misch_wi_hl = ' + str(alpha_misch_wi_hl)
print 'alpha_misch_wi_h = ' + str(alpha_misch_wi_h)
print 'alpha_misch_ph = ' + str(alpha_misch_ph)

print 'alpha_frei_aussen_vr = ' + str(alpha_frei_aussen_vr)
print 'alpha_frei_aussen_vl = ' + str(alpha_frei_aussen_vl)
print 'alpha_frei_aussen_v = ' + str(alpha_frei_aussen_v)
print 'alpha_frei_aussen_hr = ' + str(alpha_frei_aussen_hr)
print 'alpha_frei_aussen_hl = ' + str(alpha_frei_aussen_hl)
print 'alpha_frei_aussen_h = ' + str(alpha_frei_aussen_h)


###########Temp_ult = sol[:,-1:] para odeint creo. no lo se
print ' temperature after all=',np.shape(Temperature)
print 'Temp_ult kelvin =',Temp_ult 
print 'Temp_ult celsius =',Temp_ult_celsius
print 'Q_str',Q_str
print 'summe_Q_str =', summ_Q_str
print 'Q_gehaeuse=',Q_gehaeuse
print 'Summe Q_gehauese', summe_Q_gehaeuse
print 'die Berechnungszeit beträgt =',t2-t1,'sekunden'
print 'sim_date=',sim_date
print 'shape Temp_ult_celsius',np.shape(Temp_ult_celsius)
print 'sol.status= ',sol.status
print "Qkuehl_average_instat_proportional = ",Qkuehl_average_instat_proportional
#print 'mplatte (kg) = ',mplatte
#print 'mgehaeuse (kg) = ',mgehaeuse
#print 'mluft (kg) = ',mluft
#print 'Qkuehl_average',Qkuehl_average
##########################################################################
##########################################################################





################################################################################################################################
################### CSV file schreiben ##########################################################################
################################################################################################################################
#import xlwt
#from xlwt import Workbook 
#wb = Workbook()
#sheet1 =wb.add_sheet('sheet1')

##ancho de la hoja
#sheet1.col(0).width= 7000

#sheet1.write(0,0,'temperatur')
#sheet1.write(0,1,'Wert')

#sheet1.write(1,1,sol[1000][0])


##for i in range (0,sol_len):
##	sheet1.write(i,0,sol[i][0])

#wb.save('Temperatur im Schaltschrank instationär.xls')



################## CSV file schreiben
#import xlwt
#DATA = (("Temperaturen im Schaltschrank",),("",),\
#	("Temperaturen im Gehaeuse",),("Tg0",zneu[0],),("Tg1",zneu[1],),("Tg2",zneu[2],),("Tg3",zneu[3],),\
#					("Tg4",zneu[4],),("Tg5",zneu[5],),("Tg6",zneu[6],),("Tg7",zneu[7],),\
#					("Tg8",zneu[8],),("Tg9",zneu[9],),("Tg10",zneu[10],),("Tg11",zneu[11],),\
#					("Tg12",zneu[12],),("Tg13",zneu[13],),("Tg14",zneu[14],),("Tg15",zneu[15],),\
#					("Tg16",zneu[16],),("Tg17",zneu[17],),("Tg18",zneu[18],),("Tg19",zneu[19],),\
#					("Tg20",zneu[20],),("Tg21",zneu[21],),\
#	("",),("Temperaturen in der Luft",),("Tlein",zneu[25],),("Tl0",zneu[26],),("Tl1",zneu[27],),("Tl2",zneu[28],),\
#					("Tl3",zneu[29],),("Tl4",zneu[30],),("Tl5",zneu[31],),\
#	("",),("Temperaturen in der Platte",),("Tp1",zneu[22],),("Tp2",zneu[23],),("Tp3",zneu[24],))



#wb = xlwt.Workbook()
#ws = wb.add_sheet("My Sheet")
#for i, row in enumerate(DATA):
#	for j, col in enumerate(row):
#		ws.write(i,j,col)
#ws.col(0).width = 256 * max([len(row[0]) for row in DATA])
#wb.save("Ergebnisse.xls")

####################################################################################################################################################
####################################################################################################################################################

#### change directory for write results
#### Ordner wechseln um Ergebnisse zu schreiben
#os.chdir(r'/home/fabian/Dokumente/studienarbeit/results_'+str(sim_date))

os.chdir(results_directory)

#results_path = os.getcwd()+'/results_'+str(sim_date) 
#my_file = 'graph.png'

####################################################################################################################################################
################WRITE RESULTS in TXT Datei
####################################################################################################################################################

#f = file('Ergebnisse_Kaelteleistung_messungen' + repr(NumOfData) + '.csv.txt', 'w')
f = file('results_'+str(figurename)+ '.txt', 'w')
f.write('Simulation Ergebnisse \n')
f.write('\n')

f.write('die Berechnungszeit beträgt = ' + repr(t2-t1) + ' sekunden\n')
print >>f, '\n'
f.write('alpha_misch_wi_v = ' + repr(alpha_misch_wi_v) + '\n') 
f.write('alpha_misch_wi_vr = ' + repr(alpha_misch_wi_vr) + '\n')
f.write('alpha_misch_wi_vl = ' + repr(alpha_misch_wi_vl)+ '\n')
f.write('alpha_misch_pv = ' + repr(alpha_misch_pv)+ '\n')
f.write('alpha_misch_wi_hr = ' + repr(alpha_misch_wi_hr)+ '\n')
f.write('alpha_misch_wi_hl = ' + repr(alpha_misch_wi_hl)+ '\n')
f.write('alpha_misch_wi_h = ' + repr(alpha_misch_wi_h)+ '\n')
f.write('alpha_misch_ph = ' + repr(alpha_misch_ph)+ '\n')
print >>f, '\n'
f.write('alpha_frei_aussen_vr = ' + repr(alpha_frei_aussen_vr)+ '\n')
f.write('alpha_frei_aussen_vl = ' + repr(alpha_frei_aussen_vl)+ '\n')
f.write('alpha_frei_aussen_v = ' + repr(alpha_frei_aussen_v)+ '\n')
f.write('alpha_frei_aussen_hr = ' + repr(alpha_frei_aussen_hr)+ '\n')
f.write('alpha_frei_aussen_hl = ' + repr(alpha_frei_aussen_hl)+ '\n')
f.write('alpha_frei_aussen_h = ' + repr(alpha_frei_aussen_h)+ '\n')
print >>f, '\n'
f.write('Dimension Temp_ult= ' + repr(np.shape(Temp_ult))+ '\n')
f.write('Temp_ult= '+ repr(Temp_ult)+ '\n')
print >>f, '\n'
f.write('Temp_ult celsius= '+ repr(Temp_ult_celsius)+ '\n')
print >>f, '\n'
f.write('Q_str= '+ repr(Q_str)+ '\n')
print >>f, '\n'
f.write('Q_gehaeuse= '+ repr(Q_gehaeuse) + '\n')
print >>f, '\n'
f.write('Summe Q_gehauese= ' + repr(summe_Q_gehaeuse) + '\n')
f.write('Qkuehl_average_instat_proportional = ' + repr(Qkuehl_average_instat_proportional) + '\n')


f.close()

### come back to calculation directory.
os.chdir(my_path)


#############################################################################################################
############################# PLOT RESULTS#############################################################################################################################################################
####################################################################################################################################################

#plot results

# HOW TO Graphics
#https://matplotlib.org/gallery/subplots_axes_and_figures/figure_title.html#sphx-glr-gallery-subplots-axes-and-figures-figure-title-py




# # a) Schaubild 3D

# fig=plt.figure(1)
# ax = fig.add_subplot(111,projection='3d')
# ax.set_xlabel('Tiefe [m]')
# ax.set_ylabel('Breite[m]')
# ax.set_zlabel('Hoehe [m]')


# # Geometrie im Koordinatensystem

# limit = max(h,t,b)

# ax.set_xlim3d(0,limit)
# ax.set_ylim3d(0,limit)
# ax.set_zlim3d(0,limit)

# #Schaltschrank
# X, Y, Z = [0,t,t,0,0],[0,0,b,b,0],[0,0,0,0,0]
# X1, Y1, Z1 = [0,t,t,0,0],[0,0,b,b,0],[h,h,h,h,h]
# X2, Y2, Z2 = [0,0],[0,0],[0,h]
# X3, Y3, Z3 = [t,t],[0,0],[0,h]
# X4, Y4, Z4 = [t,t],[b,b],[0,h]
# X5, Y5, Z5 = [0,0],[b,b],[0,h]
# ax.plot_wireframe(X,Y,Z)
# ax.plot_wireframe(X1,Y1,Z1)
# ax.plot_wireframe(X2,Y2,Z2)
# ax.plot_wireframe(X3,Y3,Z3)
# ax.plot_wireframe(X4,Y4,Z4)
# ax.plot_wireframe(X5,Y5,Z5)

# #Platte
# X6, Y6, Z6 = [xp,xp+tp,xp+tp,xp,xp],[yp,yp,yp+bp,yp+bp,yp],[zp,zp,zp,zp,zp]
# X7, Y7, Z7 = [xp,xp+tp,xp+tp,xp,xp],[yp,yp,yp+bp,yp+bp,yp],[zp+hp,zp+hp,zp+hp,zp+hp,zp+hp]
# X8, Y8, Z8 = [xp,xp],[yp,yp],[zp,zp+hp]
# X9, Y9, Z9 = [xp,xp],[yp+bp,yp+bp],[zp,zp+hp]
# X10, Y10, Z10 = [xp+tp,xp+tp],[yp,yp],[zp,zp+hp]
# X11, Y11, Z11 = [xp+tp,xp+tp],[yp+bp,yp+bp],[zp,zp+hp]
# plt.plot(X6,Y6,Z6,'r-')
# plt.plot(X7,Y7,Z7,'r-')
# plt.plot(X8,Y8,Z8,'r-')
# plt.plot(X9,Y9,Z9,'r-')
# plt.plot(X10,Y10,Z10,'r-')
# plt.plot(X11,Y11,Z11,'r-')

# #Zonen
# plt.plot([0,t,t,0,0],[0,0,b,b,0],[hz,hz,hz,hz,hz],'k--')
# plt.plot([0,t,t,0,0],[0,0,b,b,0],[hz*2,hz*2,hz*2,hz*2,hz*2],'k--')


# #Temperaturen
# #Platte
# plt.plot([xp+tp/2,xp+tp/2,xp+tp/2],[yp+bp/2,yp+bp/2,yp+bp/2],[zp+hp3/2,3*hz/2,2*hz+hp1/2],'ro')

# #Luft
# plt.plot([s+tih/2,s+tih/2,s+tih/2,s+tih/2],[b/2,b/2,b/2,b/2],[hz/2,3*hz/2,5*hz/2,h-s/2],'bo')
# plt.plot([xp+tp+tiv/2,xp+tp+tiv/2,xp+tp+tiv/2],[b/2,b/2,b/2],[hz/2,3*hz/2,5*hz/2],'co')


# #Runden der Temperaturwerte auf die zweite Nachkommastelle
# Temp_ult=Temp_ult.transpose()

# for i in range(0,35):
	# Temp_ult[i]=Temp_ult[i].round(2)

# #Eintragen der Temperaturwerte in die Grafik
# ax.text(s+tih/2,b/2+0.1,h-s/2+0.1, "Tlein = " + str(Temp_ult[28]) + 'K', color='blue')
# ax.text(s+tih/2,b/2+0.1,5*hz/2+0.1, "Tl0 = " + str(Temp_ult[29]) + 'K', color='blue')
# ax.text(s+tih/2,b/2+0.1,3*hz/2+0.1, "Tl1 = " + str(Temp_ult[30]) + 'K', color='blue')
# ax.text(s+tih/2,b/2+0.1,hz/2+0.1, "Tl2 = " + str(Temp_ult[31]) + 'K', color='blue')

# ax.text(xp+tp+tiv/2,b/2+0.1,5*hz/2-0.1, "Tl5 = " + str(Temp_ult[34]) + 'K', color='green')
# ax.text(xp+tp+tiv/2,b/2+0.1,3*hz/2-0.1, "Tl4 = " + str(Temp_ult[33]) + 'K', color='green')
# ax.text(xp+tp+tiv/2,b/2+0.1,hz/2-0.1, "Tl3 = " + str(Temp_ult[32]) + 'K', color='green')

# ax.text(xp+tp/2,yp+bp/2+0.1,2*hz+hp1/2, "T22 = " + str(Temp_ult[22]) + 'K', color='red')
# ax.text(xp+tp/2,yp+bp/2+0.1,3*hz/2, "T23 = " + str(Temp_ult[23]) + 'K', color='red')
# ax.text(xp+tp/2,yp+bp/2+0.1,zp+hp3/2, "T24 = " + str(Temp_ult[24]) + 'K', color='red')




###### Temperatur des Gehaeuses aus der Messung bei höhe von 1630 mm
#### Berechnung liegt bei 1666 mm

Zeit_min_g= Gehaeuse_Temperatur[:,0]
Zeit_sek_g= Gehaeuse_Temperatur[:,0]*60.
LSC_37=Gehaeuse_Temperatur[:,5]
print "Zeit_sek_g",Zeit_sek_g




# b)  Schaubild Temperaturar entwicklung im Laufe der Zeit

#plt.plot(Zeit_sek,LSC_TL_134,'.',label='Theta_Amb [C]') 
#plt.plot(Zeit_sek,Tmessung_celsius,label='Theta_Messung [C]')


plt.figure(2)
plt.plot(zeit_solver,Temperature_celsius[:,19],color = 'red',label=' Gehaeuse solver') 

if Stat_Berechnung==0:
	plt.plot(Zeit_sek_g,LSC_37,label='Gehaeuse Messung')
	
plt.title('Vergleich der Gehaeusetemperatur')
plt.xlabel('Zeit [Sek]')		
plt.ylabel('Temperature [C]')
my_file2 = 'results_Validierung_gehaeuse.png'
plt.savefig(os.path.join(results_directory, my_file2), dpi=200, bbox_inches='tight')
plt.legend(loc='lower right')
plt.grid(True)


#plt.plot(zeit_solver,Temperature_celsius[:,19]) #platte Tstat =325

#fontsize = 22
#markersize = 15
#labelpad = 10
#linewidth = 2

#plt.title('Zeitliche Verlauf der T[29]')
#plt.xlabel('Zeit (sek)')		
#plt.ylabel('Temperature (C)')

#my_file2 = 'results_T28.png'
#plt.savefig(os.path.join(results_directory, my_file2), dpi=200, bbox_inches='tight')
#plt.grid(True)


plt.figure(3)
plt.plot(zeit_solver,Temperature_celsius[:,34]) #platte Tstat =325

fontsize = 22
markersize = 15
labelpad = 10
linewidth = 2

plt.title('Zeitliche Verlauf der T[34] T_Ausgang_luft')
plt.xlabel('Zeit (sek)')		
plt.ylabel('Temperature (C)')

#name = '/home/franco/Dokumente/Studienarbeit/GPM-'+str(n2)+'.pdf'\n",
# plt.savefig(str(figurename) + '.pdf', dpi=200, bbox_inches='tight')
my_file3 = 'results_T34.png'
plt.savefig(os.path.join(results_directory, my_file3), dpi=200, bbox_inches='tight')
plt.grid(True)


plt.figure(4)
plt.plot(zeit_solver,Temperature_celsius[:,23]) #platte Tstat =325

fontsize = 22
markersize = 15
labelpad = 10
linewidth = 2

plt.title('Zeitliche Verlauf der T[23] Platte')
plt.xlabel('Zeit (sek)')		
plt.ylabel('Temperature (C)')

#name = '/home/franco/Dokumente/Studienarbeit/GPM-'+str(n2)+'.pdf'\n",
# plt.savefig(str(figurename) + '.pdf', dpi=200, bbox_inches='tight')
my_file4 = 'results_T23.png'
plt.savefig(os.path.join(results_directory, my_file4), dpi=200, bbox_inches='tight')
plt.grid(True)

##### Durchschnittwerte



################# Messwerte plotten


#numzeile_csv=696
#numspalten_csv=28

#def Openluft():
#	with open('Lufttemperatur_500W.csv') as csvfile:
#		reader = csv.reader(csvfile)
#		reader.next()
#		
#		# numzeile_csv=696
#		# numspalten_csv=28
#		aux_variable=numzeile_csv-1
#		#Matrize erstellen mit der Anzahl der elemente als Dimension.
#		
#		Matrize_Lufttemperatur = np.zeros((aux_variable,numspalten_csv)) #weil die Tabelle in dem Dokument 37 zeilen hat. Python fängt von 0 an zu zahlen.

#		i=0
#		for row in reader:
#			Matrize_Lufttemperaturr = np.genfromtxt(row)
#			Matrize_Lufttemperatur[i] = Matrize_Lufttemperaturr
#			i=i+1
#        
#		for i in range (0,aux_variable):
#			for j in range(0,numspalten_csv):
#				Matrize_Lufttemperatur[i][j]=float(Matrize_Lufttemperatur[i][j])
#	return Matrize_Lufttemperatur

Matrize_Lufttemperatur=Openluft()
	
# print 'Matrize_Lufttemperatur',Matrize_Lufttemperatur
# print 'dimension Matrize_Lufttemperatur',np.shape(Matrize_Lufttemperatur)


Zeit_min= Matrize_Lufttemperatur[:,0]
Zeit_sek= Matrize_Lufttemperatur[:,0]*60.
#print 'Zeit_min',Zeit_min
#print 'zeit_sek',Zeit_sek

#1. Grafik
LSC_TL_134= Matrize_Lufttemperatur[:,1] 
LSC_TL_136= Matrize_Lufttemperatur[:,2]
LSC_TL_137= Matrize_Lufttemperatur[:,3]

#2. Grafik
LSC_TL_115= Matrize_Lufttemperatur[:,5]
LSC_TL_116= Matrize_Lufttemperatur[:,6]
LSC_TL_117= Matrize_Lufttemperatur[:,7]
LSC_TL_118= Matrize_Lufttemperatur[:,8]
LSC_TL_119= Matrize_Lufttemperatur[:,9]

#aux_variable2 = len(LSC_TL_115)
#for i in range (0,aux_variable2):
#	LSC_TL_115[i]=LSC_TL_115[i]+273.15


##Durchschnitt der Messwerte

LSC_TL_134_Durchschnitt= np.average(LSC_TL_134)
LSC_TL_136_Durchschnitt= np.average(LSC_TL_136)
LSC_TL_137_Durchschnitt= np.average(LSC_TL_137)

LSC_TL_115_Durchschnitt= np.average(LSC_TL_115)
LSC_TL_116_Durchschnitt= np.average(LSC_TL_116)
LSC_TL_117_Durchschnitt= np.average(LSC_TL_117)
LSC_TL_118_Durchschnitt= np.average(LSC_TL_118)
LSC_TL_119_Durchschnitt= np.average(LSC_TL_119)


###### rot = Messung
###### blau = Simulation
#################### EINS rot #############################
# Temp Durchschnitt berechnet aus Messung
Temp_durchschnitt=[LSC_TL_115_Durchschnitt,LSC_TL_116_Durchschnitt,LSC_TL_117_Durchschnitt,LSC_TL_118_Durchschnitt,LSC_TL_119_Durchschnitt] # Temp Durchschnitt berechnet aus Messung

## hoehe der Messgeräte 
height = [500, 800, 1100, 1400, 1700] #Position der Messfüller


############################ ZWEI blau ############################

### Temperature Durchschnitt aus Simulation
T34_Durchschnitt=np.average(Temperature[:,34])
print 'T34_Durchschnitt',T34_Durchschnitt
T4_Durchschnitt= np.average(Temperature[:,4])
T5_Durchschnitt=np.average(Temperature[:,5])
T6_Durchschnitt=np.average(Temperature[:,6])


######### transformation von K zu °C
T4_Durchschnitt=T4_Durchschnitt-273.15
T5_Durchschnitt=T5_Durchschnitt-273.15
T6_Durchschnitt=T6_Durchschnitt-273.15

T_durchschnitt_Berechnung= [T4_Durchschnitt,T5_Durchschnitt,T6_Durchschnitt]
#print 'T_durchschnitt_Berechnung',T_durchschnitt_Berechnung
#print 'type T_durchschnitt_Berechnung', type(T_durchschnitt_Berechnung)
#for i in range (3)
#T_durchschnitt_Berechnung [1][i]=T_Durchschnitt
# T_durchschnitt_Berechnung = np.zeros(1,3) 
#for i in range (3):
#	T_durchschnitt_Berechnung[1][i]=T_durchschnitt_Berechnung[1][i]-273.15

###hoehe der Temperaturknoten

h=h*1000 #h in mm hoehe des Schrankes
height_Berechnung = [h/6,h/6*3,h/6*5]
#print ' type height_Berechnung',type(height_Berechnung)




######################## DREI gruen MESSWERTE aus der MEssung, nicht berechnet #############

height = [500, 800, 1100, 1400, 1700] #Position der Messfüller
Temp_durchschnitt_messung=[32.42,32.26,32.45,35.08,34.17] ## Temp durchschnitt aus Messung
#print 'type height',type(height)
#print 'type Temp_durchschnitt_messung',type(Temp_durchschnitt_messung)





	
# T_durchschnitt_Berechnung=T_durchschnitt_Berechnung.round(2)

#print 'height_Berechnung',height_Berechnung
#print 'T_durchschnitt_Berechnung',T_durchschnitt_Berechnung

#############################################################################################

plt.figure(5)
plt.plot(zeit_solver,Temperature_celsius[:,34],label='solver') 
if Stat_Berechnung==0:
	plt.plot(Zeit_sek,LSC_TL_115,label='Messung')
#plt.plot(height,Temp_durchschnitt_messung,'g')#,Zeit_min,LSC_TL_136,'r',Zeit_min,LSC_TL_137,'b') #platte Tstat =325

plt.title('Vergleich der Lufttemperatur Ausgangbereich')
plt.xlabel('zeit [Sek]')		
plt.ylabel('Temperature [C]')
my_file5 = 'results_Validierung.png'
plt.savefig(os.path.join(results_directory, my_file5), dpi=200, bbox_inches='tight')
plt.legend(loc='lower right')
plt.grid(True)




#plt.figure(4)
##plt.subplot(2, 1, 1)
#plt.plot(Zeit_m,LSC_TL_134,'g',Zeit_min,LSC_TL_136,'r',Zeit_min,LSC_TL_137,'b') #platte Tstat =325
## plt.show
#plt.title('Messung')
#plt.xlabel('Zeit [sek]')		
#plt.ylabel('Gemessene Temperature [K]')
#plt.savefig('Figure4.pdf', dpi=200, bbox_inches='tight')

#plt.figure(5)
##plt.subplot(2, 1, 2)
#plt.plot(Zeit_min,LSC_TL_115,Zeit_min,LSC_TL_116,Zeit_min,LSC_TL_117,Zeit_min,LSC_TL_118,Zeit_min,LSC_TL_119)
#plt.title('Messung')
#plt.xlabel('Zeit [sek]')		
#plt.ylabel('Gemessene Temperature [K]')
#plt.savefig('Figure5.pdf', dpi=200, bbox_inches='tight')














##################### Resulst in excel Datei schreiben
Temp_ult_celsius=Temp_ult_celsius.transpose()
#print 'type Temp_ult_celsius',type(Temp_ult_celsius)
print 'shape Temp_ult_celsius',np.shape(Temp_ult_celsius)







#https://biokamikazi.wordpress.com/2016/07/07/export-numpy-array-to-excel-in-python/
#Nudo=np.linspace(1.,35.,36)
Nudo = np.zeros((35,1))

for i in range (0,35):
	Nudo[i][0]=Nudo[i-1][0]+1.
#print 'Nudo',Nudo
#print 'shape Nudo',np.shape(Nudo)




###  add rows or columns to an existing array.
####https://code.i-harness.com/es/q/817d96
#results_matrix=np.concatenate((Nudo,Temp_ult_celsius,masse), axis=1)
#results_matrix=np.reshape(results_matrix,(2,-1))
#results_matrix=results_matrix.transpose()
#print results_matrix
#print 'shape results_matrix', np.shape(results_matrix)
#import xlsxwriter 

#workbook = xlsxwriter.Workbook('Results01.xlsx')
#worksheet = workbook.add_worksheet()
#row=0
#col=0

### Iterate over the data and write it out row by row.

#for Nudo,T_knoten in (results_matrix):
#	worksheet.write(row,col,results_matrix)
#	worksheet.write(row,col+1,results_matrix)
#	row +=1

#workbook.close()




#https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html

#https://stackoverflow.com/questions/41870968/dumping-numpy-array-into-an-excel-file

#import pandas as pd
#d = {'col1': Nudo, 'col2': Nudo}
#df = pd.DataFrame(data=d)

### convert your array into a dataframe

##df = pd.DataFrame (Temp_ult_celsius)

### save to xlsx file

#filepath = 'my_excel_file.xlsx'

#df.to_excel(filepath, index=False)

################################# Results in excel Datei schreiben



### Transformieren von narray zu float, vllt sind nur die einzelwerte schon float nach der Operation.
#for i in range (0,len(Temp_ult_celsius)):
#		for j in range(0,35):
#			Temp_ult_celsius[i][j]=float(Temp_ult_celsius[i][j])

#Temp_ult_celsius=Temp_ult_celsius.transpose()
##Temp_ult_celsius_float=np.array(all_data,dtype=float)
#print 'Temp_ult_celsius',Temp_ult_celsius
#print 'type Temp_ult_celsius',type(Temp_ult_celsius)

#import xlwt
#DATA = (("Temperaturen im Schaltschrank",),("",),\
#	("Temperaturen im Gehaeuse",),("Tg0",Temp_ult_celsius[0],),("Tg1",Temp_ult_celsius[1],),("Tg2",Temp_ult_celsius[2],),("Tg3",Temp_ult_celsius[3],),\
#					("Tg4",Temp_ult_celsius[4],),("Tg5",Temp_ult_celsius[5],),("Tg6",Temp_ult_celsius[6],),("Tg7",Temp_ult_celsius[7],),\
#					("Tg8",Temp_ult_celsius[8],),("Tg9",Temp_ult_celsius[9],),("Tg10",Temp_ult_celsius[10],),("Tg11",Temp_ult_celsius[11],),\
#					("Tg12",Temp_ult_celsius[12],),("Tg13",Temp_ult_celsius[13],),("Tg14",Temp_ult_celsius[14],),("Tg15",Temp_ult_celsius[15],),\
#					("Tg16",Temp_ult_celsius[16],),("Tg17",Temp_ult_celsius[17],),("Tg18",Temp_ult_celsius[18],),("Tg19",Temp_ult_celsius[19],),\
#					("Tg20",Temp_ult_celsius[20],),("Tg21",Temp_ult_celsius[21],),\
#	("",),("Temp_ult_celsiusn in der Luft",),("Tlein",Temp_ult_celsius[28],),("Tl0",Temp_ult_celsius[29],),("Tl1",Temp_ult_celsius[30],),("Tl2",Temp_ult_celsius[31],),\
#					("Tl3",Temp_ult_celsius[32],),("Tl4",Temp_ult_celsius[33],),("Tl5",Temp_ult_celsius[34],),\
#	("",),("Temp_ult_celsiusn in der Platte",),("Tp1",Temp_ult_celsius[22],),("Tp2",Temp_ult_celsius[23],),("Tp3",Temp_ult_celsius[24],))

#print'prueba shape DATA',np.shape(DATA)

#wb = xlwt.Workbook()
#ws = wb.add_sheet("My Sheet")
#for i, row in enumerate(DATA):
#	for j, col in enumerate(row):
#		ws.write(i,j,col)
#ws.col(0).width = 256 * max([len(row[0]) for row in DATA])
#wb.save("Ergebnisse.xls")



plt.show()

