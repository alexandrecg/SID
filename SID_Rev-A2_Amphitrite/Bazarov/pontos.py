import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
import os

#############################################
########## Início - Input de dados ##########
#############################################
     
passo = 1.0

#lado do hexágono de cada proveta
a = 13.2
b = a*np.sin(120*np.pi/180)/np.sin(30*np.pi/180)

#Arquivo de dados do teste a ser aberto
arquivo = "III-L03 C-00.00-12 (Ox)"

#############################################
########## Fim - Input de dados #############
#############################################

#Número das provetas
ID = np.arange(1,61+1,1)
'''
#Criando o vetor com as massas medidas de cada proveta
mass = []

##########
#####Adicionando valores randômicos de massa para as provetas para teste do código
for i in ID:
    mass.append(np.sin((i/61)*np.pi)*50)
#####
'''


print("\nUtilizando dados do arquivo " + arquivo + ".txt"+"\n")

#file=open(os.path.join(os.path.dirname("ponto.py"),os.pardir)+"\\2-Dados\\"+arquivo,"r")
file=open(os.path.join(os.path.dirname("ponto.py"),os.pardir)+"/2-Dados/"+arquivo+".txt","r")

mass = []

for i in ID:
    mass.append(float(file.readline()))

file.close()    

##########
    
#Lendo arquivo .txt contendo valores de x de posição das provetas
file=open("proveta_x.txt","r")

xp = []

for i in ID:
    xp.append(float(file.readline()))

file.close()

#Lendo arquivo .txt contendo valores de y de posição das provetas
file=open("proveta_y.txt","r")

yp = []

for i in ID:
    yp.append(float(file.readline()))
    
file.close()


zp = []

for m in mass:
    zp.append(m/max(mass))

#Criando as posições x e y da plotagem

xi = np.arange(min(xp)-a,max(xp)+a,passo)

yi = np.arange(min(yp)-b,max(xp)+b,passo)

xi, yi = np.meshgrid(xi,yi)

zi = inter.griddata((xp,yp),zp,(xi,yi),method='cubic',fill_value=0)


for i in range(zi.shape[0]):
    for j in range(zi.shape[1]):
        if(zi[i][j] < 0):
            zi[i][j] = 0
        elif(zi[i][j] >1):
            zi[i][j] = 1


#fig1 = plt.figure()
#plt.plot(xp,yp,'ro')


MaZ=[]

for i in range(len(zi[0])):
    MaZ.append(max(zi[:][i]))

print("Measured maximum value:",max(zp))
print("Interpolation maximum value:",max(MaZ))

#print("Maior valor na interpolação: ",maxvalue)

fig2 = plt.figure()
ax = fig2.add_subplot(111)
#plt.contourf(xi,yi,zi,np.arange(0,maxvalue,0.01),cmap="coolwarm")
plt.contour(xi,yi,zi,np.arange(0,1+0.2,0.2),cmap="coolwarm")
ax.plot((0,0),"ro")
plt.colorbar()
#plt.plot(xp,yp,'k.')
plt.xlabel('x [mm]',fontsize=16)
plt.ylabel('y [mm]',fontsize=16)
plt.savefig('Output/'+arquivo+'_Curvas_de_Nível.png', dpi= 300, bbox_inches = 'tight')
plt.show()

fig3 = plt.figure()
ax = fig3.add_subplot(111)
plt.contourf(xi,yi,zi,np.arange(0,1+0.05,0.01),cmap="coolwarm")
ax.plot((0,0),"ro")
plt.colorbar()
plt.xlabel('x [mm]',fontsize=16)
plt.ylabel('y [mm]',fontsize=16)
plt.savefig('Output/'+arquivo+'_Colorido.png', dpi= 300, bbox_inches = 'tight')
plt.show()


fig4 = plt.figure()
ax = fig4.add_subplot(111)
plt.contourf(xi,yi,zi,np.arange(0,1+0.05,0.01),cmap="coolwarm")
ax.plot((0,0),"ro")
plt.colorbar()
plt.plot(xp,yp,'k.')

for i in range(61):
        if(zp[i] > 0.5/100): plt.annotate("%.0f %%"%(zp[i]*100),(xp[i],yp[i]))
        
plt.xlabel('x [mm]',fontsize=16)
plt.ylabel('y [mm]',fontsize=16)
plt.savefig('Output/'+arquivo+'_Colorido_com_valores.png', dpi= 300, bbox_inches = 'tight')
plt.show()


##############
##############

view_angle1 = 45
view_angle2 = 0

##############
##############

fig5 = plt.figure()
ax = Axes3D(fig5)#fig5.gca(projection='3d')
ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap="coolwarm", antialiased= True, linewidth=0, vmin=0, vmax=1, alpha = None)
plt.xlabel('x [mm]',fontsize=16)
plt.ylabel('y [mm]',fontsize=16)
#ax.scatter(xp,yp,zp, color="black")
ax.view_init(view_angle1,view_angle2)
plt.savefig('Output/'+arquivo+'_Superfície('+str(view_angle1)+';'+str(view_angle2)+').png', dpi= 300, bbox_inches = 'tight')





