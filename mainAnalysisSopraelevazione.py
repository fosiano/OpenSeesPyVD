
# COMPLETARE SMORZAMENTO RAYLEIGH ok
# COMPLETARE TIMOSHENKO (assegnare area a taglio nel main)
# CALCOLARE BILANCIO ENERGETICO
# MIGLIORARE FASE DI ANALISI (EVITARE SALVATAGGIO SU DISCO?)


import openseespy.opensees as ops
#import openseespy.postprocessing.ops_vis as opsv
#import BraineryWiz as bz
import os

# import vfo rendering module
import vfo.vfo as vfo

import opsvis as opsv
from GetModel import (GetNodes, FixNodes,rgdDphrgm, GetMasses, getBeams,getTimoshenkoBeams,getDampers,setAnalysis, setTimeHistoryAnalysis, CQC)
import copy
#from math import sqrt
#from openseespy.postprocessing.Get_Rendering import *

import openseespy.postprocessing.Get_Rendering as opsplt

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


#PARAMETRI DI VARIABILITA DISSIPATORI
#telaioNUDO=True
telaioNUDO=False
alfaK=1.0
alfaC=0.1

#******************************************************************************
#UNITA' DI MUSURA
#******************************************************************************
N=1.0
mm=1.0
Mpa=N/mm**2
sec=1.0
kN=1000.
metri=1000.
g=9.810*metri/sec**2 #mm/s2



#******************************************************************************
#INPUT DATI
#******************************************************************************

currentpath=os.getcwd()
Input=currentpath+"\\Input\\"

recordTime=Input+'SpectralPeriod.txt'
recordAcc=Input+'SpectralAcceleration.txt'

THAccfilename="AccTH-MSS-SLV"
nfileRecordAccel=1
record_dt=0.01
timeSubdivision=1 #numero di suddivisioni dello step dell'accelergramma per l'analisi
damp = 0.05
nfreq = 2

# Lx, Ly, H [m]
Lx = [4.*metri, 5.*metri, 6.*metri, 7.*metri, 8.*metri] 
Ly = [4.*metri, 4.*metri, 4.*metri, 4.*metri, 4.*metri] 
H =  [3.*metri, 3.*metri, 3.*metri, 3.*metri, 3.*metri]

teta=[]
teta.append(np.arctan(Lx[0]/H[0]))
teta.append(np.arctan(Lx[0]/H[1]))


nImp=2
nCx=3
nCy=1

x0=-sum(Lx[:nCx])/2
X=[x0] #ascisse
for i in range(0, nCx):
    X.append(X[i]+Lx[i])

y0=-sum(Ly[:nCy])/2
Y=[y0] #ordinate
for i in range(0, nCy):
    Y.append(Y[i]+Ly[i])

z0=0.0*metri
Z=[z0] #quote
for i in range(0, nImp):
    Z.append(Z[i]+H[i])


#FISSO CARATTERISTICHE MATERIALE ELEMENTI
Ec =  32308.0*Mpa #calcestruzzo
Gc =  13461.7*Mpa
gc =  25.00*kN/metri**3
#Es = 210000.#acciaio. N/mm2

Tn1 = 0.7; # sec (Natural Period)
# mass
W1 = 1000.*kN; #KN
mass1 = W1/g 

#DAMPERS
ad=0.35
C=[]
C.append(20.7452*kN/((mm/sec)**ad)*alfaC)
C.append(20.7452/20*kN/((mm/sec)**ad)*alfaC)




# Columns and Beam Properties
K1 =  ((2*np.pi/Tn1)**2)*mass1; 
I1=K1*(H[0]**3)/((nCx+1)*(nCy+1)*12*Ec) # mm4 (K=24EIc/h^3)
Ic=[I1,0.]
Ib = 1e12*I1  # mm4
Acb =  1e16 # mm2.

arrTn2=np.arange(0.72, 0.84, 0.02)
arrTn2=np.linspace(0.72, 0.82, 6)

#listTn2=[0.700+0.0200*i for i in range(0,12)]
listTn2=[0.72, 0.74, 0.76,0.77 ,0.78,0.79, 0.80,0.81, 0.82]
listTn2=[0.77]
U=[]
for Tn2 in listTn2:
    #Tn2 = 0.75; # sec (Natural Period)
    # mass
    W2 = W1/10.; #KN
    mass2 = W2/g 
    
    # Columns and Beam Properties
    K2 =  ((2*np.pi/Tn2)**2)*mass2; 
    I2 = K2*(H[1]**3)/((nCx+1)*(nCy+1)*12*Ec) # mm4 (K=24EIc/h^3)
    Ic[1]=I2
    
    # Create nodes
    P0={'x':0.0, 'y':0.0, 'z':0.0}
    vincolo={'ux':0, 'uy':0, 'uz':0,'rx':0, 'ry':0, 'rz':0}
    
    pos=[]
    cod=[]
    
    k=1
    
    z=Z[0]
    
    n=1
    for y in Y:
        for x in X:       
            cod.append(n); P0['x']=x; P0['y']=y; P0['z']=z; pos.append(copy.deepcopy(P0))
            n+=1
    
    master=[]
    slave = []
    m=100
    for z in Z[1:]:    
        cod.append(m);P0['x']=0.0; P0['y']=0.0; P0['z']=z; pos.append(copy.deepcopy(P0))

        n=m
        floorSlave=[]
        for y in Y:
            for x in X:
                n+=1
                cod.append(n); P0['x']=x; P0['y']=y; P0['z']=z; pos.append(copy.deepcopy(P0))
                floorSlave.append(n)
        slave.append(floorSlave)            
        master.append(m)    
        m+=100
    
    fixednodes=[]
    fix=[]
    
    vincolo['ux']=1;vincolo['uy']=1;vincolo['uz']=1
    vincolo['rx']=1;vincolo['ry']=1;vincolo['rz']=1
    
    n=0
    for j in range(0,nCy+1):
        for i in range(0,nCx+1):
            n+=1
            fixednodes.append(n)
            fix.append(copy.deepcopy(vincolo))

    
    vincolo['ux']=0;vincolo['uy']=0;vincolo['uz']=1
    vincolo['rx']=1;vincolo['ry']=1;vincolo['rz']=0
    for j in range(1,nImp+1):
        fixednodes.append(100*j)
        fix.append(copy.deepcopy(vincolo))
    
    
    #creo le masse
    lmassX1 = mass1
    lmassX2 = mass2
    
    lmassY = 0.
    lmassZ = 0.
    Masses=[]
    #Memorizzo le masse dei nodi master dei tre impalcati 
    Masses.append([100, lmassX1, lmassY, lmassZ, 0.00001, 0.00001, 0.00001])
    Masses.append([200, lmassX2, lmassY, lmassZ, 0.00001, 0.00001, 0.00001])
    
    
    #definizione elementi beam
    
    beamCod=[]
    Beams=[]
    
    ele={"i":0, "j":0, "A":0, "E":0, "G":0, "J":0, "Iy":0,"Iz":0, "Avy":0,"Avz":0,"gTT":0} 
            
    
    gTTagz = 1
    gTTagx = 2
    gTTagy = 3
    
    #ELEMENTI PILASTRI
    #b=bc #Wz
    #h=hc #Wy
    A,  Iy = Acb, Ib
    
    #b1=min(b,h)
    #a1=max(b,h)
    J =  Ib*10.#a1*b1**3/3.*(1-0.577*b1/a1)
    
    # Define column elements
    
    for n in range(1,nImp+1):
        q=(n-1)*100
        for p in range(0,nCy+1):
            for m in range(0,nCx+1):
                q+=1
                beamCod.append(q)
                Iz=Ic[n-1]
                ele["i"]=q; ele["j"]=q+100; ele["A"]=A; ele["E"]=Ec; ele["G"]=Gc;
                ele["J"]=J; ele["Iy"]=Iy; ele["Iz"]=Iz; ele["gTT"]=gTTagz;
                ele["Avy"]=5.0/6.0*A;ele["Avz"]=5.0/6.0*A;
                Beams.append(copy.deepcopy(ele))
        
    
    #ELEMENTI TRAVI
    #b=bb #Wz
    #h=hb #Wy
    A, Iz, Iy = Acb, Ib, Ib
    J = Ib*10.# h*b**3/3.*(1-0.63*h/b)#    
    
    # Define X-beam elment
    #p=1
    
    for n in range(1,nImp+1):
        q=n*100
        for p in range(0,nCy+1):
            q+=1
            for m in range(0,nCx):                
                i=q
                j=q+1
                bmCd=int(str(i)+str(j))
                beamCod.append(bmCd)
                ele["i"]=i; ele["j"]=i+1; ele["A"]=A; ele["E"]=Ec; ele["G"]=Gc;
                ele["J"]=J; ele["Iy"]=Iy; ele["Iz"]=Iz; ele["gTT"]=gTTagx;
                ele["Avy"]=5.0/6.0*A;ele["Avz"]=5.0/6.0*A;
                Beams.append(copy.deepcopy(ele)) 
                q+=1

        
    #Damper Properties
    
    
    
    D=[
    #  Kd                  Cd   ad
      [25.0*kN/mm*alfaK,  C[0], ad],
      #
      [25.0*kN/mm*alfaK,  C[1], ad]
      ]
    
    
    # m = 1.0
    # arr = (np.array(D))*m
    # D = arr.tolist()
    
    DamperCode=[]
    Dampers=[]
    
    damper={"i":0, "j":0, "Kd":25., "Cd":20.7452, "ad":0.35} 
    
    # Define dampers elements
    for n in range(1,nImp+1):
        d=n*1000
        for p in range(0,nCy+1):
            d+=1
            DamperCode.append(d)
            i=p*(nCx+1)+1+(n-1)*100
            j=p*(nCx+1)+2+n*100
            damper["i"]=i; damper["j"]=j; #damper["matTag"]=1
            damper["Kd"]=D[n-1][0]; damper["Cd"]=D[n-1][1]; damper["ad"]=D[n-1][2]   
            Dampers.append(copy.deepcopy(damper)) 
        
    
    #******************************************************************************
    

    #******************************************************************************
    #GENERAZIONE DEL MODELLO
    #******************************************************************************
    ops.wipe()
    # Create ModelBuilder (with 3+-dimensions and 6 DOF/node)
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    
    #Crea nodi
    GetNodes(cod,pos)
    
    #Vincola nodi
    FixNodes(fixednodes,fix)
    
    #Get rigid floor constraint
    rgdDphrgm(master,slave)
    
    #Assegnazione masse ai nodi master dei tre impalcati
    GetMasses(Masses)
    
    #Generazione pilastri e travi
    
    coordTransf = 'Linear'
    ops.geomTransf(coordTransf, gTTagz, 0., -1., 0.)
    ops.geomTransf(coordTransf, gTTagx, 0., -1., 0.)
    ops.geomTransf(coordTransf, gTTagy, 1.,  0., 0.)
    
    # getBeams(beamCod,Beams)
    
    getTimoshenkoBeams(beamCod,Beams)
    
    if not telaioNUDO:
        getDampers(DamperCode,Dampers)
    
    #******************************************************************************
    
    # ------------------------------
    # End of model generation
    # ------------------------------

    opsplt.plot_model("nodes","elements","local_axes")
    opsv.plot_model()

      
    # render the model after defining all the nodes and elements
    #vfo.plot_model(show_nodes='yes')
    
    #
    # Set time series for the response spectrum function
    # This command is used to construct a TimeSeries object which represents 
    # the relationship between the time in the domain, t,and the load factor 
    # applied to the loads , 𝜇, in the load pattern with which the TimeSeries 
    # object is associated, i.e. 𝜇 = 𝐺(t).
    ####timeSeries(tsType, tsTag, *tsArgs)
    ops.timeSeries('Path', 1, '-filePath', recordAcc, '-fileTime', recordTime, '-factor', g)
    
    setAnalysis()
    
    
    #define a recorder for the (use a higher precision otherwise the results
    # won't match with those obtained from eleResponse)
    filename = 'ele_1_sec_i.txt'
    
    #ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 1, 'section', '1', 'force')
    ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 1, 'forces')
    # perform the analysis
    #ops.analyze(1)
    
    
    #eigValues = ops.eigen('-fullGenLapack',nfreq)
    eigValues = ops.eigen("-genBandArpack",nfreq)
    #eigs = eigen("-genBandArpack", 7)
    print('MODAL Analysis DONE!')
    print()
    print('Spectral Analysis just STARTED!')
    # compute the modal properties
    ops.modalProperties("-print", "-file", "ModalReport.txt", "-unorm")
        
    # some settings for the response spectrum analysis
    tsTag = 1 # use the timeSeries 1 as response spectrum function
    direction = 1 # excited DOF = Ux
    
    # currently we use same damping for each mode
    dmp = [damp]*len(eigValues)
    
    # we don't want to scale some modes...
    scalf = [1.0]*len(eigValues)
    # The scale factor (-scale $scale) for the output modal displacements is not use
    # now, and it’s there for future implementations. When your model is linear 
    # elastic there is no need to use this option. It will be useful in future 
    # when we will allow using this command on a nonlinear model as a linear
    # perturbation about a certain nonlinear state. In that case, the scale factor 
    # can be set to a very small number, so that the computed modal displacements 
    # will be very small (linear perturbation) and will not alter the nonlinear 
    # state of your model. Then the inverse of the scale factor can be used to 
    # post-multiply any result for post-processing.
       
    #ps.responseSpectrum $tsTag $direction <-scale $scale> <-mode $mode>
    #ops.responseSpectrum(tsTag, direction)
        
    # ========================================================================
    # run a response spectrum analysis mode-by-mode.
    # grab results during the loop, not using the recorder
    # then do modal combination in post-processing.
    # ========================================================================
    #remove('recorder', 0)
    My = []
    U200x=[]
    U100x=[]
    for i in range(len(eigValues)):
        #ops.responseSpectrum(tsTag, direction, '-mode', i+1)#3.4.0.2
        
        ops.responseSpectrumAnalysis(tsTag, direction, '-mode', i+1)#3.4.0.5
        force = ops.eleResponse(1, 'force')
        My.append(force[2])
        u200x = ops.nodeDisp(200, 1)
        U200x.append(u200x)
        u100x = ops.nodeDisp(100, 1)
        U100x.append(u100x)



    T=[]
    T.append(2*np.pi/(eigValues[0]**0.5))   
    T.append(2*np.pi/(eigValues[1]**0.5))  
    den=4*np.pi*(Masses[0][1]*(U100x[0]**2)+Masses[1][1]*(U200x[0]**2))
    num=T[0]*(C[0]*(U100x[0]**2)*(np.cos(teta[0]))**2+C[1]*((U200x[0]-U100x[0])**2)*(np.cos(teta[1]))**2)
    csi1=num/den
    
    # post process the results doing the CQC modal combination for the My response (3rd column in section forces)
    MyCQC = CQC(My, eigValues, dmp, scalf)
    
    U200xCQC = CQC(U200x, eigValues, dmp, scalf)
    print("Nodo 200: UxCQC = ", U200xCQC)
    
    U100xCQC = CQC(U100x, eigValues, dmp, scalf)
    print("Nodo 100: UxCQC = ", U100xCQC)
    uu=[Tn2, U100xCQC, U200xCQC]
    U.append(uu)
    print('Spectral Analysis DONE!')
    
    
# =============================================================================
# print('\n\ANALISI 01:\nRun a Response Spectrum Analysis mode-by-mode.\nGrab results during the loop.\nDo CQC combination in post processing.\n')
# print('{0: >10}{1: >15}'.format('Mode', 'Vx [kN]'))
# for i in range(len(eigValues)):
# 	print('{0: >10}{1: >15f}'.format(i+1, My[i]/1000.))
# print('{0: >10}{1: >15f}'.format('CQC', MyCQC/1000.))
# print('{0: >10}{1: >15f}'.format('V', 4.*MyCQC/1000.) )
# print(" Tagliante complessivo alla base\n")
# 
# print("\n")
# 
# =============================================================================

opsv.plot_mode_shape(1)
plt.title('1st mode') 
opsv.plot_mode_shape(2)
plt.title('2nd mode')                     




#*****************************************************************************
#TIMEHISTORY
#*****************************************************************************

#SMORZAMENTO DI RAYLEIGH
MpropSwitch = 1.0 # M-prop. damping;  D = alphaM*M
KcurrSwitch = 1.0 # current-K;        +beatKcurr*KCurrent
KinitSwitch = 0.0 # initial-K;        +beatKinit*Kini
KcommSwitch = 0.0 # last-committed K; +betaKcomm*KlastCommitt

# Pick your modes and damping ratios
wi = eigValues[0]**0.5; zetai = damp # 5% in mode i
wj = eigValues[1]**0.5; zetaj = damp # 5% in mode j

W = np.array([[1/wi, wi],[1/wj, wj]])
zeta = np.array([zetai,zetaj])
ab = np.linalg.solve(W,2*zeta)
#ab[0]=2*wi*wj*zetai/(wi+wj)
#ab[1]=2*zetai/(wi+wj)


# =============================================================================
# period =  2*np.pi/wi
# print("Period T = ",period)
# 
# alfaM=2*damp*wi
# ops.rayleigh(alfaM, 0.,  0., 0.)    
# =============================================================================

# =============================================================================
# #                 M                   KT                 KI                 Kn
ops.rayleigh(MpropSwitch*ab[0], KcurrSwitch*ab[1], KinitSwitch*ab[1], KcommSwitch*ab[1])
# 
# =============================================================================




   #mode5dispacement=ops.nodeDisp(5,1)  
   #opsplt.plot_modeshape(1, 5000)
   #opsplt.plot_modeshape(2, 5000)
    
   #opsplt.plot_modeshape(3, 50000)
   #opsv.plot_model(node_labels=1, element_labels=1,axis_off=0)
   #plt.title('3d 3 floors: undeformed')
   
   #bz.PlotModel (plotmode=6,show_nodes_tag=True, onhover_message=True, vertical_axis=3,plot_legends=True, show_constrained=True, constrained_size=8)

# Record--------------------------------------------------------------------

max1FloordispRec=[]
maxTOPDsplcmtRec=[]

AllOutput=[]
  
def THrun(fileNumber):
    ops.reset()
    ops.wipeAnalysis() # clear previously-define analysis parameters
    ops.remove('recorders')
    ops.loadConst('-time', 0) 
    
    AccfileRecord=Input+THAccfilename+" ("+str(fileNumber)+").th"
    
    Acc=np.loadtxt(AccfileRecord,dtype='float',usecols = (0) )
    recordNumberPoint=Acc.size
   
    #timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor>
    ops.timeSeries('Path', fileNumber+1, '-dt', record_dt,'-filePath', AccfileRecord, '-factor', 1.0*g)# define acceleration vector from file (dt=0.02 is associated with the input file gm)
      
    # create data directory
    
    filespath=currentpath+"\\Output ("+str(fileNumber)+')'
    # Check whether the  
    # specified path is an 
    # existing directory or not 
    isdir = os.path.isdir(filespath) 
    if not isdir:
        os.makedirs(filespath)
    Output=filespath+"\\"
    AllOutput.append(Output)
    
    
    # # Define RECORDERS -------------------------------------------------------------
    # recorder Node -file $Output/Disp.out -time -node 4 -dof 1  disp;
    ops.recorder('Node', '-file', Output+'Node100relDisp.out', '-time', '-node', 100, '-dof', 1, 'disp')
    ops.recorder('Node', '-file', Output+'Node200relDisp.out', '-time', '-node', 200, '-dof', 1, 'disp')

    ops.recorder('Node', '-file', Output+'Node1absoluteAccel.out','-timeSeries', fileNumber+1, '-time', '-node', 1, '-dof', 1, 'accel')           
    ops.recorder('Node', '-file', Output+'Node100absoluteAccel.out','-timeSeries', fileNumber+1, '-time', '-node', 100, '-dof', 1, 'accel')
    ops.recorder('Node', '-file', Output+'Node200absoluteAccel.out','-timeSeries', fileNumber+1, '-time', '-node', 200, '-dof', 1, 'accel')    

    for dc in DamperCode:
        dfiledisp=Output+'Damper'+str(dc)+'disp.out'
        dfileforce=Output+'Damper'+str(dc)+'force.out'
        ops.recorder('Element', '-file', dfiledisp, '-time', '-ele', dc, 'deformations')
        ops.recorder('Element', '-file', dfileforce, '-time', '-ele', dc, 'localForce')
        
    
    ops.recorder('Element', '-file', Output+'Frameforce1.out', '-time', '-ele', 1, 2, '-dof', 1, 'force')
    ops.recorder('Element', '-file', Output+'Frameforce2.out', '-time', '-ele', 101, 102, '-dof', 1, 'force')
    
    # # recorder Node -file $Output/Base.out -time -node 1 2 -dof 1  reaction;  # support reaction
    # ops.recorder('Node', '-file', Output+'Base.out', '-time', '-node', 1, 2, '-dof', 1, 'reaction')
    
    # # recorder Node -file $Output/NBase.out -time -node 1 2 -dof 2  reaction;	 # support reaction
    # ops.recorder('Node', '-file', Output+'NBase.out', '-time', '-node', 1, 2, '-dof', 2, 'reaction')
    # # 		
    # # END Define RECORDERS -------------------------------------------------------------
   
    
    #   pattern  UniformExcitation  $patternTag $dir -accel      $tsTag <-vel0 $ver0>	         
    ops.pattern("UniformExcitation", fileNumber, 1, '-accel', fileNumber+1) # define where and how (pattern tag, dof) acceleration is applied
    setTimeHistoryAnalysis() 	
    err=ops.analyze(timeSubdivision*recordNumberPoint, record_dt/timeSubdivision) # apply timeSubdivision*recordNumberPoint steps for record_dt/timeSubdivision-sec time steps in analysis   	
    print("Errore = ", err)
          
    filename=Output+'Node200relDisp.out'
    disp_time=np.loadtxt(filename,dtype='float',usecols = (1,0) )
    time=disp_time[:,1]      
    dispTOP=disp_time[:,0]
    maxTOPDsplcmt=np.max(np.absolute(dispTOP))              
    tmaxdisp=time[np.argmax(np.absolute(dispTOP))]   
    
    print("max TOP displacement = ", maxTOPDsplcmt, "at time t = ", tmaxdisp, " sec")
    
    #massimi.append(massimo) / massimi.append(max(abs(max(disp)),abs(min(disp))))
    
    filename=Output+'Node100relDisp.out'
    disp1stFloor=np.loadtxt(filename,dtype='float',usecols = (1) )    
    
    max1Floordisp=np.max(np.absolute(disp1stFloor))              
    tmaxdisp=time[np.argmax(np.absolute(disp1stFloor))]  
    
    print("max First Floor displacement = ", max1Floordisp, "at time t = ", tmaxdisp, " sec")

    max1FloordispRec.append(max1Floordisp)
    maxTOPDsplcmtRec.append(maxTOPDsplcmt)

print()
for fileNumber in range(1,nfileRecordAccel+1):
    print('Time History Analysis Number',fileNumber,'just STARTED!')
    THrun(fileNumber)
    print('Time History Analysis Number',fileNumber,' DONE!')
    print()
    
print()
print('All',nfileRecordAccel,'Time History Analysis DONE!')
print()
max1FloordispRec=np.array(max1FloordispRec)
maxTOPDsplcmtRec=np.array(maxTOPDsplcmtRec)

min1F=np.min(max1FloordispRec)
mean1F=np.mean(max1FloordispRec)
max1F=np.max(max1FloordispRec)
firstFdev='{:.1f}%'.format(100*(max1F-mean1F)/mean1F)
k=np.argmin(np.absolute(np.ones(nfileRecordAccel)*mean1F-max1FloordispRec))

minTOP=np.min(maxTOPDsplcmtRec)
meanTOP=np.mean(maxTOPDsplcmtRec)
maxTOP=np.max(maxTOPDsplcmtRec)
TOPdev='{:.1f}%'.format(100*(maxTOP-meanTOP)/meanTOP)

print()
print('min of max 1stFloor Displacement = ', '{:.2f}'.format(min1F),'mm')
print('mean of max 1stFloor Displacement = ', '{:.2f}'.format(mean1F),'mm')
print('max of max 1stFloor Displacement = ', '{:.2f}'.format(max1F),'mm')
print('1stFloor Displacement Deviation between max and mean = ', firstFdev)
print()
print('min of max Top Floor Displacement = ', '{:.2f}'.format(minTOP),'mm')
print('mean of max Top Floor Displacement = ', '{:.2f}'.format(meanTOP),'mm')
print('max of max Top Floor Displacement = ', '{:.2f}'.format(maxTOP),'mm')
print('Top Floor Displacement Deviation between max and mean= ', TOPdev)

ops.wipe()
print()


print('The Nearest 1st Floor displacement to the mean is: ',max1FloordispRec[k])
print('for number ',k+1, ' time history acceleration')

#GRAFICI PER L'ULTIMO ACCELEROGRAMMA
Output=AllOutput[k]

filename=Output+'Node200relDisp.out'
disp_time=np.loadtxt(filename,dtype='float',usecols = (1,0) )
time=disp_time[:,1]      
dispTOP=disp_time[:,0]
maxTOPDsplcmt=np.max(np.absolute(dispTOP))              

filename=Output+'Node100relDisp.out'
disp1stFloor=np.loadtxt(filename,dtype='float',usecols = (1) )    
max1Floordisp=np.max(np.absolute(disp1stFloor))              

filename=Output+'Frameforce1.out'

Frameforces=np.loadtxt(filename,dtype='float',usecols = (1,2) )
V1stFloor=np.sum(Frameforces, axis=1)

filename=Output+'Frameforce2.out'
Frameforces=np.loadtxt(filename,dtype='float',usecols = (1,2) )
V2ndFloor=np.sum(Frameforces, axis=1)

Eet = np.absolute(V1stFloor*disp1stFloor/2)+np.absolute(V2ndFloor*(dispTOP-disp1stFloor)/2)

#EDamper1001=V1stFloor*0.0
#EDamper1002=V1stFloor*0.0

Damperdisp=[]
Damperforce=[]
EDamper=[]
if not telaioNUDO:
    for dc in DamperCode:
        filename=Output+'Damper'+str(dc)+'disp.out'
        Damperdisp.append(np.loadtxt(filename,dtype='float',usecols = (1) ))
        filename=Output+'Damper'+str(dc)+'force.out'   
        Damperforce.append(-np.loadtxt(filename,dtype='float',usecols = (1) ))
        EDamper.append(integrate.cumulative_trapezoid(Damperforce, Damperdisp, initial=0))
    
'''   
    filename=Output+'Damper1002disp.out'
    Damper1002disp=np.loadtxt(filename,dtype='float',usecols = (1) ) 
    
    filename=Output+'Damper1002force.out'
    Damper1002force=-np.loadtxt(filename,dtype='float',usecols = (1) )
    
    EDamper1002 = integrate.cumulative_trapezoid(Damper1002force, Damper1002disp, initial=0)
'''

filename=Output+'Node1absoluteAccel.out'
ag=np.loadtxt(filename,dtype='float',usecols = (1) )
vg=integrate.cumulative_trapezoid(ag,time, initial=0)

filename=Output+'Node100absoluteAccel.out'
a100=np.loadtxt(filename,dtype='float',usecols = (1) )
v100t=integrate.cumulative_trapezoid(a100,time, initial=0)

filename=Output+'Node200absoluteAccel.out'
a200=np.loadtxt(filename,dtype='float',usecols = (1) )
v200t=integrate.cumulative_trapezoid(a200,time, initial=0)

Eitabslt = integrate.cumulative_trapezoid((lmassX1*a100+lmassX2*a200)*vg, time, initial=0)

Ektabslt = (lmassX1*np.power(v100t,2)+lmassX2*np.power(v200t,2))/2

Eiabsltmax=np.max(Eitabslt) 
print('absolute Eimax = ',Eiabsltmax)

#vfo.plot_modeshape(modenumber=2, scale=300)

#opsplt.plot_model("nodes")

fgr, axs = plt.subplots(3)
fgr.suptitle('Time History del taglio alla base e dello spostamento ai nodo 101/102')
axs[0].plot(time,dispTOP,linewidth=1, color='r')
axs[0].set_title('TOP displacement relative to the base ground') 
axs[0].set(xlabel='t [sec]', ylabel='u [mm]')
axs[0].set_ylim([-maxTOPDsplcmt*1.10, maxTOPDsplcmt*1.10])

axs[1].plot(time,disp1stFloor,linewidth=1, color='g')
axs[1].set_title('First Floor Displacement relative to the base ground') 
axs[1].set(xlabel='t [sec]', ylabel='u [mm]')
axs[1].set_ylim([-maxTOPDsplcmt*1.10, maxTOPDsplcmt*1.10])

axs[2].plot(time,V1stFloor,linewidth=1)  
axs[2].set_title('Base Shear')      
axs[2].set(xlabel='t [sec]', ylabel='V [N]')  
plt.tight_layout()

if not telaioNUDO:
    fgr1, axs1 = plt.subplots(2)
    axs1[0].plot(Damperdisp[0],Damperforce[0],linewidth=1, color='r') 
    axs1[0].set_title('First floor Damper')
    axs1[1].plot(Damperdisp[2],Damperforce[2],linewidth=1, color='b') 
    axs1[1].set_title('Second floor Damper')
    plt.tight_layout()
# =============================================================================
#SISTEMARE LABEL PER LEGENDA
plt.figure(7)
#plt.plot(time,EDamper1001,              linewidth=1, color='red',   label='(1) first floor non-Linear Viscous Damped Energy')
#plt.plot(time,EDamper1001+EDamper1002,  linewidth=1, color='b',   label='(2)=(1) + second floor non-Linear Viscous Damped Energy')
#plt.plot(time,EDamper1001+EDamper1002+Ektabslt,    linewidth=1, color='blue',  label='(3)=(2)+ Kinetic Energy')
#plt.plot(time,EDamper1001+EDamper1002+Ektabslt+Eet,        linewidth=1, color='green', label='??(4)=(3)+ Linear Viscous Damped Energy')      
# plt.plot(time,Edt+Ecsit+Ekt,    linewidth=1, color='blue',  label='(3)=(2)+ Kinetic Energy')
# plt.plot(time,Edt+Ecsit+Ekt+Eet,linewidth=1, color='orange',label='(4)=(3)+ Elastic Energy')  
plt.plot(time,Eitabslt,              linewidth=2, color='black', label='  Total absolute Input Energy')            
plt.xlabel('t [sec]')
plt.ylabel('E [Nmm]')
# 
plt.legend(loc='lower right')
plt.xlim([0.0, time[-1]*1.05,])
# plt.ylim([0.0, Eimax*1.05])
# =============================================================================
#fg.set_title('Energia')



# fig_wi_he = 22., 14.
#fig_wi_he = 30., 20.
