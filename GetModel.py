# -*- coding: utf-8 -*-
"""
Created on Fri May 20 18:15:57 2022

@author: lucky
"""
from math import sqrt
import openseespy.opensees as ops

# Create nodes
def GetNodes(nodes,xyz):
    for i in range(0,len(nodes)):
        ops.node(nodes[i], xyz[i]['x'], xyz[i]['y'], xyz[i]['z'])

# Fix node
def FixNodes(fixednodes,fix):
    for i in range(0,len(fixednodes)):
        ops.fix(fixednodes[i], fix[i]['ux'], fix[i]['uy'], fix[i]['uz'], 
                               fix[i]['rx'], fix[i]['ry'], fix[i]['rz'])

# rigid floor constraint
#Create a multi-point constraint between nodes.
#These objects will constraint certain degrees-of-freedom 
#at the listed slave nodes to move as if in a rigid plane with the master node. 
#To enforce this constraint, Transformation constraint is recommende
def rgdDphrgm(master,slave):
    for i in range(0,len(master)):
        for j in range(0,len(slave[i])):
            ops.rigidDiaphragm(3, master[i], slave[i][j])
            
            
#Assegno le masse ai nodi 
def GetMasses(Masses):
    for i in range(0,len(Masses)):    
        ops.mass(Masses[i][0], Masses[i][1], Masses[i][2], Masses[i][3], 
                               Masses[i][4], Masses[i][5], Masses[i][6])
        
        
# Create Beams
def getBeams(beamCode,Beams):
    for i in range(0,len(beamCode)):
        ops.element('elasticBeamColumn', beamCode[i], Beams[i]["i"], Beams[i]["j"],
                                         Beams[i]["A"], Beams[i]["E"], Beams[i]["G"], 
                                         Beams[i]["J"], Beams[i]["Iy"], Beams[i]["Iz"],
                                         Beams[i]["gTT"]) 





# Create Beams
def getTimoshenkoBeams(beamCode,Beams):
    for i in range(0,len(beamCode)):
        #   element('ElasticTimoshenkoBeam',    
        ops.element('ElasticTimoshenkoBeam', beamCode[i], Beams[i]["i"], Beams[i]["j"], #eleTag,*eleNodes,
                                             Beams[i]["E"], Beams[i]["G"], Beams[i]["A"], #E_mod, G_mod, Area,                                           
                                             Beams[i]["J"], Beams[i]["Iy"], Beams[i]["Iz"], #Jxx, Iy, Iz,
                                             Beams[i]["Avy"], Beams[i]["Avz"], #Avy, Avz
                                             Beams[i]["gTT"]) #transfTag





 


# Create Dampers

def getDampers(DamperCode,Dampers):   
    # Damper Properties
    # Kd = 25.
    # Cd = 20.7452
    # ad = 0.35
    
    for i in range(0,len(DamperCode)):
        ops.uniaxialMaterial('ViscousDamper', i+1, Dampers[i]["Kd"], Dampers[i]["Cd"], Dampers[i]["ad"])
        #   element twoNodeLink $eleTag $iNode $jNode -mat $matTags -dir $dirs
        ops.element('twoNodeLink', DamperCode[i], Dampers[i]["i"], Dampers[i]["j"],'-mat', i+1, '-dir', 1)


def setAnalysis():
# DEFINE SOME ANALYSIS SETTINGS

    # This command is used to construct the ConstraintHandler object. 
    # The ConstraintHandler object determines how
    # the constraint equations are enforced in the analysis. 
    # Constraint equations enforce a specified value for a DOF,
    # or a relationship between DOFs.
    # ##constraints(constraintType, *constraintArgs)
    ops.constraints("Transformation")
    # Transformation Method
    # This command is used to construct a transformation constraint handler, 
    # which enforces the constraints using the transformation method.
  
        
    # This command is used to construct the DOF_Numberer object. 
    # The DOF_Numberer object determines the mapping between equation numbers 
    # and degrees-of-freedom – how degrees-of-freedom are numbered.
    # ##numberer(numbererType, *numbererArgs)
    ops.numberer("RCM")
    # RCM Numberer
    # This command is used to construct an RCM degree-of-freedom numbering object 
    # to provide the mapping between the degrees-of-freedom at the nodes and 
    # the equation numbers. An RCM numberer uses the reverse Cuthill-McKee scheme 
    # to order the matrix equations.
    
    
    # This command is used to construct the LinearSOE and LinearSolver objects 
    # to store and solve the system of equations in the analysis.
    # ##system(systemType, *systemArgs)
    ops.system("UmfPack")
    # This command is used to construct a sparse system of equations which uses 
    # the UmfPack solver.

    # This command is used to construct the LinearSOE and LinearSolver objects 
    #  to store and solve the test of equations in the analysis.   
    # ##test(testType, *testArgs)
    ops.test("NormUnbalance", 0.0001, 10)
    # ##test(’NormUnbalance’, tol, iter, pFlag=0, nType=2, maxIncr=maxIncr)
    # Create a NormUnbalance test, which uses the norm of the right hand side 
    # of the matrix equation to determine if convergence has been reached  
    
 
    # This command is used to construct a SolutionAlgorithm object, 
    # which determines the sequence of steps taken to solve the non-linear equation.   
    #   algorithm(algoType, *algoArgs)        
    ops.algorithm("Linear")  
    #ops.algorithm("Newton")
    #   algorithm(’Linear’, secant=False, initial=False, factorOnce=False)
    # Create a Linear algorithm which takes one iteration to solve the system 
    # of equations.
    

    # This command is used to construct the Integrator object. 
    # The Integrator object determines the meaning of the terms in the system 
    # of equation object Ax=B.
    #   integrator(intType, *intArgs)
    ops.integrator("LoadControl", 0.0)
    #   integrator(’LoadControl’, incr, numIter=1, minIncr=incr, maxIncr=incr)
    # Create a OpenSees LoadControl integrator object.



    # This command is used to construct the Analysis object, which defines
    # what type of analysis is to be performed.
    # • determine the predictive step for time t+dt
    # • specify the tangent matrix and residual vector at any iteration
    # • determine the corrective step based on the displacement increment dU    
    #   analysis(analysisType)    
    ops.analysis("Static")
    # Currently 3 valid options:
    #   1. 'Static' - for static analysis
    #   2. 'Transient' - for transient analysis constant time step
    #   3. 'VariableTransient' - for transient analysis with variable time step
    #   4. 'PFEM' - for PFEM analysis.


def setTimeHistoryAnalysis():
    # create the analysis				    
    ops.constraints('Transformation')     	 # how it handles boundary conditions
    ops.numberer('RCM')					     # renumber dof's to minimize band-width (optimization), if you want to
    ops.system('UmfPack')					 # how to store and solve the system of equations in the analysis (large model: try UmfPack)
    ops.test('EnergyIncr', 1.0e-16, 100) 	 # test Eneregy incerment 
    #ops.test('EnergyIncr', 1.0e-16, 100) 	 # test Eneregy incerment 
    ops.algorithm('KrylovNewton')            # use Kyrlow-Newton algorithm
    #ops.algorithm('RaphsonNewton')		 	  
    ops.integrator('Newmark', 0.5, 0.25)	 # determine the next time step for an analysis
    #ops.integrator('Newmark', 0.7, 0.35)	 # determine the next time step for an analysis
    ops.analysis('Transient')				 # define type of analysis: time-dependent




# CQC function
def CQC(mu, lambdas, dmp, scalf):
	u = 0.0
	ne = len(lambdas)
	for i in range(ne):
		for j in range(ne):
			di = dmp[i]
			dj = dmp[j]
			bij = lambdas[i]/lambdas[j]
			rho = ((8.0*sqrt(di*dj)*(di+bij*dj)*(bij**(3.0/2.0))) /
				((1.0-bij**2.0)**2.0 + 4.0*di*dj*bij*(1.0+bij**2.0) + 
				4.0*(di**2.0 + dj**2.0)*bij**2.0))
			u += scalf[i]*mu[i] * scalf[j]*mu[j] * rho
	return sqrt(u)