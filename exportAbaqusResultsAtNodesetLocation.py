## The goal here is to export any Abaqus results at an user specified nodeset
from odbAccess import *
import os
import sys
import csv
import numpy as np
import multiprocessing


def export_field_at_node_location_one_timestep(odb, stepName, increment, nodesetUser, nodesetName, COORD, fieldNameToOutput, dirToSave):
    print("--Reading time step "+ stepName)
    print("--Field to output ", fieldNameToOutput)
    tmp = odb.steps[stepName].frames[increment].fieldOutputs[fieldNameToOutput]

    if fieldNameToOutput == "S":
        field = tmp.getSubset(position=ELEMENT_NODAL, region=nodesetUser)
        #field = tmp.getSubset(position=CENTROID,region=nodesetUser)
    elif fieldNameToOutput == "U":
        field = tmp.getSubset(position=NODAL, region=nodesetUser)    

    nNodes = len(field.values)
    field_out = np.zeros((nNodes, 1+3))         # nodeLabel, dispX, dispY dispZ
    if fieldNameToOutput == "S":
        field_out = np.zeros((nNodes, 1+6))     # nodeLabel, S11, S22, S33, S12, S13, S23

    if stepName == "step_1" or stepName == "Geostatic":
        print("--Number of nodes to export "+str(nNodes))
        
    for iNode in range(0, nNodes):
        nodeLabel = field.values[iNode].nodeLabel
        data = field.values[iNode].data
        
        field_out[iNode, 0]=nodeLabel
        field_out[iNode, 1]= data[0]
        field_out[iNode, 2]= data[1]
        field_out[iNode, 3]= data[2]
        if fieldNameToOutput == "S":
            field_out[iNode, 4]= data[3]
            field_out[iNode, 5]= data[4]
            field_out[iNode, 6]= data[5]

    # Write the output
    print("--Writing time step " + stepName)
    fileNameToSave=dirToSave + stepName + ".dat"
    f=open(fileNameToSave, "w")
    for i in range(0, nNodes):
        if i == 0:
            if fieldNameToOutput == "U":
                f.write('NODE_LABEL, dispX, dispY, dispZ\n')
            elif fieldNameToOutput == "S":
                f.write('NODE_LABEL, Sxx, Syy, Szz, Sxy, Sxz, Syz\n')

        if fieldNameToOutput == "U":
	        f.write('%d, %0.3f, %0.3f, %0.3f\n' % (field_out[i, 0], field_out[i, 1], field_out[i, 2], field_out[i, 3]))
        elif fieldNameToOutput == "S":
            f.write('%d, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n' % (field_out[i, :]))
    f.close()

def main_fun(odb, fieldNameToOutput, nodesetName, stepNameUser, dir_odb):
    print("--Working on nodeset "+nodesetName)
    
    increment=-1 # last increment of each time step
    dirToSave= dir_odb + "extracted/nodeset_" + nodesetName + "/" + fieldNameToOutput + "/"    
    if not os.path.exists(dirToSave):
        os.makedirs(dirToSave)
    nodesetUser = odb.rootAssembly.instances['PART-1-1'].nodeSets[nodesetName.replace(" ","")]
    
    # Opening the file to save the fault tractions
    if stepNameUser == "All":
        stepNamesAll = list(odb.steps.keys()) ## get a list of all time steps names
    elif stepNameUser == "Geostatic":
        stepNamesAll = [stepNameUser]
    else:
        raise Exception("\n \t Error!! User input time step name not found")
        
    # READING COORDINATES
    a=odb.rootAssembly.instances['PART-1-1'].nodes
    alen=len(a)

    #print("--Reading COORD data")
    COORD=np.zeros((alen,4))        
    for i in range(0,alen): 
	    COORD[i,0]=a[i].label               # NODE LABEL
	    COORD[i,1]=a[i].coordinates[0]      # X COORDINATE
	    COORD[i,2]=a[i].coordinates[1]      # Y COORDINATE
	    COORD[i,3]=a[i].coordinates[2]      # Y COORDINATE        

    for stepName in stepNamesAll:
        export_field_at_node_location_one_timestep(odb, stepName, increment, nodesetUser, nodesetName, COORD, fieldNameToOutput, dirToSave)


# Inputs
fieldNameToOutput = "U"     # or U
nodesetName="TOPRANGERSURFACE"       
stepNameUser="Geostatic"    # "All" or "Geostatic"
partName="Part_1"           #

# Paths
dir_odb = "./output/elastic/"
fileNameToRead=dir_odb + "wilmington_elastic.odb"
print("\n \t Reading file "+fileNameToRead+"\n ")

# Get data
odb = openOdb(path=fileNameToRead)
main_fun(odb, fieldNameToOutput, nodesetName, stepNameUser, dir_odb)
