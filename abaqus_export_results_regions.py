## The goal here is to export any Abaqus results at an user specified nodeset
from odbAccess import *
import os
import sys
import csv
import numpy as np
import multiprocessing


def export_field_at_node_location_one_timestep(odb, stepName, increment, regionType, regionUser, node_list, fieldNameToOutput, dirToSave):
    
    tmp = odb.steps[stepName].frames[increment].fieldOutputs[fieldNameToOutput]

    # Get nodes if elementset
    if regionType == "elementset":
        if fieldNameToOutput == "S":
            field = tmp.getSubset(position=ELEMENT_NODAL)
        elif fieldNameToOutput == "U":
            field = tmp.getSubset(position=NODAL)    
    else:
        if fieldNameToOutput == "S":
            # We need ELEMENT_NODAL if we are requesting a node set.
            field = tmp.getSubset(position=ELEMENT_NODAL, region=regionUser)
        elif fieldNameToOutput == "U":
            field = tmp.getSubset(position=NODAL, region=regionUser)    
    
    nNodes = len(field.values)
    field_out = np.zeros((nNodes, 1+3))         # nodeLabel, dispX, dispY dispZ
    if fieldNameToOutput == "S":
        field_out = np.zeros((nNodes, 1+6))     # nodeLabel, S11, S22, S33, S12, S13, S23

    if stepName == "step_1" or stepName == "Geostatic":
        print("--Number of entries to export "+str(nNodes))

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
    
    # extract subset if elementset requested
    if regionType == "elementset":
        field_out_all = field_out
        field_out = np.array([])
        field_out = field_out_all[np.isin(field_out_all[:, 0], node_list), :] 
        nNodes = field_out.shape[0]
        
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
            f.write('%d, %0.0f, %0.0f, %0.0f, %0.0f, %0.0f, %0.0f\n' % (field_out[i, 0], field_out[i, 1], field_out[i, 2], field_out[i, 3], field_out[i, 4], field_out[i, 5], field_out[i, 6]))
    f.close()

def main_fun(odb, fieldNameToOutput, regionType, regName, stepNameUser, dir_odb):
    print("--Working on region " + regName)
    
    increment=-1 # last increment of each time step
    dirToSave= dir_odb + "extracted/" + fieldNameToOutput + "/" + regName + "/"
    if not os.path.exists(dirToSave):
        os.makedirs(dirToSave)
    if regionType == "elementset":
        regionUser = odb.rootAssembly.instances['PART-1-1'].elementSets[regName.replace(" ","")]
    elif regionType == "nodeset":
        regionUser = odb.rootAssembly.instances['PART-1-1'].nodeSets[regName.replace(" ","")]

    # Opening the file to save the fault tractions
    if stepNameUser == "All":
        stepNamesAll = list(odb.steps.keys()) ## get a list of all time steps names
    elif stepNameUser == "Geostatic":
        stepNamesAll = [stepNameUser]
    elif type(stepNameUser) == list:
        stepNamesAll = list(odb.steps.keys())
        stepNamesAll = [stepNamesAll[i] for i in stepNameUser]
    else:
        raise Exception("\n \t Error!! User input time step name not found")
        
    # READING COORDINATES
    #a=odb.rootAssembly.instances['PART-1-1'].nodes
    #alen=len(a)
    #print("--Reading COORD data")
    #COORD=np.zeros((alen,4))        
    #for i in range(0,alen): 
	#    COORD[i,0]=a[i].label               # NODE LABEL
	#    COORD[i,1]=a[i].coordinates[0]      # X COORDINATE
	#    COORD[i,2]=a[i].coordinates[1]      # Y COORDINATE
	#    COORD[i,3]=a[i].coordinates[2]      # Y COORDINATE        

    # Get nodes if elementset
    if regionType == "elementset":
        node_list = np.array([])
        for i in range(len(regionUser.elements)):   # This is slow
            node_list = np.append(node_list, regionUser.elements[i].connectivity)
        node_list = np.unique(node_list.astype(int))
        print("--Number of entries to export "+ str(len(node_list)))
    else:
        node_list = []

    for stepName in stepNamesAll:
        export_field_at_node_location_one_timestep(odb, stepName, increment, regionType, regionUser, node_list, fieldNameToOutput, dirToSave)


# Paths
dir_odb = "./output/elastic/"
fileNameToRead=dir_odb + "wilmington_elastic.odb"
print("\n \t Reading file "+fileNameToRead+"\n ")

# Inputs
regionType = "nodeset"                                   # elementset or nodeset
field_to_extract = "stress"                           # displacement or stress
stepNameUser=[0] + list(range(12, 301, 12))                 # "All" or "Geostatic" or list of numbers with tsteps
partName="Part_1"                                           #
if regionType == 'nodeset':
    if field_to_extract == "displacement":
        fieldNameToOutput = "U"    
        regionNames=["TOPRANGERSURFACE", "TOPBC"]   
    elif field_to_extract == "stress":
        fieldNameToOutput = "S"    
        regionNames=["COMPTONLOSALAMITOS_A", "COMPTONLOSALAMITOS_B", "COMPTONLOSALAMITOS_C", "NEWPORT", "PALOSVERDES", "TEAR_FAULT_2", "TEAR_FAULT_4", "TEAR_FAULT_5", "TEAR_FAULT_6", "TEAR_FAULT_8A", "TEAR_FAULT_8B", "TEAR_FAULT_9", "TEAR_FAULT_10A", "TEAR_FAULT_10B", "TEAR_FAULT_11", "TEAR_FAULT_13A", "TEAR_FAULT_13B", "TEAR_FAULT_15", "TEAR_FAULT_19", "TEAR_FAULT_21A", "TEAR_FAULT_21B", "TEAR_FAULT_26", "TEAR_FAULT_27"] 
elif regionType == 'elementset':
    regionNames=["RANGER"]
    if field_to_extract == "displacement":
        fieldNameToOutput = "U"          
    elif field_to_extract == "stress":
        fieldNameToOutput = "S"    

# Get data
odb = openOdb(path=fileNameToRead)
for regName in regionNames:
    main_fun(odb, fieldNameToOutput, regionType, regName, stepNameUser, dir_odb)
