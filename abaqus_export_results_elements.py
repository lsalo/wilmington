from odbAccess import *
import os

# Paths and options
dir_odb = "./output/plastic_initialization/"
fileNameToRead = dir_odb + "wilmington_plastic_initialization.odb"
field_name = "S"    
elset_name = "ENTIRE_DOMAIN"
step_name = "Geostatic"

# Open the output database
print("\n \t Reading file " + fileNameToRead + "\n ")
odb = openOdb(path = fileNameToRead)

# Get step and frame and field output
step = odb.steps[step_name]
lastFrame = step.frames[-1]
stressField = lastFrame.fieldOutputs[field_name]
# elementSet = odb.rootAssembly.elementSets[elset_name]     # Direct elsets such as ALL ELEMENTS
# instance = odb.rootAssembly.instances['INSTANCE-1']       # Elsets under an instance
# elementSet = instance.elementSets['ENTIRE_DOMAIN']
part = odb.rootAssembly.instances["PART-1-1"]               # Elsets under a part: replace 'PART-1-1' with your part name
elementSet = part.elementSets[elset_name]
elementStressField = stressField.getSubset(region=elementSet)

# Open a file for writing
dirToSave= dir_odb + "extracted/" + field_name + "/" + elset_name + "/"
if not os.path.exists(dirToSave):
    os.makedirs(dirToSave)
fileNameToSave = dirToSave + step_name + ".dat"

with open(fileNameToSave, 'w') as f:
    # Write the header
    f.write('ELEMENT_ID, Sxx, Syy, Szz, Sxy, Sxz, Syz\n')

    # Loop over each element in the element set
    for element in elementStressField.values:
        # Write the stress components to the file
        f.write('%d, %0.0f, %0.0f, %0.0f, %0.0f, %0.0f, %0.0f\n' % (element.elementLabel, element.data[0], element.data[1], element.data[2], element.data[3], element.data[4], element.data[5]))

# Close the output database
odb.close()