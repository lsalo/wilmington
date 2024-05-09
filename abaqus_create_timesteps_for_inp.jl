"""
This file writes timestep info to a .inp file for Abaqus
"""

# Setup
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using Revise
println("Finished loading packages")

# Filename and Directories
dir_p_abq = "C:/Users/lsalo/matlab/sim_data/mrst/wilmington/liquidRate_JS_LS_lowRsprops_owc1800/result_data/dp_for_abaqus/"
dir_out = "C:/Users/lsalo/abaqus/simulation_data/wilmington/"

# Write timestep info
fileNames=readdir(dir_p_abq)
fileNameToSave=dir_out*"tsteps_1936_to_1961_p_varnum_2.inp"
p_variable_num = 2
open(fileNameToSave, "w") do f
    str = ""
    t_inc_ini_s = 3600*24*2 # 2 days
    tstep_size_s = 2629746  # from mrst: year/12
    t_inc_max_s = tstep_size_s
    for (i, fileName) in enumerate(fileNames)
        stepNumber=i
        # Note: triple quotes for multiline strings, otherwise it will be interpreted as character
        str = """*step, name=step_$stepNumber, inc=1000, nlgeom=no, unsymm=yes
        *STATIC
        $t_inc_ini_s, $tstep_size_s, 0, $t_inc_max_s
        *Field, variable=$p_variable_num
        *Include, input=./nodal_dp_1936_to_1961/$fileName
        *contact controls, stabilize
        *Controls, ANALYSIS=DISCONTINUOUS
        *end step
        **
        """
        write(f, str)
    end
end

