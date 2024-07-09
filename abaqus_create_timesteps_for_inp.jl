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
dir_p_abq = "C:/Users/lsalo/matlab/sim_data/mrst/wilmington/liquidRate_JS_LS_lowRsprops_owc1800/result_data/dp_for_abaqus_1936_1975/"
dir_out = "C:/Users/lsalo/julia/output/wilmington/"

# Write timestep info
fileNames=readdir(dir_p_abq)
dirName = "nodal_dp_1936_to_1975"
fileNameToSave=dir_out*"tsteps_1936_to_1975_plastic_camclay.inp"
p_variable_num = 2
tstep_conv = 300
inc = 1000
open(fileNameToSave, "w") do f
    str = ""
    tstep_size_s = 2629746  # from mrst: year/12
    t_inc_max_s = tstep_size_s
    for (i, fileName) in enumerate(fileNames)
        stepNumber=i

        # Note: triple quotes for multiline strings, otherwise it will be interpreted as character
        if i <= tstep_conv
            t_inc_ini_s = 3600*24*5 # 5 days
            t_inc_min_s = 8640
            str = """*step, name=step_$stepNumber
            *STATIC
            $t_inc_ini_s, $tstep_size_s, $t_inc_min_s, $t_inc_max_s
            *Field, variable=$p_variable_num
            *Include, input=../$dirName/$fileName
            *end step
            **
            """
        elseif i < length(fileNames)
            t_inc_ini_s = 3600*24*1 # 1 days
            t_inc_min_s = 900
            str = """*step, name=step_$stepNumber, inc=$inc, nlgeom=no, unsymm=no
            *STATIC, stabilize
            $t_inc_ini_s, $tstep_size_s, $t_inc_min_s, $t_inc_max_s
            *Field, variable=$p_variable_num
            *Include, input=../$dirName/$fileName
            *Controls, analysis=discontinuous
            *end step
            **
            """
        else
            t_inc_ini_s = 3600*24*1 # 1 days
            t_inc_min_s = 900
            str = """*step, name=step_$stepNumber, inc=$inc, nlgeom=no, unsymm=no
            *STATIC, stabilize
            $t_inc_ini_s, $tstep_size_s, $t_inc_min_s, $t_inc_max_s
            *Field, variable=$p_variable_num
            *Include, input=../$dirName/$fileName
            *Controls, analysis=discontinuous
            *Restart, WRITE
            *end step
            **
            """
        end

        write(f, str)
    end
end

