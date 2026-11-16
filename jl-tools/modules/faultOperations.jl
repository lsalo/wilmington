module FaultOperations
    """
    This module performs fault operations, including computing fault tractions. 
    Primarily intended for use with abaqus .inp input files.
    Support is currently limited to first-order tetrahedral meshes.  
    """

    using LinearAlgebra
    using DelimitedFiles
    using Statistics

    include("./arrayOps.jl")
    using .ArrayOps
    
    export check_t0_tractions, fault_failure, fault_tractions, tractions_from_sij_and_n, face_unit_normal, face_centroid, face_normals_and_centroids

    function check_t0_tractions(tr0::AbstractMatrix, label::AbstractString=""; mu_max::Real=0.61)
        """Sanity check at t0: compressive normal tractions and |τ0| ≤ mu_max·|σ'n0| (tension positive).
        mu_max = largest friction used in the stress initialization (+ tolerance)."""
        @assert all(tr0[:, 2] .< 0) "$label: t0 normal tractions not all compressive"
        @assert all(tr0[:, 3] .<= mu_max .* abs.(tr0[:, 2])) "$label: t0 shear exceeds μ = $mu_max — check stress component order"
        return nothing
    end
    
    function fault_failure(tractions::AbstractMatrix, mu::Union{Number, Vector})
        """
        INPUTS
        tractions: Matrix with two columns: normal traction (tn), and shear traction magnitude (ts)
        mu: scalar or a vector (with length = number of rows in tractions Matrix)
        """
        cff = tractions[:, 2] .+ mu.*tractions[:, 1]

        return cff
    end

    
    function fault_tractions(s::Matrix, n::Matrix)
        """
        Compute surface tractions from stress tensor and unit normal vectors.
        Note: This computes the shear traction magnitude only. In some cases we may need the shear in a given direction on the surface (e.g., fault rake), which needs extension.
        Note: Normals should have consistent sign (e.g., upward normal vectors)

        INPUT
        s:  Matrix of size m by 6, with s11 (=sxx), s22 (=syy), s33 (=szz), s12 (=sxy), s13 (=sxz), s23 (=syz)
        n:  Matrix of size m by 3, with nx, ny, nz
        Any given row of s must correspond to the same row of n
        """

        s11 = s[:,1]
        s22 = s[:,2]
        s33 = s[:,3]
        s12 = s[:,4]
        s13 = s[:,5]
        s23 = s[:,6]

        tx = s11.*n[:,1] + s12.*n[:,2] + s13.*n[:,3]
        ty = s12.*n[:,1] + s22.*n[:,2] + s23.*n[:,3]
        tz = s13.*n[:,1] + s23.*n[:,2] + s33.*n[:,3]

        t = hcat(tx, ty, tz)
        tn = dot.(eachrow(n), eachrow(t))
        ts = sqrt.(norm.(eachrow(t)).^2 - tn.^2)

        return tn, ts
    end

    function tractions_from_sij_and_n(stress_vals::Matrix, stress_coords::AbstractMatrix, face_unit_normals::Matrix{Float64}, face_centroids::Matrix{Float64}, dir_path::String=nothing, fname::String=nothing)
        """
        Get tractions acting on a surface from stress tensor and normals at specified points.

        INPUTS
        stress_vals:            n by 7 matrix with column ids being identifier (e.g. node number), s11 (=sxx), s22 (=syy), s33 (=szz), s12 (=sxy), s13 (=sxz), s23 (=syz).
        stress_coords:          n by 4 matrix with column ids being identifier (e.g. node number), x, y, z
        face_unit_normals:      Dimensions nface x 3
        face_centroids:         Dimensions nface x 3
        fname:                  filename where fault traction output will be saved

        OUTPUT
        writes file with fault tractions to location in fname
        nodal_tractions: Matrix with three columns: id, normal traction (tn), and shear traction magnitude (ts)
        """
        
        # Sort stress matrix and check for repeated ids
        stress_vals = sortslices(stress_vals, dims=1, by=x->x[1])
        unique_ids = sort(unique(stress_vals[:, 1]))
        if length(unique_ids) < size(stress_vals)[1]
            repeated_stress_ids = true
        else
            repeated_stress_ids = false
        end
        
        # Average stresses in case of repeated nodes (extrapolation of stress to the nodes)
        if repeated_stress_ids
            stress_vals = ArrayOps.reduce_matrix_with_repeated_ids(stress_vals, "mean")
        end

        # Assign face normals to locations where stress is available (e.g., closest node to face centroids).
        # Currently, this is limited to closest point (min. Euclidean distance).
        stress_unit_normals = zeros(Union{Int64, Float64}, length(unique_ids), 4)
        for (i, id) in enumerate(unique_ids)
            stress_unit_normals[i, 1] = id
            id_coord = stress_coords[stress_coords[:,1] .== id, 2:end]
            dist = sqrt.(sum((face_centroids .- id_coord) .^2, dims=2))
            stress_unit_normals[i, 2:end] = face_unit_normals[argmin(dist)[1], :]
        end

        # Compute fault tractions
        tn, ts = fault_tractions(stress_vals[:, 2:end], stress_unit_normals[:, 2:end])    
        nodal_tractions = Matrix{Union{Int64, Float64}}(undef, length(unique_ids), 3)
        nodal_tractions[:, 1] = Int.(round.(unique_ids))
        nodal_tractions[:, 2:end] = hcat(tn, ts)

        # Export a file with node, ts, tn, 
        if dir_path !== nothing && fname !== nothing
            if !isdir(dir_path)
                mkdir(dir_path)
            end
            open(joinpath([dir_path, fname]), "w") do f
                for row in eachrow(nodal_tractions)
                    println(f, join(row, "\t"))
                end
            end
        end

        return nodal_tractions
    end

    function face_unit_normal(coords::Matrix)
        """
        Get upward-pointing unit normal vectors from the coordinates of three vertices (nodes) in the face
        """

        f = coords[2,:]-coords[1,:]
        g = coords[3,:]-coords[1,:]
        n = cross(f, g) 
        for i in 3:-1:1 # ensure that we get upward normal vector, start with z, and normalize (unit vector)
            if n[i] != 0
                sgn = n[i] >= 0 ? 1 : -1
                n = n .* sgn / norm(n) 
                break
            end
        end

        return n

    end

    function face_centroid(coords::Matrix, face_type::String)
        """
        Get centroid from the coordinates of an appropriate number of vertices (nodes) in the face
        """

        if face_type == "triangle"
            x = (coords[1,1] + coords[2,1] + coords[3,1]) / 3
            y = (coords[1,2] + coords[2,2] + coords[3,2]) / 3 
            z = (coords[1,3] + coords[2,3] + coords[3,3]) / 3
        end

        return [x, y, z]
    end

    function face_normals_and_centroids(face_connectivity::Vector{Matrix{Int64}}, nodes_inp::AbstractMatrix)
        """
        Uses face connectivity (nodes for each face) and computes normals and centroids
        """

        face_unit_normals = [zeros(Float64, 0, 3) for _ in 1:length(face_connectivity)]
        face_centroids = [zeros(Float64, 0, 3) for _ in 1:length(face_connectivity)]
        for nset in eachindex(face_connectivity)
            faces = face_connectivity[nset]

            for face_node_ids in eachrow(faces)
                # Coords of face nodes
                coords = nodes_inp[in.(nodes_inp[:,1], Ref(face_node_ids)), 2:end]

                # unit upward-pointing normal from two edges
                n = face_unit_normal(coords)
                face_unit_normals[nset] = vcat(face_unit_normals[nset], n')

                # centroids
                c = face_centroid(coords, "triangle")
                face_centroids[nset] = vcat(face_centroids[nset], c')   
            end

            print(join(["Face normals and centroids for nodeset ", nset, "/", length(face_connectivity), " completed.\n"]))
        end

        return face_unit_normals, face_centroids
    end


end