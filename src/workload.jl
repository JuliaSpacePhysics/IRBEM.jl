function workload()
    kext = T89
    t = DateTime("2015-02-02T06:12:43")
    𝐫 = GDZ(651.0, 63.0, 15.9)
    maginput = (; Kp = 40.0)
    return (
        # Compute L* and related parameters
        make_lstar(t, 𝐫, maginput; kext),
        # Find the magnetic equator
        find_magequator(t, 𝐫, maginput; kext),
        # Calculate MLT
        get_mlt(t, 𝐫),
    )
end
