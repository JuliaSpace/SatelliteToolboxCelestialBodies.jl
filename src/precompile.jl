## Description #############################################################################
#
# Precompilation workload.
#
############################################################################################

PrecompileTools.@compile_workload begin
    jd_tdb   = 2_451_545.0
    date_tdb = DateTime(2000, 1, 1, 12, 0, 0)

    # == Sun ===============================================================================

    # Sun state, position, and velocity with `Float64` Julian day and `DateTime` inputs.
    sun_state_mod(jd_tdb)
    sun_state_mod(date_tdb)
    sun_position_mod(jd_tdb)
    sun_position_mod(date_tdb)
    sun_velocity_mod(jd_tdb)
    sun_velocity_mod(date_tdb)

    # == Moon ==============================================================================

    # Moon state, position, and velocity with the default model.
    moon_state_mod(jd_tdb)
    moon_state_mod(date_tdb)
    moon_position_mod(jd_tdb)
    moon_position_mod(date_tdb)
    moon_velocity_mod(jd_tdb)
    moon_velocity_mod(date_tdb)

    # Moon state, position, and velocity with each model explicitly selected.
    for model in (Val(:Meeus), Val(:Vallado))
        moon_state_mod(jd_tdb, model)
        moon_state_mod(date_tdb, model)
        moon_position_mod(jd_tdb, model)
        moon_position_mod(date_tdb, model)
        moon_velocity_mod(jd_tdb, model)
        moon_velocity_mod(date_tdb, model)
    end
end
