project (irbem Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")


set(IRBEM_FILES
        AE8_AP8.f
        AFRL_CRRES_models.f
        Alexeev2000.f
        calcul_Lstar_o.f
        CoordTrans.f
        date_util.f
        dipole.f
        Dmag.f
        drift_bounce_orbit.f
        esa_trara.f
        field_line_tracing.f
        field_line_tracing_towards_Earth.f
        find_bm.f
        find_foot.f
        geopack_08.f
        get_bderivs.f
        get_hemi.f
        heliospheric_transformation.f
        IGE.f
        igrf_coef.f
        igrf.f
        init_nouveau.f
        init_ts07d_coeffs.f
        init_ts07d_tlpr.f
        internal_field.f
        loc_equator.f
        mead.f
        Mead_Tsyganenko.f
        meo_gnss.f
        msis86.f
        msise90.f
        myOwnMagField.f
        myOwnMagField_init.f
        nrlmsise00.f
        Olson_Pfitzer_dyn.f
        onera_desp_lib.f
        Ostapenko.f
        Pfitzer_q.f
        prot_and_elbr_bas2.f
        sgp4_ele.f
        sgp4ext.f
        sgp4io.f
        sgp4_orb.f
        sgp4_tle.f
        sgp4unit.f
        shieldose2.f
        t01_s.f
        trace_drift_shell.f
        TS07j_2015.f
        Tsy_and_Sit07_2015.f
        Tsy_and_Sit07_Jul2017.f
        Tsyganenko01.f
        Tsyganenko04.f
        Tsyganenko87l.f
        Tsyganenko87s.f
        Tsyganenko89.f
        Tsyganenko96.f
        )
add_library(irbem ${IRBEM_FILES})
