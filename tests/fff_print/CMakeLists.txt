get_filename_component(_TEST_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
add_executable(${_TEST_NAME}_tests 
	${_TEST_NAME}_tests.cpp
	test_helpers.cpp
	test_helpers.hpp
	test_adaptive_pressure_advance.cpp
	test_calib.cpp
	test_cooling.cpp
    test_pressure_equalizer.cpp
    test_adaptive_pa_processor.cpp
    test_deposition_provenance_export.cpp
    test_deposition_section_export.cpp
    test_deposition_provenance_project.cpp
    test_thermoelastic_source_ownership.cpp
    test_thermoelastic_source_diagnostics.cpp
    test_thermoelastic_conflict_ownership.cpp
	test_dock_fan.cpp
    test_thermal_airflow.cpp
    test_thermal_airflow_timing.cpp
    test_thermomechanical_gcode.cpp
    test_thermomechanical_layer_preview.cpp
    test_thermoelastic_toolpath.cpp
    test_thermoelastic_gcode_sections.cpp
    test_section_gcode_history.cpp
    test_section_airflow_scenario.cpp
    test_section_cooling_study.cpp
    test_thermoelastic_section_preview.cpp
	test_extrusion_entity.cpp
	test_extrusion_processor.cpp
	test_exclusion_volumes.cpp
	test_fan_mover.cpp
	test_fill.cpp
	test_flow.cpp
	test_gcode_processor.cpp
	test_gcode_statistics.cpp
	test_gcode_timing.cpp
	test_gcode.cpp
	test_gcodeprocessor.cpp
	test_gcodewriter.cpp
	test_imex_mode_gcode.cpp
	test_imex_slice_offset.cpp
	test_model.cpp
	test_multifilament.cpp
	test_multimaterial_segmentation.cpp
	test_perimeters.cpp
	# Seam geometry and pipeline regressions share the existing Print fixtures.
	test_precise_seam.cpp
	test_print.cpp
	test_printobject.cpp
	test_mixed_filament.cpp
	test_seam_placer.cpp
	test_skirt_brim.cpp
	test_warp_mitigation.cpp
	test_thermoelastic_gcode.cpp
    test_thermoelastic_ownership.cpp
	test_deposition_mesh.cpp
	test_deposition_geometry.cpp
	test_tool_dispatch.cpp
	test_octoprint_timing.cpp
	test_slicing_pipeline_hook.cpp
	test_spiral_vase.cpp
	test_support_material.cpp
	test_texture_displacement.cpp
	test_tree_support.cpp
	test_trianglemesh.cpp
	test_wipe.cpp
	test_wipe_tower.cpp
	)
target_link_libraries(${_TEST_NAME}_tests test_common libslic3r Catch2::Catch2WithMain)
if (SLIC3R_PCH AND NOT SLIC3R_SYNTAXONLY)
    add_precompiled_header(${_TEST_NAME}_tests ${CMAKE_SOURCE_DIR}/src/libslic3r/pchheader.hpp FORCEINCLUDE)
endif ()

set_property(TARGET ${_TEST_NAME}_tests PROPERTY FOLDER "tests")

orcaslicer_copy_test_dlls()

orcaslicer_discover_tests(${_TEST_NAME}_tests)

# Exercise the emitted-G-code audit itself, including known-invalid controls.
# A validator which silently accepts unsupported walls must fail the test run.
add_test(NAME gcode_audit_regressions
    COMMAND cmake "-DTEST_SCRIPT=${CMAKE_SOURCE_DIR}/test_audit_gcode_moves.py"
            -P "${CMAKE_SOURCE_DIR}/tests/run_python_test.cmake")
set_tests_properties(gcode_audit_regressions PROPERTIES
    LABELS "GCode;NonplanarSurface;Regression" TIMEOUT 120)
