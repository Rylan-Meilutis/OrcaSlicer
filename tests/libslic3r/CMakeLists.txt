get_filename_component(_TEST_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Observe the geometry backend's actual CGAL layout without duplicating its flags.
# This source is only added when the test suite is enabled.
target_sources(libslic3r_cgal PRIVATE texture_import_cgal.cpp)

add_executable(${_TEST_NAME}_tests
    ${_TEST_NAME}_tests.cpp
    test_3mf.cpp
    test_amf.cpp
    # Import naming policy is independent of the GUI and slicing pipeline.
    test_import_names.cpp
    # Round-trip seam metadata and active/dormant volume settings in both formats.
    test_precise_seam_3mf.cpp
    # Pure perimeter extraction is independent of Print/Layer fixtures.
    test_precise_seam.cpp
    # Contour correspondence is pure geometry; Print integration stays in fff_print.
    test_seam_contour.cpp
    test_aabbindirect.cpp
    test_appconfig.cpp
    test_arc_thermal.cpp
    test_app_update.cpp
    test_arachne_walls.cpp
    test_arrange.cpp
    test_assemble_list.cpp
    test_bambu_networking.cpp
    test_belt_brim.cpp
    test_warp_mitigation.cpp
    test_thermoelastic_solver.cpp
    test_thermoelastic_reduction.cpp
    test_thermoelastic_relaxation.cpp
    test_thermoelastic_material.cpp
    test_thermal_material_references.cpp
    test_thermoelastic_thermal_strain.cpp
    test_thermomechanical_simulation.cpp
    test_thermoelastic_mesh.cpp
    test_thermoelastic_layers.cpp
    test_thermoelastic_boundary.cpp
    test_thermoelastic_exact_sections.cpp
    test_deposition_provenance.cpp
    test_deposition_section.cpp
    test_thermoelastic_section_topology.cpp
    test_thermoelastic_sections.cpp
    test_thermoelastic_sparse_factor.cpp
    test_thermoelastic_section_study.cpp
    test_thermoelastic_local_sections.cpp
    test_thermoelastic_peeling.cpp
    test_thermoelastic_homogenization.cpp
    test_thermal_homogenization.cpp
    test_tool_mapping.cpp
    test_cooling_assistant.cpp
    test_buildvolume.cpp
    test_brim.cpp
    test_calib.cpp
    test_octoprint_spool_cache.cpp
    test_indx_profile.cpp
    test_thermal_history.cpp
    test_section_integration.cpp
    test_section_thermal_assembly.cpp
    test_section_growth.cpp
    test_section_growth_stream.cpp
    test_section_thermal_sequence.cpp
    test_section_thermal_stream.cpp
    test_section_thermal_skin.cpp
    test_thermal_deposition.cpp
    test_thermomechanical_history.cpp
    test_learned_timing.cpp
    test_clipper_offset.cpp
    test_clipper_utils.cpp
    test_cut_connectors.cpp
    test_config.cpp
    test_config_variant_expansion.cpp
    test_locales_utils.cpp
    test_drc.cpp
    test_gcodereader.cpp
    test_local_extrusion_groups.cpp
    test_toolordering_nozzle_group.cpp
    test_parallel_resolve.cpp
    test_preset_bundle_loading.cpp
    test_preset_setting_id.cpp
    test_preset_diff.cpp
    test_sequential_gantry_geometry.cpp
    test_nonplanar_surface.cpp
    test_publish_settings.cpp
    test_vendor_cache.cpp
    test_preset_options.cpp
    test_edge_grid.cpp
    test_elephant_foot_compensation.cpp
    test_exclusion_volumes.cpp
    test_fill_corner_smoothing.cpp
    test_filament_mixer.cpp
    test_fill_plane_path.cpp
    test_adaptive_infill_field.cpp
    test_geometry.cpp
    test_kdtree.cpp
    test_lightning.cpp
    test_imex_arrange.cpp
    test_imex_helpers.cpp
    test_imex_zones.cpp
    test_multimaterial_segmentation.cpp
    test_placeholder_parser.cpp
    test_png_read_write.cpp
    test_polygon.cpp
    test_mutable_polygon.cpp
    test_mutable_priority_queue.cpp
    test_minimum_spanning_tree.cpp
    test_nozzle_volume_type.cpp
    test_obj.cpp
    test_step.cpp
    test_stl.cpp
    test_triangle_mesh_deal.cpp
    test_triangle_selector.cpp
    test_meshboolean.cpp
    test_mesh_repair.cpp
    test_trianglemesh_slicer.cpp
    test_marchingsquares.cpp
    test_lay_on_face.cpp
    test_model.cpp
    test_utils.cpp
    test_miniz_extension.cpp
    test_timeutils.cpp
    test_semver.cpp
    test_voronoi.cpp
    test_wipe_tower_estimate.cpp
    test_wipe_tower.cpp
    test_wipe_path.cpp
    test_optimizers.cpp
    test_ordering_strategies.cpp
    test_shortest_path.cpp
    # test_png_io.cpp
    test_indexed_triangle_set.cpp
    test_connected_bodies.cpp
    test_texture_displacement.cpp
    test_texture_import.cpp
    test_instance_lock.cpp
    ../libnest2d/printer_parts.cpp
    )

if (SLIC3R_CAD)
    target_sources(${_TEST_NAME}_tests PRIVATE
        test_caddocument.cpp
        test_sketchconstraints.cpp
        test_sketchedit.cpp
        test_sketchprofile.cpp
        test_sketchimport.cpp
        test_sketchinference.cpp
        test_slvs_constraints.cpp)
endif ()

if (TARGET OpenVDB::openvdb)
    target_sources(${_TEST_NAME}_tests PRIVATE test_hollowing.cpp)
endif()

# libslic3r links Draco privately; test_drc.cpp encodes its own fixtures.
find_package(Draco REQUIRED)
target_link_libraries(${_TEST_NAME}_tests test_common libslic3r draco::draco Catch2::Catch2WithMain)
target_include_directories(${_TEST_NAME}_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
if (SLIC3R_PCH AND NOT SLIC3R_SYNTAXONLY)
    add_precompiled_header(${_TEST_NAME}_tests ${CMAKE_SOURCE_DIR}/src/libslic3r/pchheader.hpp FORCEINCLUDE)
endif ()

set_property(TARGET ${_TEST_NAME}_tests PROPERTY FOLDER "tests")

orcaslicer_copy_test_dlls()

orcaslicer_discover_tests(${_TEST_NAME}_tests)
