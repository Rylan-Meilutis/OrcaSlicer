# The wx fork can already contain the macOS colour fix. Accept only an exact
# reverse-check match; unrelated changes must not silently bypass this patch.
execute_process(COMMAND "${GIT_EXECUTABLE}" apply --reverse --check "${PATCH_FILE}"
    RESULT_VARIABLE already_applied OUTPUT_QUIET ERROR_QUIET)
if(already_applied EQUAL 0)
    message(STATUS "wxWidgets macOS colour fix is already present")
    return()
endif()
execute_process(COMMAND "${GIT_EXECUTABLE}" apply --check "${PATCH_FILE}"
    RESULT_VARIABLE can_apply ERROR_VARIABLE patch_error)
if(NOT can_apply EQUAL 0)
    message(FATAL_ERROR "wxWidgets macOS colour patch does not match: ${patch_error}")
endif()
execute_process(COMMAND "${GIT_EXECUTABLE}" apply "${PATCH_FILE}"
    RESULT_VARIABLE applied ERROR_VARIABLE patch_error)
if(NOT applied EQUAL 0)
    message(FATAL_ERROR "Cannot apply wxWidgets macOS colour patch: ${patch_error}")
endif()
