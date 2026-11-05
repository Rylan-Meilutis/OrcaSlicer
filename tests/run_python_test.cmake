# Test artifacts run on a different runner without the build's bundled Python.
# Resolve the host interpreter at test time instead of baking in its build path.
find_program(TEST_PYTHON NAMES python3 python REQUIRED)
execute_process(COMMAND "${TEST_PYTHON}" "${TEST_SCRIPT}"
                RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Python test failed: ${result}")
endif()
