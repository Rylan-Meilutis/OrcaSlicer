# Eigen 5.0.1's SSE2 int64 comparison treats the low word as signed even
# when the high words are equal. This reverses bounding boxes across 2^31.
# Use CMake rather than git apply: extracted dependencies can reside inside
# another repository, where git silently skips paths outside the current prefix.
set(_header "${EIGEN_SOURCE_DIR}/Eigen/src/Core/arch/SSE/PacketMath.h")
file(READ "${_header}" _content)
set(_old "  Packet2l lo_lt = Packet2l(_mm_shuffle_epi32(lt, (shuffle_mask<0, 0, 2, 2>::mask)));")
set(_new [=[  // Equal signed high words require an unsigned low-word comparison.
  const __m128i sign = _mm_set1_epi32(static_cast<int>(0x80000000u));
  Packet4i unsigned_lt = pcmp_lt<Packet4i>(Packet4i(_mm_xor_si128(a, sign)), Packet4i(_mm_xor_si128(b, sign)));
  Packet2l lo_lt = Packet2l(_mm_shuffle_epi32(unsigned_lt, (shuffle_mask<0, 0, 2, 2>::mask)));]=])
string(FIND "${_content}" "${_new}" _already_patched)
if(NOT _already_patched EQUAL -1)
    message(STATUS "Eigen SSE2 int64 comparison already patched")
    return()
endif()
string(FIND "${_content}" "${_old}" _match)
if(_match EQUAL -1)
    message(FATAL_ERROR "Cannot find the Eigen 5.0.1 SSE2 int64 comparison in ${_header}")
endif()
string(REPLACE "${_old}" "${_new}" _patched "${_content}")
file(WRITE "${_header}" "${_patched}")
message(STATUS "Patched Eigen SSE2 int64 comparison")
