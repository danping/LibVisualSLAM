project(libVisualSLAM)
cmake_minimum_required(VERSION 2.8)

#The version number
set (libVisualSLAM_VERSION_MAJOR 1)
set (libVisualSLAM_VERSION_MINOR 0)


aux_source_directory(src SRC)
aux_source_directory(src/calibration SRC_CALIBRATION)
aux_source_directory(src/extern SRC_EXTERN)
aux_source_directory(src/extern/sba-1.6 SRC_SBA)
aux_source_directory(src/extern/csparse SRC_CSPARSE)

FILE(GLOB SRC_LEVMAR "src/extern/levmar/*.[hc]")
FILE(GLOB levmar_cores "src/extern/levmar/*_core.c")
list(REMOVE_ITEM SRC_LEVMAR ${levmar_cores})

aux_source_directory(src/geometry SRC_GEOMETRY)
aux_source_directory(src/imgproc SRC_IMGPROC)
aux_source_directory(src/lap SRC_LAP)
aux_source_directory(src/matching SRC_MATCHING)
aux_source_directory(src/math SRC_MATH)
aux_source_directory(src/tools SRC_TOOLS)
aux_source_directory(src/util SRC_UTIL)


include_directories(src)

add_library(VisualSLAM ${SRC} ${SRC_CALIBRATION} ${SRC_EXTERN} ${SRC_SBA} ${SRC_CSPARSE} ${SRC_LEVMAR} ${SRC_GEOMETRY} ${SRC_IMGPROC} ${SRC_LAP} ${SRC_MATCHING} ${SRC_MATH} ${SRC_TOOLS} ${SRC_UTIL})

install (TARGETS VisualSLAM DESTINATION lib)
FILE(GLOB headers_calibration "src/calibration/*.h")
FILE(GLOB headers_extern "src/extern/*.h")
FILE(GLOB headers_sba "src/extern/sba-1.6/*.h")
FILE(GLOB headers_levmar "src/extern/levmar/*.h")
FILE(GLOB headers_geometry "src/geometry/*.h")
FILE(GLOB headers_imgproc "src/imgproc/*.h")
FILE(GLOB headers_lap "src/lap/*.h")
FILE(GLOB headers_matching "src/matching/*.h")
FILE(GLOB headers_math "src/math/*.h")
FILE(GLOB headers_tools "src/tools/*.h")
FILE(GLOB headers_util "src/util/*.h")
FILE(GLOB headers "src/*.h")

install(FILES 
    ${headers_calibration}
    DESTINATION
    include/VisualSLAM/calibration
)

install(FILES 
    ${headers_extern}
    DESTINATION
    include/VisualSLAM/extern
)

install(FILES 
    ${headers_sba}
    DESTINATION
    include/VisualSLAM/extern/sba-1.6
)

install(FILES
    ${headers_geometry} 
    DESTINATION
    include/VisualSLAM/geometry
)
install(FILES
    ${headers_imgproc} 
    DESTINATION
    include/VisualSLAM/imgproc
)
install(FILES
    ${headers_lap} 
    DESTINATION
    include/VisualSLAM/lap
)

install(FILES
    ${headers_matching} 
    DESTINATION
    include/VisualSLAM/matching
)

install(FILES
    ${headers_math}
    DESTINATION
    include/VisualSLAM/math
)

install(FILES
    ${headers_tools}
    DESTINATION
    include/VisualSLAM/tools
)

install(FILES
    ${headers_util}
    DESTINATION
    include/VisualSLAM/util
)

install(FILES
    ${headers}
    DESTINATION
    include/VisualSLAM
)


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)