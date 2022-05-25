find_path(OPENMM_INCLUDE_DIR
    NAMES
        OpenMM.h
    PATHS
        ${OPENMM_ROOT}
        /usr/openmm
        /usr/local/openmm
    PATH_SUFFIXES
        include
    # to ignore detect pyenv openmm prior to OPENMM_ROOT
    NO_DEFAULT_PATH
    )
message(STATUS "OpenMM include path : ${OPENMM_INCLUDE_DIR}")

find_library(OPENMM_LIBRARY
    NAMES
        libOpenMM.so
    PATHS
        ${OPENMM_ROOT}
        /usr/openmm
        /usr/local/openmm
    PATH_SUFFIXES
        lib
    # to ignore detect pyenv openmm prior to OPENMM_ROOT
    NO_DEFAULT_PATH
    )
message(STATUS "OpenMM library path : ${OPENMM_LIBRARY}")

mark_as_advanced(OPENMM_INCLUDE_DIR OPENMM_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMM
    REQUIRED_VARS
        OPENMM_INCLUDE_DIR
        OPENMM_LIBRARY
    )

if(OPENMM_FOUND AND NOT TARGET OpenMM::OpenMM)
    add_library(OpenMM::OpenMM UNKNOWN IMPORTED)
    set_target_properties(OpenMM::OpenMM PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION             "${OPENMM_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENMM_INCLUDE_DIR}"
        )
endif()

