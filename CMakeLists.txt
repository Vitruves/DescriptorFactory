cmake_minimum_required(VERSION 3.16)
project(desfact VERSION 0.1.0 LANGUAGES CXX C)


# PGO Configuration - now optional
option(ENABLE_PGO "Enable Profile-Guided Optimization" OFF)
if(ENABLE_PGO AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate=.")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-generate=.")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-generate=.")
    message(STATUS "PGO: Enabled profile generation with -fprofile-generate=.")
else()
    message(STATUS "PGO: Disabled")
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Option to build GUI ---
option(BUILD_GUI "Build the desfact-gui executable" ON)
# ---------------------------

if(BUILD_GUI)
    # Enable Qt features only if building GUI
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)

    # Find Qt (Try Qt6 first, fallback to Qt5) - Make required only if building GUI
    find_package(Qt6 COMPONENTS Widgets QUIET)
    if(NOT Qt6_FOUND)
        find_package(Qt5 COMPONENTS Widgets REQUIRED) # REQUIRED only if BUILD_GUI
        set(QT_VERSION 5)
        set(QT_WIDGETS_LIB Qt5::Widgets)
        message(STATUS "Using Qt 5 for GUI")
    else()
        set(QT_VERSION 6)
        set(QT_WIDGETS_LIB Qt6::Widgets)
        message(STATUS "Using Qt 6 for GUI")
    endif()
endif()

# TBB
find_package(TBB QUIET)
if(TBB_FOUND)
    add_compile_definitions(WITH_TBB)
    message(STATUS "Found TBB via find_package: ${TBB_VERSION}")
else()
    find_library(TBB_LIBRARY NAMES tbb)
    find_path(TBB_INCLUDE_DIR NAMES tbb/tbb.h)
    if(TBB_LIBRARY AND TBB_INCLUDE_DIR)
        set(TBB_FOUND TRUE)
        add_compile_definitions(WITH_TBB)
        message(STATUS "Found TBB library manually: ${TBB_LIBRARY}")
        include_directories(${TBB_INCLUDE_DIR})
    else()
        message(STATUS "TBB not found, will use std::async for parallelization")
    endif()
endif()

# RPATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/usr/local/lib:/usr/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Find RDKit package
find_package(RDKit REQUIRED)

# Set RDKit paths directly since target properties are not available
set(RDKIT_INCLUDE_DIR "/usr/include/rdkit/")
set(RDKIT_LIBRARY_DIR "/lib/")

# Add RDKit include directory
include_directories(${RDKIT_INCLUDE_DIR})

# Find RDKit libraries using find_library with proper paths
find_library(RDKIT_GRAPHMOL_LIB
    NAMES GraphMol RDKitGraphMol
    PATHS ${RDKIT_LIBRARY_DIR}
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_SMILEPARSE_LIB
    NAMES SmilesParse RDKitSmilesParse
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_FILEPARSERS_LIB
    NAMES FileParsers RDKitFileParsers
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_CHEMTRANSFORMS_LIB
    NAMES ChemTransforms RDKitChemTransforms
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_DESCRIPTORS_LIB
    NAMES Descriptors RDKitDescriptors
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_RDGENERAL_LIB
    NAMES RDGeneral RDKitRDGeneral
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_MOLSTANDARDIZE_LIB
    NAMES MolStandardize RDKitMolStandardize
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_SUBSTRUCTMATCH_LIB
    NAMES SubstructMatch RDKitSubstructMatch
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_RDGEOMETRYLIB_LIB
    NAMES RDGeometryLib RDKitRDGeometryLib
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_DATASTRUCTS_LIB
    NAMES DataStructs RDKitDataStructs
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_FINGERPRINTS_LIB
    NAMES Fingerprints RDKitFingerprints
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
find_library(RDKIT_DEPICTOR_LIB
    NAMES Depictor RDKitDepictor
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)
# Add MolStandardize library explicitly
find_library(RDKIT_MOLSTANDARDIZE_LIB
    NAMES MolStandardize RDKitMolStandardize
    PATHS ${RDKIT_LIBRARY_DIR} /lib/ /usr/local/lib
    PATH_SUFFIXES rdkit
    NO_DEFAULT_PATH
)


# Verify all required libraries were found
foreach(lib
    RDKIT_GRAPHMOL_LIB RDKIT_SMILEPARSE_LIB RDKIT_FILEPARSERS_LIB
    RDKIT_CHEMTRANSFORMS_LIB RDKIT_DESCRIPTORS_LIB RDKIT_RDGENERAL_LIB
    RDKIT_MOLSTANDARDIZE_LIB RDKIT_SUBSTRUCTMATCH_LIB RDKIT_RDGEOMETRYLIB_LIB
    RDKIT_DATASTRUCTS_LIB RDKIT_FINGERPRINTS_LIB RDKIT_DEPICTOR_LIB RDKIT_MOLSTANDARDIZE_LIB)
    if(NOT ${lib})
        message(FATAL_ERROR "Required RDKit library ${lib} not found")
    endif()
endforeach()

list(APPEND CMAKE_INSTALL_RPATH "${RDKIT_LIBRARY_DIR}")

# Add Eigen dependency with multiple search paths
find_package(Eigen3 QUIET)

if(NOT EIGEN3_FOUND)
    # Try pkg-config as a fallback
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(EIGEN3 QUIET eigen3)
        if(EIGEN3_FOUND)
            set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS})
            message(STATUS "Found Eigen3 via pkg-config: ${EIGEN3_INCLUDE_DIR}")
        endif()
    endif()
endif()

# Look in common installation directories if still not found
if(NOT EIGEN3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Core Eigen/Dense
        PATHS
        /usr/include/eigen3
        /usr/local/include/eigen3
        /opt/local/include/eigen3
        /usr/include
        /usr/local/include
        /opt/local/include
    )
    if(EIGEN3_INCLUDE_DIR)
        set(EIGEN3_FOUND TRUE)
        message(STATUS "Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
    endif()
endif()

# Final check and inclusion
if(NOT EIGEN3_FOUND)
    message(FATAL_ERROR "Eigen3 not found. Please install Eigen3 library.")
else()
    include_directories(${EIGEN3_INCLUDE_DIR})
    message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")
endif()

# CUDA (Optional, kept for future potential use)
option(WITH_CUDA "Build with CUDA support" OFF)
if(WITH_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DWITH_CUDA)
endif()

set(RDKIT_LIBRARIES
    ${RDKIT_GRAPHMOL_LIB}
    ${RDKIT_SMILEPARSE_LIB}
    ${RDKIT_FILEPARSERS_LIB}
    ${RDKIT_CHEMTRANSFORMS_LIB}
    ${RDKIT_DESCRIPTORS_LIB}
    ${RDKIT_RDGENERAL_LIB}
    ${RDKIT_MOLSTANDARDIZE_LIB}
    ${RDKIT_SUBSTRUCTMATCH_LIB}
    ${RDKIT_RDGEOMETRYLIB_LIB}
    ${RDKIT_DATASTRUCTS_LIB}
    ${RDKIT_FINGERPRINTS_LIB}
    ${RDKIT_DEPICTOR_LIB}
    ${RDKIT_MOLSTANDARDIZE_LIB} # Ensure MolStandardize is linked
)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${RDKIT_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

link_directories(
    ${RDKIT_LIBRARY_DIR}
)

# Library source list
set(LIB_SOURCES
    src/utils.cpp
    src/io.cpp
    src/descriptors.cpp
    src/descriptors/fractional.cpp
    src/descriptors/sum.cpp
    src/descriptors/strings.cpp
    src/descriptors/vague3.cpp
    src/descriptors/vague4.cpp
    src/descriptors/electronic.cpp
    src/descriptors/vague5.cpp
    src/descriptors/vague6.cpp
    src/descriptors/vague7.cpp
    src/descriptors/selfies.cpp
    src/descriptors/pka.cpp
    src/descriptors/vague8.cpp
    src/descriptors/morgan.cpp
    src/descriptors/eigen.cpp
    src/descriptors/solubility.cpp
    src/descriptors/counts.cpp
    src/descriptors/image.cpp
)

if(WITH_CUDA)
    set(CUDA_SOURCES src/cuda/descriptors_cuda.cu)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3 -arch=sm_60)
    list(APPEND LIB_SOURCES ${CUDA_SOURCES})
endif()

# Shared library
add_library(desfact SHARED ${LIB_SOURCES})
target_link_libraries(desfact
    ${RDKIT_LIBRARIES}
)

if(WITH_CUDA)
    target_link_libraries(desfact ${CUDA_LIBRARIES})
endif()

if(TBB_FOUND)
    if(TARGET TBB::tbb)
        target_link_libraries(desfact TBB::tbb)
    else()
        target_link_libraries(desfact ${TBB_LIBRARY})
    endif()
endif()

# Main CLI executable
add_executable(desfact_main src/main.cpp)
set_target_properties(desfact_main PROPERTIES
    OUTPUT_NAME "desfact" # Keep executable name as 'desfact'
    BUILD_WITH_INSTALL_RPATH FALSE
    SKIP_BUILD_RPATH FALSE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/usr/local/lib:/usr/lib:${RDKIT_LIBRARY_DIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)
target_link_libraries(desfact_main
    desfact
    ${RDKIT_LIBRARIES}
)

if(TBB_FOUND)
    if(TARGET TBB::tbb)
        target_link_libraries(desfact_main TBB::tbb)
    else()
        target_link_libraries(desfact_main ${TBB_LIBRARY})
    endif()
endif()

# --- GUI Executable (Conditional Build) ---
if(BUILD_GUI)
    message(STATUS "Building desfact-gui")
    set(GUI_SOURCES
        src/gui/gui_main.cpp
        src/gui/mainwindow.cpp
        src/gui/mainwindow.h  # Explicitly include the header
    )

    add_executable(desfact_gui WIN32 MACOSX_BUNDLE ${GUI_SOURCES})
    set_target_properties(desfact_gui PROPERTIES
        OUTPUT_NAME "desfact-gui"
        BUILD_WITH_INSTALL_RPATH FALSE
        SKIP_BUILD_RPATH FALSE
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/usr/local/lib:/usr/lib:${RDKIT_LIBRARY_DIR}"
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    target_link_libraries(desfact_gui
        desfact              # Link against our core library
        ${QT_WIDGETS_LIB}    # Link against Qt Widgets
        ${RDKIT_LIBRARIES}   # Link against RDKit
    )

    if(TBB_FOUND)
        if(TARGET TBB::tbb)
            target_link_libraries(desfact_gui TBB::tbb)
        else()
            target_link_libraries(desfact_gui ${TBB_LIBRARY})
        endif()
    endif()
else()
    message(STATUS "Skipping desfact-gui build (BUILD_GUI=OFF)")
endif() # End BUILD_GUI block
# -----------------------------------------


# Removed csvtool target assuming main.cpp covers its functionality

# Install targets - make GUI installation conditional
install(TARGETS desfact desfact_main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)

if(BUILD_GUI)
    install(TARGETS desfact_gui
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib)
endif()

install(DIRECTORY include/ DESTINATION include/desfact)

# Ensure RDKit libraries are found even if not in standard paths
list(APPEND CMAKE_PREFIX_PATH "/usr/local")