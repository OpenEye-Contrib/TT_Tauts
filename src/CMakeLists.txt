cmake_minimum_required(VERSION 2.8.11)

project(TT_TAUTS)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

find_package(OEToolkits COMPONENTS oedepict oechem oesystem oeplatform)
#find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Widgets)
message( "Qt5Widgets_INCLUDE_DIRS : ${Qt5Widgets_INCLUDE_DIRS}" )
message( "Qt5Widgets_DEFINITIONS : ${Qt5Widgets_DEFINITIONS}" )
message( "Qt5Widgets_COMPILE_DEFINITIONS : ${Qt5Widgets_COMPILE_DEFINITIONS}" )
message( "Qt5Widgets_EXECUTABLE_COMPILE_FLAGS : ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall   -Wconversion -Wno-error=conversion -Wshadow -Wno-error=shadow -D_GLIBCXX_USE_CXX11_ABI=1 -std=c++11 -O3 -ffast-math -fomit-frame-pointer -funroll-loops -fvisibility=hidden -fvisibility-inlines-hidden -fPIC -Werror" )
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall -Wshadow -Wno-error=unused-function -fomit-frame-pointer -ffast-math -fvisibility=hidden -fvisibility-inlines-hidden -fPIC  -Werror " )
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine" )
else()
  message( "Not a Centos 6 machine. Hopefully it has a decent compiler." )
endif()

#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg" )
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64 -Wl,-z,defs -Wl,-rpath,'/'" )

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS program_options filesystem system date_time REQUIRED)

find_package(Cairo REQUIRED)
message( "CAIRO_INCLUDE_DIR : ${CAIRO_INCLUDE_DIR}" )

include_directories( SYSTEM ${OEToolkits_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
      ${Boost_INCLUDE_DIRS})

set(TAUT_STAND_SRCS TautStand.cc smirks_helper_fns.cc)
set(TAUT_STAND_INCS TautStand.H taut_enum_default_standardise_smirks.H
    taut_enum_default_vector_bindings.H)

#############################################################################
## tt_tauts - the test harness with the GUI
#############################################################################

set(TT_TAUTS_SRCS tt_tauts.cc
  TTTauts.cc
  TTTautsSettings.cc
  MolGridDisplay.cc
  make_taut_skeleton.cc
  TautomerGenerator.cc)

set(TT_TAUTS_INCS TTTauts.H
  TTTautsSettings.H
  MolGridDisplay.H
  TautomerGenerator.H)

set(TT_TAUTS_DACLIB_SRCS draw_oemol_to_qimage.cc
  QT4SelectItems.cc
  QTMolDisplay2D.cc
  apply_daylight_arom_model_to_oemol.cc
  build_time.cc
  check_oechem_licence.cc
  create_cansmi.cc
  extract_smarts_from_smirks.cc
  radical_atoms.cc
  split_smiles_into_atom_bits.cc)

set(TT_TAUTS_DACLIB_INCS DACOEMolAtomIndex.H
  DACOEMolBondIndex.H
  FileExceptions.H
  SMARTSExceptions.H
  QT4SelectItems.H
  QTMolDisplay2D.H
  stddefs.H)

set(LIBS ${LIBS}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CAIRO_LIBRARY})

set(SYSLIBS z pthread m stdc++)
set(EXECUTABLE_OUTPUT_PATH ${TT_TAUTS_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

add_executable(tt_tauts ${TT_TAUTS_SRCS} ${TT_TAUTS_DACLIB_SRCS}
    ${TAUT_STAND_SRCS}
    ${TT_TAUTS_INCS} ${TT_TAUTS_DACLIB_INCS}
    ${TAUT_STAND_INCS})
target_link_libraries(tt_tauts z ${LIBS} Qt5::Widgets ${SYSLIBS})

#############################################################################
## tt_tauts_batch - the test harness without the GUI
#############################################################################

set(TT_TAUTS_BATCH_SRCS tt_tauts_batch.cc
  make_taut_skeleton.cc
  TautomerGenerator.cc
  TTTautsSettings.cc
  TTTautsBatchSettings.cpp)

set(TT_TAUTS_BATCH_INCS TautomerGenerator.H
    TTTautsSettings.H
    TTTautsBatchSettings.h)

set(TT_TAUTS_BATCH_DACLIB_INCS DACOEMolAtomIndex.H
  DACOEMolBondIndex.H
  FileExceptions.H
  stddefs.H)

set(TT_TAUTS_BATCH_DACLIB_SRCS
  apply_daylight_arom_model_to_oemol.cc
  build_time.cc
  check_oechem_licence.cc
  create_cansmi.cc
  radical_atoms.cc)

add_executable(tt_tauts_batch ${TT_TAUTS_BATCH_SRCS}
    ${TAUT_STAND_SRCS} ${TT_TAUTS_BATCH_INCS}
  ${TT_TAUTS_BATCH_DACLIB_INCS} ${TT_TAUTS_BATCH_DACLIB_SRCS}
  ${TAUT_STAND_INCS})
target_link_libraries(tt_tauts_batch z ${LIBS} ${SYSLIBS})

#############################################################################
## gen_t_skel - generate tautomer skeleton (compound key) for molecules
#############################################################################

set(GEN_T_SKEL_SRCS gen_t_skel.cc
  make_taut_skeleton.cc
  TautomerGenerator.cc
  TTTautsSettings.cc
  TTTautsBatchSettings.cpp)

set(GEN_T_SKEL_BATCH_INCS TautomerGenerator.H
    TTTautsSettings.H
    TTTautsBatchSettings.h)

set(GEN_T_SKEL_DACLIB_INCS DACOEMolAtomIndex.H
  DACOEMolBondIndex.H
  FileExceptions.H
  stddefs.H)

set(GEN_T_SKEL_DACLIB_SRCS
  apply_daylight_arom_model_to_oemol.cc
  build_time.cc
  check_oechem_licence.cc
  create_cansmi.cc
  radical_atoms.cc)

add_executable(gen_t_skel ${GEN_T_SKEL_SRCS}
    ${TAUT_STAND_SRCS} ${GEN_T_SKEL_BATCH_INCS}
  ${GEN_T_SKEL_DACLIB_INCS} ${GEN_T_SKEL_DACLIB_SRCS}
  ${TAUT_STAND_INCS})
target_link_libraries(gen_t_skel z ${LIBS} ${SYSLIBS})

#############################################################################
## test_tt_tauts: regression tests
## uses catch.hpp from https://github.com/philsquared/Catch
## see there for acknowledgement and license (boost 1.0).
#############################################################################

set(TEST_TT_TAUTS_SRCS test_tt_tauts.cpp
  make_taut_skeleton.cc
  TautomerGenerator.cc)

set(TEST_TT_TAUTS_BATCH_INCS TautomerGenerator.H)

set(TEST_TT_TAUTS_DACLIB_INCS DACOEMolAtomIndex.H
  DACOEMolBondIndex.H
  FileExceptions.H
  stddefs.H)

set(TEST_TT_TAUTS_DACLIB_SRCS
  apply_daylight_arom_model_to_oemol.cc
  build_time.cc
  check_oechem_licence.cc
  create_cansmi.cc
  radical_atoms.cc)

add_executable(test_tt_tauts ${TEST_TT_TAUTS_SRCS}
    ${TAUT_STAND_SRCS} ${TT_TAUTS_BATCH_INCS}
  ${TEST_TT_TAUTS_DACLIB_INCS} ${TEST_TT_TAUTS_DACLIB_SRCS}
  ${TAUT_STAND_INCS} catch.hpp)
target_link_libraries(test_tt_tauts z ${LIBS} ${SYSLIBS})

