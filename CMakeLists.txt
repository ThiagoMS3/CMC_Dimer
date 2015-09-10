project( CMC_Ising )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

set( CMAKE_CXX_FLAGS "-O0 -g -DBORDER_TYPE=0 -std=c++11")

if ( COMMAND cmake_policy )

cmake_policy( SET CMP0003 NEW )

endif()

# Boost and its components
find_package( Boost COMPONENTS serialization REQUIRED )

if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

# inclued for Boost directory
include_directories(${Boost_INCLUDE_DIRS})

# set directories
set ( SOURCE_DIR src)
set ( COMMON_DIR ${SOURCE_DIR}/Common_MC)
set ( CMC_MIXED_DIR ${SOURCE_DIR}/CMC_Mixed)
set ( MERGE_DIR ${SOURCE_DIR}/Merge)
set ( SET_MANUAL_DIR ${SOURCE_DIR}/CMC_Manual)
set ( SET_MOESSNER_DIR ${SOURCE_DIR}/CMC_Moessner)
set ( SET_PART_DIR ${SOURCE_DIR}/CMC_Part)

# include for local directory
include_directories( ${SOURCE_DIR} ${COMMON_DIR} ${CMC_MIXED_DIR} ${MERGE_DIR} ${SET_MANUAL_DIR} ${SET_MOESSNER_DIR} ${SET_PART_DIR})

# set source sets
file ( GLOB COMMON_SRC "${COMMON_DIR}/*.h" "${COMMON_DIR}/*.cpp") 
file ( GLOB CMC_MIXED_SRC "${CMC_MIXED_DIR}/*.h" "${CMC_MIXED_DIR}/*.cpp")
file ( GLOB MERGE_SRC "${MERGE_DIR}/*.h" "${MERGE_DIR}/*.cpp")
file ( GLOB SET_MANUAL_SRC "${MERGE_DIR}/*.h" "${MERGE_DIR}/*.cpp")
file ( GLOB SET_MOESSNER_SRC "${SET_MOESSNER_DIR}/*.h" "${SET_MOESSNER_DIR}/*.cpp")
file ( GLOB SET_PART_SRC "${SET_PART_DIR}/*.h" "${SET_PART_DIR}/*.cpp")

# Creating entries for target: intersections_2D
# ############################

add_executable( CMC_Mixed ${COMMON_SRC} ${CMC_MIXED_SRC})
add_executable( Merge ${COMMON_SRC} ${MERGE_SRC})
add_executable( Set_Manual ${COMMON_SRC} ${SET_MANUAL_SRC})
add_executable( Set_Moessner ${COMMON_SRC} ${SET_MOESSNER_SRC})
add_executable( Set_Part ${COMMON_SRC} ${SET_PART_SRC})

# set ( LIB_DIR /opt/local/lib)
# set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib ${LIB_DIR}/libboost_serialization-mt.dylib)

# Link the executable to CGAL and third-party libraries
target_link_libraries( CMC_Mixed ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries( Merge ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries( Set_Manual ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries( Set_Moessner ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries( Set_Part ${BOOST_LIBRARIES} ${Boost_LIBRARIES})

