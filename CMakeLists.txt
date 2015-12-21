project( CMC_Ising )

cmake_minimum_required(VERSION 2.8.8)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.8)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.8)
    cmake_policy(VERSION 2.8.8)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
message( ${CMAKE_CXX_COMPILER} )

set( CMAKE_CXX_FLAGS "-O3 -s -DBORDER_TYPE=0 -std=gnu++11 -fvisibility=default")
set( COMPLEX CACHE BOOLEAN false )

if( COMPLEX )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOMPLEX_PARAMS")
endif()

if ( COMMAND cmake_policy )

cmake_policy( SET CMP0003 NEW )

endif()

# Boost and its components
set ( Boost_GCC_INCLUDE_DIRS /Users/breubreubreu/Execs/boost_1_60_0)
include_directories(${Boost_GCC_INCLUDE_DIRS})
message("Using Boost found at ${Boost_GCC_INCLUDE_DIRS}")
message("Using the following linker command")

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
add_library( COMMON_OBJS OBJECT ${COMMON_SRC})
set ( LIB_DIR ${Boost_GCC_INCLUDE_DIRS}/stage/lib)
set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem.a ${LIB_DIR}/libboost_serialization.a)

if( COMPLEX )
	add_executable( CMC_Mixed_CMPLX $<TARGET_OBJECTS:COMMON_OBJS> ${CMC_MIXED_SRC})
	add_executable( Merge_CMPLX $<TARGET_OBJECTS:COMMON_OBJS> ${MERGE_SRC})
	add_executable( Set_Manual $<TARGET_OBJECTS:COMMON_OBJS> ${SET_MANUAL_SRC})
	add_executable( Set_Moessner_CMPLX $<TARGET_OBJECTS:COMMON_OBJS> ${SET_MOESSNER_SRC})
	add_executable( Set_Part $<TARGET_OBJECTS:COMMON_OBJS> ${SET_PART_SRC})
else()
	add_executable( CMC_Mixed $<TARGET_OBJECTS:COMMON_OBJS> ${CMC_MIXED_SRC})
	add_executable( Merge $<TARGET_OBJECTS:COMMON_OBJS> ${MERGE_SRC})
	add_executable( Set_Manual $<TARGET_OBJECTS:COMMON_OBJS> ${SET_MANUAL_SRC})
	add_executable( Set_Moessner $<TARGET_OBJECTS:COMMON_OBJS> ${SET_MOESSNER_SRC})
	add_executable( Set_Part $<TARGET_OBJECTS:COMMON_OBJS> ${SET_PART_SRC})
endif()

# Link the executable to CGAL and third-party libraries
if( COMPLEX )
	target_link_libraries( CMC_Mixed_CMPLX ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Merge_CMPLX ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Manual ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Moessner_CMPLX ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Part ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
else()
	target_link_libraries( CMC_Mixed ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Merge ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Manual ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Moessner ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
	target_link_libraries( Set_Part ${BOOST_LIBRARIES} ${Boost_LIBRARIES})
endif()