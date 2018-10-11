cmake_minimum_required( VERSION 2.8.9 )

project(RDKitTest)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wno-deprecated -Wno-unused-function -fno-strict-aliasing -Wno-format -Wno-logical-op-parentheses -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)

find_package(MyRDKit COMPONENTS SmilesParse RDGeneral GraphMol Subgraphs FileParsers)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)


include_directories(${RDKIT_INCLUDE_DIR})
  message( "Example : ${exnum}" )
  set( exe "test" )
  set( src "${exe}.cpp" sha1.cpp sha1.hpp )
  add_executable( ${exe} ${src} )
  target_link_libraries( ${exe} ${RDKIT_LIBRARIES} )

# add_executable( test_example test_example.cpp )
# target_link_libraries( test_example ${LIBS} )
