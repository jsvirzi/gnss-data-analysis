cmake_minimum_required(VERSION 3.6)
project(IplDataAnalysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -m64 -I/home/ubuntu/utils/root/include")

set(SOURCE_FILES main.cpp)
add_executable(IplDataAnalysis ${SOURCE_FILES})

target_link_libraries(IplDataAnalysis -L/home/ubuntu/utils/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic)