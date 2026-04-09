cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(DispFromEcho)

include_directories(K:/lib-x64-msvc2022/eigen)

add_executable(DispFromEcho MACOSX_BUNDLE 	DispFromEcho.cpp DispFromEcho.h 
											Correlation.cpp Correlation.h
											EchoSequence.cpp EchoSequence.h
											EchoPair.cpp EchoPair.h
											RigidTransform.cpp RigidTransform.h
											Calibration.cpp Calibration.h
											CImg.h)
											
target_compile_features(DispFromEcho PRIVATE cxx_std_17)

#target_link_libraries(Cylinder PRIVATE ${VTK_LIBRARIES})
