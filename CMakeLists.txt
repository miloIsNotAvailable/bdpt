cmake_minimum_required(VERSION 3.10.0)
project(shithead VERSION 0.1.0 LANGUAGES C CXX)

# include_directories( headers )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(shithead 
                main.cpp 
                ./build/_deps/glad/src/glad.c
                ./src/loadObj.cpp
                ./src/intersection.cpp )

include(CTest)
enable_testing()

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")

include(FetchContent)

FetchContent_Declare(
  tinyobjloader
  GIT_REPOSITORY https://github.com/tinyobjloader/tinyobjloader.git
  GIT_TAG release  
)

FetchContent_Declare(
  glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG 3.3-stable
)

FetchContent_Declare(
	glm
	GIT_REPOSITORY	https://github.com/g-truc/glm.git
	GIT_TAG 	bf71a834948186f4097caa076cd2663c69a10e1e #refs/tags/1.0.1
)

FetchContent_MakeAvailable(glm)

FetchContent_MakeAvailable(tinyobjloader)
FetchContent_MakeAvailable(glfw)

# https://github.com/Dav1dde/glad/issues/186
add_subdirectory(build/_deps/glad/)

# remember to add the fucking opengl32 you fuck, its obvious dumbass
# how do you not know it , no docs fucking needed ye, fuckwad dumbass piece of shit 
# I fucking hate this, there has got to be a better fucking to fucking do this
# too fucking bad I guess ba ham bag bitch 
target_link_libraries(shithead PRIVATE glfw glad opengl32 glm::glm)

# include headers file in cmake
target_include_directories( shithead PRIVATE 
                            ${CMAKE_CURRENT_SOURCE_DIR}/headers 
                            ${CMAKE_CURRENT_SOURCE_DIR}/src 
                            ${CMAKE_CURRENT_SOURCE_DIR}/build/_deps/tinyobjloader-src )

