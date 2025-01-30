cmake_minimum_required(VERSION 3.20.0)
project(Bevel_And_Emboss_M CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(${PROJECT_NAME} SHARED)

set(LUA51_INCLUDE_DIR CACHE PATH "Lua51 include directory")
set(LUA51_LIBRARY_DIR CACHE PATH "Lua51 library directory")

target_sources(${PROJECT_NAME} PRIVATE Bevel_And_Emboss_M.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${LUA51_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/aviutl_exedit_sdk)
target_link_directories(${PROJECT_NAME} PRIVATE ${LUA51_LIBRARY_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE lua51)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    set_property(TARGET ${PROJECT_NAME} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    target_compile_options(${PROJECT_NAME} PRIVATE /source-charset:utf-8 /execution-charset:shift_jis)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_link_options(${PROJECT_NAME} PRIVATE -static-libstdc++ -static)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")