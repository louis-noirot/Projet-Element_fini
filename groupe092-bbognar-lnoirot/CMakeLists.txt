cmake_minimum_required(VERSION 3.6)
project(myFemProject)

set(EXEC_NAME myFem)

# Set GMSH Version : example gmsh-4.12.2 -> GMSH_VERSION=gmsh-4.12, GMSH_MINOR=2
set(GMSH_VERSION gmsh-4.13)
set(GMSH_MINOR 1)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLFW
add_subdirectory(glfw)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/glfw/include)

# GMSH
if(MSVC OR MSYS OR MINGW)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/lib/gmsh.lib)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/include)
    #Copy .dll into Debug
    set(GMSH_DLL ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/lib/${GMSH_VERSION}.dll)
    if(MSVC)
        file(COPY ${GMSH_DLL} DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE})
    else()
        file(COPY ${GMSH_DLL} DESTINATION ${CMAKE_BINARY_DIR})
    endif()
endif(MSVC OR MSYS OR MINGW)

if(APPLE)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSARM-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSARM-sdk/include)
    else()    
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSX-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSX-sdk/include)
    endif()
endif(APPLE)

if(UNIX AND NOT APPLE)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Linux64-sdk/include)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Linux64-sdk/lib/libgmsh.so)
endif(UNIX AND NOT APPLE)

message("GMSH_INCLUDE_DIR : ${GMSH_INCLUDE_DIR}")
message("GMSH_LIBRARY : ${GMSH_LIBRARY}")
include_directories(${GMSH_INCLUDE_DIR})

add_executable(${EXEC_NAME} src/geo_fem.c src/geo_glfem.c src/geometry.c src/main.c)
target_link_libraries(${EXEC_NAME} glfw ${OPENGL_gl_LIBRARY} ${GMSH_LIBRARY})

# === Autres exécutables pour tests (non graphiques)
add_executable(test1 src/geo_fem.c src/geo_glfem.c src/geometry.c src/main.c)
add_executable(test2 src/slv_fem.c src/slv_glfem.c src/solver.c src/main.c)
add_executable(test3 src/elast_fem.c src/elast_glfem.c src/elasticity.c src/main.c)

# Liens (même dépendances)
target_link_libraries(test1 glfw ${OPENGL_gl_LIBRARY} ${GMSH_LIBRARY})
target_link_libraries(test2 glfw ${OPENGL_gl_LIBRARY} ${GMSH_LIBRARY})
target_link_libraries(test3 glfw ${OPENGL_gl_LIBRARY} ${GMSH_LIBRARY})

# Définitions de compilation spécifiques à chaque test
target_compile_definitions(test1 PRIVATE TEST1)
target_compile_definitions(test2 PRIVATE TEST2)
target_compile_definitions(test3 PRIVATE TEST3)

# Exécutable principal reste graphique
target_compile_definitions(${EXEC_NAME} PRIVATE graphic)

# Commandes pratiques
add_custom_target(run1 COMMAND ./test1 DEPENDS test1 COMMENT "Run TEST1 (mesh generation)")
add_custom_target(run2 COMMAND ./test2 DEPENDS test2 COMMENT "Run TEST2 (band solver and remesh)")
add_custom_target(run3 COMMAND ./test3 DEPENDS test3 COMMENT "Run TEST3 (elasticity solver)")

if(APPLE)
    target_compile_options(${EXEC_NAME} PRIVATE "-Wdeprecated-declarations")
    target_compile_definitions(${EXEC_NAME} PRIVATE -D GL_SILENCE_DEPRECATION)
endif()

# set exec as the startup project in visual studio
if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC_NAME})
    target_compile_options(${EXEC_NAME} PRIVATE "/W2")
endif()

target_compile_definitions(${EXEC_NAME} PRIVATE graphic)

# add a custom target so we can do make run directly
add_custom_target(run COMMAND ${EXEC_NAME} DEPENDS ${EXEC_NAME} COMMENT "./${EXEC_NAME}")

