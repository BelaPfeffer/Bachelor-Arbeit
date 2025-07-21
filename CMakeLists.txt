

cmake_minimum_required(VERSION 3.20...3.31)
project(seqan3_tutorial CXX)

# Pfad zu SeqAn3 (hier relativ zum Parent-Ordner, passe ggf. an)
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/cmake")

# SeqAn3 mit Versionsanforderung finden
find_package(seqan3 3.0 REQUIRED)

# Executable mit deinen Quellcodes
add_executable(my_app
    src/main.cpp
    src/compressV2.cpp
    src/suffix_array.cpp
)

# Include-Verzeichnis für Header-Dateien
target_include_directories(my_app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# SeqAn3 als Link-Library hinzufügen
target_link_libraries(my_app seqan3::seqan3)
