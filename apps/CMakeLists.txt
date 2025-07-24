add_executable(my_app main.cpp)

# Link your code and SeqAn3
target_link_libraries(my_app PRIVATE compressionTools seqan3::seqan3)
