#include "src/input/ReadGenesisInput.cpp"

#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>

TEST(ReadGenesisInput, PreprocessTopFile)
{
    std::cout << "current path: " << std::filesystem::current_path() << std::endl;

    const std::vector<std::string> file_contents =
        preprocess_top_file("test_preprocess_top_file_input.top", "input/");

    EXPECT_EQ("first element in test_preprocess_top_file_input.top",
              file_contents[0]);
    EXPECT_EQ("second element in test_preprocess_top_file_input.top with comment",
              file_contents[1]);
    EXPECT_EQ("included.inp was correctly included",
              file_contents[2]);
    EXPECT_EQ("included_with_comment.inp was correctly included",
              file_contents[3]);
}
