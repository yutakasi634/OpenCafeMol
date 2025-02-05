#include "src/input/ReadGenesisInput.cpp"

#include <gtest/gtest.h>

TEST(ReadGenesisInput, PreprocessTopFile)
{
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

TEST(ReadGenesisInput, ReadTopFile)
{
    std::map<std::string, std::vector<std::string>> top_contents =
        read_top_file("test_read_top_file_input.top", "input/");

    EXPECT_EQ("first elem of first section", top_contents.at("first section")[0]);
    EXPECT_EQ("second elem of second section", top_contents.at("second section")[1]);
    EXPECT_EQ("third elem of third section", top_contents.at("third section")[2]);
}
