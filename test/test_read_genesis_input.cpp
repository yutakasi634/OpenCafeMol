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

TEST(ReadGenesisInput, ParseGroTopFile)
{
    std::map<std::string, std::vector<std::string>> top_contents =
        parse_grotop_file("test_read_top_file_input.top", "input/");

    EXPECT_EQ("first elem of first section", top_contents.at("first section")[0]);
    EXPECT_EQ("second elem of second section", top_contents.at("second section")[1]);
    EXPECT_EQ("third elem of third section", top_contents.at("third section")[2]);
}

TEST(ReadGenesisInput, ReadInpFile)
{
    std::map<std::string, std::map<std::string, std::string>> inp_contents =
        read_inp_file("input/test_read_inp_file_input.inp");

    EXPECT_EQ("first_value",  inp_contents.at("FIRST_SECTION") .at("first_key"));
    EXPECT_EQ("second_value_in_second_section",
            inp_contents.at("SECOND_SECTION").at("second_key"));
    EXPECT_EQ("third_value_in_third_section",
            inp_contents.at("THIRD_SECTION") .at("third_key"));
}
