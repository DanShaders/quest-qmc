#include <gtest/gtest.h>

#include "SRC/config.h"

struct Test {
    std::string_view config;
    std::map<std::string, dqmc::ParameterValue, std::less<>> parameters = {};
    std::optional<std::string> error;
};

std::vector<Test> tests = {
    // Simple
    {
        .config = R"(a = 1
b = 1.0
)",
        .parameters = {
            { "a", { "1", 1 } },
            { "b", { "1.0", 2 } },
        },
    },
    // State::BeforeKey
    {
        .config = "     # This is a comment!",
        .parameters = {},
    },
    {
        .config = "=",
        .error = "Unexpected '=' at 1:1 in the configuration file",
    },
    {
        .config = "",
        .parameters = {},
    },
    // State::ReadingKey
    {
        .config = "a# This is not allowed",
        .error = "Unexpected '#' at 1:2 in the configuration file",
    },
    {
        .config = "abcdef = 123",
        .parameters = { { "abcdef", { "123", 1 } } },
    },
    {
        .config = "abcdef= 123",
        .parameters = { { "abcdef", { "123", 1 } } },
    },
    {
        .config = "abcdef=123",
        .parameters = { { "abcdef", { "123", 1 } } },
    },
    {
        .config = "abcdef\n",
        .error = "Unexpected end of line at 1:7 in the configuration file",
    },
    // State::BeforeEquals
    {
        .config = "abcdef      \n",
        .error = "Unexpected end of line at 1:13 in the configuration file",
    },
    {
        .config = "abcdef   #",
        .error = "Unexpected '#' at 1:10 in the configuration file",
    },
    // State::BeforeValue
    {
        .config = "abcdef =   ",
        .error = "Unexpected end of line at 1:12 in the configuration file",
    },
    {
        .config = "abcdef = \"123\"",
        .parameters = { { "abcdef", { "123", 1 } } },
    },
    {
        .config = "abcdef = #",
        .error = "Unexpected '#' at 1:10 in the configuration file",
    },
    // State::ReadingUnquotedValue
    {
        .config = "key = value # and comment",
        .parameters = { { "key", { "value", 1 } } },
    },
    {
        .config = "key = value with spaces  ",
        .parameters = { { "key", { "value with spaces", 1 } } },
    },
    {
        .config = "key =    value with\tspaces and a comment#hello",
        .parameters = { { "key", { "value with\tspaces and a comment", 1 } } },
    },
    {
        .config = "key = \\n\\r ",
        .parameters = { { "key", { "\\n\\r", 1 } } },
    },
    // State::ReadingQuotedValue
    {
        .config = "key = \"1234\"  ",
        .parameters = { { "key", { "1234", 1 } } },
    },
    {
        .config = "key = \"123 # comment",
        .error = "Unexpected end of line at 1:21 in the configuration file",
    },
    {
        .config = "key = \"123",
        .error = "Unexpected end of line at 1:11 in the configuration file",
    },
    {
        .config = R"(key = "1\n2 \r3\t 4\\ 5\" ")",
        .parameters = { { "key", { "1\n2 \r3\t 4\\ 5\" ", 1 } } },
    },
    {
        .config = R"(key = "\)",
        .error = "Unexpected end of line at 1:9 in the configuration file",
    },
    {
        .config = R"(key = "\a")",
        .error = "Unexpected escape sequence '\\a' at 1:9 in the configuration file",
    },
    // State::AfterQuotedValue
    {
        .config = R"(key = "" a)",
        .error = "Unexpected 'a' at 1:10 in the configuration file",
    },
    {
        .config = R"(key = "1" #a)",
        .parameters = { { "key", { "1", 1 } } },
    },
    // Duplicated parameters
    {
        .config = "a = 1\na=1\n",
        .error = "Parameter 'a' has duplicate definitions at lines 1 and 2 of the configuration file",
    },
    {
        .config = "a = 1\nA=1\n",
        .error = "Parameter 'A' has duplicate definitions at lines 1 and 2 of the configuration file",
    },
    // Normalization
    {
        .config = "nOrmaliZEme = BUT NOT me\nhUH? = \"HI\"\n",
        .parameters = {
            { "normalizeme", { "BUT NOT me", 1 } },
            { "huh?", { "HI", 2 } },
        },
    },
    // Real-world example
    {
        .config = R"(
# ==========================
# Output file name header
# ==========================

ofile  = small

# ==========================
# lattice dimension
# ==========================

nx     = 8
ny     = 8
)",
        .parameters = {
            { "ofile", { "small", 6 } },
            { "nx", { "8", 12 } },
            { "ny", { "8", 13 } },
        },
    },
};

TEST(ConfigParameters, create)
{
    for (auto const& test : tests) {
        std::cerr << "Checking |" << test.config << "|...\n";
        std::istringstream stream { std::string(test.config) };
        auto result = dqmc::ConfigParameters::create(stream);

        if (test.error.has_value()) {
            EXPECT_FALSE(result.has_value());
            if (!result.has_value()) {
                EXPECT_EQ(result.error(), *test.error);
            }
        } else {
            EXPECT_TRUE(result.has_value());
            if (result.has_value()) {
                EXPECT_EQ(result.value()->parameters(), test.parameters);
            }
        }
    }
}

static std::unique_ptr<dqmc::ConfigParameters> parse_config(std::string config)
{
    std::istringstream stream { std::move(config) };
    auto result = dqmc::ConfigParameters::create(stream);
    if (!result.has_value()) {
        std::println("{}", result.error());
    }
    return std::move(result.value());
}

TEST(ConfigParameters, read_string)
{
    auto config = parse_config("foo = bar");
    EXPECT_EQ(config->read_string("foo"), "bar");
    EXPECT_EQ(config->read_string("a").error(), "Parameter 'a' is missing in the configuration file");
}

TEST(ConfigParameters, read_integer)
{
    auto config = parse_config(R"(
foo = 123
notnumber = notnumber
integer_with_junk=123junk)");
    EXPECT_EQ(config->read_integer("foo"), 123);
    EXPECT_EQ(config->read_integer("a").error(), "Parameter 'a' is missing in the configuration file");
    EXPECT_EQ(config->read_integer("notnumber").error(), "Parameter 'notnumber' has an invalid integer value 'notnumber' in the configuration file");
    EXPECT_EQ(config->read_integer("integer_with_junk").error(), "Parameter 'integer_with_junk' has an invalid integer value '123junk' in the configuration file");
}

TEST(ConfigParameters, read_double)
{
    auto config = parse_config(R"(
junk = a
integer = 123
double = 8.2890466e-317
double_with_junk = 8.2890466e-317junk
negative = -1.25)");
    EXPECT_EQ(config->read_double("a").error(), "Parameter 'a' is missing in the configuration file");
    EXPECT_EQ(config->read_double("junk").error(), "Parameter 'junk' has an invalid floating-point value 'a' in the configuration file");
    EXPECT_EQ(config->read_double("integer"), 123.);
    EXPECT_EQ(config->read_double("double"), 8.2890466e-317);
    EXPECT_EQ(config->read_double("double_with_junk").error(), "Parameter 'double_with_junk' has an invalid floating-point value '8.2890466e-317junk' in the configuration file");
    EXPECT_EQ(config->read_double("negative"), -1.25);
}

TEST(ConfigParameters, read_boolean)
{
    auto config = parse_config(R"(
t1 = true
t2 = 1
f1 = false
f2 = 0
n1 = n
n2 = y
n3 = "true "
n4 = " false"
n5 = -1
n6 = 2
n7 = yes
)");
    EXPECT_EQ(config->read_boolean("a").error(), "Parameter 'a' is missing in the configuration file");
    EXPECT_EQ(config->read_boolean("t1"), true);
    EXPECT_EQ(config->read_boolean("t2"), true);
    EXPECT_EQ(config->read_boolean("f1"), false);
    EXPECT_EQ(config->read_boolean("f2"), false);
    EXPECT_EQ(config->read_boolean("n1").error(), "Parameter 'n1' has an invalid boolean value 'n' in the configuration file");
    EXPECT_EQ(config->read_boolean("n2").error(), "Parameter 'n2' has an invalid boolean value 'y' in the configuration file");
    EXPECT_EQ(config->read_boolean("n3").error(), "Parameter 'n3' has an invalid boolean value 'true ' in the configuration file");
    EXPECT_EQ(config->read_boolean("n4").error(), "Parameter 'n4' has an invalid boolean value ' false' in the configuration file");
    EXPECT_EQ(config->read_boolean("n5").error(), "Parameter 'n5' has an invalid boolean value '-1' in the configuration file");
    EXPECT_EQ(config->read_boolean("n6").error(), "Parameter 'n6' has an invalid boolean value '2' in the configuration file");
    EXPECT_EQ(config->read_boolean("n7").error(), "Parameter 'n7' has an invalid boolean value 'yes' in the configuration file");
}

TEST(ConfigParameters, read_double_array)
{
    auto config = parse_config(R"(
g1 =    1,      2,  3 
g2 = 1.25, -1,    8.2890466e-317, 0
g3 = 0
g4 = ""
g5 = 1,2,3,4
n1 = 0bad
n2 = 1,     bad
n3 = 1,a2
)");
    EXPECT_EQ(config->read_double_array("a").error(), "Parameter 'a' is missing in the configuration file");
    EXPECT_EQ(config->read_double_array("g1"), (std::vector<f64> { 1, 2, 3 }));
    EXPECT_EQ(config->read_double_array("g2"), (std::vector<f64> { 1.25, -1, 8.2890466e-317, 0 }));
    EXPECT_EQ(config->read_double_array("g3"), std::vector<f64> { 0 });
    EXPECT_EQ(config->read_double_array("g4"), std::vector<f64> {});
    EXPECT_EQ(config->read_double_array("g5"), (std::vector<f64> { 1, 2, 3, 4 }));
    EXPECT_EQ(config->read_double_array("n1").error(), "Parameter 'n1' has an invalid floating-point array value '0bad' in the configuration file");
    EXPECT_EQ(config->read_double_array("n2").error(), "Parameter 'n2' has an invalid floating-point array value 'bad' in the configuration file");
    EXPECT_EQ(config->read_double_array("n3").error(), "Parameter 'n3' has an invalid floating-point array value 'a2' in the configuration file");
}
