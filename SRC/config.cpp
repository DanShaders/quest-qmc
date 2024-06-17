#include "SRC/config.h"

#include <charconv>
#include <fast_float/fast_float.h>

namespace dqmc {

namespace {

struct KeyValuePair {
    std::string key;
    std::string original_key;
    std::string value;
};

std::expected<std::optional<KeyValuePair>, std::string> parse_line(std::string const& line, int line_count)
{
    enum class State {
        BeforeKey,
        ReadingKey,
        BeforeEquals,
        BeforeValue,
        ReadingUnquotedValue,
        ReadingQuotedValue,
        AfterQuotedValue,
    } state
        = State::BeforeKey;
    std::string key;
    std::string original_key;
    std::string value;

    for (size_t i = 0, column_count = 1; i < line.size(); ++i, ++column_count) {
        char c = line[i];
        if (state == State::BeforeKey) {
            if (std::isspace(c)) {
                continue;
            } else if (c == '#') {
                break;
            } else if (c == '=') {
                return std::unexpected { std::format(
                    "Unexpected '=' at {}:{} in the configuration file",
                    line_count, column_count) };
            } else {
                key += ascii_to_lower(c);
                original_key += c;
                state = State::ReadingKey;
            }
        } else if (state == State::ReadingKey) {
            if (std::isspace(c)) {
                state = State::BeforeEquals;
            } else if (c == '#') {
                return std::unexpected { std::format(
                    "Unexpected '#' at {}:{} in the configuration file",
                    line_count, column_count) };
            } else if (c == '=') {
                state = State::BeforeValue;
            } else {
                key += ascii_to_lower(c);
                original_key += c;
            }
        } else if (state == State::BeforeEquals) {
            if (std::isspace(c)) {
                continue;
            } else if (c == '=') {
                state = State::BeforeValue;
            } else {
                return std::unexpected { std::format(
                    "Unexpected '{}' at {}:{} in the configuration file",
                    c, line_count, column_count) };
            }
        } else if (state == State::BeforeValue) {
            if (std::isspace(c)) {
                continue;
            } else if (c == '#') {
                return std::unexpected { std::format(
                    "Unexpected '#' at {}:{} in the configuration file",
                    line_count, column_count) };
            } else if (c == '"') {
                state = State::ReadingQuotedValue;
            } else {
                value += c;
                state = State::ReadingUnquotedValue;
            }
        } else if (state == State::ReadingUnquotedValue) {
            if (c == '#') {
                break;
            } else {
                value += c;
            }
        } else if (state == State::ReadingQuotedValue) {
            if (c == '"') {
                state = State::AfterQuotedValue;
            } else if (c == '\\') {
                static std::map<char, char> const escape_characters = {
                    { 'n', '\n' },
                    { 'r', '\r' },
                    { 't', '\t' },
                    { '\\', '\\' },
                    { '"', '"' },
                };

                if (i + 1 >= line.size()) {
                    return std::unexpected { std::format(
                        "Unexpected end of line at {}:{} in the configuration file",
                        line_count, column_count + 1) };
                } else if (!escape_characters.contains(line[i + 1])) {
                    return std::unexpected { std::format(
                        "Unexpected escape sequence '\\{}' at {}:{} in the configuration file",
                        line[i + 1], line_count, column_count + 1) };
                } else {
                    value += escape_characters.at(line[++i]);
                    ++column_count;
                }
            } else {
                value += c;
            }
        } else if (state == State::AfterQuotedValue) {
            if (std::isspace(c)) {
                continue;
            } else if (c == '#') {
                break;
            } else {
                return std::unexpected { std::format(
                    "Unexpected '{}' at {}:{} in the configuration file",
                    c, line_count, column_count) };
            }
        } else {
            VERIFY_NOT_REACHED();
        }
    }

    if (state != State::BeforeKey && state != State::ReadingUnquotedValue && state != State::AfterQuotedValue) {
        return std::unexpected { std::format(
            "Unexpected end of line at {}:{} in the configuration file",
            line_count, line.size() + 1) };
    }

    if (state == State::ReadingUnquotedValue) {
        while (!value.empty() && std::isspace(value.back())) {
            value.pop_back();
        }
        VERIFY(!value.empty());
    }

    if (state == State::BeforeKey) {
        return std::nullopt;
    } else {
        return KeyValuePair { .key = key, .original_key = original_key, .value = value };
    }
}

} // namespace

std::expected<std::unique_ptr<ConfigParameters>, std::string> ConfigParameters::create(std::istream& input_stream)
{
    auto result = std::unique_ptr<ConfigParameters> { new ConfigParameters {} };
    auto& parameters = result->m_parameters;

    for (int line_count = 1;; ++line_count) {
        std::string line;
        std::getline(input_stream, line);
        if (!input_stream) {
            break;
        }

        auto parsed_line = TRY(parse_line(line, line_count));
        if (!parsed_line.has_value()) {
            continue;
        }
        auto [key, original_key, value] = std::move(parsed_line).value();

        if (auto it = parameters.find(key); it != parameters.end()) {
            return std::unexpected { std::format(
                "Parameter '{}' has duplicate definitions at lines {} and {} of the configuration file",
                original_key, it->second.line_number, line_count) };
        }
        parameters[key] = {
            .value = value,
            .line_number = line_count,
        };
    }

    if (input_stream.bad()) {
        return std::unexpected { std::format("A stream error occured while reading configuration file") };
    }

    return std::move(result);
}

// CLEANUP: Make these static
std::expected<int, Empty> parse_integer(std::string_view value);
std::expected<bool, Empty> parse_boolean(std::string_view value);
std::expected<double, Empty> parse_double(std::string_view value);
std::expected<std::vector<double>, std::string_view> parse_double_array(std::string_view value);

std::expected<int, Empty> parse_integer(std::string_view value)
{
    int result;
    auto [ptr, ec] = std::from_chars(value.begin(), value.end(), result);
    if (ec != std::errc() || ptr != value.end()) {
        return std::unexpected { Empty {} };
    }
    return result;
}

std::expected<bool, Empty> parse_boolean(std::string_view value)
{
    if (value == "true" || value == "1") {
        return true;
    } else if (value == "false" || value == "0") {
        return false;
    } else {
        return std::unexpected { Empty {} };
    }
}

std::expected<double, Empty> parse_double(std::string_view value)
{
    f64 result;
    auto [ptr, ec] = fast_float::from_chars(value.begin(), value.end(), result);
    if (ec != std::errc() || ptr != value.end()) {
        return std::unexpected { Empty {} };
    }
    return result;
}

std::expected<std::vector<double>, std::string_view> parse_double_array(std::string_view value)
{
    std::vector<f64> result;
    size_t start = 0;
    while (start < value.size()) {
        size_t end = value.find(',', start);
        if (end == std::string::npos) {
            end = value.size();
        }

        auto number_span = std::string_view { value }.substr(start, end - start);
        f64 number = TRY(parse_double(number_span).transform_error([&](auto) {
            return number_span;
        }));

        result.push_back(number);
        start = end + 1;
        while (start < value.size() && std::isspace(value[start])) {
            ++start;
        }
    }
    return result;
}

std::expected<std::string_view, std::string> ConfigParameters::read_string(std::string_view key, std::optional<std::string_view> default_value) const
{
    auto it = m_parameters.find(key);
    if (it == m_parameters.end()) {
        if (default_value.has_value()) {
            return *default_value;
        }
        return std::unexpected { std::format(
            "Parameter '{}' is missing in the configuration file",
            key) };
    }

    return it->second.value;
}

std::expected<int, std::string> ConfigParameters::read_integer(std::string_view key, std::optional<int> default_value) const
{
    auto value_or_error = read_string(key);
    if (!value_or_error.has_value() && default_value.has_value()) {
        return *default_value;
    }

    auto value = TRY(value_or_error);
    return parse_integer(value).transform_error([&](auto) {
        return std::format(
            "Parameter '{}' has an invalid integer value '{}' in the configuration file",
            key, value);
    });
}

std::expected<bool, std::string> ConfigParameters::read_boolean(std::string_view key, std::optional<bool> default_value) const
{
    auto value_or_error = read_string(key);
    if (!value_or_error.has_value() && default_value.has_value()) {
        return *default_value;
    }

    auto value = TRY(value_or_error);
    return parse_boolean(value).transform_error([&](auto) {
        return std::format(
            "Parameter '{}' has an invalid boolean value '{}' in the configuration file",
            key, value);
    });
}

std::expected<f64, std::string> ConfigParameters::read_double(std::string_view key, std::optional<f64> default_value) const
{
    auto value_or_error = read_string(key);
    if (!value_or_error.has_value() && default_value.has_value()) {
        return *default_value;
    }

    auto value = TRY(value_or_error);
    return parse_double(value).transform_error([&](auto) {
        return std::format(
            "Parameter '{}' has an invalid floating-point value '{}' in the configuration file",
            key, value);
    });
}

std::expected<std::vector<f64>, std::string> ConfigParameters::read_double_array(std::string_view key, std::optional<std::vector<f64>> default_value) const
{
    auto value_or_error = read_string(key);
    if (!value_or_error.has_value() && default_value.has_value()) {
        return *default_value;
    }

    auto value = TRY(value_or_error);
    return parse_double_array(value).transform_error([&](auto element) {
        return std::format(
            "Parameter '{}' has an invalid floating-point array value '{}' in the configuration file",
            key, element);
    });
}

} // namespace dqmc
