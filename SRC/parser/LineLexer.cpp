#include <charconv>
#include <fast_float/fast_float.h>

#include "SRC/parser/LineLexer.h"

namespace dqmc::parser {

using namespace std::literals;

static bool is_space(char c)
{
    return c == ' ' || c == '\t';
}

LineLexer::LineLexer(std::string_view data, SourceRange end_of_line, DiagnosticEngine& diag)
    : m_data(data)
    , m_end_of_line(end_of_line)
    , m_diag(diag)
{
}

bool LineLexer::skip_whitespace()
{
    while (m_position < m_data.size() && is_space(m_data[m_position])) {
        ++m_position;
    }
    return m_position != m_data.size() && m_data[m_position] != '#';
}

auto LineLexer::read_integer(std::string_view name) -> std::expected<Token<int>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(), "expected {} (an integer value)",
            name);
        return Empty::error();
    }

    int value = 0;
    auto token = maybe_read_token().value();
    auto [value_end, ec] = std::from_chars(token.begin(), token.end(), value);

    if (value_end != token.end()) {
        m_diag.error(token, "expected {} (an integer value) but found '{}'",
            name, token);
        return Empty::error();
    }

    if (ec == std::errc::result_out_of_range) {
        m_diag.error(token, "value {} for {} overflows 32-bit integral variable",
            token, name);
        return Empty::error();
    }

    return Token { value, token };
}

auto LineLexer::read_double(std::string_view name) -> std::expected<Token<f64>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(), "expected {} (a floating-point value)",
            name);
        return Empty::error();
    }

    f64 value = 0;
    auto token = maybe_read_token().value();
    auto [value_end, ec] = fast_float::from_chars(token.begin(), token.end(), value);

    if (value_end != token.end() && *value_end == 'd') {
        std::string mutated_token { token };
        mutated_token[value_end - token.begin()] = 'e';
        auto begin = mutated_token.data(), end = mutated_token.data() + mutated_token.size();
        auto [mutated_ptr, mutated_ec] = fast_float::from_chars(begin, end, value);
        if (mutated_ptr == end) {
            value_end = token.end();
            ec = mutated_ec;
        }
    }

    if (value_end != token.end()) {
        m_diag.error(token, "expected {} (a floating-point value) but found '{}'",
            name, token);
        return Empty::error();
    }

    if (ec == std::errc::result_out_of_range) {
        m_diag.error(token, "value {} for {} overflows 64-bit floating point variable",
            token, name);
        return Empty::error();
    }

    return Token { value, token };
}

auto LineLexer::read_string(std::string_view name) -> std::expected<Token<std::string>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(), "expected {} (a string)",
            name);
        return Empty::error();
    }

    if (m_data[m_position] != '"') {
        auto token = maybe_read_token().value();
        return Token { std::string { token }, token };
    }

    auto quote_range = range_for_current_position();

    std::string value;
    size_t end_position = m_position + 1;
    for (;; ++end_position) {
        if (end_position == m_data.size()) {
            m_diag.error(m_end_of_line, "expected closing '\"'");
            m_diag.note(quote_range, "to match this '\"'");
            return Empty::error();
        }

        char c = m_data[end_position];
        if (c == '"') {
            ++end_position;
            break;
        } else if (c == '\\') {
            if (end_position + 1 == m_data.size()) {
                m_diag.error(m_end_of_line, "expected escape character after '\\'");
                return Empty::error();
            }

            static std::map<char, char> const escape_characters = {
                { 'n', '\n' },
                { 'r', '\r' },
                { 't', '\t' },
                { '\\', '\\' },
                { '"', '"' },
            };
            char escape = m_data[++end_position];
            if (!escape_characters.contains(escape)) {
                m_diag.error(m_data.substr(end_position - 1, 2), "invalid escape sequence '\\{}'",
                    escape);
                return Empty::error();
            }
            value += escape_characters.at(escape);
        } else {
            value += c;
        }
    }

    auto token = m_data.substr(m_position, end_position - m_position);
    m_position = end_position;
    return Token { value, token };
}

std::expected<void, Empty> LineLexer::read_comma()
{
    auto token = maybe_read_token();
    if (!token.has_value()) {
        m_diag.error(range_for_current_position(), "expected ','");
        return Empty::error();
    }
    if (token.value() != ",") {
        m_diag.error(token.value(), "expected ',' but found '{}'", token.value());
        return Empty::error();
    }
    return {};
}

std::expected<void, Empty> LineLexer::read_equals()
{
    auto token = maybe_read_token();
    if (!token.has_value()) {
        m_diag.error(range_for_current_position(), "expected '='");
        return Empty::error();
    }
    if (token.value() != "=") {
        m_diag.error(token.value(), "expected '=' but found '{}'", token.value());
        return Empty::error();
    }
    return {};
}

std::expected<void, Empty> LineLexer::expect_eof()
{
    if (skip_whitespace()) {
        m_diag.error(range_for_current_position(), "expected end of line");
        return Empty::error();
    }
    return {};
}

SourceRange LineLexer::range_for_current_position()
{
    if (m_position == m_data.size()) {
        return m_end_of_line;
    }
    return m_data.substr(m_position, 1);
}

std::optional<std::string_view> LineLexer::maybe_read_token()
{
    if (!skip_whitespace()) {
        return std::nullopt;
    }
    size_t start_position = m_position;
    if (m_data[start_position] != ',' && m_data[start_position] != '=') {
        auto is_valid_token_continuation = [&](char c) {
            return !is_space(c) && c != '#' && c != ',' && c != '=';
        };
        while (m_position < m_data.size() && is_valid_token_continuation(m_data[m_position])) {
            ++m_position;
        }
    } else {
        ++m_position;
    }
    return m_data.substr(start_position, m_position - start_position);
}

} // namespace dqmc::parser
