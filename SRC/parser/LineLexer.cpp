#include <charconv>
#include <fast_float/fast_float.h>

#include "SRC/parser/LineLexer.h"

namespace dqmc::parser {

using namespace std::literals;

static bool is_space(char c)
{
    return c == ' ' || c == '\n' || c == '\t';
}

LineLexer::LineLexer(std::shared_ptr<FileView> file, std::string_view data, DiagnosticEngine& diag)
    : m_file(std::move(file))
    , m_data(data)
    , m_diag(diag)
{
    VERIFY(data.back() == '\n');
}

bool LineLexer::skip_whitespace()
{
    while (m_position < m_data.size() - 1 && is_space(m_data[m_position])) {
        ++m_position;
    }
    return is_valid_token_character(m_data[m_position]);
}

auto LineLexer::read_named_integer(std::string_view name) -> std::expected<Token<int>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(),
            "expected {} (an integer value) but found line end", name);
        return Empty::error();
    }

    int value = 0;
    auto [ptr, ec] = std::from_chars(m_data.data() + m_position, m_data.data() + m_data.size(), value);
    auto token = maybe_read_token().value();

    if (is_valid_token_character(*ptr)) {
        m_diag.error({ m_file, token },
            "expected {} (an integer value) but found '{}'", name, token);
        return Empty::error();
    }
    size_t new_position = ptr - m_data.data();
    VERIFY(new_position == m_position);

    if (ec == std::errc::result_out_of_range) {
        m_diag.error({ m_file, token },
            "value {} for {} overflows 32-bit integral variable", token, name);
        return Empty::error();
    }

    return Token { value, token };
}

auto LineLexer::read_named_double(std::string_view name) -> std::expected<Token<f64>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(),
            "expected {} (a floating point value) but found line end", name);
        return Empty::error();
    }

    f64 value = 0;
    auto data_start = m_data.data() + m_position, data_end = m_data.data() + m_data.size();
    auto [ptr, ec] = fast_float::from_chars(data_start, data_end, value);
    auto token = maybe_read_token().value();

    if (*ptr == 'd') {
        std::string mutated_token { token };
        mutated_token[ptr - data_start] = 'e';
        auto begin = mutated_token.data(), end = mutated_token.data() + mutated_token.size();
        auto [mutated_ptr, mutated_ec] = fast_float::from_chars(begin, end, value);
        if (mutated_ptr == end) {
            ptr = token.end();
            ec = mutated_ec;
        }
    }

    if (is_valid_token_character(*ptr)) {
        m_diag.error({ m_file, token },
            "expected {} (a floating point value) but found '{}'", name, token);
        return Empty::error();
    }
    size_t new_position = ptr - m_data.data();
    VERIFY(new_position == m_position);

    if (ec == std::errc::result_out_of_range) {
        m_diag.error({ m_file, token },
            "value {} for {} overflows 64-bit floating point variable", token, name);
        return Empty::error();
    }

    return Token { value, token };
}

auto LineLexer::read_named_string(std::string_view name) -> std::expected<Token<std::string>, Empty>
{
    if (!skip_whitespace()) {
        m_diag.error(range_for_current_position(),
            "expected {} (a string) but found line end", name);
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
        if (end_position + 1 == m_data.size()) {
            m_diag.error({ m_file, m_data.substr(end_position, 1) },
                "expected closing '\"'");
            m_diag.note(quote_range, "to match this '\"'");
            return Empty::error();
        }

        char c = m_data[end_position];
        if (c == '"') {
            break;
        } else if (c == '\\') {
            static std::map<char, char> const escape_characters = {
                { 'n', '\n' },
                { 'r', '\r' },
                { 't', '\t' },
                { '\\', '\\' },
                { '"', '"' },
            };
            char escape = m_data[++end_position];
            if (!escape_characters.contains(escape)) {
                m_diag.error({ m_file, m_data.substr(end_position - 1, 2) },
                    "invalid escape character {:?}", escape);
                return Empty::error();
            }
            value += escape_characters.at(escape);
        } else {
            value += c;
        }
    }

    auto token = m_data.substr(m_position, end_position - m_position + 1);
    m_position = end_position + 1;
    return Token { value, token };
}

auto LineLexer::expect_eof() -> std::expected<void, Empty>
{
    if (skip_whitespace()) {
        m_diag.error(range_for_current_position(),
            "expected end of line");
        return Empty::error();
    }
    return {};
}

SourceRange LineLexer::range_for_current_position()
{
    return { m_file, m_data.substr(m_position, 1) };
}

bool LineLexer::is_valid_token_character(char c)
{
    return !is_space(c) && c != '#';
}

std::optional<std::string_view> LineLexer::maybe_read_token()
{
    if (!skip_whitespace()) {
        return std::nullopt;
    }
    size_t start_position = m_position;
    while (!is_space(m_data[m_position]) && m_data[m_position] != '#') {
        ++m_position;
    }
    return m_data.substr(start_position, m_position - start_position);
}

} // namespace dqmc::parser
