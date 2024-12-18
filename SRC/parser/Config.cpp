#include "SRC/common/Ascii.h"
#include "SRC/parser/Config.h"

namespace dqmc::parser {

ConfigParser::ConfigParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag)
    : m_file(std::move(file))
    , m_diag(diag)
{
    m_diag.register_file(m_file);
}

DiagnosticOr<void> ConfigParser::parse(ReportUnusedKeys report_unused_keys)
{
    Lexer lexer { m_file->content(), SourceRange::at_end_of_file(*m_file), m_diag };
    while (true) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            break;
        }

        auto key = TRY(line->read_string("key"));
        std::ranges::transform(key.value, key.value.begin(), ascii_to_lower);
        TRY(line->read_equals());

        m_parameters.emplace(key.value,
            Parameter {
                .location = key.location,
                .lexer = std::move(*line),
            });
    }

    for (auto& parser : m_parsers) {
        parser->parse(*this, m_diag);
    }

    if (report_unused_keys == ReportUnusedKeys::Yes) {
        for (auto& [key, parameter] : m_parameters) {
            if (!parameter.is_claimed) {
                m_diag.warning(parameter.location, "key '{}' is not used", key);
            }
        }
    }

    return {};
}

DiagnosticOr<Token<f64>> ConfigParser::claim_double(
    std::string_view key,
    std::optional<f64> default_value)
{
    if (default_value.has_value() && !m_parameters.contains(key)) {
        return Token { default_value.value() };
    }

    auto& lexer = *TRY(claim_lexer(key));
    auto result = TRY(lexer.read_double(std::format("{} value", key)));
    TRY(lexer.expect_eof());
    return result;
}

template<std::integral T>
DiagnosticOr<Token<T>> ConfigParser::claim_integer(
    std::string_view key,
    std::optional<T> default_value)
{
    if (default_value.has_value() && !m_parameters.contains(key)) {
        return Token { default_value.value() };
    }

    auto& lexer = *TRY(claim_lexer(key));
    auto result = TRY(lexer.read_integer<T>(std::format("{} value", key)));
    TRY(lexer.expect_eof());
    return result;
}

template DiagnosticOr<Token<i32>> ConfigParser::claim_integer<i32>(std::string_view, std::optional<i32>);
template DiagnosticOr<Token<u32>> ConfigParser::claim_integer<u32>(std::string_view, std::optional<u32>);
template DiagnosticOr<Token<i64>> ConfigParser::claim_integer<i64>(std::string_view, std::optional<i64>);
template DiagnosticOr<Token<u64>> ConfigParser::claim_integer<u64>(std::string_view, std::optional<u64>);

DiagnosticOr<ArrayWithSourceLocation<f64>> ConfigParser::claim_double_array(
    std::string_view key,
    std::optional<std::vector<f64>> default_value)
{
    if (default_value.has_value() && !m_parameters.contains(key)) {
        return ArrayWithSourceLocation {
            .value = default_value.value(),
            .locations = std::vector<SourceRange>(default_value->size()),
            .end_of_array = {},
        };
    }

    auto& lexer = *TRY(claim_lexer(key));
    ArrayWithSourceLocation<f64> result;
    while (true) {
        auto value = TRY(lexer.read_double(std::format("{} array element", key)));
        result.value.push_back(value.value);
        result.locations.push_back(value.location);

        bool has_more_tokens = lexer.skip_whitespace();
        if (!has_more_tokens) {
            break;
        }
        TRY(lexer.read_comma());
    }
    TRY(lexer.expect_eof());
    result.end_of_array = lexer.range_for_current_position();
    return result;
}

DiagnosticOr<LineLexer*> ConfigParser::claim_lexer(std::string_view key)
{
    auto it = m_parameters.find(key);
    if (it == m_parameters.end()) {
        return m_diag.error(*m_file, "required key '{}' is not defined",
            key);
    }

    auto& parameter = it->second;
    VERIFY(!parameter.is_claimed);
    parameter.is_claimed = true;

    return &parameter.lexer;
}

} // namespace dqmc::parser
