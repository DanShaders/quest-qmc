#pragma once

#include <expected>
#include <memory>

#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"

namespace dqmc::parser {

class LineLexer {
public:
    LineLexer(std::string_view line, SourceRange end_of_line, DiagnosticEngine& diag);

    bool skip_whitespace();

    std::expected<Token<int>, Empty> read_integer(std::string_view name);
    std::expected<Token<f64>, Empty> read_double(std::string_view name);
    std::expected<Token<std::string>, Empty> read_string(std::string_view name);

    std::expected<void, Empty> read_comma();
    std::expected<void, Empty> read_equals();

    std::expected<void, Empty> expect_eof();

    SourceRange range_for_current_position();

private:
    std::optional<std::string_view> maybe_read_token();

    std::string_view m_data;
    SourceRange m_end_of_line;
    DiagnosticEngine& m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
