#pragma once

#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"

namespace dqmc::parser {

class LineLexer {
public:
    LineLexer(std::string_view line, SourceRange end_of_line, DiagnosticEngine& diag);

    bool skip_whitespace();

    template<std::integral T = i32>
    DiagnosticOr<Token<T>> read_integer(std::string_view name);
    DiagnosticOr<Token<f64>> read_double(std::string_view name);
    DiagnosticOr<Token<std::string>> read_string(std::string_view name);

    DiagnosticOr<void> read_comma();
    DiagnosticOr<void> read_equals();

    DiagnosticOr<void> expect_eof();

    SourceRange range_for_current_position();

private:
    std::optional<std::string_view> maybe_read_token();

    std::string_view m_data;
    SourceRange m_end_of_line;
    DiagnosticEngine& m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
