#pragma once

#include "SRC/parser/LineLexer.h"

namespace dqmc::parser {

using namespace std::literals;

class Lexer {
public:
    Lexer(std::string_view data, SourceRange end_of_data, DiagnosticEngine& diag);

    std::expected<LineLexer, SourceRange> nonempty_line();
    DiagnosticOr<void> expect_section_end();

private:
    std::string_view m_data;
    SourceRange m_end_of_data;
    std::reference_wrapper<DiagnosticEngine> m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
