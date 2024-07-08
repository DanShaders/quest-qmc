#include "SRC/parser/Lexer.h"

namespace dqmc::parser {

Lexer::Lexer(std::string_view data, SourceRange end_of_data, DiagnosticEngine& diag)
    : m_data(data)
    , m_end_of_data(end_of_data)
    , m_diag(diag)
{
}

DiagnosticOr<void> Lexer::expect_section_end()
{
    auto line = nonempty_line();
    if (line.has_value()) {
        return m_diag.get().error(line->range_for_current_position(), "expected section end");
    }
    return {};
}

std::expected<LineLexer, SourceRange> Lexer::nonempty_line()
{
    while (m_position < m_data.size()) {
        size_t end_of_line_position = m_data.find('\n', m_position);
        if (end_of_line_position == std::string_view::npos) {
            end_of_line_position = m_data.size();
        }

        SourceRange end_of_line_range;
        if (end_of_line_position == m_data.size()) {
            end_of_line_range = m_end_of_data;
        } else {
            end_of_line_range = m_data.substr(end_of_line_position, 1);
        }

        auto line = m_data.substr(m_position, end_of_line_position - m_position);
        LineLexer lexer { line, end_of_line_range, m_diag };

        m_position = end_of_line_position + 1;

        if (lexer.skip_whitespace()) {
            return lexer;
        }
    }
    return std::unexpected { m_end_of_data };
}

} // namespace dqmc::parser
