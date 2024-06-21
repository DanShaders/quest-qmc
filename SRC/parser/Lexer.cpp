#include "SRC/parser/Lexer.h"

namespace dqmc::parser {

Lexer::Lexer(std::shared_ptr<FileView> file, std::string_view data, std::string_view past_end_token, DiagnosticEngine& diag)
    : m_file(std::move(file))
    , m_data(data)
    , m_past_end_token(past_end_token)
    , m_diag(diag)
{
    VERIFY(data.back() == '\n');
}

std::expected<LineLexer, std::string_view> Lexer::nonempty_line()
{
    while (m_position != m_data.size()) {
        size_t end = m_data.find('\n', m_position);
        VERIFY(end != std::string_view::npos);
        LineLexer lexer { m_file, m_data.substr(m_position, end - m_position + 1), m_diag };
        m_position = end + 1;
        if (lexer.skip_whitespace()) {
            return lexer;
        }
    }
    return std::unexpected { m_past_end_token };
}

std::expected<void, Empty> Lexer::expect_section_end()
{
    auto line = nonempty_line();
    if (line.has_value()) {
        m_diag.get().error(line->range_for_current_position(),
            "expected section end");
        return Empty::error();
    }
    return {};
}

} // namespace dqmc::parser
