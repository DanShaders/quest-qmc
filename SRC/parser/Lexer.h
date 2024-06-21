#pragma once

#include "SRC/parser/LineLexer.h"

namespace dqmc::parser {

using namespace std::literals;

class Lexer {
public:
    Lexer(std::shared_ptr<FileView> file, std::string_view data, std::string_view past_end_token, DiagnosticEngine& diag);

    std::expected<LineLexer, std::string_view> nonempty_line();
    std::expected<void, Empty> expect_section_end();

    std::string_view data() const { return m_data; }

private:
    std::shared_ptr<FileView> m_file;
    std::string_view m_data;
    std::string_view m_past_end_token;
    std::reference_wrapper<DiagnosticEngine> m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
