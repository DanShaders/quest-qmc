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

    std::string_view combine(std::string_view token1, std::string_view token2)
    {
        if (token1.empty()) {
            return token2;
        }
        if (token2.empty()) {
            return token1;
        }
        size_t start1 = token1.data() - m_data.data();
        size_t end1 = start1 + token1.size();
        size_t start2 = token2.data() - m_data.data();
        size_t end2 = start2 + token2.size();
        size_t start = std::min(start1, start2);
        return m_data.substr(start, std::max(end1, end2) - start);
    }

private:
    std::shared_ptr<FileView> m_file;
    std::string_view m_data;
    std::string_view m_past_end_token;
    std::reference_wrapper<DiagnosticEngine> m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
