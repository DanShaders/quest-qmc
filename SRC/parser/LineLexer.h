#pragma once

#include <expected>
#include <memory>

#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"

namespace dqmc::parser {

class LineLexer {
public:
    LineLexer(std::shared_ptr<FileView> file, std::string_view data, DiagnosticEngine& diag);

    template<typename T>
    struct Token {
        T value;
        std::string_view token;
    };

    bool skip_whitespace();

    std::expected<Token<int>, Empty> read_named_integer(std::string_view name);
    std::expected<Token<f64>, Empty> read_named_double(std::string_view name);
    std::expected<Token<std::string>, Empty> read_named_string(std::string_view name);

    std::expected<void, Empty> expect_eof();

    SourceRange range_for_current_position();

private:
    bool is_valid_token_character(char c);
    std::optional<std::string_view> maybe_read_token();

    std::shared_ptr<FileView> m_file;
    std::string_view m_data;
    DiagnosticEngine& m_diag;

    size_t m_position = 0;
};

} // namespace dqmc::parser
