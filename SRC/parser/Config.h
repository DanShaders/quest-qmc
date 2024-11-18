#pragma once

#include "SRC/common/Error.h"
#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"
#include "SRC/parser/Lexer.h"
#include "SRC/parser/LineLexer.h"

namespace dqmc::parser {

class ConfigParser;

class ParametersParser {
    DQMC_MAKE_NONCOPYABLE(ParametersParser);
    DQMC_MAKE_NONMOVABLE(ParametersParser);

public:
    ParametersParser() { }
    virtual ~ParametersParser() = default;

    virtual DiagnosticOr<void> parse(ConfigParser& parser, DiagnosticEngine& diag) = 0;
};

template<typename T>
struct ArrayWithSourceLocation {
    std::vector<T> value;
    std::vector<SourceRange> locations;
    SourceRange end_of_array;
};

class ConfigParser {
public:
    enum class ReportUnusedKeys {
        Yes,
        No,
    };

    ConfigParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag);

    template<std::derived_from<ParametersParser> T, typename... Args>
    requires(std::is_constructible_v<T, Args...>)
    T& register_parser(Args&&... args)
    {
        m_parsers.push_back(std::make_unique<T>(std::forward<Args>(args)...));
        return *static_cast<T*>(m_parsers.back().get());
    }

    DiagnosticOr<void> parse(ReportUnusedKeys = ReportUnusedKeys::Yes);

    DiagnosticOr<Token<f64>> claim_double(
        std::string_view key,
        std::optional<f64> default_value = std::nullopt);

    template<std::integral T = i32>
    DiagnosticOr<Token<T>> claim_integer(
        std::string_view key,
        std::optional<T> default_value = std::nullopt);

    DiagnosticOr<ArrayWithSourceLocation<f64>> claim_double_array(
        std::string_view key,
        std::optional<std::vector<f64>> default_value = std::nullopt);

private:
    DiagnosticOr<LineLexer*> claim_lexer(std::string_view key);

    struct Parameter {
        SourceRange location;
        LineLexer lexer;
        bool is_claimed = false;
    };

    std::shared_ptr<FileView> m_file;
    DiagnosticEngine& m_diag;

    std::vector<std::unique_ptr<ParametersParser>> m_parsers;

    std::map<std::string, Parameter, std::less<>> m_parameters;
};

} // namespace dqmc::parser
