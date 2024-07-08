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

    virtual std::expected<void, Empty> parse(ConfigParser& parser, DiagnosticEngine& diag) = 0;
};

template<typename T>
struct ArrayWithSourceLocation {
    std::vector<T> value;
    std::vector<SourceRange> locations;
    SourceRange end_of_array;
};

class ConfigParser {
public:
    ConfigParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag);

    template<std::derived_from<ParametersParser> T, typename... Args>
    requires(std::is_constructible_v<T, Args...>)
    T& register_parser(Args&&... args)
    {
        m_parsers.push_back(std::make_unique<T>(std::forward<Args>(args)...));
        return *static_cast<T*>(m_parsers.back().get());
    }

    std::expected<void, Empty> parse();

    std::expected<Token<f64>, Empty> claim_double(
        std::string_view key,
        std::optional<f64> default_value = std::nullopt);

    std::expected<ArrayWithSourceLocation<f64>, Empty> claim_double_array(
        std::string_view key,
        std::optional<std::vector<f64>> default_value = std::nullopt);

private:
    std::expected<LineLexer*, Empty> claim_lexer(std::string_view key);

    struct Parameter {
        LineLexer lexer;
        bool is_claimed = false;
    };

    std::shared_ptr<FileView> m_file;
    DiagnosticEngine& m_diag;

    std::vector<std::unique_ptr<ParametersParser>> m_parsers;

    std::map<std::string, Parameter, std::less<>> m_parameters;
};

} // namespace dqmc::parser
