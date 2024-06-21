#pragma once

#include <map>

#include "SRC/common/FileView.h"

namespace dqmc::parser {

struct SourceRange {
    std::shared_ptr<FileView> file;
    std::string_view range;
};

class DiagnosticEngine {
public:
    template<typename... Args>
    void error(SourceRange location, std::format_string<Args...> fmt, Args&&... args)
    {
        m_messages.push_back({
            .severity = Severity::Error,
            .location = std::move(location),
            .message = std::format(fmt, std::forward<Args>(args)...),
        });
        m_has_errors = true;
    }

    template<typename... Args>
    void note(SourceRange location, std::format_string<Args...> fmt, Args&&... args)
    {
        m_messages.push_back({
            .severity = Severity::Note,
            .location = std::move(location),
            .message = std::format(fmt, std::forward<Args>(args)...),
        });
    }

    ErrorOr<void> with(std::ostream& output, auto&& func)
    {
        m_has_errors = false;
        m_messages.clear();
        std::expected<void, Empty> result = func();
        format_diagnostics(output);
        VERIFY(m_has_errors != result.has_value());
        if (result.has_value()) {
            return {};
        } else {
            return Error::formatted("Errors have been encountered");
        }
    }

    void format_diagnostics(std::ostream& output);

    bool has_errors() const { return m_has_errors; }

private:
    enum class Severity {
        Error,
        Warning,
        Note,
    };

    struct Message {
        Severity severity;
        SourceRange location;
        std::string message;
    };

    std::vector<size_t> const& line_positions_for(FileView& file_view);
    std::string_view severity_to_string_view(Severity severity);

    std::map<FileView*, std::unique_ptr<std::vector<size_t>>> m_line_positions;
    std::vector<Message> m_messages;
    bool m_has_errors = false;
};

} // namespace dqmc::parser
