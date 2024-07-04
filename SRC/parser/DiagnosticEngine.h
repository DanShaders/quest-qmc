#pragma once

#include <map>
#include <utility>

#include "SRC/common/FileView.h"

namespace dqmc::parser {

class DiagnosticEngine;

class SourceRange {
    DQMC_MAKE_DEFAULT_COPYABLE(SourceRange);
    DQMC_MAKE_DEFAULT_MOVABLE(SourceRange);

public:
    SourceRange() = default;

    SourceRange(FileView const& file)
        : m_start(file.content().data())
        , m_length(0)
    {
    }

    SourceRange(std::string_view substring)
        : m_start(substring.data())
        , m_length(substring.size())
    {
    }

    static SourceRange at_end_of_file(FileView const& file)
    {
        return { file.content().data(), eof_marker };
    }

    bool is_null() const { return m_start == nullptr && m_length == 0; }
    bool is_file() const { return m_start != nullptr && (m_length == 0 || m_length == eof_marker); }
    bool is_source_range() const { return m_length != 0 && m_length != eof_marker; }

    SourceRange combined_with(SourceRange other);

private:
    static constexpr size_t eof_marker = std::numeric_limits<size_t>::max();

    friend class DiagnosticEngine;

    SourceRange(char const* start, size_t length)
        : m_start(start)
        , m_length(length)
    {
    }

    char const* m_start = nullptr;
    size_t m_length = 0;
};

template<typename T>
struct Token {
    T value;
    SourceRange location;
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
    void warning(SourceRange location, std::format_string<Args...> fmt, Args&&... args)
    {
        m_messages.push_back({
            .severity = Severity::Warning,
            .location = std::move(location),
            .message = std::format(fmt, std::forward<Args>(args)...),
        });
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

    void register_file(std::shared_ptr<FileView> file);
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

    struct FileData {
        std::shared_ptr<FileView> file;
        std::string_view filename;
        std::vector<size_t> line_positions;
    };

    std::string_view severity_to_string_view(Severity severity);
    std::vector<size_t> compute_line_positions(FileView& file_view);

    std::map<std::pair<uintptr_t, uintptr_t>, FileData> m_file_data;
    std::vector<Message> m_messages;
    bool m_has_errors = false;
};

} // namespace dqmc::parser
