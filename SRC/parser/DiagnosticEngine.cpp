#include <algorithm>
#include <utility>

#include "SRC/parser/DiagnosticEngine.h"

namespace dqmc::parser {

SourceRange SourceRange::combined_with(SourceRange other)
{
    VERIFY(!is_file() && !other.is_file());
    if (is_null()) {
        return other;
    }
    if (other.is_null()) {
        return *this;
    }

    uintptr_t start1 = std::bit_cast<uintptr_t>(m_start);
    uintptr_t start2 = std::bit_cast<uintptr_t>(other.m_start);
    uintptr_t end1 = std::bit_cast<uintptr_t>(m_start + m_length);
    uintptr_t end2 = std::bit_cast<uintptr_t>(other.m_start + other.m_length);
    char const* start = std::bit_cast<char const*>(std::min(start1, start2));
    char const* end = std::bit_cast<char const*>(std::max(end1, end2));
    return { start, static_cast<size_t>(end - start) };
}

void DiagnosticEngine::format_diagnostics(std::ostream& output)
{
    for (auto const& message : m_messages) {
        auto const& location = message.location;

        auto severity = severity_to_string_view(message.severity);

        if (location.is_null()) {
            std::println("{}: {}",
                severity, message.message);
            continue;
        }

        auto start_pointer = std::bit_cast<uintptr_t>(location.m_start);
        auto file_data_iterator = m_file_data.upper_bound({ start_pointer, std::numeric_limits<uintptr_t>::max() });
        VERIFY(m_file_data.begin() != file_data_iterator);
        --file_data_iterator;
        auto const& [file_data_key, file_data] = *file_data_iterator;
        VERIFY(start_pointer >= file_data_key.first);
        VERIFY(start_pointer < file_data_key.second);

        auto filename = file_data.filename;
        auto const& content = file_data.file->content();

        size_t start_offset;
        size_t end_offset;

        if (location.is_file()) {
            if (location.m_length == 0) {
                std::println("{}:1:1: {}: {}",
                    filename, severity, message.message);
                continue;
            } else {
                VERIFY(location.m_length == location.eof_marker);
                start_offset = content.size();
                end_offset = content.size();
            }
        } else {
            VERIFY(location.is_source_range());
            start_offset = start_pointer - file_data_iterator->first.first;
            end_offset = start_offset + location.m_length - 1;
        }

        auto const& line_offsets = file_data.line_positions;

        auto find_position = [&](size_t offset) -> std::pair<size_t, size_t> {
            size_t line = static_cast<size_t>(std::ranges::lower_bound(line_offsets, offset) - line_offsets.begin());
            size_t line_start_offset = line ? line_offsets[line - 1] + 1 : 0;
            size_t column = static_cast<size_t>(offset) - line_start_offset;
            return { line, column };
        };

        auto [start_line, start_column] = find_position(start_offset);
        auto [end_line, end_column] = find_position(end_offset);

        std::println(output, "{}:{}:{}: {}: {}",
            filename, start_line + 1, start_column + 1, severity, message.message);

        for (size_t line = start_line; line <= end_line; ++line) {
            size_t line_start_offset = line ? line_offsets[line - 1] + 1 : 0;
            size_t line_end_offset = line < line_offsets.size() ? line_offsets[line] : content.size();
            std::println(output, "{:5} | {}", line + 1, content.substr(line_start_offset, line_end_offset - line_start_offset));

            std::string squiggle(line_end_offset - line_start_offset + 1, ' ');
            std::fill(
                squiggle.begin() + static_cast<ssize_t>(std::max(line_start_offset, start_offset) - line_start_offset),
                squiggle.begin() + static_cast<ssize_t>(std::min(line_end_offset, end_offset) - line_start_offset + 1),
                '~');
            if (line == start_line) {
                squiggle[start_column] = '^';
            }
            std::println(output, "      | {}", squiggle);
        }
    }
}

std::vector<size_t> DiagnosticEngine::compute_line_positions(FileView& file_view)
{
    auto content = file_view.content();
    std::vector<size_t> result;
    for (size_t i = 0; i < content.size(); ++i) {
        if (content[i] == '\n') {
            result.push_back(i);
        }
    }
    result.push_back(content.size());
    return result;
}

std::string_view DiagnosticEngine::severity_to_string_view(Severity severity)
{
    using namespace std::literals;

    if (severity == Severity::Error) {
        return "error"sv;
    } else if (severity == Severity::Warning) {
        return "warning"sv;
    } else {
        return "note"sv;
    }
}

void DiagnosticEngine::register_file(std::shared_ptr<FileView> file)
{
    auto content = file->content();
    auto filename = file->filename();

    std::pair file_data_key = {
        std::bit_cast<uintptr_t>(content.data()),
        std::bit_cast<uintptr_t>(content.data()) + content.size(),
    };

    if (!m_file_data.contains(file_data_key)) {
        m_file_data[file_data_key] = {
            .file = file,
            .filename = filename,
            .line_positions = compute_line_positions(*file),
        };
    }
}

} // namespace dqmc::parser
