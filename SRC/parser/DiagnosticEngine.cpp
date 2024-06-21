#include <algorithm>
#include <utility>

#include "SRC/parser/DiagnosticEngine.h"

namespace dqmc::parser {

void DiagnosticEngine::format_diagnostics(std::ostream& output)
{
    for (auto const& message : m_messages) {
        auto const& location = message.location;
        auto const& positions = line_positions_for(*location.file);
        auto const& content = location.file->content();

        if (location.range.empty()) {
            std::println(output, "{}:1:1: {}: {}",
                message.location.file->filename(),
                severity_to_string_view(message.severity), message.message);
            continue;
        }

        auto find_position = [&](char const* pointer) -> std::tuple<size_t, size_t, size_t> {
            ptrdiff_t offset = pointer - content.data();
            VERIFY(std::cmp_less_equal(0, offset) && std::cmp_less(offset, content.size()));

            size_t line = static_cast<size_t>(std::ranges::lower_bound(positions, offset) - positions.begin());
            size_t line_start_offset = line ? positions[line - 1] + 1 : 0;
            size_t column = static_cast<size_t>(offset) - line_start_offset;
            return { offset, line, column };
        };

        auto [start_offset, start_line, start_column] = find_position(location.range.data());
        auto [end_offset, end_line, end_column] = find_position(location.range.data() + location.range.size() - 1);

        std::println(output, "{}:{}:{}: {}: {}",
            message.location.file->filename(), start_line + 1, start_column + 1,
            severity_to_string_view(message.severity), message.message);

        for (size_t line = start_line; line <= end_line; ++line) {
            size_t line_start_offset = line ? positions[line - 1] + 1 : 0;
            size_t line_end_offset = positions[line];
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

std::vector<size_t> const& DiagnosticEngine::line_positions_for(FileView& file_view)
{
    if (auto it = m_line_positions.find(&file_view); it != m_line_positions.end()) {
        return *it->second;
    }

    auto content = file_view.content();
    auto& positions = m_line_positions[&file_view];
    positions = std::make_unique<std::vector<size_t>>();
    for (size_t i = 0; i < content.size(); ++i) {
        if (content[i] == '\n') {
            positions->push_back(i);
        }
    }
    return *positions;
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

} // namespace dqmc::parser
