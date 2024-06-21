#include "SRC/parser/FreeformGeometryParser.h"

namespace dqmc::parser {

static constexpr auto ordinal_names = std::to_array({ "1st"sv, "2nd"sv, "3rd"sv });
static constexpr auto component_names = std::to_array({ "x"sv, "y"sv, "z"sv });

FreeformGeometryParser::FreeformGeometryParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag)
    : m_file(std::move(file))
    , m_diag(diag)
{
}

std::expected<ParsedFreeformGeometry, Empty> FreeformGeometryParser::parse()
{
    split_into_sections();

    (void)parse_number_of_dimensions();
    if (!m_geometry.dimensions) {
        return Empty::error();
    }
    (void)parse_lattice_basis();
    (void)parse_supercell_basis();
    (void)parse_primitive_cell_sites();

    if (m_diag.has_errors()) {
        return Empty::error();
    }

    return m_geometry;
}

void FreeformGeometryParser::split_into_sections()
{
    using namespace std::literals;

    std::string_view contents = m_file->content();

    std::optional<Section> section_header;
    size_t section_start = 0;

    std::map<Section, SourceRange> read_headers;

    auto finalize_current_section = [&](size_t section_end, std::string_view past_end_token) {
        if (section_start == section_end) {
            return;
        }
        auto current_section_data = contents.substr(section_start, section_end - section_start);
        Lexer current_section { m_file, current_section_data, past_end_token, m_diag };
        if (section_header.has_value()) {
            auto location = read_headers.at(*section_header);
            m_sections.emplace(*section_header, std::pair { std::move(current_section), location });
        } else {
            parse_preamble(current_section);
        }
        section_start = section_end;
    };

    for (size_t i = 0; i < contents.size(); ++i) {
        size_t line_end = contents.find('\n', i);
        VERIFY(line_end != std::string_view::npos);

        auto line = contents.substr(i, line_end - i);

        for (size_t section_header_index = 0; section_header_index < section_headers.size(); ++section_header_index) {
            auto current_header_string = section_headers[section_header_index];
            if (line.starts_with(current_header_string)) {
                auto current_header = static_cast<Section>(section_header_index);
                std::string_view current_header_token = line.substr(0, current_header_string.size());
                SourceRange current_header_location = { m_file, current_header_token };

                finalize_current_section(i, current_header_token);

                if (read_headers.contains(current_header)) {
                    m_diag.error(current_header_location,
                        "duplicate section '{}'", current_header_string.substr(1));
                    m_diag.note(read_headers[current_header],
                        "previously defined here");
                }

                section_header = current_header;
                read_headers[current_header] = current_header_location;
                break;
            }
        }
        i = line_end;
    }

    finalize_current_section(contents.size(), contents.substr(contents.size() - 1));
}

void FreeformGeometryParser::parse_preamble(Lexer section)
{
    auto line = section.nonempty_line();
    if (line.has_value()) {
        m_diag.error(line->range_for_current_position(),
            "only comments are allowed before first section header");
    }
}

std::expected<void, Empty> FreeformGeometryParser::parse_number_of_dimensions()
{
    if (!m_sections.contains(Section::NumberOfDimensions)) {
        m_diag.error({ m_file },
            "missing required 'NDIM' section");
        return Empty::error();
    }

    auto& [lexer, _] = m_sections.at(Section::NumberOfDimensions);
    auto maybe_line = lexer.nonempty_line();
    if (!maybe_line.has_value()) {
        m_diag.error({ m_file, maybe_line.error() },
            "expected number of dimensions (an integer value) but found section end");
        return Empty::error();
    }

    auto [number_of_dimensions, token] = TRY(maybe_line->read_named_integer("number of dimensions"));
    if (number_of_dimensions < 1 || number_of_dimensions > 3) {
        m_diag.error({ m_file, token },
            "number of dimensions must be 1, 2, or 3");
        return Empty::error();
    }
    m_geometry.dimensions = number_of_dimensions;

    TRY(maybe_line->expect_eof());
    TRY(lexer.expect_section_end());
    return {};
}

std::expected<void, Empty> FreeformGeometryParser::parse_lattice_basis()
{
    if (!m_sections.contains(Section::LatticeBasis)) {
        m_diag.error({ m_file },
            "missing required 'PRIM' section");
        return Empty::error();
    }

    auto& [lexer, location] = m_sections.at(Section::LatticeBasis);
    m_geometry.lattice_basis_location = location;

    for (size_t i = 0; i < m_geometry.dimensions; ++i) {
        auto maybe_line = lexer.nonempty_line();
        if (!maybe_line.has_value()) {
            m_diag.error({ m_file, maybe_line.error() },
                "expected {} lattice basis vector but found section end", ordinal_names[i]);
            return Empty::error();
        }

        for (size_t j = 0; j < m_geometry.dimensions; ++j) {
            auto component_name = std::format("{} lattice basis vector {}-component", ordinal_names[i], component_names[j]);
            auto [value, token] = TRY(maybe_line->read_named_double(component_name));
            m_geometry.lattice_basis[i][j] = value;
        }

        TRY(maybe_line->expect_eof());
    }

    TRY(lexer.expect_section_end());
    return {};
}

std::expected<void, Empty> FreeformGeometryParser::parse_supercell_basis()
{
    if (!m_sections.contains(Section::SupercellBasis)) {
        m_diag.error({ m_file },
            "missing required 'SUPER' section");
        return Empty::error();
    }

    auto& [lexer, location] = m_sections.at(Section::SupercellBasis);
    m_geometry.supercell_basis_location = location;

    for (size_t i = 0; i < m_geometry.dimensions; ++i) {
        auto maybe_line = lexer.nonempty_line();
        if (!maybe_line.has_value()) {
            m_diag.error({ m_file, maybe_line.error() },
                "expected {} supercell basis vector but found section end", ordinal_names[i]);
            return Empty::error();
        }

        for (size_t j = 0; j < m_geometry.dimensions; ++j) {
            auto component_name = std::format("{} supercell basis vector {} component", ordinal_names[i], ordinal_names[j]);
            auto [value, token] = TRY(maybe_line->read_named_integer(component_name));
            m_geometry.supercell_basis[i][j] = value;
        }

        TRY(maybe_line->expect_eof());
    }

    TRY(lexer.expect_section_end());
    return {};
}

std::expected<void, Empty> FreeformGeometryParser::parse_primitive_cell_sites()
{
    if (!m_sections.contains(Section::PrimitiveCellSites)) {
        m_diag.error({ m_file },
            "missing required 'ORB' section");
        return Empty::error();
    }

    auto& [lexer, _] = m_sections.at(Section::PrimitiveCellSites);
    for (int i = 0;; ++i) {
        auto maybe_line = lexer.nonempty_line();
        if (!maybe_line.has_value()) {
            if (i == 0) {
                m_diag.error({ m_file, maybe_line.error() },
                    "at least 1 primitive cell site is required");
            }
            break;
        }

        auto [label, label_token] = TRY(maybe_line->read_named_string("primitive cell site label"));
        ParsedFreeformGeometry::PrimitiveCellSite site {
            .label = label,
        };
        for (size_t i = 0; i < 3; ++i) {
            auto name = std::format("{} site displacement {}-component", label, component_names[i]);
            auto [value, token] = TRY(maybe_line->read_named_double(name));
            site.displacement[i] = value;
        }
        m_geometry.primitive_cell_sites.emplace_back(std::move(site));

        TRY(maybe_line->expect_eof());
    }
    return {};
}

} // namespace dqmc::parser
