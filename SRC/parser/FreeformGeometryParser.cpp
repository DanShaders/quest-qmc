#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/common/Ascii.h"

namespace dqmc::parser {

static constexpr auto ordinal_names = std::to_array({ "1st"sv, "2nd"sv, "3rd"sv });
static constexpr auto component_names = std::to_array({ "x"sv, "y"sv, "z"sv });

FreeformGeometryParser::FreeformGeometryParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag)
    : m_file(std::move(file))
    , m_data(m_file->content())
    , m_diag(diag)
{
    m_diag.register_file(m_file);
}

DiagnosticOr<ParsedFreeformGeometry> FreeformGeometryParser::parse()
{
    split_into_sections();

    (void)parse_number_of_dimensions();
    if (!m_geometry.dimensions) {
        return std::unexpected { DiagnosedError {} };
    }
    (void)parse_lattice_basis();
    (void)parse_supercell_basis();
    (void)parse_primitive_cell_sites();
    (void)parse_hamiltonian();
    (void)parse_symmetries();

    if (m_diag.has_errors()) {
        return std::unexpected { DiagnosedError {} };
    }

    return m_geometry;
}

void FreeformGeometryParser::split_into_sections()
{
    using namespace std::literals;

    std::optional<Section> section_header;
    size_t section_start = 0;

    std::map<Section, SourceRange> read_headers;

    auto finalize_current_section = [&](size_t section_end, SourceRange past_end_token) {
        if (section_start == section_end) {
            return;
        }
        auto current_section_data = m_data.substr(section_start, section_end - section_start);
        Lexer current_section { current_section_data, past_end_token, m_diag };
        if (section_header.has_value()) {
            auto location = read_headers.at(*section_header);
            m_sections.emplace(*section_header, FreeformGeometryParser::SplitSection { std::move(current_section), location });
        } else {
            parse_preamble(current_section);
        }
        section_start = section_end;
    };

    for (size_t i = 0; i < m_data.size(); ++i) {
        size_t line_end = m_data.find('\n', i);
        if (line_end == std::string_view::npos) {
            line_end = m_data.size();
        }

        auto line = m_data.substr(i, line_end - i);

        for (size_t section_header_index = 0; section_header_index < section_headers.size(); ++section_header_index) {
            auto current_header_string = section_headers[section_header_index];
            if (line.starts_with(current_header_string)) {
                auto current_header = static_cast<Section>(section_header_index);
                SourceRange current_header_location = line.substr(0, current_header_string.size());

                finalize_current_section(i, current_header_location);

                if (read_headers.contains(current_header)) {
                    m_diag.error(current_header_location, "duplicate section '{}'",
                        current_header_string.substr(1));
                    m_diag.note(read_headers[current_header], "previously defined here");
                }

                section_header = current_header;
                read_headers[current_header] = current_header_location;
                break;
            }
        }
        i = line_end;
    }

    finalize_current_section(m_data.size(), SourceRange::at_end_of_file(*m_file));
}

void FreeformGeometryParser::parse_preamble(Lexer section)
{
    auto line = section.nonempty_line();
    if (line.has_value()) {
        m_diag.error(line->range_for_current_position(), "only comments are allowed before first section header");
    }
}

DiagnosticOr<void> FreeformGeometryParser::parse_number_of_dimensions()
{
    if (!m_sections.contains(Section::NumberOfDimensions)) {
        return m_diag.error(*m_file, "missing required 'NDIM' section");
    }

    auto& [lexer, _] = m_sections.at(Section::NumberOfDimensions);
    auto line = lexer.nonempty_line();
    if (!line.has_value()) {
        return m_diag.error(line.error(), "expected number of dimensions (an integer value)");
    }

    auto [number_of_dimensions, token] = TRY(line->read_integer("number of dimensions"));
    if (number_of_dimensions < 1 || number_of_dimensions > 3) {
        return m_diag.error(token, "number of dimensions must be 1, 2, or 3");
    }
    m_geometry.dimensions = number_of_dimensions;

    TRY(line->expect_eof());
    TRY(lexer.expect_section_end());
    return {};
}

DiagnosticOr<void> FreeformGeometryParser::parse_lattice_basis()
{
    if (!m_sections.contains(Section::LatticeBasis)) {
        return m_diag.error(*m_file, "missing required 'PRIM' section");
    }

    auto& [lexer, location] = m_sections.at(Section::LatticeBasis);
    m_geometry.lattice_basis_location = location;

    for (size_t i = 0; i < m_geometry.dimensions; ++i) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            return m_diag.error(line.error(), "expected {} lattice basis vector",
                ordinal_names[i]);
        }

        for (size_t j = 0; j < m_geometry.dimensions; ++j) {
            auto component_name = std::format("{} lattice basis vector {}-component", ordinal_names[i], component_names[j]);
            auto [value, token] = TRY(line->read_double(component_name));
            m_geometry.lattice_basis[i][j] = value;
        }

        TRY(line->expect_eof());
    }

    TRY(lexer.expect_section_end());
    return {};
}

DiagnosticOr<void> FreeformGeometryParser::parse_supercell_basis()
{
    if (!m_sections.contains(Section::SupercellBasis)) {
        return m_diag.error(*m_file, "missing required 'SUPER' section");
    }

    auto& [lexer, location] = m_sections.at(Section::SupercellBasis);
    m_geometry.supercell_basis_location = location;

    for (size_t i = 0; i < m_geometry.dimensions; ++i) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            return m_diag.error(line.error(), "expected {} supercell basis vector",
                ordinal_names[i]);
        }

        for (size_t j = 0; j < m_geometry.dimensions; ++j) {
            auto component_name = std::format("{} supercell basis vector {} component", ordinal_names[i], ordinal_names[j]);
            auto [value, token] = TRY(line->read_integer(component_name));
            m_geometry.supercell_basis[i][j] = value;
        }

        TRY(line->expect_eof());
    }

    TRY(lexer.expect_section_end());
    return {};
}

DiagnosticOr<void> FreeformGeometryParser::parse_primitive_cell_sites()
{
    if (!m_sections.contains(Section::PrimitiveCellSites)) {
        return m_diag.error(*m_file, "missing required 'ORB' section");
    }

    auto& [lexer, _] = m_sections.at(Section::PrimitiveCellSites);
    for (int i = 0;; ++i) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            if (i == 0) {
                m_diag.error(line.error(), "at least 1 primitive cell site is required");
            }
            break;
        }

        ParsedFreeformGeometry::PrimitiveCellSite site;

        auto [label, label_token] = TRY(line->read_string("primitive cell site label"));
        site.label = std::move(label);
        for (size_t i = 0; i < 3; ++i) {
            auto name = std::format("{} site displacement {}-component", label, component_names[i]);
            auto [value, token] = TRY(line->read_double(name));
            site.displacement[i] = value;
        }
        m_geometry.primitive_cell_sites.emplace_back(std::move(site));

        TRY(line->expect_eof());
    }
    return {};
}

DiagnosticOr<void> FreeformGeometryParser::parse_hamiltonian()
{
    if (!m_sections.contains(Section::Hamiltonian)) {
        return m_diag.error(*m_file, "missing required 'HAMILT' section");
    }

    auto& [lexer, _] = m_sections.at(Section::Hamiltonian);
    while (true) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            break;
        }

        auto [from, from_token] = TRY(line->read_integer("site index"));
        auto [to, to_token] = TRY(line->read_integer("site index"));

        for (auto [index, token] : { std::tuple { from, from_token }, { to, to_token } }) {
            if (index < 0 || index > m_geometry.primitive_cell_sites.size()) {
                return m_diag.error(token, "site index must be in the range [0, {})",
                    m_geometry.primitive_cell_sites.size());
            }
        }

        f64 delta[3] {};
        SourceRange delta_token;
        for (size_t i = 0; i < 3; ++i) {
            auto [value, token] = TRY(line->read_double(std::format("site coordinate delta {}-component", component_names[i])));
            delta[i] = value;
            delta_token = delta_token.combined_with(token);
        }

        bool is_interaction = from == to && delta[0] == 0 && delta[1] == 0 && delta[2] == 0;

        if (is_interaction) {
            auto [mu_up, _1] = TRY(line->read_double("on-site chemical potential shift for spin-up electron"));
            auto [mu_down, _2] = TRY(line->read_double("on-site chemical potential shift for spin-down electron"));
            auto [u, _3] = TRY(line->read_double("on-site interaction energy"));

            m_geometry.hamiltonian.emplace_back(ParsedFreeformGeometry::OnSiteInteraction {
                .site = from,
                .mu_up_offset = mu_up,
                .mu_down_offset = mu_down,
                .u = u,
            });
        } else {
            auto [t_up, _1] = TRY(line->read_double("hopping integral for spin-up electron"));
            auto [t_down, _2] = TRY(line->read_double("hopping integral for spin-down electron"));
            auto [u, u_token] = TRY(line->read_double("on-site interaction energy"));

            if (u != 0) {
                return m_diag.error(u_token, "value of on-site interaction energy parameter must be zero for hopping terms");
            }

            m_geometry.hamiltonian.emplace_back(ParsedFreeformGeometry::Hopping {
                .location = to_token,
                .from = from,
                .to = to,
                .coordinate_delta = { delta[0], delta[1], delta[2] },
                .t_up = t_up,
                .t_down = t_down,
            });
        }
    }
    return {};
}

DiagnosticOr<void> FreeformGeometryParser::parse_symmetries()
{
    if (!m_sections.contains(Section::Symmetry)) {
        return {};
    }

    auto& [lexer, _] = m_sections.at(Section::Symmetry);
    while (true) {
        auto line = lexer.nonempty_line();
        if (!line.has_value()) {
            break;
        }

        auto [type, type_token] = TRY(line->read_string("symmetry type"));
        VERIFY(!type.empty());

        bool type_understood = TRY([&] -> DiagnosticOr<bool> {
            char c = ascii_to_lower(type[0]);
            if (c == 'c' || c == 's') {
                if (type.size() != 2 || !std::isdigit(type[1])) {
                    return false;
                }
                int order = type[1] - '0';
                if (order != 2 && order != 3 && order != 4 && order != 6) {
                    return m_diag.error(type_token.substr(1),
                        "rotational symmetry order must be 2, 3, 4, or 6");
                }
                auto x = TRY(line->read_double("x-coordinate of a point on the axis")).value;
                auto y = TRY(line->read_double("y-coordinate of a point on the axis")).value;
                auto z = TRY(line->read_double("z-coordinate of a point on the axis")).value;
                auto dx = TRY(line->read_double("x-coordinate of axis direction vector")).value;
                auto dy = TRY(line->read_double("y-coordinate of axis direction vector")).value;
                auto dz = TRY(line->read_double("z-coordinate of axis direction vector")).value;
                m_geometry.symmetries.emplace_back(ParsedFreeformGeometry::Rotation {
                    .order = order,
                    .center = { x, y, z },
                    .direction = { dx, dy, dz },
                    .is_rotoreflection = (c == 's'),
                });
                return true;
            } else if (c == 'd') {
                auto x = TRY(line->read_double("x-coordinate of a point on the reflection plane")).value;
                auto y = TRY(line->read_double("y-coordinate of a point on the reflection plane")).value;
                auto z = TRY(line->read_double("z-coordinate of a point on the reflection plane")).value;
                auto dx = TRY(line->read_double("x-coordinate of reflection plane normal vector")).value;
                auto dy = TRY(line->read_double("y-coordinate of reflection plane normal vector")).value;
                auto dz = TRY(line->read_double("z-coordinate of reflection plane normal vector")).value;
                m_geometry.symmetries.emplace_back(ParsedFreeformGeometry::Reflection {
                    .point = { x, y, z },
                    .normal = { dx, dy, dz },
                });
                return true;
            } else if (c == 'i') {
                auto x = TRY(line->read_double("x-coordinate of inversion point")).value;
                auto y = TRY(line->read_double("y-coordinate of inversion point")).value;
                auto z = TRY(line->read_double("z-coordinate of inversion point")).value;
                m_geometry.symmetries.emplace_back(ParsedFreeformGeometry::Inversion {
                    .center = { x, y, z },
                });
                return true;
            }
            return false;
        }());

        if (!type_understood) {
            return m_diag.error(type_token, "unrecognized symmetry type '{}'", type);
        }

        TRY(line->expect_eof());
    }

    return {};
}

} // namespace dqmc::parser
