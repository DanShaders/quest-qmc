#pragma once

#include "SRC/parser/Lexer.h"

namespace dqmc::parser {

struct ParsedFreeformGeometry {
    struct PrimitiveCellSite {
        std::string label;
        f64 displacement[3] = { 0, 0, 0 };
    };

    int dimensions = 0;

    f64 lattice_basis[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };
    SourceRange lattice_basis_location;

    int supercell_basis[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };
    SourceRange supercell_basis_location;

    std::vector<PrimitiveCellSite> primitive_cell_sites;
};

class FreeformGeometryParser {
public:
    FreeformGeometryParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag);

    std::optional<ParsedFreeformGeometry> parse();

private:
    enum class Section {
        NumberOfDimensions,
        LatticeBasis,
        SupercellBasis,
        PrimitiveCellSites,
        HAMILT,
        SYMM,
        PHASE,
        BONDS,
        PAIR,
        DILUT,
    };

    static constexpr auto section_headers = std::to_array({
        "#NDIM"sv,
        "#PRIM"sv,
        "#SUPER"sv,
        "#ORB"sv,
        "#HAMILT"sv,
        "#SYMM"sv,
        "#PHASE"sv,
        "#BONDS"sv,
        "#PAIR"sv,
        "#DILUT"sv,
    });

    void split_into_sections();
    void parse_preamble(Lexer section);
    std::expected<void, Empty> parse_number_of_dimensions();
    std::expected<void, Empty> parse_lattice_basis();
    std::expected<void, Empty> parse_supercell_basis();
    std::expected<void, Empty> parse_primitive_cell_sites();

    std::shared_ptr<FileView> m_file;
    std::map<Section, std::pair<Lexer, SourceRange>> m_sections;
    ParsedFreeformGeometry m_geometry;
    DiagnosticEngine& m_diag;
};

} // namespace dqmc::parser
