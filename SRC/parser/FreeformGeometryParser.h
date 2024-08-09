#pragma once

#include "SRC/common/Variant.h"
#include "SRC/parser/Lexer.h"

namespace dqmc::parser {

struct ParsedFreeformGeometry {
    struct PrimitiveCellSite {
        std::string label;
        f64 displacement[3] = {};
    };

    struct Hopping {
        int from;
        int to;
        f64 coordinate_delta[3];
        f64 t_up;
        f64 t_down;
    };

    struct OnSiteInteraction {
        int site;
        f64 mu_up_offset;
        f64 mu_down_offset;
        f64 u;
    };

    int dimensions;

    // TODO: Check if we can just use 1 for basis and fake 1e3 in
    //       FreeformGeometry::legacy_compatible_format_into.
    f64 lattice_basis[3][3] = {
        { 1e3, 0, 0 },
        { 0, 1e3, 0 },
        { 0, 0, 1e3 },
    };
    SourceRange lattice_basis_location;

    int supercell_basis[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };
    SourceRange supercell_basis_location;

    std::vector<PrimitiveCellSite> primitive_cell_sites;

    std::vector<std::variant<Hopping, OnSiteInteraction>> hamiltonian;
};

class FreeformGeometryParser {
public:
    FreeformGeometryParser(std::shared_ptr<FileView> file, DiagnosticEngine& diag);

    DiagnosticOr<ParsedFreeformGeometry> parse();

private:
    enum class Section {
        NumberOfDimensions,
        LatticeBasis,
        SupercellBasis,
        PrimitiveCellSites,
        Hamiltonian,
        SYMM,
        PHASE,
        BONDS,
        PAIR,
        DILUT,
    };

    struct SplitSection {
        Lexer lexer;
        SourceRange header_location;
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
    DiagnosticOr<void> parse_number_of_dimensions();
    DiagnosticOr<void> parse_lattice_basis();
    DiagnosticOr<void> parse_supercell_basis();
    DiagnosticOr<void> parse_primitive_cell_sites();
    DiagnosticOr<void> parse_hamiltonian();

    std::shared_ptr<FileView> m_file;
    std::string_view m_data;
    std::map<Section, SplitSection> m_sections;
    ParsedFreeformGeometry m_geometry;
    DiagnosticEngine& m_diag;
};

} // namespace dqmc::parser
