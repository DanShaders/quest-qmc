#include "SRC/parser/FreeformGeometryParametersParser.h"

namespace dqmc::parser {

std::expected<void, Empty> FreeformGeometryParametersParser::parse(ConfigParser& parser, DiagnosticEngine& diag)
{
    auto mu_up_token = TRY(parser.claim_double("mu_up", 0));
    m_parameters.mu_up = mu_up_token.value;

    auto mu_down_token = TRY(parser.claim_double("mu_dn", 0));
    m_parameters.mu_down = mu_down_token.value;

    auto phase_shift = TRY(parser.claim_double_array("bcond", std::vector<double> { 0, 0, 0 }));
    if (phase_shift.value.size() != 3) {
        auto location = phase_shift.value.size() < 3 ? phase_shift.end_of_array : phase_shift.locations[3];
        diag.error(location, "expected an array with 3 elements for bcond");
        return Empty::error();
    }
    for (int i = 0; i < 3; ++i) {
        m_parameters.phase_shift_for_supercell_vector[i] = phase_shift.value[i];
    }
    return {};
}

} // namespace dqmc::parser
