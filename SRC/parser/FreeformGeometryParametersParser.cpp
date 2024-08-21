#include "SRC/parser/FreeformGeometryParametersParser.h"
#include "SRC/common/Enumerate.h"

namespace dqmc::parser {

DiagnosticOr<void> FreeformGeometryParametersParser::parse(ConfigParser& parser, DiagnosticEngine& diag)
{
    auto mu_up_token = TRY(parser.claim_double("mu_up", 0));
    m_parameters.mu_up = mu_up_token.value;

    auto mu_down_token = TRY(parser.claim_double("mu_dn", 0));
    m_parameters.mu_down = mu_down_token.value;

    auto phase_shift = TRY(parser.claim_double_array("bcond", std::vector<double> { 0, 0, 0 }));
    if (phase_shift.value.size() != 3) {
        auto location = phase_shift.value.size() < 3 ? phase_shift.end_of_array : phase_shift.locations[3];
        return diag.error(location, "expected an array with 3 elements for bcond");
    }
    for (auto [i, element] : enumerate(phase_shift.value)) {
        if (std::round(element) != element) {
            return diag.error(phase_shift.locations[i],
                "only integer multiples of pi for phase shifts on boundaries are allowed");
        }
        bool is_odd = std::fmod(element, 2) != 0;
        m_parameters.should_negate_phase[i] = is_odd;
    }
    return {};
}

} // namespace dqmc::parser
