#include "SRC/parser/SimulationParametersParser.h"

namespace dqmc::parser {

DiagnosticOr<void> HubbardModelParametersParser::parse(ConfigParser& parser, DiagnosticEngine& diag)
{
    m_parameters.time_slices = TRY(parser.claim_integer<u32>("l")).value;
    m_parameters.dtau = TRY(parser.claim_double("dtau")).value;
    m_parameters.warm_up_steps = TRY(parser.claim_integer<u64>("nwarm")).value;
    m_parameters.measurement_steps = TRY(parser.claim_integer<u64>("npass")).value;
    return {};
}

} // namespace dqmc::parser
