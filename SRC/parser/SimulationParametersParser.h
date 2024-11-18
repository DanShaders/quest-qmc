#pragma once

#include "SRC/parser/Config.h"
#include "SRC/simulation/SimulationParameters.h"

namespace dqmc::parser {

class HubbardModelParametersParser final : public ParametersParser {
public:
    virtual DiagnosticOr<void> parse(ConfigParser& parser, DiagnosticEngine& diag) override;

    simulation::HubbardModelParameters const& parameters() const& { return m_parameters; }

private:
    simulation::HubbardModelParameters m_parameters {};
};

} // namespace dqmc::parser
