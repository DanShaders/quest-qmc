#pragma once

#include "SRC/parser/Config.h"

namespace dqmc::parser {

struct ParsedFreeformGeometryParameters {
    f64 mu_up;
    f64 mu_down;
    bool should_negate_phase[3];
};

class FreeformGeometryParametersParser final : public ParametersParser {
public:
    virtual DiagnosticOr<void> parse(ConfigParser& parser, DiagnosticEngine& diag) override;

    ParsedFreeformGeometryParameters const& parameters() const& { return m_parameters; }

private:
    ParsedFreeformGeometryParameters m_parameters {};
};

} // namespace dqmc::parser
