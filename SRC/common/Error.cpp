#include <print>

#include "SRC/common/Error.h"

namespace dqmc::detail {

void verify_assertion_failed(char const* expression, char const* file, int line)
{
    std::println(stderr, "Assertion failed: `{}` evaluated to false at {}:{}!", expression, file, line);
    abort();
}

void reachability_assertion_failed(char const* file, int line)
{
    std::println(stderr, "Reached VERIFY_NOT_REACHED() statement at {}:{}!", file, line);
    abort();
}

void must_assertion_failed(char const* expression, char const* file, int line)
{
    std::println(stderr, "Assertion failed: `{}` did not have a value at {}:{}!", expression, file, line);
    abort();
}

} // namespace dqmc::detail
