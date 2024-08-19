#pragma once

namespace dqmc {

inline char ascii_to_lower(char c)
{
    if (c >= 'A' && c <= 'Z') {
        return c + ('a' - 'A');
    } else {
        return c;
    }
}

} // namespace dqmc
