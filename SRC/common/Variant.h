#pragma once

#include <variant>

namespace dqmc {

namespace detail {

template<typename... Ts>
struct OverloadSet : Ts... {
    using Ts::operator()...;
};

template<typename... Ts>
OverloadSet(Ts...) -> OverloadSet<Ts...>;

} // namespace detail

template<typename Variant, typename... Visitors>
decltype(auto) visit(Variant&& variant, Visitors&&... visitors)
{
    return std::visit(detail::OverloadSet { std::forward<Visitors>(visitors)... }, std::forward<Variant>(variant));
}

} // namespace dqmc
