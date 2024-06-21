#pragma once

#include <expected>
#include <format>
#include <string.h>
#include <string>

#include "SRC/common/Noncopyable.h"

namespace dqmc {

class Error {
    DQMC_MAKE_DEFAULT_COPYABLE(Error);
    DQMC_MAKE_DEFAULT_MOVABLE(Error);

public:
    template<typename... Args>
    static std::unexpected<Error> formatted(std::format_string<Args...> fmt, Args&&... args)
    {
        return std::unexpected { Error { std::format(fmt, std::forward<Args>(args)...) } };
    }

    static std::unexpected<Error> from_errno_with_context(int errno_, std::string const& context)
    {
        auto message = std::format("{}: {} ({}, {})", context, strerrordesc_np(errno_), strerrorname_np(errno_), errno_);
        return std::unexpected { Error { message } };
    }

    template<typename... Args>
    static std::unexpected<Error> from_errno_with_context(int errno_, std::format_string<Args...> fmt, Args&&... args)
    {
        auto context = std::format(fmt, std::forward<Args>(args)...);
        return from_errno_with_context(errno_, context);
    }

    std::string const& message() const& { return m_message; }
    std::string&& message() && { return std::move(m_message); }

private:
    Error(std::string message)
        : m_message(std::move(message))
    {
    }

    std::string m_message;
};

struct Empty {
    static std::unexpected<Empty> error() { return std::unexpected<Empty> { {} }; }
};

template<typename T, typename E = Error>
using ErrorOr = std::expected<T, E>;

namespace detail {

void verify_assertion_failed(char const* expression, char const* file, int line);
void reachability_assertion_failed(char const* file, int line);
void must_assertion_failed(char const* expression, char const* file, int line);

} // namespace detail

#define VERIFY(expr)                                                            \
    do {                                                                        \
        if (!(expr)) [[unlikely]] {                                             \
            ::dqmc::detail::verify_assertion_failed(#expr, __FILE__, __LINE__); \
        }                                                                       \
    } while (false)

#define VERIFY_NOT_REACHED() ::dqmc::detail::reachability_assertion_failed(__FILE__, __LINE__)

#define TRY(expr)                                                 \
    ({                                                            \
        auto&& _value = (expr);                                   \
        if (!_value.has_value())                                  \
            return std::unexpected { std::move(_value).error() }; \
        std::move(_value).value();                                \
    })

#define MUST(expr)                                                            \
    ({                                                                        \
        auto&& _value = (expr);                                               \
        if (!_value.has_value()) [[unlikely]] {                               \
            std::println(stderr, "{}", _value.error().message());             \
            ::dqmc::detail::must_assertion_failed(#expr, __FILE__, __LINE__); \
        }                                                                     \
        std::move(_value).value();                                            \
    })

} // namespace dqmc
