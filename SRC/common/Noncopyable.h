#pragma once

#define DQMC_MAKE_NONCOPYABLE(c) \
private:                         \
    c(c const&) = delete;        \
    c& operator=(c const&) = delete

#define DQMC_MAKE_NONMOVABLE(c) \
private:                        \
    c(c&&) = delete;            \
    c& operator=(c&&) = delete

#define DQMC_MAKE_DEFAULT_MOVABLE(c) \
public:                              \
    c(c&&) = default;                \
    c& operator=(c&&) = default

#define DQMC_MAKE_CONSTEXPR_MOVABLE(c) \
public:                                \
    constexpr c(c&&) = default;        \
    constexpr c& operator=(c&&) = default

#define DQMC_MAKE_DEFAULT_COPYABLE(c) \
public:                               \
    c(c const&) = default;            \
    c& operator=(c const&) = default

#define DQMC_MAKE_CONSTEXPR_COPYABLE(c) \
public:                                 \
    constexpr c(c const&) = default;    \
    constexpr c& operator=(c const&) = default
