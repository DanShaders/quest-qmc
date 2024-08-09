#include <complex>
#include <cstddef>
#include <cstdint>
#include <print>

#include "SRC/common/Variant.h"

using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using f32 = float;
using f64 = double;

using CFI_index_t = i64;
using CFI_rank_t = i8;
using CFI_attribute_t = i8;
using CFI_type_t = i16;

struct CFI_dim_t {
    CFI_index_t lower_bound;
    CFI_index_t extent;
    CFI_index_t sm;
};

struct CFI_cdesc_t {
    void* base_addr;
    size_t elem_len;
    int version;
    CFI_rank_t rank;
    CFI_attribute_t attribute;
    CFI_type_t type;
    CFI_dim_t dim[];
};

namespace dqmc {

inline void assertion_failed(char const* expr, char const* file, int line)
{
    std::println(stderr, "Assertion failed: `{}` evaluated to false at {}:{}!\n", expr, file, line);
    abort();
}

inline std::string_view as_string_view(CFI_cdesc_t* desc)
{
    return { reinterpret_cast<char*>(desc->base_addr), desc->elem_len };
}

inline char ascii_to_lower(char c)
{
    if ('A' <= c && c <= 'Z') {
        return c - 'A' + 'a';
    }
    return c;
}

struct Empty { };

} // namespace dqmc

#define VERIFY(expr)                                             \
    do {                                                         \
        if (!(expr)) [[unlikely]] {                              \
            ::dqmc::assertion_failed(#expr, __FILE__, __LINE__); \
        }                                                        \
    } while (false)

#define VERIFY_NOT_REACHED() VERIFY(false)

#define TRY(expr)                                                 \
    ({                                                            \
        auto&& _value = (expr);                                   \
        if (!_value.has_value())                                  \
            return std::unexpected { std::move(_value).error() }; \
        std::move(_value).value();                                \
    })

#define MUST(expr)                                     \
    ({                                                 \
        auto&& _value = (expr);                        \
        if (!_value.has_value()) [[unlikely]] {        \
            std::println("Error: {}", _value.error()); \
            VERIFY_NOT_REACHED();                      \
        }                                              \
        std::move(_value).value();                     \
    })

using namespace std::literals;
