#include <complex>
#include <cstddef>
#include <cstdint>
#include <print>

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

inline void assertion_failed(char const* expr, char const* file, int line) {
    std::println(stderr, "Assertion failed: `{}` evaluated to false at {}:{}!\n", expr, file, line);
    abort();
}

#define VERIFY(expr)                                     \
    do {                                                 \
        if (!(expr)) [[unlikely]] {                      \
            assertion_failed(#expr, __FILE__, __LINE__); \
        }                                                \
    } while (false)

#define VERIFY_NOT_REACHED() VERIFY(false)
