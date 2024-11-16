#include <fast_float/fast_float.h>
#include <fstream>
#include <ranges>
#include <utility>

#include "SRC/config.h"

namespace dqmc {

std::expected<int, Empty> parse_integer(std::string_view value);
std::expected<bool, Empty> parse_boolean(std::string_view value);
std::expected<double, Empty> parse_double(std::string_view value);
std::expected<std::vector<double>, std::string_view> parse_double_array(std::string_view value);

} // namespace dqmc

namespace {

struct config {
    std::unique_ptr<dqmc::ConfigParameters> config_parameters;
};
static_assert(sizeof(config) == 8);

std::map<std::string_view, std::string_view, std::less<>> const config_default_values = {
    { "hsf"sv, "-1"sv },
    { "hsfin"sv, "HSF.in"sv },
    { "hsfout"sv, "HSF.out"sv },
    { "hsftype"sv, "0"sv },
    { "l"sv, "12"sv },
    { "u"sv, "0"sv },
    { "accept"sv, "0"sv },
    { "bcond"sv, "0,0,0"sv },
    { "debug"sv, "0"sv },
    { "delta1"sv, "1"sv },
    { "delta2"sv, "1"sv },
    { "difflim"sv, "0.001"sv },
    { "dmu"sv, "0"sv },
    { "dtau"sv, "0.125"sv },
    { "errrate"sv, "0.001"sv },
    { "fixwrap"sv, "0"sv },
    { "gamma"sv, "0"sv },
    { "gfile"sv, "geom.def"sv },
    { "mu_dn"sv, "0"sv },
    { "mu_up"sv, "0"sv },
    { "n"sv, "16"sv },
    { "nbin"sv, "10"sv },
    { "nhist"sv, "0"sv },
    { "nitvl"sv, "4"sv },
    { "north"sv, "12"sv },
    { "npass"sv, "5000"sv },
    { "ntry"sv, "0"sv },
    { "nwarm"sv, "1000"sv },
    { "nwrap"sv, "12"sv },
    { "nx"sv, "4"sv },
    { "ny"sv, "4"sv },
    { "nz"sv, "2"sv },
    { "ofile"sv, "quest"sv },
    { "reject"sv, "0"sv },
    { "seed"sv, "0"sv },
    { "ssxx"sv, "0"sv },
    { "t_dn"sv, "1"sv },
    { "t_up"sv, "1"sv },
    { "tausk"sv, "10"sv },
    { "tdm"sv, "0"sv },
};

} // namespace

extern "C" void detour_dqmc_cfg_dqmc_readln(CFI_cdesc_t* str, i32* ipt, i32* status);
extern "C" void fortran_dqmc_cfg_dqmc_readln(CFI_cdesc_t* str, i32* ipt, i32* status)
{
    return detour_dqmc_cfg_dqmc_readln(str, ipt, status);
}

extern "C" void fortran_dqmc_cfg_dqmc_default_def(config* cfg)
{
}

extern "C" void fortran_dqmc_cfg_dqmc_config_free(config* cfg)
{
    cfg->config_parameters.~unique_ptr();
}

extern "C" void binding_config_filename(CFI_cdesc_t* filename);
extern "C" void binding_free_config_filename(CFI_cdesc_t* filename);

extern "C" void fortran_dqmc_cfg_dqmc_read_config(config* cfg)
{
    CFI_cdesc_t filename_descriptor {};
    binding_config_filename(&filename_descriptor);

    std::string_view config_filename = dqmc::as_string_view(&filename_descriptor);

    std::ifstream config_reader { config_filename };
    VERIFY(config_reader.good());

    auto configuration_or_error = dqmc::ConfigParameters::create(config_reader);
    if (!configuration_or_error.has_value()) {
        std::println("ConfigurationFile::create: {}", configuration_or_error.error());
        VERIFY_NOT_REACHED();
    }
    cfg->config_parameters = std::move(configuration_or_error.value());

    binding_free_config_filename(&filename_descriptor);
}

static std::string key_from_descriptor(CFI_cdesc_t* name)
{
    return dqmc::as_string_view(name)
        | std::views::transform(dqmc::ascii_to_lower)
        | std::ranges::to<std::string>();
}

extern "C" void fortran_dqmc_cfg_dqmc_config_seti(config* cfg, CFI_cdesc_t* name, i32* value)
{
    cfg->config_parameters->parameters()[key_from_descriptor(name)] = {
        .value = std::to_string(*value),
    };
}

extern "C" void fortran_dqmc_cfg_dqmc_config_setr(config* cfg, CFI_cdesc_t* name, f64* value)
{
    cfg->config_parameters->parameters()[key_from_descriptor(name)] = {
        .value = std::format("{}", *value),
    };
}

extern "C" void fortran_dqmc_cfg_dqmc_config_setpr(config* cfg, CFI_cdesc_t* name, i32* n, f64* value)
{
    cfg->config_parameters->parameters()[key_from_descriptor(name)] = {
        .value = std::format("{:n}", std::span<f64> { value, value + *n }),
    };
}

extern "C" void fortran_dqmc_cfg_dqmc_config_sets(config* cfg, CFI_cdesc_t* name, CFI_cdesc_t* value)
{
    cfg->config_parameters->parameters()[key_from_descriptor(name)] = {
        .value = std::string { dqmc::as_string_view(value) },
    };
}

extern "C" i32 fortran_dqmc_cfg_dqmc_config_isset(config* cfg, CFI_cdesc_t* name)
{
    return cfg->config_parameters->parameters().contains(key_from_descriptor(name));
}

extern "C" void fortran_dqmc_cfg_dqmc_config_geti(config* cfg, CFI_cdesc_t* name, i32* value)
{
    auto key = key_from_descriptor(name);
    auto const& parameters = cfg->config_parameters->parameters();

    if (auto it = parameters.find(key); it != parameters.end()) {
        *value = MUST(cfg->config_parameters->read_integer(key_from_descriptor(name)));
    } else {
        auto default_value = config_default_values.at(key);
        std::println(stderr, "Using default value '{}' for an uninitialized key '{}'.", default_value, key);
        *value = dqmc::parse_integer(default_value).value();
    }
}

extern "C" void fortran_dqmc_cfg_dqmc_config_getr(config* cfg, CFI_cdesc_t* name, f64* value)
{
    auto key = key_from_descriptor(name);
    auto const& parameters = cfg->config_parameters->parameters();

    if (auto it = parameters.find(key); it != parameters.end()) {
        *value = MUST(cfg->config_parameters->read_double(key_from_descriptor(name)));
    } else {
        auto default_value = config_default_values.at(key);
        std::println(stderr, "Using default value '{}' for an uninitialized key '{}'.", default_value, key);
        *value = dqmc::parse_double(default_value).value();
    }
}

extern "C" void binding_allocate_double_array(i64 n, CFI_cdesc_t* descriptor);

extern "C" void fortran_dqmc_cfg_dqmc_config_getpr(config* cfg, CFI_cdesc_t* name, i32* n, CFI_cdesc_t* value)
{
    auto key = key_from_descriptor(name);
    auto const& parameters = cfg->config_parameters->parameters();
    std::vector<f64> result;

    if (auto it = parameters.find(key); it != parameters.end()) {
        result = MUST(cfg->config_parameters->read_double_array(key_from_descriptor(name)));
    } else {
        auto default_value = config_default_values.at(key);
        std::println(stderr, "Using default value '{}' for an uninitialized key '{}'.", default_value, key);
        result = dqmc::parse_double_array(default_value).value();
    }

    VERIFY(std::in_range<i32>(result.size()));
    *n = static_cast<i32>(result.size());
    binding_allocate_double_array(*n, value);
    std::copy(result.begin(), result.end(), reinterpret_cast<f64*>(value->base_addr));
}

extern "C" void fortran_dqmc_cfg_dqmc_config_gets(config* cfg, CFI_cdesc_t* name, CFI_cdesc_t* value)
{
    auto key = key_from_descriptor(name);
    auto const& parameters = cfg->config_parameters->parameters();

    std::string_view result;

    if (auto it = parameters.find(key); it != parameters.end()) {
        result = MUST(cfg->config_parameters->read_string(key_from_descriptor(name)));
    } else {
        result = config_default_values.at(key);
        std::println(stderr, "Using default value '{}' for an uninitialized key '{}'.", result, key);
    }

    if (!std::cmp_less(result.size(), value->elem_len)) {
        std::println(stderr, "Fortran string is too small to store the value of key '{}'.", key);
        VERIFY_NOT_REACHED();
    }

    char* data_start = reinterpret_cast<char*>(value->base_addr);
    std::copy(result.begin(), result.end(), data_start);
    std::fill(data_start + result.size(), data_start + value->elem_len, ' ');
}
