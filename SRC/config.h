#pragma once

#include "SRC/common.h"

#include <expected>
#include <map>

namespace dqmc {

struct ParameterValue {
    std::string value;
    int line_number = 0;
    bool is_used = false;

    bool operator==(ParameterValue const&) const = default;
};

class ConfigParameters {
public:
    static std::expected<std::unique_ptr<ConfigParameters>, std::string> create(std::istream&);

    std::expected<int, std::string> read_integer(std::string_view key, std::optional<int> default_value = std::nullopt) const;
    std::expected<std::string_view, std::string> read_string(std::string_view key, std::optional<std::string_view> default_value = std::nullopt) const;
    std::expected<bool, std::string> read_boolean(std::string_view key, std::optional<bool> default_value = std::nullopt) const;
    std::expected<f64, std::string> read_double(std::string_view key, std::optional<f64> default_value = std::nullopt) const;
    std::expected<std::vector<f64>, std::string> read_double_array(std::string_view key, std::optional<std::vector<f64>> default_value = std::nullopt) const;

    decltype(auto) parameters(this auto&& self) { return std::forward_like<decltype(self)>(self.m_parameters); }

private:
    ConfigParameters() = default;

    std::map<std::string, ParameterValue, std::less<>> m_parameters;
};

} // namespace dqmc
