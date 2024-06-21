#pragma once

#include "SRC/common/FileView.h"

namespace dqmc {

class Filesystem final : public FilesystemView {
public:
    ErrorOr<std::shared_ptr<FileView>> read(std::filesystem::path const& path);
};

} // namespace dqmc
