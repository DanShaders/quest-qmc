#pragma once

#include <filesystem>

#include "SRC/common/Error.h"

namespace dqmc {

class FileView {
    DQMC_MAKE_NONCOPYABLE(FileView);
    DQMC_MAKE_NONMOVABLE(FileView);

public:
    virtual ~FileView() = default;

    std::string_view filename() const { return m_filename; }
    std::string_view content() const { return m_content; }

protected:
    FileView() = default;

    std::string_view m_filename;
    std::string_view m_content;
};

class FilesystemView {
public:
    virtual ~FilesystemView() = default;

    virtual ErrorOr<std::shared_ptr<FileView>> read(std::filesystem::path const& path) = 0;
};

} // namespace dqmc
