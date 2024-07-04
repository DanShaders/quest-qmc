#include <fcntl.h>
#include <unistd.h>

#include "SRC/common/FileView.h"
#include "SRC/common/Filesystem.h"

namespace dqmc {

namespace {

class File final : public FileView {
public:
    File(std::string&& filename_owner, std::string&& content_owner)
        : m_filename_owner(std::move(filename_owner))
        , m_content_owner(std::move(content_owner))
    {
        m_filename = m_filename_owner;
        m_content = m_content_owner;
    }

private:
    std::string m_filename_owner;
    std::string m_content_owner;
};

} // namespace

ErrorOr<std::shared_ptr<FileView>> Filesystem::read(std::filesystem::path const& path)
{
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd == -1) {
        return Error::from_errno_with_context(errno, "Failed to open '{}' for reading", path.string());
    }

    std::string data;
    int saved_errno = 0;

    while (true) {
        constexpr int min_read_size = 4096;
        size_t old_size = data.size();
        data.resize_and_overwrite(
            std::max(data.size() + min_read_size, data.capacity()),
            [&](char* buffer, size_t buffer_size) {
                ssize_t nread = ::read(fd, buffer + old_size, buffer_size - old_size);
                if (nread == -1) {
                    saved_errno = errno;
                    VERIFY(saved_errno != 0);
                    nread = 0;
                }
                VERIFY(nread >= 0);
                return old_size + static_cast<size_t>(nread);
            });
        if (old_size == data.size()) {
            break;
        }
    }

    if (saved_errno != 0) {
        return Error::from_errno_with_context(saved_errno, "Failed to read '{}'", path.string());
    }

    return std::make_shared<File>(path.string(), std::move(data));
}

} // namespace dqmc
