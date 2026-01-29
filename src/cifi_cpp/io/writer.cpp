#include "writer.hpp"
#include <stdexcept>
#include <algorithm>

namespace cifi {

// PlainFastqWriter

PlainFastqWriter::PlainFastqWriter(const std::string& path) : out_(path) {
    if (!out_) {
        throw std::runtime_error("Cannot open for writing: " + path);
    }
}

PlainFastqWriter::~PlainFastqWriter() {
    close();
}

void PlainFastqWriter::write(const std::string& name,
                              const std::string& seq,
                              const std::string& qual) {
    out_ << '@' << name << '\n' << seq << "\n+\n" << qual << '\n';
}

void PlainFastqWriter::close() {
    if (out_.is_open()) {
        out_.close();
    }
}

// GzipFastqWriter

GzipFastqWriter::GzipFastqWriter(const std::string& path) {
    gz_ = gzopen(path.c_str(), "wb");
    if (!gz_) {
        throw std::runtime_error("Cannot open for gzip writing: " + path);
    }
}

GzipFastqWriter::~GzipFastqWriter() {
    close();
}

void GzipFastqWriter::write(const std::string& name,
                             const std::string& seq,
                             const std::string& qual) {
    gzprintf(gz_, "@%s\n%s\n+\n%s\n",
             name.c_str(), seq.c_str(), qual.c_str());
}

void GzipFastqWriter::close() {
    if (gz_) {
        gzclose(gz_);
        gz_ = nullptr;
    }
}

// Factory

bool ends_with_gz(const std::string& path) {
    if (path.size() < 3) return false;
    std::string suffix = path.substr(path.size() - 3);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    return suffix == ".gz";
}

std::unique_ptr<FastqWriter> make_writer(const std::string& path, bool force_gzip) {
    bool use_gzip = force_gzip || ends_with_gz(path);

    if (use_gzip) {
        std::string gz_path = ends_with_gz(path) ? path : path + ".gz";
        return std::make_unique<GzipFastqWriter>(gz_path);
    }
    return std::make_unique<PlainFastqWriter>(path);
}

} // namespace cifi
