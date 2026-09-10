#include "writer.hpp"
#include <stdexcept>
#include <algorithm>
#include <limits>

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
    // gzprintf() formats into a fixed internal buffer (8KB by default) and
    // silently writes nothing when the record does not fit. CiFi fragments run
    // to tens of kb, so build the record ourselves and hand it to gzwrite,
    // which has no length limit.
    buf_.clear();
    buf_.reserve(name.size() + seq.size() + qual.size() + 6);
    buf_ += '@';
    buf_ += name;
    buf_ += '\n';
    buf_ += seq;
    buf_ += "\n+\n";
    buf_ += qual;
    buf_ += '\n';

    if (buf_.size() > static_cast<size_t>(std::numeric_limits<unsigned>::max())) {
        throw std::runtime_error("FASTQ record too large to write: " + name);
    }

    int written = gzwrite(gz_, buf_.data(), static_cast<unsigned>(buf_.size()));
    if (written != static_cast<int>(buf_.size())) {
        int err = 0;
        const char* msg = gzerror(gz_, &err);
        throw std::runtime_error("Failed writing FASTQ record " + name + ": " +
                                 (msg ? msg : "short write"));
    }
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
