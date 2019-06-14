#ifndef UTILS_LOGGING_H
#define UTILS_LOGGING_H

#include "system.h"
#include "timer.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

namespace options {
class OptionParser;
}

namespace utils {
/*
  Simple logger that prepends time and peak memory info to messages.
  Logs are written to stdout.

  Usage:
        utils::g_log << "States: " << num_states << endl;
*/
struct Log {
    template<typename T>
    std::ostream &operator<<(const T &elem) const {
        return std::cout << "[t=" << g_timer << ", "
                         << get_peak_memory_in_kb() << " KB] " << elem;
    }
};

extern Log g_log;

// See add_verbosity_option_to_parser for documentation.
enum class Verbosity {
    SILENT,
    NORMAL,
    VERBOSE,
    DEBUG
};

extern void add_verbosity_option_to_parser(options::OptionParser &parser);

class TraceBlock {
    std::string block_name;
public:
    explicit TraceBlock(const std::string &block_name);
    ~TraceBlock();
};

extern void trace(const std::string &msg = "");
}

namespace std {
template<class T>
ostream &operator<<(ostream &stream, const vector<T> &vec) {
    stream << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i != 0)
            stream << ", ";
        stream << vec[i];
    }
    stream << "]";
    return stream;
}

template<class T>
ostream &operator<<(ostream &stream, const set<T> &set_) {
    stream << "{";
    for (typename set<T>::iterator it = set_.begin(); it != set_.end(); ++it) {
        stream << *it << ", ";
    }
    stream << "}";
    return stream;
}

template<class T>
ostream &operator<<(ostream &stream, const pair<T, T> &pair_) {
    stream << "(" << pair_.first << ", " << pair_.second << ")";
    return stream;
}
}

#endif
