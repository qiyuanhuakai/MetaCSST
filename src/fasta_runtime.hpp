#ifndef FASTA_RUNTIME_HPP
#define FASTA_RUNTIME_HPP

#include <algorithm>
#include <iostream>
#include <string>

#include "fun_modern.hpp"

namespace metacsst::app {

inline std::string extract_fasta_name(const std::string& header, bool stop_at_bracket) {
  if (header.empty() || header[0] != '>') {
    return "";
  }

  std::size_t end_pos = header.find_first_of(" \n");
  if (stop_at_bracket) {
    const std::size_t bracket_pos = header.find('[');
    end_pos = std::min(end_pos, bracket_pos);
  }

  if (end_pos == std::string::npos || end_pos <= 1) {
    return header.substr(1);
  }
  return header.substr(1, end_pos - 1);
}

template <typename SequenceHandler>
inline int stream_fasta_sequences(const std::string& input_path, bool stop_at_bracket, SequenceHandler handler) {
  metacsst::LineReader in(input_path);
  if (!in.is_open()) {
    std::cerr << "Error: Cannot open input file: " << input_path << std::endl;
    return 1;
  }

  std::string line;
  std::string name;
  while (in.getline(line)) {
    if (!line.empty() && line[0] == '>') {
      name = extract_fasta_name(line, stop_at_bracket);
      continue;
    }
    if (!line.empty()) {
      handler(name, line);
    }
  }
  return 0;
}

}

#endif
