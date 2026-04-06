#ifndef SUB_FORMATTER_HPP
#define SUB_FORMATTER_HPP

#include <iomanip>
#include <ostream>
#include <string>

#include "fun_modern.hpp"

namespace metacsst::formatter {

inline void write_sub_header(std::ostream& out,
                             const std::string& name,
                             const std::string& sequence) {
  out << name << "\n";
  out << sequence << "\n";
}

inline void write_sub_score(std::ostream& out,
                            float score,
                            int strand,
                            int start,
                            int end,
                            const std::string& match_seq) {
  if (strand == 1) {
    out << "Score:" << std::fixed << std::setprecision(2) << score
        << "\t+\tmatchSeq(" << start << "-" << end << "):" << match_seq << "\n";
  } else {
    const std::string complementary_seq = metacsst::complementary(match_seq);
    out << "Score:" << std::fixed << std::setprecision(2) << score
        << "\t-\tmatchSeq(" << start << "-" << end << "):" << complementary_seq << "\n";
  }
}

}

#endif
