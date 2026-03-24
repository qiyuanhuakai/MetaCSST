#ifndef MAIN_FORMATTER_HPP
#define MAIN_FORMATTER_HPP

#include <iomanip>
#include <ostream>
#include <string>

#include "fun_modern.hpp"

namespace metacsst::formatter {

inline const char* main_feature_label(int feature_type) {
  if (feature_type == 1) {
    return "TR";
  }
  if (feature_type == 3) {
    return "RT";
  }
  return "";
}

inline void write_main_feature(std::ostream& out,
                               const std::string& name,
                               int feature_type,
                               float score,
                               int strand,
                               int start,
                               int end,
                               const std::string& match_seq) {
  out << name << "\t" << main_feature_label(feature_type) << "\t";
  if (strand == 1) {
    out << std::fixed << std::setprecision(2) << score << "\t+\t" << start << "\t" << end << "\t" << match_seq << "\n";
  } else {
    const std::string complementary_seq = metacsst::complementary(match_seq);
    out << std::fixed << std::setprecision(2) << score << "\t-\t" << start << "\t" << end << "\t" << complementary_seq << "\n";
  }
}

inline void write_main_dgr(std::ostream& out,
                           const std::string& name,
                           float score,
                           int start,
                           int end,
                           const std::string& match_seq) {
  out << name << "\tDGR\t" << std::fixed << std::setprecision(2) << score << "\t*\t" << start << "\t" << end << "\t" << match_seq << "\n";
}

}

#endif
