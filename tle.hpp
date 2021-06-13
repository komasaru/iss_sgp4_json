#ifndef ISS_SGP4_JSON_TLE_HPP_
#define ISS_SGP4_JSON_TLE_HPP_

#include "time.hpp"

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace iss_sgp4_json {

class Tle {
public:
  std::vector<std::string> tle;        // TLE
  Tle(struct timespec);                // コンストラクタ
  std::vector<std::string> get_tle();  // TLE 読み込み

private:
  struct timespec ut1;  // UT1
};

}  // namespace iss_sgp4_json

#endif

