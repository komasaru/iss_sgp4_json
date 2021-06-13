#ifndef ISS_SGP4_JSON_EOP_HPP_
#define ISS_SGP4_JSON_EOP_HPP_

#include "time.hpp"

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace iss_sgp4_json {

class Eop {

public:
  Eop(struct timespec);  // コンストラクタ
  std::string eop;       // EOP データ
  double pm_x;           // 極運動(x)
  double pm_y;           // 極運動(y)
  double dut1;           // DUT1
  double lod;            // LOD

private:
  std::string get_eop(struct timespec);  // EOP データ取得
};

}  // namespace iss_sgp4_json

#endif

