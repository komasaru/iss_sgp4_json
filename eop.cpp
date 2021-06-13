#include "eop.hpp"

namespace iss_sgp4_json {

// 定数
static constexpr char         kFEop[]    = "eop.txt";
//static constexpr unsigned int kSecDay =  86400;  // Seconds in a day

/*
 * @brief      コンストラクタ
 *
 * @param[in]  UTC (timespec)
 */
Eop::Eop(struct timespec utc) {
  std::string lod_str;
  eop = get_eop(utc);
  pm_x    = stod(eop.substr(22,  9));
  pm_y    = stod(eop.substr(41,  9));
  dut1    = stod(eop.substr(62, 10));
  lod_str =      eop.substr(83,  7);
  lod = 0.0;
  if (lod_str != "       ") { lod = stod(lod_str); }
}

/*
 * @brief   EOP データ取得
 *
 * @param[in]  UTC (timespec)
 * @return     EOP データ (string)
 */
std::string Eop::get_eop(struct timespec utc) {
  std::string f(kFEop);  // ファイル名
  std::string str_utc;   // 対象の UTC 文字列
  std::string eop = "";  // EOP データ

  try {
    // 対象の UTC 年月日
    str_utc = gen_time_str(utc).substr(0, 10);

    // ファイル OPEN
    std::ifstream ifs(f);
    if (!ifs) throw;  // 読み込み失敗

    // ファイル READ
    while (getline(ifs, eop)) {
      if (eop.substr(0, 10) == str_utc) { break; }
    }
  } catch (...) {
    throw;
  }

  return eop;
}

}  // namespace iss_sgp4_json

