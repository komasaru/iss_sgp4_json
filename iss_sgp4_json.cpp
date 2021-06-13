/***********************************************************
  指定日時から48時間分の ISS 位置／速度を計算し、 JSON 出力
  : 但し、座標系は BLH(Beta(Latitude), Lambda(Longitude), Height)

    DATE        AUTHOR       VERSION
    2021.06.10  mk-mode.com  1.00 新規作成

  Copyright(C) 2021 mk-mode.com All Rights Reserved.
  ---
  引数 : JST（日本標準時）
           書式：最大23桁の数字
                 （先頭から、西暦年(4), 月(2), 日(2), 時(2), 分(2), 秒(2),
                             1秒未満(9)（小数点以下9桁（ナノ秒）まで））
                 無指定なら現在(システム日時)と判断。
  ---
  MEMO:
    TEME: True Equator, Mean Equinox; 真赤道面平均春分点
     PEF: Pseudo Earth Fixed; 擬地球固定座標系
    ECEF: Earth Centered, Earth Fixed; 地球中心・地球固定直交座標系
***********************************************************/
#include "blh.hpp"
#include "eop.hpp"
#include "sgp4.hpp"
#include "time.hpp"
#include "tle.hpp"

#include <cstdlib>   // for EXIT_XXXX
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace iss_sgp4_json {

static constexpr char         kFOut[] = "iss.json";  // 書き込みファイル
static constexpr unsigned int kDay    = 2;           // 計算日数(日)
static constexpr unsigned int kSec    = 10;          // 計算間隔(秒)
static constexpr double       kSecD   = 86400.0;     // 秒数(1日分)
}

int main(int argc, char* argv[]) {
  namespace ns = iss_sgp4_json;
  std::string     f(ns::kFOut);  // 書き込みファイル
  std::string     tm_str;        // time string
  unsigned int    s_tm;          // size of time string
  unsigned int    i;             // loop index
  unsigned int    j;             // loop index
  int             s_nsec;        // size of nsec string
  int             ret;           // return of functions
  struct          tm t = {};     // for work
  struct timespec jst;           // JST
  struct timespec utc;           // UTC
  struct timespec ut1;           // UT1
  struct timespec tai;           // TAI
  struct timespec jst_wk;        // JST(作業用)
  struct timespec utc_wk;        // UTC(作業用)
  struct timespec ut1_wk;        // UT1(作業用)
  struct timespec tai_wk;        // TAI(作業用)
  double          pm_x;          // 極運動(x)
  double          pm_y;          // 極運動(y)
  double          lod;           // LOD
  std::string     lod_str;       // LOD 文字列
  std::vector<std::string> tle;  // TLE
  ns::Satellite   sat;           // 衛星情報
  ns::PvTeme      teme;          // 位置・速度(TEME)
  ns::PvBlh       blh;           // 位置・速度(BLH)

  try {
    // 現在日時(UT1) 取得
    if (argc > 1) {
      // コマンドライン引数より取得
      tm_str = argv[1];
      s_tm = tm_str.size();
      if (s_tm > 23) { 
        std::cout << "[ERROR] Over 23-digits!" << std::endl;
        return EXIT_FAILURE;
      }
      s_nsec = s_tm - 14;
      std::istringstream is(tm_str);
      is >> std::get_time(&t, "%Y%m%d%H%M%S");
      jst.tv_sec  = mktime(&t);
      jst.tv_nsec = 0;
      if (s_tm > 14) {
        jst.tv_nsec = std::stod(
            tm_str.substr(14, s_nsec) + std::string(9 - s_nsec, '0'));
      }
    } else {
      // 現在日時の取得
      ret = std::timespec_get(&jst, TIME_UTC);
      if (ret != 1) {
        std::cout << "[ERROR] Could not get now time!" << std::endl;
        return EXIT_FAILURE;
      }
    }

    // 書き込みファイル open
    std::ofstream ofs(f);
    if (!ofs) return 0;

    // LOOP (日)
    ofs << "{" << std::endl;
    ofs << "  \"counts\": " << ns::kSecD * ns::kDay / ns::kSec
        << "," << std::endl;
    ofs << "  \"data\": [" << std::endl;
    ofs << std::setprecision(12);
    for (i = 0; i < ns::kDay; ++i) {
      jst = ns::ts_add(jst, i * ns::kSecD);

      // EOP データ取得
      utc = ns::jst2utc(jst);
      ns::Eop o_e(utc);
      pm_x = o_e.pm_x;
      pm_y = o_e.pm_y;
      lod  = o_e.lod;
      ut1  = ns::utc2ut1(utc);
      tai  = ns::utc2tai(utc);

      // LOOP (指定秒間隔)
      for (j = 0; j < int(ns::kSecD); j += ns::kSec) {
        jst_wk = ns::ts_add(jst, j);
        utc_wk = ns::ts_add(utc, j);
        ut1_wk = ns::ts_add(ut1, j);
        tai_wk = ns::ts_add(tai, j);
        //std::cout << ns::gen_time_str(jst_wk) << " JST" << std::endl;

        // TLE 読み込み, gravconst 取得
        ns::Tle o_t(ut1_wk);
        tle = o_t.get_tle();

        // ISS 初期位置・速度の取得
        ns::Sgp4 o_s(ut1_wk, tle);
        sat = o_s.twoline2rv();

        // 指定 UT1 の ISS 位置・速度の取得
        teme = o_s.propagate(sat);

        // TEME -> BLH 変換
        ns::Blh o_b(ut1_wk, tai_wk, pm_x, pm_y, lod);
        blh = o_b.teme2blh(teme);

        // 結果出力
        ofs << "    {" << std::endl;
        ofs << "      \"jst\": \"" << ns::gen_time_str(jst_wk) << "\","
            << std::endl;
        ofs << "      \"utc\": \"" << ns::gen_time_str(utc_wk) << "\","
            << std::endl;
        ofs << "      \"latitude\": "  << blh.r.b << "," << std::endl;
        ofs << "      \"longitude\": " << blh.r.l << "," << std::endl;
        ofs << "      \"height\": "    << blh.r.h << "," << std::endl;
        ofs << "      \"velocity\": "  << blh.v   << std::endl;
        if (i == 1 && j == int(ns::kSecD) - ns::kSec) {
          ofs << "    }"  << std::endl;
        } else {
          ofs << "    }," << std::endl;
        }
      }
    }
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;

    // 書き込みファイル close
    ofs.close();
  } catch (...) {
      std::cerr << "EXCEPTION!" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

