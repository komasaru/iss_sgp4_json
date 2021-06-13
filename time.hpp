#ifndef ISS_SGP4_JSON_TIME_HPP_
#define ISS_SGP4_JSON_TIME_HPP_

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace iss_sgp4_json {

struct DateTime {
  unsigned int year;
  unsigned int month;
  unsigned int day;
  unsigned int hour;
  unsigned int minute;
  double       second ;
};

std::string gen_time_str(struct timespec ts);     // 日時文字列生成
struct timespec ts_add(struct timespec, double);  // 時刻 + 秒数
DateTime days2ymdhms(unsigned int, double);       // 年+経過日数 => 年月日時分秒
double jday(DateTime);                            // 年月日時分秒 => ユリウス日
double gstime(double);                            // Greenwich sidereal time calculation
double get_dut1(struct timespec);                 // DUT1 取得(EOP 読み込み)
unsigned int get_dat(struct timespec);            // DAT (= TAI - UTC)（うるう秒総和）取得
struct timespec jst2utc(struct timespec);         // JST -> UTC
struct timespec utc2ut1(struct timespec);         // UTC -> UT1
struct timespec utc2tai(struct timespec);         // UTC -> TAI
struct timespec tai2tt(struct timespec);          // TAI -> TT
double gc2jd(struct timespec);                    // Gregorian Calendar -> Julian Day
double jd2jcn(double);                            // Julian Day -> Julian Century Number

}  // namespace iss_sgp4_json

#endif

