#include "time.hpp"

namespace iss_sgp4_json {

static constexpr char         kFEop[]    = "eop.txt";
static constexpr char         kFDat[]    = "Leap_Second.dat";
static constexpr double       kE9        = 1.0e9;
static constexpr double       kEm6       = 1.0e-6;
static constexpr double       kPi        = atan(1.0) * 4.0;  // 円周率
static constexpr double       kPi2       = kPi * 2.0;        // 円周率 * 2
static constexpr double       kDeg2Rad   = kPi / 180.0;      // 0.0174532925199433
static constexpr double       kJstOffset = 32400.0;          //  9 * 60 * 60
static constexpr double       kTtTai     = 32.184;           // TT - TAI
static constexpr unsigned int kJ2k       = 2451545;          // Julian Day of 2000-01-01 12:00:00
static constexpr unsigned int kDayJc     = 36525;            // Days per Julian century

/*
 * @brief      日時文字列生成
 *
 * @param[in]  日時 (timespec)
 * @return     日時文字列 (string)
 */
std::string gen_time_str(struct timespec ts) {
  struct tm t;
  std::stringstream ss;
  std::string str_tm;

  try {
    localtime_r(&ts.tv_sec, &t);
    ss << std::setfill('0')
       << std::setw(4) << t.tm_year + 1900 << "-"
       << std::setw(2) << t.tm_mon + 1     << "-"
       << std::setw(2) << t.tm_mday        << " "
       << std::setw(2) << t.tm_hour        << ":"
       << std::setw(2) << t.tm_min         << ":"
       << std::setw(2) << t.tm_sec         << "."
       << std::setw(3) << round(ts.tv_nsec * kEm6);
    return ss.str();
  } catch (...) {
    throw;
  }
}

/*
 * @brief      時刻 + 秒数
 *
 * @param[in]  時刻 (timespec)
 * @param[in]  秒数 (double)
 * @return     時刻 (timespec)
 */
struct timespec ts_add(struct timespec ts_src, double s) {
  struct timespec ts;

  try {
    ts.tv_sec  = ts_src.tv_sec + int(s);
    ts.tv_nsec = ts_src.tv_nsec + (s - int(s)) * kE9;
    while (ts.tv_nsec > kE9) {
      ++ts.tv_sec;
      ts.tv_nsec -= kE9;
    }
    while (ts.tv_nsec < 0) {
      --ts.tv_sec;
      ts.tv_nsec += kE9;
    }
  } catch (...) {
    throw;
  }

  return ts;
}

/*
 * @brief       年+経過日数 => 年月日時分秒
 *              * this procedure converts the day of the year, days, to the 
 *                equivalent month, day, hour, minute and second.
 *              * algorithm: set up array for the number of days per month 
 *                           find leap year - use 1900 because 2000 is a leap 
 *                           year loop through a temp value while the value 
 *                           is < the days perform int conversions to the 
 *                           correct day and month convert remainder into 
 *                           h m s using type conversions
 *
 * @param[in]  年 (unsigned int)
 * @param[in]  日数 (double)
 * @return     日時 (DateTime)
 */
DateTime days2ymdhms(unsigned int year, double days) {
  unsigned int lmonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  unsigned int dayofyr;
  unsigned int i;
  unsigned int inttemp;
  double       temp;
  DateTime     dt;

  try {
    dayofyr = int(days);

    // ----------------- find month and day of month ----------------
    if (year % 4 == 0) { lmonth[1] = 29; }
    i = 0;
    inttemp = 0;
    while (dayofyr > inttemp + lmonth[i] && i < 12) {
      inttemp += lmonth[i];
      ++i;
    }
    dt.year  = year;
    dt.month = i + 1;
    dt.day   = dayofyr - inttemp;

    // ----------------- find hours minutes and seconds -------------
    temp   = (days - dayofyr) * 24.0;
    dt.hour   = int(temp);
    temp   = (temp - dt.hour) * 60.0;
    dt.minute = int(temp);
    dt.second = (temp - dt.minute) * 60.0;
  } catch (...) {
    throw;
  }

  return dt;
}

/*
 * @brief      年月日時分秒 => ユリウス日
 *             * this procedure finds the julian date given the year, month, day, and time.
 *               the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
 *             * algorithm: calculate the answer in one step for efficiency
 *
 * @param[in]  日時 (DateTime)
 * @return     ユリウス日 (double)
 */
double jday(DateTime dt) {
  double jd;

  try {
    jd = (367.0 * dt.year
       - int(7.0 * (dt.year + int((dt.month + 9.0) / 12.0)) * 0.25)
       + int(275.0 * dt.month / 9.0) + dt.day + 1721013.5
       + ((dt.second / 60.0 + dt.minute) / 60.0 + dt.hour) / 24.0);
  } catch (...) {
    throw;
  }

  return jd;
}

/*
 * @brief      Greenwich sidereal time calculation
 *
 * @param[in]  ユリウス日(UT1) (double)
 * @return     GST (double
 */
double gstime(double jdut1) {
  double tut1;
  double gst;

  try {
    tut1 = (jdut1 - 2451545.0) / 36525.0;
    gst = 67310.54841
        + ((876600.0 * 3600.0 + 8640184.812866)
        + (0.093104
        - 6.2e-6
        * tut1) * tut1) * tut1;
    gst = fmod(gst * kDeg2Rad / 240.0, kPi2);
    if (gst < 0.0) { gst += kPi2; }
  } catch (...) {
    throw;
  }

  return gst;
}

/*
 * @brief      DUT1 取得 (EOP 読み込み)
 *
 * @param<in>  UTC (timespec)
 * @return     DUT1 (double)
 */
double get_dut1(struct timespec ts) {
  std::string f(kFEop);    // ファイル名
  std::string str_utc;     // 対象の UTC 文字列
  std::string eop;         // 1行分バッファ
  double      dut1 = 0.0;  // DUT1

  try {
    // 対象の UTC 年月日
    str_utc = gen_time_str(ts).substr(0, 10);

    // ファイル OPEN
    std::ifstream ifs(f);
    if (!ifs) throw;  // 読み込み失敗

    // ファイル READ
    while (getline(ifs, eop)) {
      if (eop.substr(0, 10) == str_utc) { break; }
    }

    // DUT1 取得
    dut1 = stod(eop.substr(62, 10));
  } catch (...) {
    throw;
  }

  return dut1;
}

/*
 * @brief      DAT (= TAI - UTC)（うるう秒の総和）取得
 *
 * @param<in>  UTC (timespec)
 * @return     DAT (int)
 */
unsigned int get_dat(struct timespec ts) {
  std::string  f(kFDat);   // ファイル名
  std::string  str_utc;    // 対象の UTC 文字列
  std::string  str_wk;     // UTC 文字列（作業用）
  std::string  buf;        // 1行分バッファ
  std::string  ls;         // 該当行
  unsigned int dat = 0.0;  // DUT1

  try {
    // 対象の UTC 年月日
    str_utc = gen_time_str(ts);

    // ファイル OPEN
    std::ifstream ifs(f);
    if (!ifs) throw;  // 読み込み失敗

    // ファイル READ
    while (getline(ifs, buf)) {
      if (buf.substr(0, 1) == "#") { continue; }
      std::stringstream ss;
      ss << std::setfill('0')
         << std::setw(4) << stoi(buf.substr(20, 4)) << "-"
         << std::setw(2) << stoi(buf.substr(17, 2)) << "-"
         << std::setw(2) << stoi(buf.substr(14, 2));
      str_wk = ss.str();
      if (str_utc < str_wk) { break; }
      ls = buf;
    }

    // DAT 取得
    dat = stoi(ls.substr(31, 2));
  } catch (...) {
    throw;
  }

  return dat;
}

/*
 * @brief      JST -> UTC
 *
 * @param[in]  JST (timespec)
 * @return     UTC (timespec)
 */
struct timespec jst2utc(struct timespec jst) {
  struct timespec utc;

  try {
    utc = ts_add(jst, -kJstOffset);
  } catch (...) {
    throw;
  }

  return utc;
}

/*
 * @brief      UTC -> UT1
 *
 * @param[in]  UTC (timespec)
 * @return     UT1 (timespec)
 */
struct timespec utc2ut1(struct timespec utc) {
  double dut1;
  struct timespec ut1;

  try {
    dut1 = get_dut1(utc);
    ut1 = ts_add(utc, dut1);
  } catch (...) {
    throw;
  }

  return ut1;
}

/*
 * @brief      UTC -> TAI
 *
 * @param[in]  UTC (timespec)
 * @return     TAI (timespec)
 */
struct timespec utc2tai(struct timespec utc) {
  int dat;
  struct timespec tai;

  try {
    dat = get_dat(utc);
    tai = ts_add(utc, dat);
  } catch (...) {
    throw;
  }

  return tai;
}

/*
 * @brief      TAI -> TT
 *
 * @param[in]  TAI (timespec)
 * @return     TT  (timespec)
 */
struct timespec tai2tt(struct timespec tai) {
  struct timespec tt;

  try {
    tt = ts_add(tai, kTtTai);
  } catch (...) {
    throw;
  }

  return tt;
}

/*
 * @brief      Gregrorian Calendar -> Julian Day
 *
 * @param[in]  GC (timespec)
 * @return     JD (double)
 */
double gc2jd(struct timespec ts) {
  struct tm t;
  unsigned int year;
  unsigned int month;
  unsigned int day;
  unsigned int hour;
  unsigned int min;
  unsigned int sec;
  double jd;

  try {
    localtime_r(&ts.tv_sec, &t);
    year  = t.tm_year + 1900;
    month = t.tm_mon + 1;
    day   = t.tm_mday;
    hour  = t.tm_hour;
    min   = t.tm_min;
    sec   = t.tm_sec;
    // 1月,2月は前年の13月,14月とする
    if (month < 3) {
      --year;
      month += 12;
    }
    // 日付(整数)部分
    jd = static_cast<int>(365.25 * year)
       + static_cast<int>(year / 400.0)
       - static_cast<int>(year / 100.0)
       + static_cast<int>(30.59 * (month - 2))
       + day
       + 1721088.5;
    // 時間(小数)部分
    jd += (sec / 3600.0 + min / 60.0 + hour) / 24.0;
    // 時間(ナノ秒)部分
    jd += ts.tv_nsec / 1000000000.0 / 3600.0 / 24.0;
  } catch (...) {
    throw;
  }

  return jd;
}

/*
 * @brief      Julian Day -> Julian Century Number
 *
 * @param[in]  JD  (double)
 * @return     JCN (double)
 */
double jd2jcn(double jd) {
  double jcn;

  try {
    jcn = (jd - kJ2k) / kDayJc;
  } catch (...) {
    throw;
  }

  return jcn;
}

}  // namespace iss_sgp4_json

