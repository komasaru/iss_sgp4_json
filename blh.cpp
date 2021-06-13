#include "blh.hpp"

namespace iss_sgp4_json {

// 定数
static constexpr unsigned int kJ2k   = 2451545;          // Julian Day of 2000-01-01 12:00:00
static constexpr unsigned int kDayJc = 36525;            // Days per Julian century
static constexpr double       kPi    = atan(1.0) * 4.0;  // 円周率
static constexpr double       kPi2   = kPi * 2.0;        // 円周率 * 2
static constexpr double       kPi180 = kPi / 180.0;      // 円周率 / 180.0
static constexpr double       kSecD  = 86400.0;          // Seconds per day
// [ WGS84 座標パラメータ ]
// a(地球楕円体長半径(赤道面平均半径))
static constexpr double kA   = 6378137.0;
// 1 / f(地球楕円体扁平率=(a - b) / a)
static constexpr double k1F  = 298.257223563;
// b(地球楕円体短半径)
static constexpr double kB   = kA * (1.0 - 1.0 / k1F);
// e^2 = 2 * f - f * f
//     = (a^2 - b^2) / a^2
static constexpr double kE2  = (1.0 / k1F) * (2.0 - (1.0 / k1F));
// e'^2= (a^2 - b^2) / b^2
static constexpr double kEd2 = kE2 * kA * kA / (kB * kB);

/*
 * @brief      コンストラクタ
 *
 * @param[in]  UTC       (timespec)
 * @param[in]  DUT1      (double)
 * @param[in]  極運動(x) (double)
 * @param[in]  極運動(y) (double)
 * @param[in]  LOD       (double)
 */
Blh::Blh(
    struct timespec ut1, struct timespec tai,
    double pm_x, double pm_y, double lod) {
  struct timespec tt;
  double jd_tt;
  this->pm_x = pm_x;
  this->pm_y = pm_y;
  this->lod  = lod;
  tt         = tai2tt(tai);
  jd_ut1     = gc2jd(ut1);
  jd_tt      = gc2jd(tt);
  jcn_tt     = jd2jcn(jd_tt);
}

/*
 * @brief   TEME -> BLH
 *
 * @param   TEME (PvTeme)
 * @return  BLH  (PvBlh)
 */
PvBlh Blh::teme2blh(PvTeme teme) {
  double gmst;
  double om;
  double gmst_g;
  std::vector<std::vector<double>> mtx_z( 3, std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> mtx_pm(3, std::vector<double>(3, 0.0));
  Coord    r_pef  = {0.0, 0.0, 0.0};
  //Coord    v_pef  = {0.0, 0.0, 0.0};
  Coord    r_ecef = {0.0, 0.0, 0.0};
  //Coord    v_wk   = {0.0, 0.0, 0.0};
  Coord    om_e   = {0.0, 0.0, 0.0};
  CoordBlh blh_wk;
  PvBlh    blh;

  try {
    // GMST（グリニッジ平均恒星時）計算
    gmst   = calc_gmst();
    // Ω（月の平均昇交点黄経）計算（IAU1980章動理論）
    om     = calc_om();
    // GMST に運動項を適用（1997年より新しい場合）
    gmst_g = apply_kinematic(gmst, om);
    // GMST 回転行列（z軸を中心とした回転）
    mtx_z  = gen_mtx_rz(gmst_g);
    // 極運動(Polar Motion)回転行列
    mtx_pm = gen_mtx_rpm();
    // PEF 座標の計算（GMST 回転行列の適用）
    r_pef = apply_mtx(mtx_z, teme.r);
    // ECEF 座標（位置）の計算（極運動(Polar Motion)の適用）
    r_ecef = apply_mtx(mtx_pm, r_pef);
    // Ω_earth値の計算
    om_e = calc_om_e();
    // PEF 座標（速度）の計算（GMST 回転行列の適用）
    //v_pef = apply_mtx(mtx_z, teme.v);
    //v_wk  = v_cross(om_e, r_pef);
    //v_pef.x -= v_wk.x;
    //v_pef.y -= v_wk.y;
    //v_pef.z -= v_wk.z;
    // ECEF 座標 => BLH(Beta, Lambda, Height) 変換
    blh_wk = ecef2blh(r_ecef);
    blh.r.b = blh_wk.b;
    blh.r.l = blh_wk.l;
    blh.r.h = blh_wk.h / 1000.0;
    // 速度は BLH 変換しない
    blh.v = sqrt(teme.v.x * teme.v.x
               + teme.v.y * teme.v.y
               + teme.v.z * teme.v.z);
  } catch (...) {
    throw;
  }

  return blh;
}  // teme2blh

/********************************************
 **** 以下、 private function/procedures ****
 ********************************************/

/*
 * @brief   GMST（グリニッジ平均恒星時）計算
 *          * IAU1982理論(by David Vallado)によるもの
 *              GMST = 18h 41m 50.54841s
 *                   + 8640184.812866s T + 0.093104s T^2 - 0.0000062s T^3 
 *              (但し、 T = 2000年1月1日12時(UT1)からのユリウス世紀単位)
 *
 * @param   <none>
 * @return  gmst (double)
 */
double Blh::calc_gmst() {
  double t_ut1;
  double gmst;

  try {
    t_ut1 = (jd_ut1 - kJ2k) / kDayJc;
    gmst = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866
         + (0.093104 -  6.2e-6 * t_ut1) * t_ut1) * t_ut1;
    gmst = fmod(gmst * kPi180 / 240.0, kPi2);
    while (gmst < 0.0 ) { gmst += kPi2; }
    while (gmst > kPi2) { gmst -= kPi2; }
  } catch (...) {
    throw;
  }

  return gmst;
}

/*
 * @brief   Ω（月の平均昇交点黄経）計算（IAU1980章動理論）
 *          * Ω = 125°02′40″.280
 *               - ((5 * 360)° + 134°08′10″.539) * T
 *               + 7″.455 * T2
 *               + 0″.008 * T3
 *
 * @param   <none>
 * @return  Ω (double)
 */
double Blh::calc_om() {
  double om;

  try {
    om = 125.04452222 + ((-6962890.5390 + (7.455 + 0.008
       * jcn_tt) * jcn_tt) * jcn_tt) / 3600.0;
    om = fmod(om, 360.0);
    while (om <   0.0) { om += 360.0; }
    while (om > 360.0) { om -= 360.0; }
    om *= kPi180;
  } catch (...) {
    throw;
  }

  return om;
}

/*
 * @brief   GMST に運動項を適用（1997年より新しい場合）
            * gmst_g = gmst \
                     + 0.00264 * PI / (3600 * 180) * sin(om) \
                     + 0.000063 * PI / (3600 * 180) * sin(2.0 * om)
 *
 * @param[in]  GMST (適用前) (double)
 * @param[in]  Ω (double)
 * @return     GMST (適用後) (double)
 */
double Blh::apply_kinematic(double gmst, double om) {
  double gmst_g;

  try {
    if (jd_ut1 > 2450449.5) {
      gmst_g = gmst
             + 0.00264  * kPi / (3600.0 * 180.0) * sin(om)
             + 0.000063 * kPi / (3600.0 * 180.0) * sin(om * 2.0);
    } else {
      gmst_g = gmst;
    }
    gmst_g = fmod(gmst_g, kPi2);
  } catch (...) {
    throw;
  }

  return gmst_g;
}

/*
 * @brief      z 軸を軸とした座標軸回転行列
 *
 * @param[in]  回転量(rad) (double)
 * @return     回転行列(3x3) (vector<vector<double>>)
 */
std::vector<std::vector<double>> Blh::gen_mtx_rz(double ang) {
  double c;
  double s;
  std::vector<std::vector<double>> mtx(3, std::vector<double>(3, 0.0));

  try {
    c = cos(ang);
    s = sin(ang);
    mtx[0][0] =   c;
    mtx[0][1] =   s;
    mtx[1][0] =  -s;
    mtx[1][1] =   c;
    mtx[2][2] = 1.0;
  } catch (...) {
    throw;
  }

  return mtx;
}

/*
 * @brief   極運動の座標軸回転行列
 *
 * @param   <none>
 * @return  回転行列(3x3) (vector<vector<double>>)
 */
std::vector<std::vector<double>> Blh::gen_mtx_rpm() {
  double pm_x_r;
  double pm_y_r;
  double conv;
  double c_xp;
  double s_xp;
  double c_yp;
  double s_yp;
  double sp;
  double s_sp;
  double c_sp;
  std::vector<std::vector<double>> mtx(3, std::vector<double>(3, 0.0));

  try {
    pm_x_r = pm_x * kPi / (180.0 * 60.0 * 60.0 * 1000.0);
    pm_y_r = pm_y * kPi / (180.0 * 60.0 * 60.0 * 1000.0);
    conv = kPi / (3600.0 * 180.0);
    c_xp = cos(pm_x_r);
    s_xp = sin(pm_x_r);
    c_yp = cos(pm_y_r);
    s_yp = sin(pm_y_r);
    // approximate sp value in rad
    sp = -47.0e-6 * jcn_tt * conv;
    s_sp = sin(sp);
    c_sp = cos(sp);
    mtx[0][0] =  c_xp * c_sp;
    mtx[0][1] =  c_xp * s_sp;
    mtx[0][2] =  s_xp;
    mtx[1][0] = -c_yp * s_sp + s_yp * s_xp * c_sp;
    mtx[1][1] =  c_yp * c_sp + s_yp * s_xp * s_sp;
    mtx[1][2] = -s_yp * c_xp;
    mtx[2][0] = -s_yp * s_sp - c_yp * s_xp * c_sp;
    mtx[2][1] =  s_yp * c_sp - c_yp * s_xp * s_sp;
    mtx[2][2] =  c_yp * c_xp;
  } catch (...) {
    throw;
  }

  return mtx;
}

/*
 * @brief       回転行列適用
 *
 * @param[in]  回転行列(3x3) (vector<vector<double>>)
 * @param[in]  座標（回転前）(Coord)
 * @return     座標（回転後）(Coord)
 */
Coord Blh::apply_mtx(std::vector<std::vector<double>> mtx_r, Coord cd_src) {
  Coord cd;

  try {
    cd.x = mtx_r[0][0] * cd_src.x
         + mtx_r[0][1] * cd_src.y
         + mtx_r[0][2] * cd_src.z;
    cd.y = mtx_r[1][0] * cd_src.x
         + mtx_r[1][1] * cd_src.y
         + mtx_r[1][2] * cd_src.z;
    cd.z = mtx_r[2][0] * cd_src.x
         + mtx_r[2][1] * cd_src.y
         + mtx_r[2][2] * cd_src.z;
  } catch (...) {
    throw;
  }

  return cd;
}

/*
 * @brief   Ω_earch 計算
 *
 * @param   <none>
 * @return  Ω_earth (Coord)
 */
Coord Blh::calc_om_e() {
  Coord om_e = {0.0, 0.0, 0.0};

  try {
    om_e.z = 7.29211514670698e-05 * (1.0 - lod / kSecD);
  } catch (...) {
    throw;
  }

  return om_e;
}

/*
 * @brief   ベクトルの外積計算
 *
 * @param[in]  ベクトル (Coord)
 * @param[in]  ベクトル (Coord)
 * @return     ベクトルの外積 (Coord)
 */
Coord Blh::v_cross(Coord cd_a, Coord cd_b) {
  Coord cd = {0.0, 0.0, 0.0};

  try {
    cd.x = cd_a.y * cd_b.z - cd_a.z * cd_b.y;
    cd.y = cd_a.z * cd_b.x - cd_a.x * cd_b.z;
    cd.z = cd_a.x * cd_b.y - cd_a.y * cd_b.x;
  } catch (...) {
    throw;
  }

  return cd;
}

/*
 * @brief      関数 N (ECEF -> BLH 変換処理用)
 *
 * @param[in]  X (double)
 * @return     計算結果 (double)
 */
double Blh::n(double x) {
  double res;

  try {
    res = kA / sqrt(1.0 - kE2 * pow(sin(x * kPi180), 2));
  } catch (...) {
    throw;
  }

  return res;
}

/*
 * @brief      ECEF -> BLH
 *
 * @param[in]  ECEF 座標 (Coord)
 * @return     BLH  座標 (CoordBlh)
 *             (x -> b, y -> l, z -> h と読み替え)
 */
CoordBlh Blh::ecef2blh(Coord ecef) {
  double   x;
  double   y;
  double   z;
  double   p;
  double   theta;
  CoordBlh blh;

  try {
    x = ecef.x * 1.0e3;
    y = ecef.y * 1.0e3;
    z = ecef.z * 1.0e3;
    p = sqrt(x * x + y * y);
    theta = atan2(z * kA, p * kB) / kPi180;
    blh.b = atan2(
      z + kEd2 * kB * pow(sin(theta * kPi180), 3),
      p - kE2  * kA * pow(cos(theta * kPi180), 3)
    ) / kPi180;                                 // Beta(Latitude)
    blh.l = atan2(y, x) / kPi180;                // Lambda(Longitude)
    blh.h = (p / cos(blh.b * kPi180)) - n(blh.b);  // Height
  } catch (...) {
    throw;
  }

  return blh;
}

}  // namespace iss_sgp4_json

