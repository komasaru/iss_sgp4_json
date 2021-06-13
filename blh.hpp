#ifndef ISS_SGP4_JSON_BLH_HPP_
#define ISS_SGP4_JSON_BLH_HPP_

#include "sgp4.hpp"
#include "time.hpp"

#include <iomanip>
#include <string>
#include <vector>

namespace iss_sgp4_json{

// 座標構造体
struct CoordBlh {
  double b;  // Beta(Latitude)
  double l;  // Lambda(Longitude)
  double h;  // Height
};
// 位置・速度構造体(BLH)
struct PvBlh {
  CoordBlh r;  // 位置
  double   v;  // 速度
};

class Blh{
  double pm_x;    // 極運動(x)
  double pm_y;    // 極運動(y)
  double lod;     // LOD
  double jd_ut1;  // JD(UT1)
  double jcn_tt;  // JCN(TT)

public:
  Blh(struct timespec, struct timespec, double, double, double);  // コンストラクタ
  PvBlh teme2blh(PvTeme);                  // TEME -> BLH

private:
  double calc_gmst();                      // GMST (グリニッジ平均恒星時) 計算
  double calc_om();                        // Ω (月の平均昇交点黄経; IAU1980章動理論) 計算
  double apply_kinematic(double, double);  // GMST に運動項を適用(1997年より新しい場合)
  std::vector<std::vector<double>> gen_mtx_rz(double);   // z軸を中心とした座標軸回転行列生成
  std::vector<std::vector<double>> gen_mtx_rpm();
                                           // 極運動の座標軸回転行列生成
  Coord apply_mtx(std::vector<std::vector<double>>, Coord);  // 回転行列適用
  Coord calc_om_e();                       // Ω_earch 計算
  Coord v_cross(Coord, Coord);             // ベクトルの外積計算
  double n(double);                        // 関数 N (ECEF -> BLH 変換処理用)
  CoordBlh ecef2blh(Coord);                // ECEF -> BLH
};

}  // namespace iss_sgp4_json

#endif

