# iss_sgp4_json

概要
====

* 指定日時から48時間分の ISS 位置／速度を計算し、 JSON 出力  
  （但し、座標系は BLH(Beta(Latitude), Lambda(Longitude), Height)）

ビルド方法
==========

`make`

（やり直す場合は、 `make clean` をしてから）

準備
====

* TLE（Two-line elements; 2行軌道要素形式）データを [こちら](https://celestrak.com/NORAD/elements/supplemental/iss.txt") からダウンロード。実行プログラムと同じディレクトリ内に配置する。（ファイル名: `tle.txt`）
* EOP（Earth Orientation Parameter; 地球姿勢（回転）パラメータ） データを [こちら](https://gist.github.com/komasaru/7b8fa4da835920d78c2c5f4742acaee9 "GitHub Gist") を参考にダウンロード＆整形。実行プログラムと同じディレクトリ内に配置する。（ファイル名: `eop.txt`）
* うるう秒データを [こちら](https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat file/Leap_Second.dat) からダウンロード。実行プログラムと同じディレクトリ内に配置する。（ファイル名: `Leap_Second.dat`）
* 当 GitHub リポジトリにアップされている `tle.txt`, `eop.txt` は古いものなので、注意。（場合によっては、 `Leap_Second.dat` も古い可能性がある）

実行方法
========

`./iss_sgp4_json [YYYYMMDDHHMMSSMMMMMMMMM]`

* コマンドライン引数には JST（日本標準時） を指定する。
* JST（日本標準時）は「年・月・日・時・分・秒・ナノ秒」を最大23桁で指定する。
* JST（日本標準時）を指定しない場合は、システム日時を JST とみなす。
* JST（日本標準時）を先頭から部分的に指定した場合は、指定していない部分を 0 とみなす。
* 正常に終了すれば、実行プログラムと同じディレクトリ内に `iss.json` が生成される。

