# pymatgen 使い方メモ
 
 pymatgenは基本的にVASPを使った様々な計算を支援してくれるツールである。結構便利なのが、Band計算のときのパスの生成、
 POTCARの生成である。POTCARを作らせるためにはpymatgenにPOTCARの在り処を知らせないといけない。
 この方法を忘れがちなのでメモしておく  

 設定にはpymatgenのコマンドラインツールpmgをつかう。

VASPのPOTCARを解凍済みのディレクトリを<VASP_PP_DIR>とする。pymatgenが処理できる形式にpotcar をパースしたものの置き場所を<pmg_pp_dir>とする。まず
 ```
 pmg config -p <VASP_PP_DIR> <pmg_pp_dir>
 ```
 でpymatgenが読み込める形に変換する。
 その後、
 ```
  pmg config --add PMG_VASP_PSP_DIR <pmg_pp_dir>
 ```
 で環境変数を設定する。


```
from pymatgen.io.vasp.inputs import Potcar

```