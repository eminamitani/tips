# ASE メモ

Atomistic Simulation Environment(ASE)は各種原子シミュレーションするためのセットアップなどを支援してくれる神pythonモジュールである。もはや無いと仕事にならない。神機能が色々あるのだけれども、色々と癖があって、なれないと「これどうするの…」ということが多いので、自分が直面した問題をベースにメモを作っておく。

## VASP関係
自分はVASPで「表面に何かつけて計算するやで〜」という仕事が多い。この時に使う機能あれこれについてメモしておく。
VASPのPOTCARなどを自動で作成するには、ASEのCalculatorクラスにあるvasp.Vaspモジュールを使うのだけれども、このモジュールがPOTCARの場所を参照できるように、
```
export VASP_PP_PATH=$HOME/vasp_paw
```
みたいに、VASPのポテンシャルのファイルがある場所を環境変数で指定する必要がある。

### 分子を吸着させる系の構造作成
このときはだいたい  
- 基板の構造
- 分子の構造
  
を別個に作って構造緩和させておいて、それを合体するパターンが多い。分子は中心をユニットセル原点に来るように作ることが多いのだが、それだと回転してくっつけたりするときにとんでもない操作が行われて崩壊することが多い。一回真ん中に寄せておくのが正解。その時の操作がこんな感じ。

```
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms

#read structure from CONTCAR
mol = read('CONTCAR', format='vasp')
#move to center
mol.translate([mol.cell[0][0]/2,mol.cell[1][1]/2,mol.cell[2][2]/2])
mol.wrap(pbc=[1,1,1])
```

基板についても読み込ませる。
```
slab=read('CONTCAR_slab', format='vasp')
```

分子の構造をコピーして、ユニットセルの形をslabのものに合わせる。（このときに変形させないようにscaleのオプションは切っておく）。さらに、指定した場所を中心にｚ軸回りに回転させたり、
吸着サイトに合うように平行移動したりする。最後にextendの機能をつかって合体。centerとvacuumオプションで真空層もつけておく
```
adsorbates = Atoms(mol.symbols, mol.positions)
adsorbates.set_cell(NbSe2.get_cell(), scale_atoms=False)
#set center to rotate the molecule
#for example origin
c = [0.0,0.0,0.0]
adsorbates.rotate(angle, 'z', center=c)
#translate vector
vector=[0.2,0.2,0.0]
adsorbates.translate(vector)
#move upward
dist=2.0
adsorbates.translate([0.0, 0.0, dist])

slab.extend(adsorbates)
slab.center(vacuum=vacuum, axis=2)
```

場合によっては、Slabの何層かを固定したいときがある。これにはFixAtomsを使う。自分はよくある条件より下の原子にはFixの条件をつけるみたいなので実装している。
```
bottom=min(slab.positions[:,2])
fix=FixAtoms(indices=[atom.index for atom in slab if abs(atom.position[2]-bottom) < 0.5])
slab.set_constraint(fix)
```

作ったものをVASP POSCAR形式として書き出す
```
write("POSCAR.vasp", slab, sort=True, format='vasp', vasp5=True, direct=True)
```

potcar もついでに作りたい。上でPOSCAR作る時にソートかけてるので、順番がもとのデータと違っている可能性があるから、作ったファイル読みこんで作ったほうがいい。どの元素でどのポテンシャル使えばいいかはVASP公式おすすめセットがすでにVaspモジュール内に準備されていて`setups='recommended'`のオプションを使えばそれを参照できる。
```
from ase.calculators.vasp import Vasp
calc = Vasp(xc='pbe',setups='recommended')
atoms=read('POSCAR.vasp',format='vasp')
calc.initialize(atoms)
calc.write_potcar()
```


