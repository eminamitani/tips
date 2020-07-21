# Mac上で各種シミュレーションをするには

## macで使えるFortran compiler
無料：gfortran
```
brew cask install gfortran
```
で入るが、openmpi+gfortranを想定している場合は
```
brew install open-mpi
```
でgccも含めて一気に入れてしまうのが便利

有料：intel ifort for Mac
高い。あとopenmpiを使えるようにするための設定がめんどくさい。
brewでも--withとかのオプションで入るかもしれないがよくわからんので、ソースコンパイルした

```
./configure --prefix=/usr/local/openmpi CC=gcc CXX=g++ F77=ifort FC=ifort
make all
sudo make install
```