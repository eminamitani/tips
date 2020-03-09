# gnuplot関係のメモ  

## Macでの設定
Macのgnuplotは
```
brew install gnuplot
```
でもインストール出来るのだが、この場合、フォントの文字化けやepsへの出力でのエラーなどに悩まされる。
なので、aquatermを入れてそれを介して出力するのが一つの手。aquatermはベクトルグラフィックもレンダーできる描画ソフトなので、eps周りもこれで解決する。

gnuplot+aquatermの場合はX11も必要になる。XQuartzを
https://www.xquartz.org/
から入手してインストールしておくと良い。

aquatermのインストール
```
brew cask install aquaterm
```

その後の作業は、brewのgnuplotインストール用の設定ファイル`gnuplot.rb`を書き換えてオプションを付けてインストールする。
手順はこれそのまま
https://hatarakutaitsu.hatenablog.com/entry/2019/09/11/141842