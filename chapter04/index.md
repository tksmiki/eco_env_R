# 第4章

## 準備1: フォルダの準備
この本の内容に沿ってRを勉強していくには、そろそろちゃんとデータファイルやサンプルコード・自分で書いたコードをフォルダに整理していく必要があります。以下のようにフォルダを作っていきましょう。

### [1.1] フォルダ名の付け方のルールを復習しよう
「3.2.3 ファイル名の付けかた」で説明したファイル名の付けかたのルールと同様に、フォルダ名を付けるときには以下のルールを厳守しましょう。<br>
　　・ルール1：半角英数字及び半角記号の一部(=, -, _)だけ使う<br>
　　・ルール2：半角でも全角でもスペースをフォルダ名の一部に混ぜない<br>

### [1.2] この本全体用の親フォルダをつくる
<b>practice_R</b>という名前のフォルダを作るとよいでしょう。

### [1.3] 各章ごとにフォルダを作る
まずは、第4章・第5章両方に使うフォルダとして, <b>chapter04_05</b>フォルダを親フォルダ(<b>practice_R</b>)の下に作りましょう。

### [1.4] chapter04_05フォルダの下の設定
<b>chapter04_05</b>フォルダの下に、さらにいくつかフォルダを作りましょう。<br>
　　・<b>chpter04_05</b>フォルダの直下に、二つのフォルダ(<b>test_phytoplankton</b>と<b>test_ecoplate</b>)を作る。<br>
　　・<b>test_ecoplate</b>フォルダの直下に、もう一つのフォルダ <b>text_file</b>を作る。<br>

## 準備2: データファイル・メタデータファイルのダウンロード
第4章・第5章では一番たくさんファイルを使います。必要なファイルはすべてここからダウンロードできます。必ず以下の手順に従って指定したフォルダにファイルを保存しましょう。別のところに間違えて保存してしまうと、Rからそのファイルを見つけられなくなってしまいます。

### [2.1] 植物プランクトンの種組成に関するメタデータ・観測データのダウンロード
新しいタブで以下のリンクを開いてファイルを<b>test_phytoplankon</b>フォルダの中に保存しましょう。
<a href="./test_phytoplankton/" target="_blank" rel="noopener noreferrer"><b>test_phytoplankton</b></a>
### [2.2] 細菌群集の炭素代謝に関するメタデータ・観測データのダウンロード
新しいタブで以下のリンクを開いてファイルを<b>test_ecoplate</b>フォルダの中に保存しましょう。
<a href="./test_ecoplate/" target="_blank" rel="noopener noreferrer"><b>test_ecoplate</b></a>

## 4.2 エクセルを使ってメタデータをまとめる
以下のサンプルスクリプトを<b>test_phytoplankton</b>フォルダの中に保存しましょう。<b>chapter04_05</b>の直下に保存してはいけません！
[phyto_test.R](./test_phytoplankton/phyto_test.R)

## 4.3 プログラミングの基礎(if, loop, 関数)
以下のサンプルスクリプトを<b>chapter04_05</b>フォルダの直下に保存しましょう。
[chapter04_4-3.R](./chapter04_4-3.R)
このRスクリプトの日本語解説込みのHTMLファイルは、次のリンクをクリックすればそのままウェブブラウザの新規タブで表示されます（Edgeを推奨、Safariも可）：<br>
<a href="./chapter04_4-3.nb.html" target="_blank" rel="noopener noreferrer">chapter04_4-3.nb.html</a><br>

## 4.4 メタデータを活用した大量のデータの自動読み込み
### 4.4.1 自動読み込み例1: 植物プランクトンの種構成データ
4.2のところでダウンロードしたサンプルスクリプトphyto_test.Rを引き続き使います。<br>
このRスクリプトの日本語解説込みのHTMLファイルは、次のリンクをクリックすればそのままウェブブラウザの新規タブで表示されます（Edgeを推奨、Safariも可）：<br>
<a href="./test_phytoplankton/phyto_test.nb.html" target="_blank" rel="noopener noreferrer">phyto_test.nb.html</a><br>

### 4.4.2 自動読み込み例2: 細菌群集の炭素代謝データ（エコプレート）
以下のサンプルスクリプトを<b>test_ecoplate</b>フォルダの中に保存しましょう。<b>chapter04_05</b>の直下や<b>test_phytoplankton</b>の下に保存してはいけません！
[ecoplate_test.R](./test_ecoplate/ecoplate_test.R)
このRスクリプトの日本語解説込みのHTMLファイルは、次のリンクをクリックすればそのままウェブブラウザの新規タブで表示されます（Edgeを推奨、Safariも可）：<br>
<a href="./test_ecoplate/ecoplate_test.nb.html" target="_blank" rel="noopener noreferrer">ecoplate_test.nb.html</a><br>


