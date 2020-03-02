# ImpDet_macllo
  - 改良型検出器のデータ解析マクロ

## sumop.cpp  
  ファイルの先頭に  
  BG⇒バックグラウンド  
  OB⇒オブジェクト  
  とついているものをそれぞれ足し合わせ以下のファイルに出力する。  
  **1_______________Background_output.txt**  
  **2__Object_output._______txt**  
  一行目に「100 100 <計測時間>」を書いておく必要がある  

## atten.cpp  
  1_Background_output.txt  
  2_Object_output.txt  
  を読み込んできてそれぞれの強度分布と減衰率を計算する

## ana00.cpp  
  00ベクトルのデータのみを解析するプログラム  
  計測時間はプログラムの中を変更する必要がある。
  だいぶ初期に作ったプログラムなので改善の余地多々あり。
  BG_00_data.txt, OB__00_data.txtをoutput.txtに変換する部分と  
  BG__00_output.txt, OB__00_output.txtからヒストグラムを作る部分に分かれている  
