# Bevel_And_Emboss_M
gometh兄貴の[【AviUtl】ベベルとエンボス スクリプト](https://www.nicovideo.jp/watch/sm39006767)を高速化する試み

元スクリプト作成者のgometh兄貴に感謝申し上げます

## Version
+ r1  
  とりあえず、そのままC++で実装  
  rikky_moduleが導入されていない環境でも動作するようにした  
  どこか実装間違えてるかもしれない...  
  汚いコードだが、疲れてしまったので、誰かがプルリクをくれると信じて公開

+ r2  
  `std::pow(hoge, 2);`  
  を  
  `hoge * hoge;`  
  に変更(こっちのほうが早い)