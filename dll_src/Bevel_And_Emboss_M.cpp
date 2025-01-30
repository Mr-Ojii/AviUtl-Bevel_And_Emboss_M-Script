#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <tuple>
#include <memory>
#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <windows.h>
#include "lua.hpp"

struct Pixel_BGRA {
    uint8_t b, g, r, a;
};

struct Pixel_Info {
    int alpha;
    double dis;
    double x;
    double y;
    double line;
    double gray;
};

struct Point {
    double x, y;
};

struct ObjField {
    double ox, oy, oz;
    double rx, ry, rz;
    double cx, cy, cz;
    double zoom;
    double alpha;
    double aspect;
};

inline ObjField getObjField(lua_State* L) {
    ObjField obj;

    lua_getglobal(L, "obj");
    lua_getfield(L, -1, "ox");
    obj.ox = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "oy");
    obj.oy = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "oz");
    obj.oz = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "rx");
    obj.rx = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "ry");
    obj.ry = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "rz");
    obj.rz = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "cx");
    obj.cx = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "cy");
    obj.cy = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "cz");
    obj.cz = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "zoom");
    obj.zoom = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "alpha");
    obj.alpha = lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, -1, "aspect");
    obj.aspect = lua_tonumber(L, -1);
    lua_pop(L, 1);

    return obj;
}

inline void setObjField(lua_State* L, ObjField obj) {
    lua_getglobal(L, "obj");
    lua_pushnumber(L, obj.ox);
    lua_setfield(L, -2, "ox");
    lua_pushnumber(L, obj.oy);
    lua_setfield(L, -2, "oy");
    lua_pushnumber(L, obj.oz);
    lua_setfield(L, -2, "oz");
    lua_pushnumber(L, obj.rx);
    lua_setfield(L, -2, "rx");
    lua_pushnumber(L, obj.ry);
    lua_setfield(L, -2, "ry");
    lua_pushnumber(L, obj.rz);
    lua_setfield(L, -2, "rz");
    lua_pushnumber(L, obj.cx);
    lua_setfield(L, -2, "cx");
    lua_pushnumber(L, obj.cy);
    lua_setfield(L, -2, "cy");
    lua_pushnumber(L, obj.cz);
    lua_setfield(L, -2, "cz");
    lua_pushnumber(L, obj.zoom);
    lua_setfield(L, -2, "zoom");
    lua_pushnumber(L, obj.alpha);
    lua_setfield(L, -2, "alpha");
    lua_pushnumber(L, obj.aspect);
    lua_setfield(L, -2, "aspect");
}

inline double rad(double deg) {
    return deg * 3.1415926535897932384626433832795028841971 / 180.0;
}

inline double squared(double num) {
    //std::pow(x, 2)は遅い
    return num * num;
}

inline void utl_blur(lua_State *L, double n_blur) {
    lua_getfield(L, -1, "effect");
    lua_pushstring(L, "ぼかし");
    lua_pushstring(L, "範囲");
    lua_pushnumber(L, n_blur);
    lua_pushstring(L, "サイズ固定");
    lua_pushnumber(L, 1);
    lua_call(L, 5, 0);
}

inline void utl_putpixeldata(lua_State *L, void* pix) {
    lua_getfield(L, -1, "putpixeldata");
    lua_pushlightuserdata(L, pix);
    lua_call(L, 1, 0);
}

inline void* utl_getpixeldata(lua_State *L) {
    lua_getfield(L, -1, "getpixeldata");
    lua_call(L, 0, 3);
    void* ptr = lua_touserdata(L, -3);
    lua_pop(L, 3);
    return ptr;
}

inline void utl_blend(lua_State *L, int ble) {
    lua_getfield(L, -1, "setoption");
    lua_pushstring(L, "blend");
    lua_pushnumber(L, ble);
    lua_call(L, 2, 0);
}

inline void utl_copybuffer(lua_State *L, const char* dst, const char* src) {
    lua_getfield(L, -1, "copybuffer");
    lua_pushstring(L, dst);
    lua_pushstring(L, src);
    lua_call(L, 2, 0);
}

inline void utl_draw(lua_State *L, double alp) {
    lua_getfield(L, -1, "draw");
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 1);
    lua_pushnumber(L, alp / 100.0);
    lua_call(L, 5, 0);
}

inline void utl_drawtarget(lua_State *L, const char* dst) {
    lua_getfield(L, -1, "setoption");
    lua_pushstring(L, "drawtarget");
    lua_pushstring(L, dst);
    lua_call(L, 2, 0);
}

inline void utl_temptarget(lua_State *L, int w, int h) {
    lua_getfield(L, -1, "setoption");
    lua_pushstring(L, "drawtarget");
    lua_pushstring(L, "tempbuffer");
    lua_pushnumber(L, w);
    lua_pushnumber(L, h);
    lua_call(L, 4, 0);
}

//点[x,y]から最も近い線分[x0,y0]:[x0+dx,y0+dy]上の点の、距離の2乗と座標を返す関数
inline std::tuple<double, double, double> distance_line(double x, double y, double x0, double y0, double dx, double dy)
{
    //t<=0については、次のループで上書きされるので計算の必要無し
    double t = std::min(1.0, (dx * (x - x0) + dy * (y - y0)) / (dx * dx + dy * dy));
    //そもそもt<=0となるピクセルはifで弾いてる
    return std::make_tuple(std::abs(squared(x0 + dx * t - x) + squared(y0 + dy * t - y)), x0 + dx * t, y0 + dy * t);
}


int bevel_and_emboss(lua_State *L) {
    int a_th = lua_tointeger(L, 1);
    double bev_w = lua_tonumber(L, 2);
    double rot = lua_tonumber(L, 3);
    double high = lua_tonumber(L, 4);
    unsigned int style = lua_tointeger(L, 5);
    double blur = lua_tonumber(L, 6);
    int col1 = lua_tointeger(L, 7);
    int ble1 = lua_tointeger(L, 8);
    double alp1 = lua_tonumber(L, 9);
    int col2 = lua_tointeger(L, 10);
    int ble2 = lua_tointeger(L, 11);
    double alp2 = lua_tonumber(L, 12);
    int bgcol = lua_tointeger(L, 13);
    double preblur = lua_tonumber(L, 14);
    double isline = lua_tonumber(L, 15);

    // 事前計算できるものは計算しておく
    high = std::cos(rad(high));
    rot = rad(rot);
    double sin_rot = std::sin(rot);
    double cos_rot = std::cos(rot);
    uint8_t bgcol_r = (bgcol >> 16) & 0xff;
    uint8_t bgcol_g = (bgcol >> 8) & 0xff;
    uint8_t bgcol_b = (bgcol) & 0xff;
    uint8_t col1_r = (col1 >> 16) & 0xff;
    uint8_t col1_g = (col1 >> 8) & 0xff;
    uint8_t col1_b = (col1) & 0xff;
    uint8_t col2_r = (col2 >> 16) & 0xff;
    uint8_t col2_g = (col2 >> 8) & 0xff;
    uint8_t col2_b = (col2) & 0xff;

    if(style % 4 == 2 || style % 4 == 3)
        bev_w /= 2.0;

    if(bev_w <= 0 || style > 15)
        return 0;

    //4近傍を取得する時にOutOfIndexにしないためのマージン&ぼかし時に輪郭が端で途切れるのを回避
    lua_getglobal(L, "obj");

    double initblur = preblur + 1;
    if (style % 4 != 1)
        initblur += std::ceil(bev_w);

    lua_getfield(L, -1, "effect");
    lua_pushstring(L, "領域拡張");
    lua_pushstring(L, "上");
    lua_pushnumber(L, initblur);
    lua_pushstring(L, "下");
    lua_pushnumber(L, initblur);
    lua_pushstring(L, "左");
    lua_pushnumber(L, initblur);
    lua_pushstring(L, "右");
    lua_pushnumber(L, initblur);
    lua_call(L, 9, 0);

    //"getpixeldata"でw,hとれるじゃん！
    //と思ったが、exedit.auf側で例外が発生する例があるので、getpixelを呼ぶ
    lua_getfield(L, -1, "getpixel");
    lua_call(L, 0, 2);
    int w = lua_tointeger(L, -2);
    int h = lua_tointeger(L, -1);
    lua_pop(L, 2);

    //キャッシュテキストスクリプト用の例外処理（サイズ0の画像でも処理が走ってエラーを起こすため）
    if(w * h == 0)
        return 0;

    std::unique_ptr<Pixel_BGRA[]> bevel_buffer = std::make_unique_for_overwrite<Pixel_BGRA[]>(w * h);
    memcpy(bevel_buffer.get(), utl_getpixeldata(L), sizeof(Pixel_BGRA) * w * h);

    //辺は滑らかになるが頂点が丸くなるので一長一短
    utl_blur(L, preblur);

    //ピクセルのアルファ値,輪郭までの最短距離,最短輪郭の座標X,Y,輪郭の内外判定値,グレスケ(-1～1)
    std::unique_ptr<Pixel_Info[]> pix_info = std::make_unique_for_overwrite<Pixel_Info[]>(w * h);

    Pixel_BGRA* pixel = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
    
    std::unique_ptr<bool[]> border_ans = std::make_unique_for_overwrite<bool[]>(w * h);
    
    for(int j = 0; j < h; j++) {
        for(int i = 0; i < w; i++) {
            pix_info[j * w + i] = {pixel[j * w + i].a, HUGE_VAL, static_cast<double>(i), static_cast<double>(j), 1, 0};
            
            //rikky_module.borderingも一緒に
            //端か、端ではなかった場合上下左右がth以下だった場合
            if(a_th < pixel[j * w + i].a && ((i == 0 || j == 0 || i == w-1 || j == h-1) || (pixel[j * w + i - 1].a <= a_th || pixel[j * w + i + 1].a <= a_th || pixel[j * w + i - w].a <= a_th || pixel[j * w + i + w].a <= a_th)))
                border_ans[j * w + i] = true;
            else
                border_ans[j * w + i] = false;
        }
    }

    //rikky_module.bordering
    std::vector< std::vector<Point> > p;
    int border_i = 0;
    while(1)
    {
        std::vector<Point> pb;
        Point tmp_point = {-1, -1};
        for(; border_i < w * h; border_i++) {
            if(border_ans[border_i]) {
                tmp_point.x = border_i % w;
                tmp_point.y = border_i / w;
                break;
            }
        }
        if(tmp_point.x == -1 || tmp_point.y == -1)
            break;

        //外側判定
        bool ccw = (tmp_point.y == 0 || (pixel[static_cast<int>((tmp_point.y - 1) * w + tmp_point.x)].a <= a_th));
        
        while(true) {
            pb.push_back(tmp_point);
            border_ans[tmp_point.y * w + tmp_point.x] = false;
            Point next_point = {-1, -1};

            if(ccw) {
                //反時計回り 左から反時計回りで探索
                if(tmp_point.x != 0 && border_ans[tmp_point.y * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y};
                else if(tmp_point.y != h - 1 && tmp_point.x != 0 &&  border_ans[(tmp_point.y + 1) * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y + 1};
                else if(tmp_point.y != h - 1 && border_ans[(tmp_point.y + 1) * w + tmp_point.x])
                    next_point = {tmp_point.x, tmp_point.y + 1};
                else if(tmp_point.y != h - 1 && tmp_point.x != w - 1 && border_ans[(tmp_point.y + 1) * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y + 1};
                else if(tmp_point.x != w - 1 && border_ans[tmp_point.y * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y};
                else if(tmp_point.y != 0 && tmp_point.x != w - 1 && border_ans[(tmp_point.y - 1) * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y - 1};
                else if(tmp_point.y != 0 && border_ans[(tmp_point.y - 1) * w + tmp_point.x] == 1)
                    next_point = {tmp_point.x, tmp_point.y - 1};
                else if(tmp_point.y != 0 && tmp_point.x != 0 &&  border_ans[(tmp_point.y - 1) * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y - 1};
            } else {
                //時計回り 右から時計回りで探索
                if(tmp_point.x != w - 1 && border_ans[tmp_point.y * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y};
                else if(tmp_point.y != h - 1 && tmp_point.x != w - 1 && border_ans[(tmp_point.y + 1) * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y + 1};
                else if(tmp_point.y != h - 1 && border_ans[(tmp_point.y + 1) * w + tmp_point.x])
                    next_point = {tmp_point.x, tmp_point.y + 1};
                else if(tmp_point.y != h - 1 && tmp_point.x != 0 &&  border_ans[(tmp_point.y + 1) * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y + 1};
                else if(tmp_point.x != 0 && border_ans[tmp_point.y * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y};
                else if(tmp_point.y != 0 && tmp_point.x != 0 &&  border_ans[(tmp_point.y - 1) * w + tmp_point.x - 1])
                    next_point = {tmp_point.x - 1, tmp_point.y - 1};
                else if(tmp_point.y != 0 && border_ans[(tmp_point.y - 1) * w + tmp_point.x])
                    next_point = {tmp_point.x, tmp_point.y - 1};
                else if(tmp_point.y != 0 && tmp_point.x != w - 1 && border_ans[(tmp_point.y - 1) * w + tmp_point.x + 1])
                    next_point = {tmp_point.x + 1, tmp_point.y - 1};
            }

            if(next_point.x == -1 || next_point.y == -1)
                break;
            
            tmp_point = next_point;
        }
        p.push_back(std::move(pb));
    }

    utl_putpixeldata(L, bevel_buffer.get());

    //近傍の透明度を参照して、輪郭を小数点以下で補正
    for(int i = 0; i < p.size(); i++) {
        auto& pi = p[i];
        for(int j = pi.size() - 1; j >= 0; j--) {
            const Point& point = pi[j];
            int index = point.y * w + point.x;
            double a5 = pix_info[index].alpha, a6 = pix_info[index + 1].alpha, a4 = pix_info[index - 1].alpha, a2 = pix_info[index + w].alpha, a8 = pix_info[index - w].alpha;
            //何故か全周囲が閾値以上のピクセルが時々あるので削除
            if(std::min({a5, a6, a4, a2, a8}) >= a_th) {
                // 削除用としてマーク
                pi[j].x = HUGE_VAL;
            } else {
                double dx = std::clamp(a6 < a_th ? (a5 - a_th) / (a5 - a6) : a4 < a_th ? (a5 - a_th) / (a4 - a5) : 0.0, -1.0, 1.0);
                double dy = std::clamp(a2 < a_th ? (a5 - a_th) / (a5 - a2) : a8 < a_th ? (a5 - a_th) / (a8 - a5) : 0.0, -1.0, 1.0);
                pi[j].x += dx / (1 + std::abs(dy / dx)); 
                pi[j].y += dy / (1 + std::abs(dx / dy));
            }
        }
        // マークされたのを削除
        std::erase_if(pi, [](const auto& pij){ return pij.x == HUGE_VAL; });
    }

    //輪郭平滑化処理
    //謎の定数（根拠無し)
    double hoge = 0.8;
    for(int i = 0; i < p.size(); i++) {
        double bufx = p[i][0].x, bufy = p[i][0].y;
        for(int j = 1; j < p[i].size(); j++) {
            double x0 = p[i][j].x, y0 = p[i][j].y, x1 = p[i][(j + 1) % p[i].size()].x, y1 = p[i][(j + 1) % p[i].size()].y;
            if(squared(x0 * 2 - bufx - x1) + squared(y0 * 2 - bufy - y1) < hoge) {
                p[i][j].x = x0 / 2 + bufx / 4 + x1 / 4;
                p[i][j].y = y0 / 2 + bufy / 4 + y1 / 4;
            }
            bufx = x0;
            bufy = y0;
        }
    }
    
    for(int i = 0; i < p.size(); i++) {
        auto& pi = p[i];
        size_t s[2] = {pi.size() - 2, pi.size() - 1};
        for(int j = pi.size() - 3; j >= 0; j--) {
            double x = pi[j].x, y = pi[j].y;
            double xn = pi[s[0]].x, yn = pi[s[0]].y;
            double xnn = pi[s[1]].x, ynn = pi[s[1]].y;
            //外積が閾値以下（同一直線上と見做せる）なら中間点を削除
            if(std::abs((x-xn)*(yn-ynn)-(y-yn)*(xn-xnn)) <= isline) {
                // 削除用としてマーク
                pi[j + 1].x = HUGE_VAL;
            } else {
                s[1] = s[0];
            }
            s[0] = j;
        }
        // マークされたのを削除
        std::erase_if(pi, [](const auto& pij) { return pij.x == HUGE_VAL; });
    }

    //ここまでで補正した輪郭線を使って立体化を作れば、側面がそこそこ滑らかな立体が出来ると思う

    for(int i = 0; i < p.size(); i++) {
        for(int j = 0; j < p[i].size(); j++) {
            double x0 = p[i][j].x, y0 = p[i][j].y, x1 = p[i][(j + 1) % p[i].size()].x, y1  = p[i][(j + 1) % p[i].size()].y;
            double dx = x1 - x0, dy = y1 - y0;
            int xmin = std::floor(std::max(1.0,std::min(x0,x1)-bev_w));
            int xmax = std::ceil(std::min(w-1.0,std::max(x0,x1)+bev_w));
            int ymin = std::floor(std::max(1.0,std::min(y0,y1)-bev_w));
            int ymax = std::ceil(std::min(h-1.0,std::max(y0,y1)+bev_w));
            for(int y = ymin; y <= ymax; y++) {
                for(int x = xmin; x <= xmax; x++) {
                    double x_ul = x - 0.5, y_ul = y - 0.5;
                    //ポインタ渡し
                    Pixel_Info* info = &pix_info.get()[y*w+x];
                    //[x0,y0]より向こう側については計算の必要無し（次のループで計算される）
                    if(dx * (x_ul - x0) > dy * (y0 - y_ul)) {
                        std::tuple<double, double, double> ret = distance_line(x_ul, y_ul, x0, y0, dx, dy);
                        //比較段階ではmath.sqrt()を使う必要はない
                        if(info->dis > std::get<0>(ret)) {
                            info->dis = std::get<0>(ret);
                            info->x = std::get<1>(ret);
                            info->y = std::get<2>(ret);
                            //輪郭線の内外判定（真理値で保持するより減算した方が速かった、当然か）
                            info->line = dy * (x_ul - x0) - dx * (y_ul - y0);
                        }
                    }
                }
            }
        }
    }

    //いくつかの変数を使い易い形に変換
    for(int j = 0; j < h; j++) {
        for(int i = 0; i < w; i++) {
            Pixel_Info* info = &pix_info[j * w + i];
            info->dis = std::sqrt(info->dis);
            info->x = (info->x - (i - 0.5)) / info->dis;
            info->y = (info->y - (j - 0.5)) / info->dis;

            // cos(atan2(y, x) + r)
            // = cos(atan2(y, x)) * cos(r) - sin(atan2(y, x)) * sin(r)
            // (x, y)は上で正規化したので
            // = x * cos(r) - y * sin(r)
            info->gray = (info->x * cos_rot - info->y * sin_rot) * high * (info->line < 0 ? -1 : 1);
        }
    }

    //四隅の色から中心の色を求めるアンチエイリアス処理
    for(int j = 1; j < h - 1; j++) {
        for(int i = 1; i < w - 1; i++) {
            Pixel_Info* info7 = &pix_info[j * w + i];
            Pixel_Info* info9 = &pix_info[j * w + i + 1];
            Pixel_Info* info1 = &pix_info[j * w + i + w];
            Pixel_Info* info3 = &pix_info[j * w + i + 1 + w];
            double p8 = info7->x == info9->x ? 0.5 : std::clamp((info9->dis - info7->dis + info9->x) / (info9->x - info7->x), 0.0, 1.0);
            double p4 = info7->y == info1->y ? 0.5 : std::clamp((info1->dis - info7->dis + info1->y) / (info1->y - info7->y), 0.0, 1.0);
            double p6 = info9->y == info3->y ? 0.5 : std::clamp((info3->dis - info9->dis + info3->y) / (info3->y - info9->y), 0.0, 1.0);
            double p2 = info1->x == info3->x ? 0.5 : std::clamp((info3->dis - info1->dis + info3->x) / (info3->x - info1->x), 0.0, 1.0);
            //4頂点の平均値なら低コスト
            double x = (p8 + p2 + 1) / 4, y = (p4 + p6 + 1) / 4;
            info7->gray = (info7->gray * (p8 * y + p4 * x) + info9->gray * ((1 - p8) * y + p6 * (1 - x)) + info1->gray * ((1 - p4) * x + p2 * (1 - y)) + info3->gray * ((1 - p6) * (1 - x) + (1 - p2) * (1 - y))) / 2;
            //このままだとピクセル左上隅の座標で計算されるので、雑に補正する
            info7->dis = (info7->dis + info9->dis + info1->dis + info3->dis) / 4.0;
        }
    }

    switch(style) {
        case 0: //ベベル(外側)
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = bgcol_r;
                pix[i].g = bgcol_g;
                pix[i].b = bgcol_b;
                pix[i].a = 0xff;
            }
            utl_putpixeldata(L, pix);
            utl_copybuffer(L, "tmp", "obj");
            utl_drawtarget(L, "tempbuffer");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = col1_r;
                pix[i].g = col1_g;
                pix[i].b = col1_b;
                pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, pix_info[i].gray));
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);
            utl_draw(L, alp1);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = col2_r;
                pix[i].g = col2_g;
                pix[i].b = col2_b;
                pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -pix_info[i].gray));
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);
            utl_draw(L, alp2);
            utl_copybuffer(L, "obj", "tmp");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = static_cast<uint8_t>(bevel_buffer[i].r * (bevel_buffer[i].a / 255.0) + pix[i].r * (1.0 - bevel_buffer[i].a / 255.0));
                pix[i].g = static_cast<uint8_t>(bevel_buffer[i].g * (bevel_buffer[i].a / 255.0) + pix[i].g * (1.0 - bevel_buffer[i].a / 255.0));
                pix[i].b = static_cast<uint8_t>(bevel_buffer[i].b * (bevel_buffer[i].a / 255.0) + pix[i].b * (1.0 - bevel_buffer[i].a / 255.0));
                pix[i].a = std::max(bevel_buffer[i].a, static_cast<uint8_t>(std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0) * 255));
            }
            utl_putpixeldata(L, pix);
            break;
        }
        case 1: //ベベル(内側)
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].a = 0xff;
            }
            utl_putpixeldata(L, pix);
            utl_copybuffer(L, "tmp", "obj");
            utl_drawtarget(L, "tempbuffer");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);
            utl_draw(L, alp1);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);
            utl_draw(L, alp2);
            utl_copybuffer(L, "obj", "tmp");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].a = bevel_buffer[i].a;
            }
            
            utl_putpixeldata(L, pix);
            break;
        }
        case 2: //エンボス
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = static_cast<uint8_t>(pix[i].r * (pix[i].a / 255.0) + bgcol_r * (1 - pix[i].a / 255.0));
                pix[i].g = static_cast<uint8_t>(pix[i].g * (pix[i].a / 255.0) + bgcol_g * (1 - pix[i].a / 255.0));
                pix[i].b = static_cast<uint8_t>(pix[i].b * (pix[i].a / 255.0) + bgcol_b * (1 - pix[i].a / 255.0));
                pix[i].a = 0xff;
            }
            utl_putpixeldata(L, pix);
            utl_copybuffer(L, "tmp", "obj");
            utl_drawtarget(L, "tempbuffer");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, pix_info[i].gray) * (pix_info[i].line < 0 ? 1 : std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);
            utl_draw(L, alp1);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -pix_info[i].gray) * (pix_info[i].line < 0 ? 1 : std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);
            utl_draw(L, alp2);
            utl_copybuffer(L, "obj", "tmp");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].a = std::max(bevel_buffer[i].a, static_cast<uint8_t>(std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0) * 255.0));
            }
            
            utl_putpixeldata(L, pix);
            break;
        }
        case 3: //ピローエンボス
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].r = static_cast<uint8_t>(pix[i].r * (pix[i].a / 255.0) + bgcol_r * (1 - pix[i].a / 255.0));
                pix[i].g = static_cast<uint8_t>(pix[i].g * (pix[i].a / 255.0) + bgcol_g * (1 - pix[i].a / 255.0));
                pix[i].b = static_cast<uint8_t>(pix[i].b * (pix[i].a / 255.0) + bgcol_b * (1 - pix[i].a / 255.0));
                pix[i].a = 0xff;
            }
            utl_putpixeldata(L, pix);
            utl_copybuffer(L, "tmp", "obj");
            utl_drawtarget(L, "tempbuffer");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, (bevel_buffer[i].a * 2.0 / 255.0 - 1.0) * pix_info[i].gray) * (pix_info[i].line < 0 ? 1 : std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);
            utl_draw(L, alp1);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -(bevel_buffer[i].a * 2.0 / 255.0 - 1.0) * pix_info[i].gray) * (pix_info[i].line < 0 ? 1 : std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);
            utl_draw(L, alp2);
            utl_copybuffer(L, "obj", "tmp");

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                pix[i].a = std::max(bevel_buffer[i].a, static_cast<uint8_t>(std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0) * 255.0));
            }

            utl_putpixeldata(L, pix);

            break;
        }
        case 4: //ベベル(外側)直接描画
        {
            utl_drawtarget(L, "framebuffer");

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            ObjField obj = getObjField(L);
            lua_getfield(L, -1, "draw");
            lua_call(L, 0, 0);

            //obj.effect()で解像度が変わる可能性があるので再読み込み
            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(std::min(255.0 - bevel_buffer[i].a, 255.0 * std::max(0.0, pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);

            utl_draw(L, alp1);
            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(std::min(255.0 - bevel_buffer[i].a, 255.0 * std::max(0.0, -pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);

            utl_draw(L, alp2);
            break;
        }
        case 5: //ベベル(内側)直接描画
        {
            utl_drawtarget(L, "framebuffer");

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            ObjField obj = getObjField(L);
            lua_getfield(L, -1, "draw");
            lua_call(L, 0, 0);

            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(std::min(static_cast<double>(bevel_buffer[i].a), 255.0 * std::max(0.0, pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);

            utl_draw(L, alp1);
            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(std::min(static_cast<double>(bevel_buffer[i].a), 255.0 * std::max(0.0, -pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            
            utl_draw(L, alp2);
            break;
        }
        case 6: //エンボス 直接描画
        {
            utl_drawtarget(L, "framebuffer");

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            ObjField obj = getObjField(L);
            lua_getfield(L, -1, "draw");
            lua_call(L, 0, 0);

            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);

            utl_draw(L, alp1);
            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);

            utl_draw(L, alp2);
            break;
        }
        case 7: //ピローエンボス 直接描画
        {
            utl_drawtarget(L, "framebuffer");

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            ObjField obj = getObjField(L);
            lua_getfield(L, -1, "draw");
            lua_call(L, 0, 0);

            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");

            setObjField(L, obj);

            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, (bevel_buffer[i].a * 2 / 255.0 - 1) * pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble1);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            
            utl_draw(L, alp1);
            utl_temptarget(L, w, h);
            utl_copybuffer(L, "obj", "tmp");
            utl_drawtarget(L, "framebuffer");
            
            setObjField(L, obj);

            pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, (bevel_buffer[i].a * 2 / 255.0 - 1) * (-pix_info[i].gray)) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            utl_blend(L, ble2);

            lua_getfield(L, -1, "effect");
            lua_call(L, 0, 0);
            
            utl_draw(L, alp2);
            break;
        }
        case 8: //ベベル(外側) ハイライトのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(std::min(255.0 - bevel_buffer[i].a, 255.0 * std::max(0.0, pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 9: //ベベル(内側) ハイライトのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(std::min(static_cast<double>(bevel_buffer[i].a), 255.0 * std::max(0.0, pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 10: //エンボス ハイライトのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 11: //ピローエンボス ハイライトのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col1_r;
                    pix[i].g = col1_g;
                    pix[i].b = col1_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, (bevel_buffer[i].a * 2 / 255.0 - 1) * pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 12: //ベベル(外側) シャドウのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(std::min(255.0 - bevel_buffer[i].a, 255.0 * std::max(0.0, -pix_info[i].gray) * std::min(1.0, bev_w-pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 13: //ベベル(内側) シャドウのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(std::min(static_cast<double>(bevel_buffer[i].a), 255.0 * std::max(0.0, -pix_info[i].gray) * std::min(1.0, bev_w - pix_info[i].dis)));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 14: //エンボス シャドウのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, -pix_info[i].gray) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
        case 15: //ピローエンボス シャドウのみ
        {
            Pixel_BGRA* pix = reinterpret_cast<Pixel_BGRA*>(utl_getpixeldata(L));
            for(int i = 0; i < w * h; i++) {
                if(pix_info[i].dis <= bev_w) {
                    pix[i].r = col2_r;
                    pix[i].g = col2_g;
                    pix[i].b = col2_b;
                    pix[i].a = static_cast<uint8_t>(255.0 * std::max(0.0, (bevel_buffer[i].a * 2 / 255.0 - 1) * (-pix_info[i].gray)) * std::clamp(bev_w - pix_info[i].dis, 0.0, 1.0));
                } else {
                    pix[i].r = 0;
                    pix[i].g = 0;
                    pix[i].b = 0;
                    pix[i].a = 0;
                }
            }
            utl_putpixeldata(L, pix);
            utl_blur(L, blur);
            break;
        }
    }

    utl_blend(L, 0);
    utl_drawtarget(L, "framebuffer");

    return 0;
}

static luaL_Reg functions[] = {
    { "Bevel_And_Emboss", bevel_and_emboss },
    { nullptr, nullptr }
};

EXTERN_C __declspec(dllexport) int luaopen_Bevel_And_Emboss_M(lua_State* L) {
    luaL_register(L, "Bevel_And_Emboss_M", functions);
    return 1;
}
