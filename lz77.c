#include "lz77.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 在 [pos - window, pos) 这段区间里，找一条和 buf[pos...] 开头最像的串。
// 命中就把距离写进 *out_dist 并返回匹配长度；没命中返回 0。
//
// 找法是最朴素的：窗口里每个起点都试一遍，从头逐字节往下比，留最长的那条。
// 复杂度 O(窗口长度 x 匹配长度)，大文件上会慢 —— 换成哈希链是后面的活。
int find_match(const unsigned char* buf, long buf_len, int pos, int window, int* out_dist) {
    int start = pos - window < 0 ? 0 : (pos - window);
    int best_len = 0;
    int best_dist = 0;
    for (int i = start; i < pos; i++) {
        int len = 0;
        while (len < MAX_MATCH && pos + len < buf_len && buf[i + len] == buf[pos + len]) len++;
        if (len >= best_len) {
            best_len = len;
            best_dist = pos - i;
        }
    }
    *out_dist = best_dist;
    return best_len;
}

// 把原文扫成 (符号流, 距离流)。返回 0 成功。
// 符号流里：字面量直接写字节值；匹配写 257 + 长度 - 3，距离另存进 dist。
// 两条流是配对的 —— 符号流里每出现一个 >= 257 的值，dist 里就有一个距离对应它。
int lz77_encode(const unsigned char* data, long data_len, unsigned short* sym,
                unsigned short* dist, int* out_sym_n, int* out_dist_n) {

    int pos = 0, mlen = 0, num = 0, n_s = 0, n_d = 0;
    while (pos < data_len) {
        // 窗口 32767 = 2^15 - 1，跟 Deflate 取齐
        mlen = find_match(data, data_len, pos, 32767, &num);
        // 最短匹配取 3：再短的匹配，光「长度符号 + 距离」的代价就盖过直接存字面量了。
        // 这个阈值跟距离的编码代价挂钩，等距离压得更便宜之后要重算。
        if (mlen >= 3) {
            sym[n_s++] = 257 + mlen - 3;
            dist[n_d++] = num;
            pos += mlen;
        } else {
            sym[n_s++] = data[pos];
            pos += 1;
        }
    }
    *out_sym_n = n_s;
    *out_dist_n = n_d;
    return 0;
}

// 把两条流还原成原文。返回 0 成功；符号流和距离流对不上号时返回 -1。
// 字面量直接吐出来；匹配则从「已经吐出来的部分」往前翻 dist 个字节，抄 len 个。
// 边抄边写，所以源和目标重叠也没关系 —— 重复的模式会自己展开（比如 "aaaa"）。
//
// *out_len 是带进带出的游标：函数从它当前的位置接着往外写，
// 所以调用方必须先把它初始化成 0。
int lz77_decode(const unsigned short* sym, int sym_n, const unsigned short* dist,
                int dist_n, unsigned char* out, long* out_len) {
    int di = 0;
    for (int i = 0; i < sym_n; i++) {
        if (sym[i] < 256)
            out[(*out_len)++] = sym[i];
        else {
            size_t len = sym[i] - 257 + 3;
            long src = *out_len - dist[di];
            for (size_t k = 0; k < len; k++) out[(*out_len)++] = out[src++];
            di++;
        }
    }
    if (di != dist_n) return -1;
    return 0;
}
#ifdef LZ77_TEST
int main(void) {
    char file_name[] = "demo.txt";
    FILE* in = fopen(file_name, "rb");
    if (!in) {
        printf("cannot open %s\n", file_name);
        return 1;
    }
    fseek(in, 0, SEEK_END);
    long size = ftell(in);
    fseek(in, 0, SEEK_SET);

    unsigned char* data = malloc((size_t)size + 1);
    if (!data) {
        fclose(in);
        return 1;
    }
    if (size > 0 && fread(data, 1, (size_t)size, in) != (size_t)size) {
        fclose(in);
        free(data);
        return 1;
    }

    fclose(in);
    unsigned short sym[8092], dist[8092];
    int sn = 0, dn = 0;

    lz77_encode(data, size, sym, dist, &sn, &dn);
    printf("%d %d\n", sn, dn);
    unsigned char* back = malloc(size + 1);
    long back_len = 0;
    lz77_decode(sym, sn, dist, dn, back, &back_len);
    printf("identical=%s (back_len=%ld)\n",
           (back_len == size && memcmp(back, data, size) == 0) ? "YES" : "NO", back_len);
    free(data);
    return 0;
}
#endif