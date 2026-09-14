#ifndef LZ77_H
#define LZ77_H

// 匹配长度上限。
// 为什么必须封顶：长度符号的值是 257 + 长度 - 3，字母表定成 513 就是照 258 反推的
//（257 + 258 - 3 = 512）。不封顶的话，长重复段会让符号值越界。
#define MAX_MATCH 258

// 把原文扫成两条流：
//   sym  —— 符号流。字面量直接写字节值；匹配写 257 + 长度 - 3
//   dist —— 距离流。每个匹配对应的回退距离，和 sym 里的匹配按顺序一一对应
// sym / dist 由调用方分配，大小至少各 data_len 个。返回 0 成功。
int lz77_encode(const unsigned char* data, long data_len,
                unsigned short* sym, unsigned short* dist,
                int* out_sym_n, int* out_dist_n);

// 把两条流还原成原文。返回 0 成功。
// 字面量直接吐出来；匹配则从「已经吐出来的部分」往前翻 dist 个字节，抄 len 个。
//
// out_len 是带进带出的游标：函数从它当前的位置接着往外写。
// 所以调用方必须先把它初始化成 0，否则会从垃圾位置开始写。
int lz77_decode(const unsigned short* sym, int sym_n, const unsigned short* dist,
                int dist_n, unsigned char* out, long* out_len);
#endif