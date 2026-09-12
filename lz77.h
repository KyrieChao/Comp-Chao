#ifndef LZ77_H
#define LZ77_H

// buf 全文，pos 当前位置，window 窗口大小
// 找到了返回匹配长度，并把距离写进 *out_dist；没找到返回 0
int find_match(const unsigned char* buf, int buf_len, int pos, int window, int* out_dist);

#endif