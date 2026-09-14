#ifndef HUFFMAN_ALGO_H
#define HUFFMAN_ALGO_H

#include <stddef.h>

#define ALPHABET_SIZE 513 // 符号有多少种取值
#define MAX_CODE_BITS 256 // 码最长能有几位

typedef struct node Node;

// ---------------- 频次统计 ----------------
// 数 data 里每个符号出现了多少次，结果写进 freq[]（下标就是符号值）。
// 返回值 = 出现过多少种不同字节，也就是叶子数。
int counter(const unsigned short* data, size_t len, unsigned int* freq);

// ---------------- 建树（编码端）----------------
// 按频次摆好一片森林：每个出现过的字节是一个叶子。
// 返回 0 正常，1 表示输入是空的，-1 表示内存不够。
int HuFF_Init(const unsigned int* freq, size_t len);

// 反复挑出权最小的两棵树合并，直到只剩一棵。根在 forest[0]。
int HuFF_Build();

// 释放整片森林，顺便把码表和全局状态清零。
void HuFF_Destroy();

// ---------------- 编码（写端）----------------
// 算出编码后的总位数：Σ(叶子频次 × 叶子深度)。depth 传 0。
void HuFF_Wpl(int* total_bits);

// 从树里读出每个字节的码，填进码表。
void HuFF_Code_Table();

// 把 data 的每个字节换成码，依次拼成一条 '0'/'1' 串写进 out。
// out 至少要能放 total_bits 个字符。
void HuFF_Get(const unsigned short* data, unsigned char* out, size_t len);

// 取某个符号的码串（'0'/'1'）。没出现过的符号返回空串。
const char* HuFF_Code(unsigned short sym);

// 把树的形状写成「结构串」：前序遍历，内部节点写 '1'，叶子写 '0'。
//
// 同一趟遍历里，按「叶子被访问到的顺序」把每个叶子的字节收集进 symbols。
// 这个顺序必须和结构串里 '0' 出现的顺序完全一致 ——
// 解压端就是照这个顺序往树上挂叶子的，错位一个就全乱。
// 注意：它不是按字节值排的，两者不一定相同。
//
// struct_out 至少要能放 2*叶子数 个字节（2*叶子数-1 个标记 + 结尾的 '\0'）。
// symbols    至少要能放 叶子数 个 unsigned short。
void HuFF_Struct(char* struct_out, int capacity, unsigned short* symbols);

// ---------------- 解码（读端）----------------
// 沿树走 bit_str 里的 '0'/'1'：'0' 走左、'1' 走右，
// 每走到一个叶子就吐出一个字节，一共吐 len 个。
void HuFF_Decode(const unsigned char* bit_str, unsigned short* out, size_t len);

// 从位串里解出一个符号，游标推进它用掉的位数。
unsigned short HuFF_Decode_One(const unsigned char* bit_str, size_t* pos);

// 从「结构串 + 符号表」把树重建出来。
// 建完之后 forest[0] 是根、count = 1，后面的 Code_Table / Decode 可以直接用。
int HuFF_Rebuild(const unsigned char* struct_bits, const unsigned short* symbols);

#endif