// main.c —— 命令行入口 + 文件读写 + 比特打包。
//
// 这里是「字节世界」和「符号世界」的边界：
//   huffman_algo 只认符号（unsigned short），不知道字节是什么
//   lz77 把字节流拆成 (符号流, 距离流)
//   这个文件负责把两条流拼成 .chao 文件，再原样拆回来
//
// .chao 文件的布局（按写入顺序）：
//
//   [data_len     : 4]   原文件字节数
//   [leaf_count   : 4]   哈夫曼叶子数
//   [struct_bytes : 4]   结构串压完占几个字节
//   [dist_bytes   : 4]   距离段占几个字节
//   [sym_n        : 4]   符号流有多少个符号
//   [dist_n       : 4]   距离流有多少个距离
//   [结构串       : struct_bytes]   前序遍历，内部 '1'、叶子 '0'，按位打包
//   [符号表       : leaf_count * 2] 每个叶子一个 unsigned short，按叶子前序遍历顺序
//   [距离段       : dist_bytes]     每条 5 位档位号 + 若干 extra，按位打包
//   [载荷         : 到文件尾]       符号流的哈夫曼码流，末字节补 0 对齐
//
// 载荷放最后，是为了它的长度能直接用 EOF - 当前位置 算出来，省掉一个长度字段。

#include "huffman_algo.h"
#include "lz77.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_STRUCT_STR (2 * ALPHABET_SIZE)
#define MAX_STRUCT_BYTES ((MAX_STRUCT_STR + 7) / 8)

static const int d_base[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
static const int d_extra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

int dist_code(int n) {
    if (n >= d_base[29]) return 29;
    for (int i = 0; i < 29; i++) {
        if ((d_base[i] <= n && n < d_base[i + 1])) return i;
    }
    return -1;
}

// 把 value 的低 n 位写进位串，高位在前，游标推进 n
static void put_bits(unsigned char* bits, size_t* pos, unsigned int value, int n) {
    for (int i = 0; i < n; i++) {
        bits[(*pos)++] = '0' + ((value >> (n - i - 1)) & 1);
    }
    bits[*pos] = '\0';
}
// 从位串读 n 位，高位在前，游标推进 n
static unsigned int get_bits(const unsigned char* bits, size_t* pos, int n) {
    int number = 0;
    for (int i = 0; i < n; i++) {
        number = (number << 1) | (bits[((*pos)++)] - '0');
    }
    return number;
}
// 把 '0'/'1' 串压成字节。
// bit_count 是位数，返回压出来的字节数 = ceil(bit_count/8)。
// 最后一个字节不够 8 位时，低位用 0 补满。
int bit_to_byte(const unsigned char* bit_str, unsigned char* packed, size_t bit_count) {
    size_t byte_count = (bit_count + 7) >> 3;
    for (size_t i = 0; i < byte_count; i++) {
        packed[i] = 0;
        size_t pos = i * 8;
        int k = 8;
        while (k-- > 0) {
            if (pos >= bit_count) break;
            packed[i] = (packed[i] << 1) | (bit_str[pos++] & 1);
        }
        packed[i] <<= 8 - (pos - i * 8);
    }
    return (int)byte_count;
}

// bit_to_byte 的逆操作：把 packed 摊回 bit_count 个 '0'/'1' 字符。
// 从每个字节的最高位开始读，跟 bit_to_byte 的写入口径一致。
// 末尾会补一个 '\0'，所以 bit_str 至少要能放 bit_count + 1 个字节。
int byte_to_bit(const unsigned char* packed, unsigned char* bit_str, int bit_count) {
    for (int i = 0; i < bit_count; i++) {
        int byte_pos = i / 8;
        int bit_pos = 7 - (i % 8);
        bit_str[i] = (packed[byte_pos] & (1 << bit_pos)) ? '1' : '0';
    }
    bit_str[bit_count] = '\0';
    return bit_count;
}

// 把 (符号流, 距离流) 编码成 .chao 文件。成功返回 0。
int compress(const char* out_name,
             const unsigned short* sym, int sym_n,
             const unsigned short* dist, int dist_n,
             int data_len) {
    unsigned int freq[ALPHABET_SIZE] = {0};
    int leaf_count = counter(sym, (size_t)sym_n, freq);
    // 空文件
    if (leaf_count == 0) return -1;

    // 建树，并算出每个字节的码
    if (HuFF_Init(freq, (size_t)sym_n) != 0) return -1;
    HuFF_Build();
    HuFF_Code_Table();

    // 编码后的总位数，后面所有长度都从它推
    int total_bits = 0;
    HuFF_Wpl(&total_bits);
    size_t payload_bytes = (total_bits + 7) / 8;
    int struct_bytes = (2 * leaf_count - 1 + 7) / 8;

    // 码串：一个字符一位，先在这儿拼好，再整块压成字节
    unsigned char* bit_str = malloc((size_t)total_bits + 1);
    unsigned char* payload = malloc(payload_bytes + 1);
    char struct_str[MAX_STRUCT_STR];
    unsigned char struct_packed[MAX_STRUCT_BYTES];
    unsigned short symbols[ALPHABET_SIZE];
    if (!bit_str || !payload) {
        free(bit_str);
        free(payload);
        HuFF_Destroy();
        return -1;
    }

    // 结构串 + 符号表：同一趟前序遍历里一起产出，顺序天然对齐
    HuFF_Struct(struct_str, MAX_STRUCT_STR, symbols);
    bit_to_byte((unsigned char*)struct_str, struct_packed, (size_t)(2 * leaf_count - 1));

    // 原文 -> 码串 -> 载荷
    HuFF_Get(sym, bit_str, (size_t)sym_n);
    bit_to_byte(bit_str, payload, (size_t)total_bits);

    HuFF_Destroy();

    // --- todo ---
    unsigned int d_freq[ALPHABET_SIZE] = {0};
    for (int i = 0; i < dist_n; i++) d_freq[dist_code(dist[i])]++;
    if (HuFF_Init(d_freq, (size_t)dist_n) != 0) {
        free(bit_str);
        free(payload);
        return -1;
    }
    HuFF_Build();
    HuFF_Code_Table();
    char d_struct_str[MAX_STRUCT_STR];
    unsigned short d_symbols[ALPHABET_SIZE];
    unsigned char d_struct_packed[MAX_STRUCT_BYTES];
    HuFF_Struct(d_struct_str, MAX_STRUCT_STR, d_symbols);
    int d_leaf_count = 0;
    for (size_t i = 0; i < ALPHABET_SIZE; i++) {
        if (d_freq[i] > 0) d_leaf_count++;
    }
    int d_struct_bytes = bit_to_byte((unsigned char*)d_struct_str, d_struct_packed, (size_t)(2 * d_leaf_count - 1));
    // -----------
    FILE* fp = fopen(out_name, "wb");
    if (!fp) {
        free(bit_str);
        free(payload);
        return -1;
    }
    // ---- 距离段：每条写 5 位档位号，后面跟若干位 extra ----
    // 宽度随档位变（近的距离 extra 短），所以先算总位数，才知道位串开多大。
    size_t dist_bit_count = 0;
    for (int i = 0; i < dist_n; i++) {
        int k = dist_code(dist[i]);
        dist_bit_count += strlen(HuFF_Code(k)) + d_extra[k];
    }

    unsigned char* dist_bit_str = malloc(dist_bit_count + 1);
    if (!dist_bit_str) {
        free(bit_str);
        free(payload);
        fclose(fp);
        return -1;
    }
    size_t pos = 0;
    for (int i = 0; i < dist_n; i++) {
        int k = dist_code(dist[i]);
        const char* code = HuFF_Code(k);
        for (int j = 0; code[j]; j++) dist_bit_str[pos++] = (unsigned char)code[j];
        put_bits(dist_bit_str, &pos, dist[i] - d_base[k], d_extra[k]);
    }
    HuFF_Destroy();

    size_t dist_bytes = (dist_bit_count + 7) / 8;
    unsigned char* dist_packed = malloc(1 + dist_bytes);
    if (!dist_packed) {
        free(bit_str);
        free(payload);
        free(dist_bit_str);
        fclose(fp);
        return -1;
    }
    int dist_packed_count = bit_to_byte((unsigned char*)dist_bit_str, dist_packed, dist_bit_count);
    fwrite(&data_len, sizeof(int), 1, fp);
    fwrite(&leaf_count, sizeof(int), 1, fp);
    fwrite(&struct_bytes, sizeof(int), 1, fp);
    fwrite(&dist_bytes, sizeof(int), 1, fp);
    fwrite(&sym_n, sizeof(int), 1, fp);
    fwrite(&dist_n, sizeof(int), 1, fp);
    fwrite(&d_leaf_count, sizeof(int), 1, fp);
    fwrite(&d_struct_bytes, sizeof(int), 1, fp);

    fwrite(struct_packed, 1, (size_t)struct_bytes, fp);
    fwrite(symbols, 1, sizeof(unsigned short) * leaf_count, fp);
    fwrite(d_struct_packed, 1, (size_t)d_struct_bytes, fp);
    fwrite(d_symbols, 1, sizeof(unsigned short) * d_leaf_count, fp);
    fwrite(dist_packed, 1, dist_packed_count, fp);
    fwrite(payload, 1, payload_bytes, fp);

    long file_size = ftell(fp);
    printf("%d byte -> %ld byte (%.1f%%)\n", data_len, file_size, 100.0 * file_size / data_len);

    fclose(fp);
    free(bit_str);
    free(payload);
    free(dist_packed);
    free(dist_bit_str);
    return 0;
}

// 把 .chao 文件还原成原文件。成功返回 0。
int decompress(const char* in_name, const char* out_name) {
    FILE* fp = fopen(in_name, "rb");
    if (!fp) return -1;

    int data_len, leaf_count, struct_bytes, sym_n, dist_n, dist_bytes, d_leaf_count, d_struct_bytes;
    if (fread(&data_len, sizeof(int), 1, fp) != 1 ||
        fread(&leaf_count, sizeof(int), 1, fp) != 1 ||
        fread(&struct_bytes, sizeof(int), 1, fp) != 1 ||
        fread(&dist_bytes, sizeof(int), 1, fp) != 1 ||
        fread(&sym_n, sizeof(int), 1, fp) != 1 ||
        fread(&dist_n, sizeof(int), 1, fp) != 1 ||
        fread(&d_leaf_count, sizeof(int), 1, fp) != 1 ||
        fread(&d_struct_bytes, sizeof(int), 1, fp) != 1) {
        fclose(fp);
        return -1;
    }

    int struct_bits = 2 * leaf_count - 1; // 结构串有多少位

    unsigned char struct_packed[MAX_STRUCT_BYTES] = {0};
    char struct_str[MAX_STRUCT_STR];
    for (int i = 0; i < struct_bytes; i++) {
        if (fread(&struct_packed[i], 1, 1, fp) != 1) {
            fclose(fp);
            return -1;
        }
    }

    unsigned short symbols[ALPHABET_SIZE] = {0};
    if (fread(symbols, sizeof(unsigned short), (size_t)leaf_count, fp) != (size_t)leaf_count) {
        fclose(fp);
        return -1;
    }

    byte_to_bit(struct_packed, (unsigned char*)struct_str, struct_bits);

    if (HuFF_Rebuild((unsigned char*)struct_str, symbols) != 0) {
        fclose(fp);
        return -1;
    }
    unsigned char* dist_bit_str = malloc((size_t)dist_bytes * 8 + 1);
    if (!dist_bit_str) {
        fclose(fp);
        return -1;
    }
    unsigned char* dist_packed = malloc(dist_bytes + 1);
    if (!dist_packed) {
        free(dist_bit_str);
        fclose(fp);
        return -1;
    }
    if (fread(dist_packed, 1, (size_t)dist_bytes, fp) != (size_t)dist_bytes) {
        HuFF_Destroy();
        free(dist_bit_str);
        free(dist_packed);
        fclose(fp);
        return -1;
    }
    byte_to_bit(dist_packed, dist_bit_str, (int)((size_t)dist_bytes * 8));

    unsigned short* dist = malloc(dist_n * sizeof(unsigned short));
    if (!dist) {
        HuFF_Destroy();
        free(dist_bit_str);
        free(dist_packed);
        fclose(fp);
        return -1;
    }
    size_t pos = 0;
    for (int i = 0; i < dist_n; i++) {
        int k = get_bits(dist_bit_str, &pos, 5);
        dist[i] = d_base[k] + get_bits(dist_bit_str, &pos, d_extra[k]);
    }
    // 4. 载荷 = 文件剩下的全部字节。
    //    载荷是文件的最后一段，所以直接拿 EOF 兜底，不必再存一个长度字段。
    long payload_start = ftell(fp);
    fseek(fp, 0, SEEK_END);
    long payload_bytes = ftell(fp) - payload_start;
    fseek(fp, payload_start, SEEK_SET);

    unsigned char* payload = malloc((size_t)payload_bytes + 1);
    unsigned char* bit_str = malloc((size_t)payload_bytes * 8 + 1);
    unsigned char* final = malloc(((size_t)data_len + 1));
    unsigned short* out_buf = malloc(((size_t)sym_n + 1) * sizeof(unsigned short));
    if (!payload || !bit_str || !out_buf || !final) {
        free(payload);
        free(bit_str);
        free(out_buf);
        HuFF_Destroy();
        fclose(fp);
        return -1;
    }
    if (payload_bytes > 0 &&
        fread(payload, 1, (size_t)payload_bytes, fp) != (size_t)payload_bytes) {
        free(payload);
        free(bit_str);
        free(out_buf);
        HuFF_Destroy();
        fclose(fp);
        return -1;
    }
    fclose(fp);

    // 5. 摊成位串再解码。
    //    载荷最后一个字节的低位可能带着补的 0，但 HuFF_Decode 只要凑够
    //    data_len 个字符就停，那些补出来的位不会被碰到。
    byte_to_bit(payload, bit_str, (int)payload_bytes * 8);
    HuFF_Decode(bit_str, out_buf, (size_t)sym_n);
    HuFF_Destroy();

    FILE* out = fopen(out_name, "wb");
    if (!out) {
        free(payload);
        free(bit_str);
        free(out_buf);
        return -1;
    }
    long final_len = 0;
    lz77_decode(out_buf, sym_n, dist, dist_n, final, &final_len);
    fwrite(final, 1, (size_t)final_len, out);
    fclose(out);

    free(payload);
    free(bit_str);
    free(out_buf);
    free(final);
    free(dist);
    free(dist_packed);
    free(dist_bit_str);
    return 0;
}

// ==================== 文件名处理 ====================

// 文件存不存在（能打开就算存在）
static int file_exists(const char* path) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

// 找最后一个 '.' 的位置，找不到就返回 len。
// 只看文件名那一段，不越过路径分隔符 —— "dir.v2/file" 里的点不算扩展名。
static size_t dot_pos(const char* path) {
    size_t len = strlen(path);
    for (size_t i = len; i > 0; i--) {
        char c = path[i - 1];
        if (c == '/' || c == '\\') break;
        if (c == '.') return i - 1;
    }
    return len;
}

// 把扩展名换成 ext："back2.txt" + ".chao" -> "back2.chao"
// 名字以 '.' 开头（隐藏文件）时，当作没有扩展名。
static void replace_ext(const char* path, const char* ext, char* out, size_t out_size) {
    size_t base = dot_pos(path);
    if (base == 0) base = strlen(path);
    size_t ext_len = strlen(ext);
    if (base + ext_len + 1 > out_size) {
        snprintf(out, out_size, "%s", path); // 装不下就原样返回
        return;
    }
    memcpy(out, path, base);
    memcpy(out + base, ext, ext_len + 1); // 连结尾的 '\0' 一起拷
}

// 在扩展名前面插一段后缀："demo.txt" + "_5675" -> "demo_5675.txt"
static void insert_suffix(const char* path, const char* suffix, char* out, size_t out_size) {
    size_t base = dot_pos(path);
    if (base == 0) base = strlen(path);
    snprintf(out, out_size, "%.*s%s%s", (int)base, path, suffix, path + base);
}

// 输入路径的短 hash（FNV-1a），用来让撞名的输出文件有个能认出来的名字
static unsigned int short_hash(const char* s) {
    unsigned int h = 2166136261u;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 16777619u;
    }
    return h;
}

// 挑一个还没被占用的输出名：
//   demo.txt 被占了 -> demo_5675.txt -> 还被占 -> demo_5675_2.txt -> ...
static void pick_output_name(const char* wanted, unsigned int hash, char* out, size_t out_size) {
    if (!file_exists(wanted)) {
        snprintf(out, out_size, "%s", wanted);
        return;
    }

    char suffix[32];
    snprintf(suffix, sizeof(suffix), "_%04x", hash & 0xFFFFu);
    insert_suffix(wanted, suffix, out, out_size);
    if (!file_exists(out)) return;

    for (int n = 2; n < 10000; n++) {
        snprintf(suffix, sizeof(suffix), "_%04x_%d", hash & 0xFFFFu, n);
        insert_suffix(wanted, suffix, out, out_size);
        if (!file_exists(out)) return;
    }
}

// ==================== 入口 ====================
// 用法：
//   huffman c <原文> [压缩包]      省略输出名 -> 把扩展名换成 .chao
//   huffman d <压缩包> [还原文件]   省略输出名 -> 把扩展名换成 .txt，撞名就插一段短 hash
int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("usage:\n");
        printf("  %s c <input> [output]    compress;   default output: <input>.chao\n", argv[0]);
        printf("  %s d <input> [output]    decompress; default output: <input>.txt\n", argv[0]);
        return 1;
    }

    char out_name[1024];

    if (argv[1][0] == 'c') {
        // 输出名：给了就用给的，没给就把原文的扩展名换成 .chao
        if (argc >= 4) {
            snprintf(out_name, sizeof(out_name), "%s", argv[3]);
        } else {
            replace_ext(argv[2], ".chao", out_name, sizeof(out_name));
        }

        FILE* in = fopen(argv[2], "rb");
        if (!in) {
            printf("cannot open %s\n", argv[2]);
            return 1;
        }
        fseek(in, 0, SEEK_END);
        long size = ftell(in);
        fseek(in, 0, SEEK_SET);
        unsigned char* raw = malloc(size);
        if (!raw) {
            fclose(in);
            return 1;
        }
        if (size > 0 && fread(raw, 1, (size_t)size, in) != (size_t)size) {
            fclose(in);
            free(raw);
            return 1;
        }
        fclose(in);
        unsigned short* sym = malloc(size * sizeof(unsigned short));
        unsigned short* dist = malloc(size * sizeof(unsigned short));
        if (!sym || !dist) {
            free(sym);
            free(dist);
            return 1;
        }
        int sn = 0, dn = 0;
        lz77_encode(raw, size, sym, dist, &sn, &dn);
        if (compress(out_name, sym, sn, dist, dn, (int)size) != 0) {
            printf("compress failed\n");
            return 1;
        }
        printf("-> %s\n", out_name);
        free(sym);
        free(dist);
        free(raw);

    } else if (argv[1][0] == 'd') {
        // 输出名：给了就用给的（明确指定就照写，允许覆盖）；
        //         没给就把扩展名换成 .txt，如果已经存在就换个带短 hash 的名字。
        if (argc >= 4) {
            snprintf(out_name, sizeof(out_name), "%s", argv[3]);
        } else {
            char wanted[1024];
            replace_ext(argv[2], ".txt", wanted, sizeof(wanted));
            pick_output_name(wanted, short_hash(argv[2]), out_name, sizeof(out_name));
        }

        if (decompress(argv[2], out_name) != 0) {
            printf("decompress failed\n");
            return 1;
        }
        printf("-> %s\n", out_name);

    } else {
        printf("first arg must be c or d\n");
        return 1;
    }
    return 0;
}