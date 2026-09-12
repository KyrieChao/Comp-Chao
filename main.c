#include "huffman_algo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_STRUCT_STR (2 * MAX_CODE_LEN)
#define MAX_STRUCT_BYTES ((MAX_STRUCT_STR + 7) / 8)

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

// bit_to_byte 的逆操作：把字节摊回 '0'/'1'。
// bit_count 是要摊多少「位」（不是字节数）。
// 末尾会补 '\0'，所以 out 至少要 bit_count+1 个字节。
int byte_to_bit(const unsigned char* packed, unsigned char* bit_str, int bit_count) {
    for (int i = 0; i < bit_count; i++) {
        int byte_pos = i / 8;
        int bit_pos = 7 - (i % 8);
        bit_str[i] = (packed[byte_pos] & (1 << bit_pos)) ? '1' : '0';
    }
    bit_str[bit_count] = '\0';
    return bit_count;
}

// 格式：
//   [data_len     : 4 字节] 原文件多少字节（解压时靠它知道要吐多少个字符）
//   [leaf_count   : 4 字节] 叶子数 = 不同字节数 = 结构串里 '0' 的个数
//   [struct_bytes : 4 字节] 结构串压完占几个字节
//   [符号表   : leaf_count 字节]   按叶子被前序遍历访问到的顺序排列（不是按字节值排的！）
//   [结构串   : struct_bytes 字节] 树的形状，已经压成字节
//   [载荷     : 一直到文件尾]       原文编码后的码流
int compress(const char* out_name, const unsigned char* data, int data_len) {
    unsigned int freq[MAX_CODE_LEN] = {0};
    int leaf_count = counter(data, (size_t)data_len, freq);
    // 空文件
    if (leaf_count == 0) return -1;

    // 建树，并算出每个字节的码
    if (HuFF_Init(freq, (size_t)data_len) != 0) return -1;
    HuFF_Build();
    HuFF_Code_Table();

    // 编码后的总位数，后面所有长度都从它推
    int total_bits = 0;
    HuFF_Wpl(&total_bits);
    int payload_bytes = (total_bits + 7) / 8;
    int struct_bytes = (2 * leaf_count - 1 + 7) / 8;

    // 码串：一个字符一位，先在这儿拼好，再整块压成字节
    unsigned char* bit_str = malloc((size_t)total_bits + 1);
    unsigned char* payload = malloc((size_t)payload_bytes + 1);
    char struct_str[MAX_STRUCT_STR];
    unsigned char struct_packed[MAX_STRUCT_BYTES];
    unsigned char symbols[MAX_CODE_LEN]; // 符号表：顺序由 HuFF_Struct 一趟产出来
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
    HuFF_Get(data, bit_str, (size_t)data_len);
    bit_to_byte(bit_str, payload, (size_t)total_bits);

    HuFF_Destroy();

    FILE* fp = fopen(out_name, "wb");
    if (!fp) {
        free(bit_str);
        free(payload);
        return -1;
    }

    fwrite(&data_len, sizeof(int), 1, fp);
    fwrite(&leaf_count, sizeof(int), 1, fp);
    fwrite(&struct_bytes, sizeof(int), 1, fp);

    // 符号表：按「叶子被前序遍历访问到的顺序」，不是按字节值排的。
    // 这个顺序必须和结构串里 '0' 出现的顺序一致，否则解压端挂叶子会错位。
    fwrite(symbols, 1, (size_t)leaf_count, fp);

    fwrite(struct_packed, 1, (size_t)struct_bytes, fp);
    fwrite(payload, 1, (size_t)payload_bytes, fp);

    long file_size = ftell(fp);
    printf("%d byte -> %ld byte (%.1f%%)\n", data_len, file_size, 100.0 * file_size / data_len);

    fclose(fp);
    free(bit_str);
    free(payload);
    return 0;
}

int decompress(const char* in_name, const char* out_name) {
    FILE* fp = fopen(in_name, "rb");
    if (!fp) return -1;

    int data_len, leaf_count, struct_bytes;
    if (fread(&data_len, sizeof(int), 1, fp) != 1 ||
        fread(&leaf_count, sizeof(int), 1, fp) != 1 ||
        fread(&struct_bytes, sizeof(int), 1, fp) != 1) {
        fclose(fp);
        return -1;
    }

    int struct_bits = 2 * leaf_count - 1; // 结构串有多少位

    // 1. 符号表：leaf_count 个单字节
    unsigned char symbols[MAX_CODE_LEN] = {0};
    for (int i = 0; i < leaf_count; i++) {
        if (fread(&symbols[i], 1, 1, fp) != 1) {
            fclose(fp);
            return -1;
        }
    }

    // 2. 结构串：读 struct_bytes 个字节，摊成 struct_bits 个 '0'/'1'
    unsigned char struct_packed[MAX_STRUCT_BYTES] = {0};
    char struct_str[MAX_STRUCT_STR];
    for (int i = 0; i < struct_bytes; i++) {
        if (fread(&struct_packed[i], 1, 1, fp) != 1) {
            fclose(fp);
            return -1;
        }
    }
    byte_to_bit(struct_packed, (unsigned char*)struct_str, struct_bits);

    // 3. 重建树。
    //    这里不调 HuFF_Code_Table（解码是沿树走，用不到码表），
    //    也不调 HuFF_Wpl（总位数是靠频次算的，解码端没有频次）。
    if (HuFF_Rebuild((unsigned char*)struct_str, symbols) != 0) {
        fclose(fp);
        return -1;
    }

    // 4. 载荷 = 文件剩下的全部字节。
    //    载荷是文件的最后一段，所以直接拿 EOF 兜底，不必再存一个长度字段。
    long payload_start = ftell(fp);
    fseek(fp, 0, SEEK_END);
    long payload_bytes = ftell(fp) - payload_start;
    fseek(fp, payload_start, SEEK_SET);

    unsigned char* payload = malloc((size_t)payload_bytes + 1);
    unsigned char* bit_str = malloc((size_t)payload_bytes * 8 + 1);
    unsigned char* out_buf = malloc((size_t)data_len + 1);
    if (!payload || !bit_str || !out_buf) {
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
    HuFF_Decode(bit_str, out_buf, (size_t)data_len);
    HuFF_Destroy();

    FILE* out = fopen(out_name, "wb");
    if (!out) {
        free(payload);
        free(bit_str);
        free(out_buf);
        return -1;
    }
    fwrite(out_buf, 1, (size_t)data_len, out);
    fclose(out);

    free(payload);
    free(bit_str);
    free(out_buf);
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
        snprintf(out, out_size, "%s", path);        // 装不下就原样返回
        return;
    }
    memcpy(out, path, base);
    memcpy(out + base, ext, ext_len + 1);           // 连结尾的 '\0' 一起拷
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

        unsigned char* data = malloc((size_t)size + 1);
        if (!data) { fclose(in); return 1; }
        if (size > 0 && fread(data, 1, (size_t)size, in) != (size_t)size) {
            fclose(in);
            free(data);
            return 1;
        }
        fclose(in);

        if (compress(out_name, data, (int)size) != 0) {
            printf("compress failed\n");
            free(data);
            return 1;
        }
        printf("-> %s\n", out_name);
        free(data);

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