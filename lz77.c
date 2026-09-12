#include "lz77.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int lz77_encode(const unsigned char* data, long data_len, unsigned short* sym,
                unsigned short* dist, int* out_sym_n, int* out_dist_n) {

    int pos = 0, mlen = 0, num = 0, n_s = 0, n_d = 0;
    while (pos < data_len) {
        mlen = find_match(data, data_len, pos, 32767, &num);
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