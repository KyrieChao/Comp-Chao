#ifndef LZ77_H
#define LZ77_H

#define MAX_MATCH 258

int lz77_encode(const unsigned char* data, long data_len,
                unsigned short* sym, unsigned short* dist,
                int* out_sym_n, int* out_dist_n);

int lz77_decode(const unsigned short* sym, int sym_n, const unsigned short* dist,
                int dist_n, unsigned char* out, long* out_len);
#endif