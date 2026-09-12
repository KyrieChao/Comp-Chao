#include <stdio.h>

int find_match(const unsigned char* buf, int buf_len, int pos, int window, int* out_dist) {
    int start = pos - window < 0 ? 0 : (pos - window);
    int best_len = 0;
    int best_dist = 0;
    for (int i = start; i < pos; i++) {
        int len = 0;
        while (pos + len < buf_len && buf[i + len] == buf[pos + len]) len++;
        if (len >= best_len) {
            best_len = len;
            best_dist = pos - i;
        }
    }
    *out_dist = best_dist;
    return best_len;
}

int main(void) {
    int num = 0;
    int len = find_match((const unsigned char*)"ABABABABC", 9, 2, 2, &num);
    printf("len=%d, num=%d\n", len, num);
    return 0;
}