#include "huffman_algo.h"
#include <stdlib.h>
#include <string.h>

// 树上的一个节点。
//   叶子：leaf = 1，c 是它代表的字节，weight 是频次（解码端建出来的叶子 weight = 0）
//   内部：leaf = 0，c 没有意义，weight 是左右孩子权之和
struct node {
    struct node *left, *right;
    unsigned short c;
    int weight;
    int leaf;
};

// 森林：建树过程中「还没被合并的树」都排在这里
Node* forest[ALPHABET_SIZE] = {0};

// 森林里有几棵树。建树时一路减少，最后剩 1 棵，那就是根
int count = 0;

// 码表：code_tab[字节] = 这个字节的 '0'/'1' 串
char code_tab[ALPHABET_SIZE][MAX_CODE_BITS];

// 后序释放一棵树（先放孩子，再放自己）
static void free_tree(Node* node) {
    if (!node) return;
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

// 从 root 往下走，把每个叶子的路径抄进码表。
// path 是「当前这条路径」，depth 是已经走了几位。
// 往左补 '0'，往右补 '1'。
static void fill_code_tab(Node* root, int depth, unsigned char* path) {
    if (root->leaf == 1) {
        path[depth] = '\0'; // 路径到此结束
        strcpy(code_tab[root->c], (char*)path);
        return;
    }
    path[depth] = '0'; // 这一步向左
    fill_code_tab(root->left, depth + 1, path);

    // 同一个格子改成 '1' 再走右边。
    // 左边递归已经跑完、码也抄走了，所以这里覆盖掉没有关系。
    path[depth] = '1';
    fill_code_tab(root->right, depth + 1, path);
}

// 累加 Σ(叶子频次 × 叶子深度)，走完整棵树就是编码后的总位数
static void accumulate_wpl(Node* root, int* total_bits, int depth) {
    if (root->leaf == 1) {
        (*total_bits) += root->weight * depth;
        return;
    }
    accumulate_wpl(root->left, total_bits, depth + 1);
    accumulate_wpl(root->right, total_bits, depth + 1);
}

// ==================== 频次统计 ====================
int counter(const unsigned short* data, size_t len, unsigned int* freq) {
    int leaf_count = 0;
    for (size_t i = 0; i < len; i++) {
        freq[(unsigned short)data[i]]++;
    }
    for (size_t i = 0; i < ALPHABET_SIZE; i++) {
        if (freq[i] > 0) leaf_count++;
    }
    return leaf_count;
}

// ==================== 建树 ====================
int HuFF_Init(const unsigned int* freq, size_t len) {
    if (len <= 0) return 1; // 空输入，连叶子都没有
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (freq[i] <= 0) continue; // 没出现过的字节不建叶子
        Node* node = malloc(sizeof(Node));
        if (!node) { // 半路没内存了，把已经建好的清掉
            for (int j = 0; j < count; j++) {
                free_tree(forest[j]);
                forest[j] = NULL;
            }
            count = 0;
            return -1;
        }
        node->c = (unsigned short)i;
        node->weight = freq[i];
        node->leaf = 1;
        node->left = node->right = NULL;
        forest[count] = node;
        count++;
    }
    return 0;
}

int HuFF_Build() {
    while (1) {
        if (count <= 1) return 0; // 只剩一棵了，它就是根
        Node* node = malloc(sizeof(Node));
        if (!node) return -1;

        // 找出权最小的两棵树：min_1 最小，min_2 次小
        int min_1 = -1, min_2 = -1;
        for (int i = 0; i < count; i++) {
            int weight = forest[i]->weight;
            if (min_1 < 0 || forest[min_1]->weight > weight) {
                min_2 = min_1;
                min_1 = i;
            } else if (min_2 < 0 || forest[min_2]->weight > weight) {
                min_2 = i;
            }
        }

        node->weight = forest[min_1]->weight + forest[min_2]->weight;
        node->c = 0;
        node->leaf = 0;
        node->left = forest[min_1];
        node->right = forest[min_2];

        // 把新树放回森林：占住下标小的那个位置；
        // 被合掉的另一棵如果不在末尾，就拿末尾那棵来填它的坑。
        int lo = min_1 < min_2 ? min_1 : min_2;
        int hi = min_1 < min_2 ? min_2 : min_1;
        forest[lo] = node;
        if (hi != count - 1) forest[hi] = forest[count - 1];
        forest[count - 1] = NULL;
        count--;
    }
    return 0;
}

void HuFF_Destroy() {
    for (int i = count - 1; i >= 0; i--) free_tree(forest[i]);
    memset(forest, 0, sizeof(forest));
    memset(code_tab, 0, sizeof(code_tab));
    count = 0;
}

// ==================== 编码 ====================
void HuFF_Wpl(int* total_bits) {
    accumulate_wpl(forest[0], total_bits, 0);
}

void HuFF_Code_Table() {
    unsigned char path[MAX_CODE_BITS] = {0};
    fill_code_tab(forest[0], 0, path);
}

void HuFF_Get(const unsigned short* data, unsigned char* out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        size_t code_len = strlen(code_tab[data[i]]);
        memcpy(out, code_tab[data[i]], code_len);
        out += code_len;
    }
}

// ==================== 解码 ====================
void HuFF_Decode(const unsigned char* bit_str, unsigned short* out, size_t len) {
    Node* root = forest[0];
    size_t bit_pos = 0;
    for (size_t i = 0; i < len; i++) {
        Node* p = root;
        while (!p->leaf) { // 还没落到叶子，就继续往下走
            if (bit_str[bit_pos] == '0')
                p = p->left;
            else
                p = p->right;
            bit_pos++;
        }
        *out++ = p->c; // 落到叶子，吐出一个字节
    }
}

// ==================== 结构串（树的形状）====================

// 前序遍历，把形状写进 struct_out：内部节点 '1'，叶子 '0'。
// 碰到叶子的时候，顺手把它的字节按「被访问到的顺序」追加进 symbols。
//
// 两个游标都用指针传，是为了让左右两棵子树共用同一个 ——
// 左子树写完之后游标停在它末尾，右子树从那儿接着写，两边不会互相覆盖。
static void dump_struct(Node* root, char* struct_out, int* pos, unsigned short* symbols, size_t* sym_pos) {
    if (root->leaf == 1) {
        struct_out[(*pos)++] = '0';
        symbols[(*sym_pos)++] = root->c;
        return;
    }
    struct_out[(*pos)++] = '1';
    dump_struct(root->left, struct_out, pos, symbols, sym_pos);
    dump_struct(root->right, struct_out, pos, symbols, sym_pos);
}

void HuFF_Struct(char* struct_out, int capacity, unsigned short* symbols) {
    int pos = 0;
    size_t sym_pos = 0;
    dump_struct(forest[0], struct_out, &pos, symbols, &sym_pos);
    struct_out[pos < capacity ? pos : capacity - 1] = '\0';
}

// ==================== 从结构串重建树（解码端）====================
static Node* rebuild(const unsigned char* struct_bits, size_t* bit_pos, const unsigned short* symbols, size_t* sym_pos) {
    unsigned char mark = struct_bits[(*bit_pos)++]; // 先吃掉当前这一位

    Node* node = malloc(sizeof(Node));
    if (!node) return NULL;
    node->left = node->right = NULL;
    // node->weight = 0; // 解码端没有频次，weight 没有意义

    if (mark == '0') { // 叶子：从符号表取下一个字节
        node->c = symbols[*sym_pos];
        (*sym_pos)++;
        node->leaf = 1;
    } else { // 内部节点：左右各递归一次
        node->leaf = 0;
        node->left = rebuild(struct_bits, bit_pos, symbols, sym_pos);
        node->right = rebuild(struct_bits, bit_pos, symbols, sym_pos);
    }
    return node;
}

int HuFF_Rebuild(const unsigned char* struct_bits, const unsigned short* symbols) {
    HuFF_Destroy();

    size_t bit_pos = 0; // 位游标：每读一位前进一格
    size_t sym_pos = 0; // 符号游标：只在遇到叶子时前进一格

    Node* root = rebuild(struct_bits, &bit_pos, symbols, &sym_pos);
    if (root == NULL) return -1;

    forest[0] = root;
    count = 1;
    return 0;
}