#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <string.h>

typedef int Element_Type;

void gene_nums(Element_Type *a, int N)
{
    srand(time(NULL));
    for (int i = 0; i < N; i++)
        *(a + i) = rand();
}

void disp_nums(Element_Type *a, int N)
{
    for (int i = 0; i < N; i++)
        printf("%-7d", *(a + i));
    printf("\n");
}

void write2file(char a[], Element_Type A[], int N)
{
    FILE *fp = NULL;
    fp = fopen(a, "w+");
    for (int i = 0; i < N; i++)
        fprintf(fp, "%d\n", A[i]);
    fclose(fp);
}

//selection sort
void Selection_sort(Element_Type A[], int N)
{
    Element_Type min_val;
    int min_index;
    for (int i = 0; i < N; i++)
    {
        min_val = A[i];
        min_index = i;
        for (int j = i; j < N; j++)
        {
            if (A[j] < min_val)
            {
                min_val = A[j];
                min_index = j;
            }
        }
        A[min_index] = A[i];
        A[i] = min_val;
    }
}

//insertion sort
void Insertion_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    int j;
    for (int i = 0; i < N; i++)
    {
        tmp = A[i];
        for (j = i; j > 0 && A[j - 1] > tmp; j--)
            A[j] = A[j - 1];
        A[j] = tmp;
    }
}

//quick sort
void Quick_sort(Element_Type A[], int left, int right)
{
    int i, j;
    Element_Type base, tmp;
    base = A[left];
    if (left >= right)
        return;
    i = left;
    j = right;
    while (i < j)
    {
        while (A[j] >= base && i < j)
            j--;
        while (A[i] <= base && i < j)
            i++;
        if (i < j)
        {
            tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;
        }
    }
    A[left] = A[j];
    A[j] = base;
    Quick_sort(A, left, i - 1);
    Quick_sort(A, i + 1, right);
}

//merge sort
void merge(Element_Type A[], int left, int mid, int right)
{
    int x = mid - left + 1, y = right - mid;
    int i = 0, j = 0, k = left;
    Element_Type *B = (Element_Type *)malloc(x * sizeof(Element_Type));
    Element_Type *C = (Element_Type *)malloc(y * sizeof(Element_Type));

    for (; i < x; i++)
        B[i] = A[left + i];
    for (; j < y; j++)
        C[j] = A[mid + j + 1];
    i = j = 0;
    while (i < x && j < y)
    {
        if (B[i] <= C[j])
            A[k++] = B[i++];
        else
            A[k++] = C[j++];
    }
    while (i < x)
        A[k++] = B[i++];
    while (j < y)
        A[k++] = C[j++];

    free(B);
    free(C);
}
void Merge_sort(Element_Type A[], int left, int right)
{
    int mid;
    if (left < right)
    {

        mid = (left + right) / 2;
        Merge_sort(A, left, mid);
        Merge_sort(A, mid + 1, right);
        merge(A, left, mid, right);
    }
}

//heap sort
void swap(Element_Type *a, Element_Type *b)
{
    Element_Type tmp = *a;
    *a = *b;
    *b = tmp;
}

void max_heap(Element_Type A[], int i, int N)
{
    int parent;

    for (int j = i + N - 1; j > i; j--)
    {
        parent = j / 2 - 1;
        if (A[parent] < A[j])
            swap(&A[parent], &A[j]);
    }
}

void Heap_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    int i;

    max_heap(A, 0, N);

    for (int i = N - 1; i > 0; i--)
    {
        swap(&A[0], &A[i]);
        max_heap(A, 0, i);
    }
}

//radix sort(LSD)
Element_Type Get_max(Element_Type A[], int N)
{
    Element_Type tmp = A[0];
    for (int i = 0; i < N; i++)
        if (A[i] > tmp)
            tmp = A[i];

    return tmp;
}
struct id_val
{
    Element_Type a;
    int index;
    struct id_val *next;
};
int Getmodel(Element_Type val, int n)
{
    while (--n)
    {
        val /= 10;
    }
    return val % 10;
}
void divide(Element_Type A[], struct id_val *P, int N, int base)
{
    struct id_val *curr[10], *node;
    int flag;
    for (int i = 0; i < 10; i++)
        curr[i] = &P[i];

    for (int i = 0; i < N; i++)
    {
        node = (struct id_val *)malloc(sizeof(struct id_val));
        node->a = A[i];
        node->index = i;
        node->next = NULL;

        flag = Getmodel(A[i], base);
        curr[flag]->next = node;
        curr[flag] = node;
    }
}
void conquer(Element_Type A[], struct id_val *P)
{
    struct id_val *node;
    for (int i = 0, j = 0; i < 10; i++)
    {
        node = &P[i];
        while (node->next != NULL)
        {
            node = node->next;
            A[j] = node->a;
            j++;
        }
    }
}

void Radix_sort(Element_Type A[], int N)
{
    struct id_val buckets[10];
    Element_Type max_val;
    int digits = 0, base = 0;

    max_val = Get_max(A, N);
    while (max_val != 0)
    {
        digits++;
        max_val /= 10;
    }

    for (int i = 0; i < digits; i++)
    {
        base += 1;
        divide(A, buckets, N, base);
        conquer(A, buckets);
    }
}

//comb sort
void Comb_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    for (int step = N / 1.3; step > 0; step /= 1.3)
    {
        for (int i = 0; i + step < N; i++)
            if (A[i] > A[i + step])
            {
                tmp = A[i];
                A[i] = A[i + step];
                A[i + step] = tmp;
            }
    }
}

//shell sort
void Shell_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    int i, P, D, j;
    /*希尔增量序列 1-> Hibbard; 2-> Sedgewick*/
    do
    {
        j++;
        D = 9 * (int)(pow(4, j)) - 9 * (int)(pow(2, j)) + 1;
    } while (D < N);

    for (j--; j >= 0; j--)
    {
        D = 9 * (int)(pow(4, j)) - 9 * (int)(pow(2, j)) + 1;
        for (P = D; P < N; P++) /*插入排序*/
        {
            tmp = A[P];
            for (i = P; i >= D && A[i - D] > tmp; i -= D)
                A[i] = A[i - D];
            A[i] = tmp;
        }
    }
}

//bubble sort
void Bubble_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    for (int i = 0; i < N - 1; i++)
        for (int j = 0; j < N-1-i; j++)
            if (A[j - 1] > A[j])
            {
                tmp = A[j - 1];
                A[j - 1] = A[j];
                A[j] = tmp;
            }
}

//cocktail sort
void Cocktail_sort(Element_Type A[], int left, int right)
{
    Element_Type tmp;
    while (left < right)
    {
        for (int i = left; i < right; i++)
        {
            if (A[i] > A[i + 1])
            {
                tmp = A[i + 1];
                A[i + 1] = A[i];
                A[i] = tmp;
            }
        }
        right--;
        for (int i = right; i > left; i--)
        {
            if (A[i] < A[i - 1])
            {
                tmp = A[i];
                A[i] = A[i - 1];
                A[i - 1] = tmp;
            }
        }
        left++;
    }
}
//odd-even sort
void Odd_even_sort(Element_Type A[], int N)
{
    int change_flag = 1, odd_flag = 1;
    Element_Type tmp;

    while (change_flag)
    {
        change_flag = 0;

        for (int i = 1 - odd_flag; i + 1 < N; i += 2)
            if (A[i] > A[i + 1])
            {
                tmp = A[i];
                A[i] = A[i + 1];
                A[i + 1] = tmp;
                change_flag = 1;
            }

        odd_flag = 1 - odd_flag;
    }
}

//gnome sort
void Gnome_sort(Element_Type A[], int N)
{
    Element_Type tmp;
    int i = 0;

    while (i < N)
    {
        if (i == 0 || A[i - 1] <= A[i])
            i++;
        else
        {
            tmp = A[i];
            A[i] = A[i - 1];
            A[--i] = tmp;
        }
    }
}

//cycle sort
void Cycle_sort(Element_Type A[], int N)
{
    Element_Type item, tmp;
    int pos;

    for (int i = 0; i < N - 1; i++)
    {
        item = A[i];
        pos = i;

        for (int j = i + 1; j < N; j++)
            if (A[j] < item)
                pos += 1;

        if (pos == i)
            continue;

        while (item == A[pos])
            pos += 1;

        tmp = A[pos];
        A[pos] = item;
        item = tmp;

        while (pos != i)
        {
            pos = i;
            for (int j = i + 1; j < N; j++)
                if (A[j] < item)
                    pos += 1;
            while (item == A[pos])
                pos += 1;
            tmp = A[pos];
            A[pos] = item;
            item = tmp;
        }
    }
}


int main(void)
{
    int amount = 10000;
    clock_t start, end;
    Element_Type *A = (Element_Type *)malloc(amount * sizeof(Element_Type));
    Element_Type *B = (Element_Type *)malloc(amount * sizeof(Element_Type));

    gene_nums(A, amount);
    printf("原始数据(%d个)\n\n", amount);
    // disp_nums(A,amount);

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Shell_sort(B, amount);
    end = clock();
    printf("希尔排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Selection_sort(B, amount);
    end = clock();
    printf("选择排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Insertion_sort(B, amount);
    end = clock();
    printf("插入排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Quick_sort(B, 0, amount - 1);
    end = clock();
    printf("快速排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Merge_sort(B, 0, amount - 1);
    end = clock();
    printf("归并排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Bubble_sort(B, amount);
    end = clock();
    printf("冒泡排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Gnome_sort(B, amount);
    end = clock();
    printf("地精排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Cocktail_sort(B, 0, amount - 1);
    end = clock();
    printf("鸡尾酒排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Heap_sort(B, amount);
    end = clock();
    printf("堆排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Comb_sort(B, amount);
    end = clock();
    printf("梳排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Odd_even_sort(B, amount);
    end = clock();
    printf("奇偶排序耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Radix_sort(B, amount);
    end = clock();
    printf("基数排序(LSD)耗时 %lf ms\n\n", difftime(end, start));

    memcpy(B, A, sizeof(Element_Type) * amount);
    start = clock();
    Cycle_sort(B, amount);
    end = clock();
    printf("圈排序(LSD)耗时 %lf ms\n\n", difftime(end, start));

    // write2file("d:/2.txt", B, amount);

    free(B);
    free(A);
    return 0;
}
