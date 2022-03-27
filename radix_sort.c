#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

const int MAX_SIZE_NAME = 1024;

struct tablo {
    int * tab;
    int size;
};

// Allocate an array of size size
struct tablo * allocate_tablo(int size) {
    struct tablo * tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->tab = malloc(size*sizeof(int));
    return tmp;
}

// Free the array tab_to_free
void free_tablo(struct tablo* tab_to_free){
    free(tab_to_free->tab);
    free(tab_to_free);
}

// Write in the array s all the number in the file file_name
void generate_array_with_input(const char* file_name, struct tablo * s){
    int fd = open(file_name, O_RDONLY);
    if (fd < 0){
        printf("Error opening the input file : %s\n Exiting\n",file_name);
        exit(1);
    }
    char c;
    char* str_nb = malloc(MAX_SIZE_NAME * sizeof(char));
    int index_tab = 0;
    int index = 0;
    int last_one = 0;
    while(read(fd, &c, 1)){
        if (index >= MAX_SIZE_NAME){
            printf("Number too big : %d, max_size = %d, try with something smaller !\n Exiting\n",index,MAX_SIZE_NAME);
            exit(1);
        }
        if (c != '\n' && c != ' '){
            str_nb[index] = c;
            last_one = 1;
            index++;
        }
        else {
            str_nb[index] = '\0';
            s->tab[index_tab] = atoi(str_nb);
            index_tab++;
            index = 0;
            last_one = 0;
        }
    }

    if(last_one){
        str_nb[index] = '\0';
        s->tab[index_tab] = atoi(str_nb);
    }
    close(fd);
    free(str_nb);
}

// Write in the file with the file descriptor fd all the number of the array s
void array_to_output(int fd, struct tablo * s, int nb_digits_max){
    char nb_str[nb_digits_max];
    for(int i = 0; i < s->size; i++){
        sprintf(nb_str, "%d ", s->tab[i]);
        write(fd,nb_str,strlen(nb_str));
    }
    nb_str[0] = '\n';
    nb_str[1] = '\0';
    write(fd,nb_str,strlen(nb_str));

}

// Write in the file file_name the first line (N, size, nb_threads, time), the input array and the output array
void generate_output(const char* file_name, int N, int array_size, int nb_threads, struct tablo* source, struct tablo* result, int time, int nb_digits_max){
    int fd = open(file_name, O_CREAT | O_TRUNC | O_WRONLY, S_IRWXU | S_IRWXG | S_IRWXO);
    if (fd < 0){
        printf("Error opening the output file : %s\n Exiting\n",file_name);
        exit(1);
    }
    char str[MAX_SIZE_NAME];
    sprintf(str,"N=%d Size=%d Threads=%d Time_in_microsec=%d\n",N,array_size,nb_threads,time);
    write(fd,str,strlen(str));
    array_to_output(fd, source, nb_digits_max);
    array_to_output(fd, result, nb_digits_max);
    close(fd);
}

// Generate an array with a size given and random value between 0 and max_value
void generate_array(struct tablo * s, int size, int max_value) {
    s->size=size;
    s->tab=malloc(s->size*sizeof(int));
    for(int i = 0; i < size; i++){
        s->tab[i] = rand() % (max_value + 1);
    }
}

// Function used in the prefix and suffix
int add(int a, int b){
    return a + b;
}

// Copy elements of an array source (size = N) in an array destination (size = 2 * N)
void copy_array(struct tablo* source, struct tablo* destination){
    int const size = source->size;
    for(int i = 0; i < source->size;i++){
        destination->tab[size+i] = source->tab[i];
    }
}

// Copy elements of an array source (size = 2 * N) in an array destination (size = N)
void copy_half_array(struct tablo* source, struct tablo* destination){
    int const size = destination->size;
    for(int i = 0; i < size;i++){
        destination->tab[i] = source->tab[size+i];
    }
}

// Copy elements of an array source (size = N) in an array destination (size = N)
void copy_array_same_size(struct tablo* source, struct tablo* destination){
    int const size = destination->size;
    for(int i = 0; i < size;i++){
        destination->tab[i] = source->tab[i];
    }
}

// Reverse all the bit in the array source
void reverse_bit(struct tablo* source){
    for(int i = 0; i < source->size; i++){
        source->tab[i] = !source->tab[i];
    }
}

// Permute the elements of the array source by changing there index by the new corresponding index in Index
void permute(struct tablo* source, struct tablo* Index){
    struct tablo* res = allocate_tablo(source->size);

    #pragma omp parallel for
    for(int i = 0; i < source->size; i++){
        res->tab[Index->tab[i]] = source->tab[i];
    }

    #pragma omp parallel for
    for(int i = 0; i < source->size; i++){
        source->tab[i] = res->tab[i];
    }
    free_tablo(res);
}

// Ascending phase for prefix and suffix for an array of size 2^n
void montee(struct tablo* A, int (*funct)(int, int)){
    int m =(int) (log(A->size/2) /log(2));
    for(int l = m-1; l >= 0; l--){
        #pragma omp parallel for
        for(int j = 1 << l; j < 1 << (l + 1);j++){
            A->tab[j] = funct(A->tab[2*j],A->tab[2*j+1]);
        }
    }
}

// Descending phase for prefix for an array of size 2^n
void descente_prefix(struct tablo* A, int (*funct)(int, int)){
    int m =(int) (log(A->size/2) /log(2));
    for(int l = 0; l < m; l++){
        const int step_max = 1 << (l + 1);
        #pragma omp parallel for
        for(int j = 1 << l; j < step_max; j++){
            A->tab[2*j+1] = funct(A->tab[j],A->tab[2*j]);
            A->tab[2*j] = A->tab[j];
        }
    }
}

// Descending phase for suffix for an array of size 2^n
void descente_suffix(struct tablo* B, int (*funct)(int, int)){
    int m = (int) log2(B->size/2) - 1;
    for(int l = 0; l<=m; l++) {
        const int step_max = 1 << (l + 1);
        #pragma omp parallel for
        for(int j = 1 << l; j < step_max; j++) {
            B->tab[2*j] = funct(B->tab[2*j + 1], B->tab[j]);
            B->tab[2*j + 1] = B->tab[j];
        }
    }
}

// Ascending and descending phase for prefix for an array of size 2^n
struct tablo* scan_prefix(struct tablo* source, int (*funct)(int, int), int eltneutre){
    struct tablo* A = allocate_tablo(source->size*2);
    struct tablo* result = allocate_tablo(source->size);

    copy_array(source,A);

    montee(A,funct);

    A->tab[1] = eltneutre;
    descente_prefix(A,funct);

    copy_half_array(A,result);
    free_tablo(A);
    return result;
}

// Ascending and descending phase for suffix for an array of size 2^n
struct tablo* scan_suffix(struct tablo* source, int (*funct)(int, int), int eltneutre){
    struct tablo* A = allocate_tablo(source->size*2);
    struct tablo* result = allocate_tablo(source->size);

    copy_array(source,A);

    montee(A,funct);

    A->tab[1] = eltneutre;
    descente_suffix(A,funct);

    copy_half_array(A,result);
    free_tablo(A);
    return result;
}

// Final step : scan to final result
void scan_to_final(struct tablo* source, struct tablo* result){
    int size_min = source->size;
    #pragma omp parallel for
    for(int i = 0; i <= size_min;i++){
        result->tab[i] = result->tab[i] + source->tab[i];
    }
}

// Return an array of the binary writing of a number in base 10
int* digit_to_binary(int digit) {
    int b;
    int k = 0;
    int nb_bit = (int) (log(digit) /log(2))+1;
    int *bits = malloc(nb_bit * sizeof(int));
    while (digit) {
        b = digit % 2;
        digit = digit / 2;
        bits[k] = b;
        k++;
    }
    return bits;
}

// Prefix scan for an array of any size
struct tablo* scan_prefix_final(struct tablo* source, int (*funct)(int, int), int eltneutre){
    int size = source->size;
    int* bin_tab = digit_to_binary(size);
    // Number max of bit
    int size_bit_tab = (int) (log(size) /log(2))+1;
    struct tablo * result = allocate_tablo(size);
    int already_used = 0;
    int last_nb = 0;
    for(int i = 0;i<size_bit_tab;i++){
        if(bin_tab[i]) {
            struct tablo *tmp = allocate_tablo(1 << i);

            #pragma omp parallel for
            for (int j = 0; j < tmp->size; j++) {
                tmp->tab[j] = source->tab[j + already_used];
            }

            struct tablo* B = scan_prefix(tmp, funct, eltneutre);

            #pragma omp parallel for
            for (int j = 0; j < tmp->size; j++) {
                result->tab[j + already_used] = last_nb + B->tab[j];
            }

            already_used += 1 << i;
            last_nb += B->tab[B->size-1] + source->tab[already_used-1];
            free_tablo(tmp);
            free_tablo(B);
        }
    }
    free(bin_tab);

    return result;
}

// Suffix for an array of any size
struct tablo* suffix_final(struct tablo* source, int (*funct)(int, int), int eltneutre){
    int size = source->size;
    int* bin_tab = digit_to_binary(size);
    // Number max of bit
    int size_bit_tab = (int) (log(size) /log(2))+1;
    struct tablo * result = allocate_tablo(size);
    int already_used = 0;
    int last_nb = 0;
    for(int i = 0;i<size_bit_tab;i++){
        if(bin_tab[i]) {
            struct tablo *tmp = allocate_tablo(1 << i);

            #pragma omp parallel for
            for (int j = 0; j < tmp->size; j++) {
                tmp->tab[j] = source->tab[size + j - already_used - tmp->size];
            }

            struct tablo* B = scan_suffix(tmp, funct, eltneutre);


            scan_to_final(tmp,B);


            #pragma omp parallel for
            for (int j = 0; j < tmp->size; j++) {
                result->tab[result->size -1 -already_used - j] = last_nb + B->tab[tmp->size - 1-j];
            }


            already_used += 1 << i;
            last_nb += B->tab[0];
            free_tablo(tmp);
            free_tablo(B);
        }
    }
    free(bin_tab);
    return result;
}

// Array where the i element is size - Suffix[i]
struct tablo* size_suffix(int size, struct tablo* flags, int (*funct)(int, int), int eltneutre){
    struct tablo* suf = suffix_final(flags, funct, eltneutre);
    struct tablo* A = allocate_tablo(size);

    for(int i = 0; i < size; i++){
        A->tab[i] = size - suf->tab[i];
    }
    free_tablo(suf);
    return A;
}

// Generate the new index array with the array of the bit number i for each element and permute
void split_bit( struct tablo* source, struct tablo* flags, int (*funct)(int, int), int eltneutre, int size){

    reverse_bit(flags);
    struct tablo* Idown = scan_prefix_final(flags,funct,eltneutre);

    reverse_bit(flags);
    struct tablo* Iup = size_suffix(size, flags, funct, eltneutre);


    struct tablo* Index = allocate_tablo(size);

    #pragma omp parallel for
    for (int i = 0; i < size; i++){
        if (flags->tab[i]){
            Index->tab[i] = Iup->tab[i];
        }
        else {
            Index->tab[i] = Idown->tab[i];
        }
    }

    permute(source, Index);

    free_tablo(Idown);
    free_tablo(Iup);
    free_tablo(Index);
}

// Sort the array source bit by bit
void sort_by_digit(struct tablo* source, int nb_bit, int (*funct)(int, int), int eltneutre, int size){
    struct tablo* tab_bit = allocate_tablo(size);

    // Iteration on all the bit
    for (int i = 0; i <= nb_bit; i++){
        // Iteration on each element
        for (int j=0; j < source->size;j++) {
            tab_bit->tab[j] = (source->tab[j] >> i) & 1;
        }
        // Sort element in relation to the bit i
        split_bit(source,tab_bit,funct, eltneutre, size);
    }
    free_tablo(tab_bit);
}

int main(int argc, char **argv) {
    struct tablo* tmp = malloc(sizeof(struct tablo));
    struct tablo* copy = malloc(sizeof(struct tablo));

    // Argument management
    if (argc<4) {
        fprintf(stdout, "Usage: %s N array_size [input_file_name] output_file_name [number_of_threads, otherwise the by default number] \n", argv[0]);
        return -1;
    }
    short nb_threads=1;
    int N = atoi(argv[1]);
    int array_size = atoi(argv[2]);
    int is_input = 0;
    char* input_file_name = "";
    char* output_file_name;
    if (argc==5){
        // CAUTION: do not put a number as the first character of the output file name!!
        // If it is a number, update the number of threads
        if (argv[4][0] >= '0' && argv[4][0] <= '9'){
            output_file_name = argv[3];
            nb_threads=atoi(argv[4]);
        }
        // Otherwise we update the name of the input file
        else {
            input_file_name = argv[3];
            output_file_name = argv[4];
            is_input = 1;
        }
    }
    else if (argc>=6){
        is_input = 1;
        input_file_name = argv[3];
        output_file_name = argv[4];
        nb_threads=atoi(argv[5]);
    }
    else {
        output_file_name = argv[3];
    }

    // To remove the output file if it already exist
    remove(output_file_name);

    tmp->size = array_size;
    copy->size = array_size;
    tmp->tab = malloc(tmp->size * sizeof(int));
    copy->tab = malloc(copy->size * sizeof(int));

    int nb_bit = (int) (log(N) /log(2))+1;
    int nb_digit_max = (int)log10(N)+1;

    #ifdef _OPENMP
    // Define the number of threads we will use
    omp_set_num_threads(nb_threads);
    // We don't need nested parallelism in this program
    omp_set_nested(0);
    // Allow parallelism
    omp_set_max_active_levels(1);
    #endif

    // Array filling
    if (is_input){
        generate_array_with_input(input_file_name,tmp);
    }
    else{
        generate_array(tmp,array_size,N);
    }

    copy_array_same_size(tmp,copy);

    // Time beginning
    struct timeval start, end;
    gettimeofday(&start, NULL);

    sort_by_digit(tmp, nb_bit, add, 0, array_size);

    //Time end
    gettimeofday(&end, NULL);

    generate_output(
            output_file_name,
            N,
            array_size,
            nb_threads,
            copy,
            tmp,
            (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec),
            nb_digit_max
            );

    free_tablo(tmp);
    free_tablo(copy);
}
