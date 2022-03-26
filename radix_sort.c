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

struct tablo * allocate_tablo(int size) {
    struct tablo * tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->tab = malloc(size*sizeof(int));
    return tmp;
}

void free_tablo(struct tablo* tab_to_free){
    free(tab_to_free->tab);
    free(tab_to_free);
}

void print_array(struct tablo * tmp) {
    printf("---- Array of size %i ---- \n", tmp->size);
    int size = tmp->size;
    int i;
    for (i = 0; i < size; ++i) {
        printf("%i ", tmp->tab[i]);
    }
    printf("\n");
}

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
            printf("File too big : %d, max_size = %d, try with something smaller !\n Exiting\n",index,MAX_SIZE_NAME);
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


// Generate an 12 size array
void generate_example_array(struct tablo * s) {
    s->size=12;
    s->tab=malloc(s->size*sizeof(int));
    s->tab[0]=5;
    s->tab[1]=1021;
    s->tab[2]=2;
    s->tab[3]=9;
    s->tab[4]=0;
    s->tab[5]=23;
    s->tab[6]=9;
    s->tab[7]=512;
    s->tab[8]=511;
    s->tab[9]=8;
    s->tab[10]=9;
    s->tab[11]=5;
}

// Generate an 8 size array
void generate_example_8_array(struct tablo * s) {
    s->size=8;
    s->tab=malloc(s->size*sizeof(int));
    s->tab[0]=5;
    s->tab[1]=1021;
    s->tab[2]=2;
    s->tab[3]=9;
    s->tab[4]=0;
    s->tab[5]=23;
    s->tab[6]=9;
    s->tab[7]=512;
}

void generate_simple_13_array(struct tablo * s) {
    s->size=13;
    s->tab=malloc(s->size*sizeof(int));
    s->tab[0]=1;
    s->tab[1]=2;
    s->tab[2]=3;
    s->tab[3]=4;
    s->tab[4]=5;
    s->tab[5]=6;
    s->tab[6]=7;
    s->tab[7]=8;
    s->tab[8]=9;
    s->tab[9]=10;
    s->tab[10]=11;
    s->tab[11]=12;
    s->tab[12]=13;
}

// Generate an array with a size given and random value between 0 and max_value
void generate_array(struct tablo * s, int size, int max_value) {
    s->size=size;
    s->tab=malloc(s->size*sizeof(int));
    for(int i = 0; i < size; i++){
        s->tab[i] = rand() % (max_value + 1);
    }
}

int add(int a, int b){
    return a + b;
}

int my_pow(int a, int p){
    int res = a;
    for(int i = 1; i < p-1;i++){
        res *= a;
    }
    return res;
}

void copy_array(struct tablo* source, struct tablo* destination){
    int const size = source->size;
    for(int i = 0; i < source->size;i++){
        destination->tab[size+i] = source->tab[i];
    }
}

void copy_half_array(struct tablo* source, struct tablo* destination){
    int const size = destination->size;
    for(int i = 0; i < size;i++){
        destination->tab[i] = source->tab[size+i];
    }
}

void copy_array_same_size(struct tablo* source, struct tablo* destination){
    int const size = destination->size;
    for(int i = 0; i < size;i++){
        destination->tab[i] = source->tab[i];
    }
}

void reverse_bit(struct tablo* source){
    for(int i = 0; i < source->size; i++){
        source->tab[i] = !source->tab[i];
    }
}

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

void montee(struct tablo* A, int (*funct)(int, int)){
    int m =(int) (log(A->size/2) /log(2));
    for(int l = m-1; l >= 0; l--){
        #pragma omp parallel for
        for(int j = 1 << l; j < 1 << (l + 1);j++){
            A->tab[j] = funct(A->tab[2*j],A->tab[2*j+1]);
        }
    }
}

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

void scan_to_final(struct tablo* source, struct tablo* result){
    int size_min = source->size;
    #pragma omp parallel for
    for(int i = 0; i <= size_min;i++){
        result->tab[i] = result->tab[i] + source->tab[i];
    }
}

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

struct tablo* scan_prefix_final(struct tablo* source, int (*funct)(int, int), int eltneutre){
    int size = source->size;
    //Convertir la taille du tableau pour faire les sous-tableaux : 17 = 10001 (tableau de taille 2^4 et 2^0)
    int* bin_tab = digit_to_binary(size);
    //Nombre de bit max
    int size_bit_tab = (int) (log(size) /log(2))+1;
    struct tablo * result = allocate_tablo(size);
    int already_used = 0;
    int last_nb = 0;
    for(int i = 0;i<size_bit_tab;i++){
        //Pour chaque sous-tableau de taile 2^i
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

struct tablo* suffix_final(struct tablo* source, int (*funct)(int, int), int eltneutre){
    int size = source->size;
    //Convertir la taille du tableau pour faire les sous-tableaux : 17 = 10001 (tableau de taille 2^4 et 2^0)
    int* bin_tab = digit_to_binary(size);
    //Nombre de bit max
    int size_bit_tab = (int) (log(size) /log(2))+1;
    struct tablo * result = allocate_tablo(size);
    int already_used = 0;
    int last_nb = 0;
    for(int i = 0;i<size_bit_tab;i++){
        //Pour chaque sous-tableau de taile 2^i
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

struct tablo* size_suffix(int size, struct tablo* flags, int (*funct)(int, int), int eltneutre){
    struct tablo* suf = suffix_final(flags, funct, eltneutre);
    struct tablo* A = allocate_tablo(size);

    for(int i = 0; i < size; i++){
        A->tab[i] = size - suf->tab[i];
    }
    free_tablo(suf);
    return A;
}

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

void sort_by_digit(struct tablo* source, int nb_bit, int (*funct)(int, int), int eltneutre, int size){
    struct tablo* tab_bit = allocate_tablo(size);

    //Iteration sur les éléments du tableau
    for (int i = 0; i <= nb_bit; i++){
        //Itération sur les bits de chaque élément
        for (int j=0; j < source->size;j++) {
            // tab[j] va contenir 0 ou 1 (valeur
            tab_bit->tab[j] = (source->tab[j] >> i) & 1;
        }
        //On va trier les éléments du tableau par rapport au bit n°i
        split_bit(source,tab_bit,funct, eltneutre, size);
    }
    free_tablo(tab_bit);
}

int main(int argc, char **argv) {
    struct tablo* tmp = malloc(sizeof(struct tablo));
    struct tablo* copy = malloc(sizeof(struct tablo));

    // GESTION DES ARGUMENTS
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
        //ATTENTION : ne pas mettre de chiffre comme premier caractère du nom du fichier de sortie !!
        //Si c'est un nombre ==> on met à jour le nombre de threads
        if (argv[4][0] >= '0' && argv[4][0] <= '9'){
            output_file_name = argv[3];
            nb_threads=atoi(argv[4]);
        }
        //Sinon on met à jour le nom du fichier d'entrée
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

    //Pour enlever l'entrer
    remove(output_file_name);

    // REMPLISSAGE DU TABLEAU

    tmp->size = array_size;
    copy->size = array_size;
    tmp->tab = malloc(tmp->size * sizeof(int));
    copy->tab = malloc(copy->size * sizeof(int));
    //nb digit : (int)log10(1)+1
    //nb bit : (int) (log(185) /log(2))+1

    int nb_bit = (int) (log(N) /log(2))+1;
    int nb_digit_max = (int)log10(N)+1;

    #ifdef _OPENMP
    omp_set_num_threads(nb_threads);
    omp_set_nested(1); // autorise l'imbrication
    omp_set_max_active_levels(2); // d'un niveau max de 2
    #endif

    //printf("N = %i, array_size = %i, input = %s; output = %s, nb_threads = %i\n",N,array_size,input_file_name,output_file_name,nb_threads);

    if (is_input){
        generate_array_with_input(input_file_name,tmp);
    }
    else{
        generate_array(tmp,array_size,N);
    }

    copy_array_same_size(tmp,copy);

    //Time beginning
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //CHANGE THE FUNCTION "add" AND THE NEUTRAL ELEMENT HERE TO USE AN OTHER OPERATION LIKE "sub" OR "mult"
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
