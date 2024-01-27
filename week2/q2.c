#include <stdio.h>
#include <stdlib.h>

int arr_dot_product(int *arr, int *arr2, int len){
    int dot = 0;
    for(int i = 0; i < len; i++){
        dot += arr[i] * arr2[i];
    }

    return dot;
}

int main(){
    int arr[] = {1,3};
    int arr2[] = {2,4};
    int len = 2;

    int result = arr_dot_product(arr, arr2, len);

    printf("%d", result);

    return 0;    
}