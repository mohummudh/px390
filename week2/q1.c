#include <stdio.h>
#include <stdlib.h>

int *move_to_front(int *arr, int arr_length){
    int* copy = (int*)malloc(arr_length * sizeof(int));

    copy[0] = arr[arr_length - 1];

    for(int i = 0; i < arr_length - 1; i++) {
        copy[i+1] = arr[i];
    }

    return copy;

}

int main(){
    int arr[] = {3, 4, 1, -9, 7, 8};
    int arr_length = 6;

    int* result = move_to_front(arr, arr_length);

    for (int i = 0; i < arr_length; i++) {
            printf("%d ", result[i]);
        }

    return 0;    
}