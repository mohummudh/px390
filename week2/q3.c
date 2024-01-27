#include <stdio.h>
#include <stdlib.h>

int  *cumsum(int *arr, int arr_length){
    int* cSum = (int*)malloc((arr_length+1) * sizeof(int));

    cSum[0] = 0;

    for(int i = 0; i < arr_length; i++) {
        cSum[i+1] = arr[i] + cSum[i];
    }

    return cSum;
}

int main(){
    int arr[] = {3,4,-6};
    int arr_length = 3;

    int* result = cumsum(arr, arr_length);

    for (int i = 0; i <= arr_length; i++) {
            printf("%d ", result[i]);
        }

    return 0;    
}