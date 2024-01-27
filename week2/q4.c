#include <stdio.h>
#include <stdlib.h>

long* mirror(long *arr, long n){
    long* twice = (long*)malloc((2*n) * sizeof(long));

    for(int i = 0; i < n; i++){
        twice[i] = arr[i];
    }

    for(int i = 0; i < n; i++){
        twice[i + n] = twice[n - i - 1];
    } 

    return twice;

}

int main(){
    long arr[] = {2,5,1};
    long n = 3;

    long* result = mirror(arr, n);

    for (int i = 0; i < 2*n; i++) {
            printf("%ld ", result[i]);
        }

    return 0;    
}