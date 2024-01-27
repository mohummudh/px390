#include <stdio.h>
#include <stdlib.h>

int main()
{
    long int A;
    long int B;
    long int result = 1;

    scanf("%ld", &A);
    scanf("%ld", &B);
    
    /*for loop to raise A to the power of B*/
    for(int i = 1; i <= B; i++) {
        result *= A;
    }

    /*printing final result*/
    printf("%ld", result);

    return 0;
}