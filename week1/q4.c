#include <stdio.h>
#include <stdlib.h>

int isprime(long k){
    /*remainder = dividend % divisor*/
    
    if(k <= 1){
        return 0;
    }

    /*found that you only need to check upto k's 2-root*/
    for(int i = 2; i*i <= k; i++){
        if(k % i == 0){
            return 0;
        }
    }
    return 1;
}

int main()
{
    printf("%d", isprime(6));
    return 0;
}