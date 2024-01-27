#include <stdio.h>
#include <stdlib.h>

int num_conseq_digits(long k){
    int result = 1;

    long temp = k;
    int digits = 0;

    /*finding the length of the number*/

    while (temp != 0){
        temp /= 10;
        digits++;
    }

    int array[digits];

    temp = k;

    /*turning it into an array of numbers*/

    for(int i = digits - 1; i >= 0; i--){
        array[i] = temp % 10;
        /*printf("%d", array[i]);*/
        temp /= 10;
    }

    /*looping one by one to get the max length*/

    int current = 1;
    int max = 1;

    for(int a = 0; a < digits - 1; a++){
        if (array[a+1] == array[a]){
            current++;
        } else {
            if(current > max){
                max = current;
            } 
            current = 1;
        }
    }

    if (current > max){
        max = current;
    }
    
    return max;
}


int main()
{
    printf("%d", num_conseq_digits(67));
    return 0;
}