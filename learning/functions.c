#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// have to give return type
void sayHi(char name[], int age){
    printf("Hello %s, you are %d\n", name, age);
}

int main()
{
    //main function is the one that gets executed when we run the file.
    printf("Top\n");
    sayHi("Mohammed", 20);
    return 0;
}

