#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    char name[20];
    printf("Enter your name:");
    fgets(name, 20, stdin);
    printf("You name is %s!", name);

    return 0;
}