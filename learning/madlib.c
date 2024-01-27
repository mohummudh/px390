#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    char colour[10];
    char pluralNoun[20];
    char celebrityF[20];
    char celebrityL[20];

    printf("Enter a colour: ");
    scanf("%s", colour);
    printf("Enter a plural noun: ");
    scanf("%s", pluralNoun);
    printf("Enter a celebrity name: ");
    scanf("%s%s", celebrityF, celebrityL);


    printf("Roses are %s\n", colour);
    printf("%s are blue\n", pluralNoun);
    printf("I love %s %s\n", celebrityF, celebrityL);

    return 0;
}