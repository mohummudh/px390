#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    /*arrays*/
    int luckyNumbers[] = {4, 8, 15, 16, 200, 2, 4, 6};
    luckyNumbers[3] = 7;
    printf("%d\n", luckyNumbers[3]);

    int newLuck[10];
    newLuck[6] = 80;
    printf("%d\n", newLuck[6]);
    return 0;
}