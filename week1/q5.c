#include <stdio.h>
#include <stdlib.h>

long minbox(long area){
    int length = 1;
    int width = 1;
    int perimeter;
    while(length * width < area){
        if(width > length){
            length++;
        } else{
            width++;
        }
    }
    perimeter = 2 * (length + width);
    return perimeter;
}

int main()
{
    long a=4;
    long c;
    c = minbox(a);
    printf("%ld\n", c);
    return 0;
}