#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

int main(int argc, char* argv[])
{
    if((argc < 6) || (strcmp(argv[1],"-h") == 0))
    {
        printf("usage: [x_1] [y_1] [x_2] [y_2] [x]\n");
        printf("output: y(x) as linearly interpolated between provided points\n");
        return 0;
    }
    double x_1 = atof(argv[1]);
    double y_1 = atof(argv[2]);

    double x_2 = atof(argv[3]);
    double y_2 = atof(argv[4]);

    double x = atof(argv[5]);

    double y = (y_2 - y_1)/(x_2 - x_1)*(x - x_1) + y_1;
    
    if(argc == 7)
    {
        printf("%.*e %s\n", DECIMAL_DIG, y, argv[6]);
    }else
    {
        printf("%.*e\n", DECIMAL_DIG, y);
    }
    return 0;
}
