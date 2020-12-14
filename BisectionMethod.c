#include<stdio.h>

//function used is x^3-2x^2+3
double func(double x)
{
    return x*x*x - 2*x*x + 3;
}

double e=0.01;
double c;

void bisection(double a,double b)
{
    if(func(a) * func(b) >= 0)
    {
        printf("Incorrect a and b");
        return;
    }

    c = a;

    while ((b-a) >= e)
    {
        c = (a+b)/2;
        if (func(c) == 0.0){
            printf("Root = %lf\n",c);
            break;
        }
        else if (func(c)*func(a) < 0){
                printf("Root = %lf\n",c);
                b = c;
        }
        else{
                printf("Root = %lf\n",c);
                a = c;
        }
    }
}

int main()
{
    double a,b;
    a=-10;
    b=20;

    printf("The function used is x^3-2x^2+3\n");
    printf("a = %lf\n",a);
    printf("b = %lf\n",b);
    bisection(a,b);
    printf("\n");
    printf("Accurate Root calculated is = %lf\n",c);

    return 0;
}
