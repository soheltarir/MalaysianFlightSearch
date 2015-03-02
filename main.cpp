#include "my_class.h"
#include <iostream>
#include <vector>

int main()
{
    vector <CartesianPts> vectPts;
    CartesianPts tmp;
    int n, x, y;
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        cin >> x >> y;
        tmp.X_ = x; tmp.Y_ = y;
        vectPts.push_back(tmp);
    }
    Grid my_grid(vectPts);
    my_grid.PrintGridPts();
    return 0;
}
