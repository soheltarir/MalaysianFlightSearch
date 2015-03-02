#ifndef _MY_CLASS_H_
#define _MY_CLASS_H_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;
#define MAX_POINTS 1000000
struct CartesianPts {
    int X_;
    int Y_;
};

struct NameEquals {
    bool operator() (const CartesianPts & rhs) const
    {
        return (rhs.X_ == name.X_ && rhs.Y_ == name.Y_ );
    }
    NameEquals(CartesianPts n) : name(n) {}

    private:
    CartesianPts name;
};
class valueComp {
    public:
        bool operator()(const CartesianPts& A, const CartesianPts& B) const
        {
            return ((!(A.X_<B.X_) && !(B.X_<A.X_)) && (!(A.Y_<B.Y_) && !(B.Y_<A.Y_))) ;
        }
};
/*****Quick Sort Implementation******/
void swap(int &a, int &b)
{
    int temp;
    temp = a;
    a = b;
    b = temp;
}
int SplitArray(int* array, int pivot, int startIndex, int endIndex)
{
    int leftBoundary = startIndex;
    int rightBoundary = endIndex;
    while(leftBoundary < rightBoundary)
    {
        while(pivot < array[rightBoundary] && rightBoundary > leftBoundary)
        {
            rightBoundary--;
        }
        swap(array[leftBoundary], array[rightBoundary]);
        while(pivot >= array[leftBoundary] && leftBoundary < rightBoundary)
        {
            leftBoundary++;
        }
        swap(array[rightBoundary], array[leftBoundary]);
    }
    return leftBoundary;
}
void QuickSort(int* array, int startIndex, int endIndex)
{
    int pivot = array[startIndex];
    int splitPoint;
    if (endIndex > startIndex)
    {
        splitPoint = SplitArray(array, pivot, startIndex, endIndex);
        array[splitPoint] = pivot;
        QuickSort(array, startIndex, splitPoint-1);
        QuickSort(array, splitPoint+1, endIndex);
    }
}
/******************************************/

class Plane {
    public:
        Plane(int id_, int range_, int economy_)
        {
            m_Id = id_;
            m_Range = range_;
            m_Economy = economy_;
        }
        int Id()    { return m_Id; }
        int Range() { return m_Range; }
        int Economy()  { return m_Economy; }
        int Cost(int x1_, int y1_, int x2_ = 0, int y2_ = 0)
        {
            return (m_Economy*(abs(x1_ - x2_) + abs(y1_ - y2_)));
        }
    private:
        int m_Id;
        int m_Range;
        int m_Economy;
};

class Grid {
    public:
        Grid(vector <CartesianPts> vectPts_)
        {
            vector <CartesianPts>::iterator it;
            m_vectBoundary = vectPts_;
            SortByX(m_vectBoundary);
            InitXY();
            InitGrid();
        }
        void InitGrid()
        {
            int Xmax = 0, Ymax = 0, Xmid = 0, Xmin = 0, Ymid, Ymin, y = 0;
            rectGrid(Xmax, Xmin, Ymax, Ymin, Xmid, Ymid);
            int intPts;
            for (int i = 0; i < m_vectBoundary.size(); i++)
            {
                CartesianPts Pts;
                Pts.Y_ = m_vectBoundary.at(i).Y_;
                if (!checkNextBoundary(true, m_vectBoundary.at(i).X_, m_vectBoundary.at(i).Y_))
                    continue;

                if (m_vectBoundary.at(i).X_ < Xmid)
                {
                    for (int x = m_vectBoundary.at(i).X_+1; x < Xmax; x++)
                    {
                        vector <CartesianPts>::iterator it;
                        vector <int>::iterator it2;
                        Pts.X_ = x;
                        intPts = x*10 + Pts.Y_; 
                        it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                        it2 = find(m_vectGridX.begin(), m_vectGridX.end(), intPts);
                        if (it != m_vectBoundary.end()) //&& x == m_vectBoundary.at(i).X_)
                        {
                            x = Xmax;
                        }
                        else if (it2 != m_vectGridX.end())
                        {
                            continue;
                        }
                        else
                        {
                            m_vectGridX.push_back(intPts);
                        }
                    }
                }
                else if (m_vectBoundary.at(i).X_ > Xmid)
                {
                    for (int x = m_vectBoundary.at(i).X_-1; x > Xmin; x--)
                    {
                        vector <CartesianPts>::iterator it;
                        vector <int>::iterator it2;
                        Pts.X_ = x;
                        intPts = x*10+Pts.Y_;
                        it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                        it2 = find(m_vectGridX.begin(), m_vectGridX.end(), intPts);
                        if ((it != m_vectBoundary.end())) //&& x == m_vectBoundary.at(i).X_))
                        {
                            x = Xmin;
                        }
                        else if (it2 != m_vectGridX.end())
                        {
                            continue;
                        }
                        else
                        {
                            m_vectGridX.push_back(intPts);
                        }
                    }
                }
            }
            for (int i = 0; i < m_vectBoundary.size(); i++)
            {
                CartesianPts Pts;
                Pts.X_ = m_vectBoundary.at(i).X_;
                if (!checkNextBoundary(false, m_vectBoundary.at(i).X_, m_vectBoundary.at(i).Y_))
                    continue;

                if (m_vectBoundary.at(i).Y_ < Ymid)
                {
                    for (int y = m_vectBoundary.at(i).Y_ + 1; y < Ymax; y++)
                    {
                        vector <CartesianPts>::iterator it;
                        vector <int>::iterator it2;
                        Pts.Y_ = y;
                        intPts = Pts.X_*10 + Pts.Y_;
                        it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                        it2 = find(m_vectGridY.begin(), m_vectGridY.end(), intPts);
                        if (it != m_vectBoundary.end()) //&& x == m_vectBoundary.at(i).X_)
                        {
                            y = Ymax;
                        }
                        else if (it2 != m_vectGridY.end())
                        {
                            continue;
                        }
                        else
                        {
                            m_vectGridY.push_back(intPts);
                        }
                    }
                }
                else if (m_vectBoundary.at(i).Y_ > Ymid)
                {
                    for (int y = m_vectBoundary.at(i).Y_-1; y > Ymin; y--)
                    {
                        vector <CartesianPts>::iterator it;
                        vector <int>::iterator it2;
                        Pts.Y_ = y;
                        intPts = Pts.X_*10 + Pts.Y_;
                        it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                        it2 = find(m_vectGridY.begin(), m_vectGridY.end(), intPts);
                        if ((it != m_vectBoundary.end())) //&& x == m_vectBoundary.at(i).X_))
                        {
                            y = Ymin;
                        }
                        else if (it2 != m_vectGridY.end())
                        {
                            continue;
                        }
                        else
                        {
                            m_vectGridY.push_back(intPts);
                        }
                    }
                }
            }
            sort(m_vectGridX.begin(), m_vectGridX.end());
            sort(m_vectGridY.begin(), m_vectGridY.end());
            set_intersection(m_vectGridX.begin(), m_vectGridX.end(), m_vectGridY.begin(), m_vectGridY.end(), back_inserter(m_vectGrid));
            sort(m_vectGrid.begin(), m_vectGrid.end());
        }
        bool checkNextBoundary(bool isAlongX, int Xstart, int Ystart)
        {
            int Xmax, Xmin, Ymax, Ymin, Xmid, Ymid;
            rectGrid(Xmax, Xmin, Ymax, Ymin, Xmid, Ymid);
            vector <CartesianPts>::iterator it;
            CartesianPts Pts; bool bReturn = true;
            if (isAlongX && Xstart < Xmid)
            {
                for (int x = Xstart; x <= Xmax; x++)
                {
                    Pts.X_ = x; Pts.Y_ = Ystart;
                    it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                    if (it == m_vectBoundary.end() && x == Xmax)    { bReturn = false;}
                }
            }
            else if (isAlongX && Xstart > Xmid)
            {
                for (int x = Xstart; x >= Xmin; x--)
                {
                    Pts.X_ = x; Pts.Y_ = Ystart;
                    it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                    if (it == m_vectBoundary.end() && x == Xmin)    { bReturn = false;}
                }
            }
            else if (!isAlongX && Ystart < Ymid)
            {
                for (int y = Ystart; y <= Ymax; y++)
                {
                    Pts.X_ = Xstart; Pts.Y_ = y;
                    it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                    if (it == m_vectBoundary.end() && y == Ymax)    { bReturn = false;}
                }
            }
            else if (!isAlongX && Ystart > Ymid)
            {
                for (int y = Ystart; y >= Ymin; y--)
                {
                    Pts.X_ = Xstart; Pts.Y_ = y;
                    it = find_if(m_vectBoundary.begin(), m_vectBoundary.end(), NameEquals(Pts));
                    if (it == m_vectBoundary.end() && y == Ymin)    { bReturn = false;}
                }
            }
            return bReturn;
        }
        void rectGrid(int &Xmax, int &Xmin, int &Ymax, int &Ymin, int &Xmid, int &Ymid)
        {
            CartesianPts V1, V2, V3, V4; //Circumscribing Rectangle;
            /*****Initialise Circumscribing Rectangle*****
             *   V4*********V3
             *   *           *
             *   *           *
             *   *           *
             *   V1*********V2
             */
            V1.X_ = m_vectX.front(); 
            V2.X_ = m_vectX.back(); 
            V1.Y_ = m_vectY.front();
            V2.Y_ = m_vectY.back();
            Xmid = (int)((V1.X_ + V2.X_)/2);
            Xmax = V2.X_; Xmin = V1.X_;
            Ymax = V2.Y_; Ymin = V1.Y_;
            Ymid = (int)((Ymax + Ymin)/2);
            /***********************************************/
            CartesianPts Pts;
            /*for (int x = V1.X_+1; x < V2.X_; x++)
            {
                for (int y = V1.Y_+1; y < V2.Y_; y++)
                {
                    Pts.X_ = x; Pts.Y_ = y;
                    m_vectGrid.push_back(Pts);
                }
            }
            SortByX(m_vectGrid);*/
        }
        void InitXY()
        {
            for (int i = 0; i < m_vectBoundary.size(); i++)
            {
                m_vectX.push_back(m_vectBoundary.at(i).X_);
                m_vectY.push_back(m_vectBoundary.at(i).Y_);
            }
            Sort(m_vectX); Sort(m_vectY);
        }
        void PrintGridPts()
        {
            for (int i = 0; i < m_vectGridX.size(); i++)
                cout << "X:" << (m_vectGridX.at(i))/10 << " " << "Y:" << (m_vectGridX.at(i))%10 << endl;
            cout << "--------------" << endl;
            for (int i = 0; i < m_vectGridY.size(); i++)
                cout << "X:" << (m_vectGridY.at(i))/10 << " " << "Y:" << (m_vectGridY.at(i))%10 << endl;
        }
    protected:
        static bool compareCoord(const CartesianPts &a, const CartesianPts &b)
        {
            if (a.X_ == b.X_)
                return a.Y_ < b.Y_;
            else
                return a.X_ < b.X_;
        }
        static bool compareCoordY(const CartesianPts &a, const CartesianPts &b)
        {
            if (a.Y_ == b.Y_)
                return a.X_ < b.X_;
            else
                return a.Y_ < b.Y_;
        }

        static bool compare(const int &a, const int &b)
        {
            return a < b;
        }
        void SortByX(vector <CartesianPts> &vectPts_)
        {
            sort(vectPts_.begin(), vectPts_.end(), compareCoord);
        }
        void SortByY(vector <CartesianPts> &vectPts_)
        {
            sort(vectPts_.begin(), vectPts_.end(), compareCoordY);
        }

        void Sort(vector <int> &vectPts)
        {
            sort(vectPts.begin(), vectPts.end());
        }

    private:
        vector <CartesianPts> m_vectBoundary;
        vector <int> m_vectGridX;
        vector <int> m_vectGridY;
        vector <int> m_vectGrid;
        vector <int> m_vectX;
        vector <int> m_vectY;
        map <CartesianPts, bool, valueComp> m_mapPts;
        bool m_isInGrid, m_isOutGrid;
};
#endif
