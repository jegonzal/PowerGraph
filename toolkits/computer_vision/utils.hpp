
// utils.hpp - miscellaneous utilities 
// Originally from Nicol N. Schraudolph's isinf package
// Later expanded by Dhruv Batra

#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <cstdlib>
#include <string>

// row-major array access
#define ARR_RM(arr, r_ind, c_ind, ncols) (*(arr + r_ind*ncols + c_ind))
// col-major array access
#define ARR_CM(arr, r_ind, c_ind, nrows) (*(arr + c_ind*nrows + r_ind))

// row-major ind2sub
#define IND2SUB_RM(ind,r,c,ncols) \
        r = floor(ind/ncols);     \
        c = ind % ncols;
// column-major ind2sub
#define IND2SUB_CM(ind,r,c,nrows) \
        c = floor(ind/nrows);     \
        r = ind % nrows;     

// row-major sub2ind
#define SUB2IND_RM(r,c,ncols) r*ncols + c

// col-major sub2ind
#define SUB2IND_CM(r,c,nrows) c*nrows + r

// operators & formatted I/O for vectors
// inner product
template <class T>
inline T operator*(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    T sum(0);
    for (size_t i = 0; i < a.size(); ++i)
        sum += a[i]*b[i];
    return sum;
}

// element-wise sum
template <class T>
inline std::vector<T>& operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    std::vector<T> sum(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        sum[i] = a[i]+b[i];
    return sum;
}

template <class T>
inline std::vector<T>& operator+=(std::vector<T>& a, const T& b)
{
    typename std::vector<T>::iterator i(a.begin());
    while (i != a.end()) *i++ += b;
    return a;
}

template <class T>
inline std::vector<T>& operator-=(std::vector<T>& a, const T& b)
{
    typename std::vector<T>::iterator i(a.begin());
    while (i != a.end()) *i++ -= b;
    return a;
}

template <class T>
inline std::vector<T>& operator*=(std::vector<T>& a, const T& b)
{
    typename std::vector<T>::iterator i(a.begin());
    while (i != a.end()) *i++ *= b;
    return a;
}

template <class T>
inline std::vector<T>& operator/=(std::vector<T>& a, const T& b)
{
    typename std::vector<T>::iterator i(a.begin());
    while (i != a.end()) *i++ /= b;
    return a;
}

template <class T>
inline std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    typename std::vector<T>::iterator i(a.begin());
    typename std::vector<T>::const_iterator j(b.begin());
    while (i != a.end()) *i++ += *j++;
    return a;
}

template <class T>
inline std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    typename std::vector<T>::iterator i(a.begin());
    typename std::vector<T>::const_iterator j(b.begin());
    while (i != a.end()) *i++ -= *j++;
    return a;
}

template <class T>
inline std::vector<T>& operator*=(std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    typename std::vector<T>::iterator i(a.begin());
    typename std::vector<T>::const_iterator j(b.begin());
    while (i != a.end()) *i++ *= *j++;
    return a;
}

template <class T>
inline std::vector<T>& operator/=(std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    typename std::vector<T>::iterator i(a.begin());
    typename std::vector<T>::const_iterator j(b.begin());
    while (i != a.end()) *i++ /= *j++;
    return a;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& x)
{
    typename std::vector<T>::const_iterator i(x.begin());
    while(i != x.end()) os << *i++ << ' ';
    return os;
}

template <class T>
inline std::istream& operator>>(std::istream& is, std::vector<T>& x)
{
    std::string s;
    const size_t n = x.size();

    while (x.size() == n)
    {
        getline(is, s);
        if (is.fail()) break;

        std::istringstream iss(s);
        T item;

        iss >> item;
        while (iss.good())
        {
            x.push_back(item);
            iss >> item;
        }
        if (!iss.fail())
            x.push_back(item);
    }

    return is;
}

// Function to write a vector to file (Assumes << is defined for type T)
// CHECK_NULL is provided by Danny Tarlow's Nymph Utils
template <typename T> void WriteToFile(std::string fname, std::vector<T> vecx)
{
    std::ofstream fout; 
    fout.open(fname.c_str());
    
    //CHECK_NULL(fout.fail(),"Could not open file for writing results\n");
    
    fout << vecx;
    
    fout.close();
}


#endif

