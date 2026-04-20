#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <cstdlib>
#include "rand_norm.h"
#include <ctime>

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) 
	{ 
		return v[index];
	}

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		return v[index]; 
	}

	int size() const { return v.size(); } // number of elements

	// add data to your vector
	void push_back(double x){v.push_back(x);}

	// infinity norm
	double LInfNorm() const
	{
		double maxAbs = 0;
		std::size_t s = size();
		for (int i=0; i<s; i++)
		{
			maxAbs = std::max(std::abs(v[i]), maxAbs);
		}
		return maxAbs;
	}

	// L2 norm
	double L2Norm() const
	{
		double sum = 0;
		std::size_t s = size();
		for (int i=0; i<s; i++)
		{
			sum+=v[i]*v[i];
		}
		return std::sqrt(sum);
	}

	auto begin() { return v.begin(); }

	auto end()   { return v.end(); }

	auto begin() const { return v.begin(); }

	auto end() const   { return v.end(); }

	// Threshold
	void Threshold(int k)
	{
		MVector w(v.size(),0);
        auto itmin = std::min_element(v.begin(), v.end());
        int indexmin = std::distance(v.begin(), itmin);
        auto itmax = std::max_element(v.begin(), v.end());
        int indexmax = std::distance(v.begin(), itmax);
        for(int i=0;i<k;i++)
        {
            if(std::abs(*itmax)>std::abs(*itmin))
            {
                w[indexmax]=*itmax;
                v[indexmax]=0;
                itmax = std::max_element(v.begin(), v.end());
                indexmax = std::distance(v.begin(), itmax);
            }
            else
            {
                w[indexmin]=*itmin;
                v[indexmin]=0;
                itmin = std::min_element(v.begin(), v.end());
                indexmin = std::distance(v.begin(), itmin);
            }
        }
		*this=w;
	}

	// sparse initialiser
	void initialise_sparse_normal(int k) 
	{
		if (k > size())
		{
			std::cout<<"invalid argument"<<std::endl;
			return;
		}
		// Step 2: generate indices [0..size()-1]
		std::vector<int> indices(size());
		std::iota(indices.begin(), indices.end(), 0);
		// Step 3: shuffle them
		std::mt19937 gen(static_cast<unsigned>(std::time(nullptr)));
		std::shuffle(indices.begin(), indices.end(), gen);
		// Step 4: set k random positions using rand_normal()
		for (int i = 0; i < k; ++i) {
			v[indices[i]] = rand_normal();
		}
	}

private:
	std::vector<double> v;
};

//1.2.1

inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]*=lhs;
	return temp;
}

inline MVector operator*(const MVector& rhs, const double& lhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]*=lhs;
	return temp;
}

inline MVector operator/(const MVector& rhs, const double& lhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]/=lhs;
	return temp;
}

inline MVector operator+(const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()!=rhs.size()) 
	{
		std::cout<<"incompatible vector size"<<std::endl;
		return MVector();
	}
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]=lhs[i]+rhs[i];
	return temp;
}

inline MVector operator-(const MVector& lhs, const MVector& rhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]=lhs[i]-rhs[i];
	return temp;
}

inline double dot(const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()!=rhs.size())
	{
		std::cout<<"incompatible vectors"<<std::endl; 
		return 1e6;
	}
	double sum=0;
	for(int i=0;i<lhs.size();i++)
	{
		sum+=lhs[i]*rhs[i];
	}
	return sum;
}

std::ostream& operator<<(std::ostream& os, const MVector& v)
{
	int n = v.size();
	os << "("<<v[0];
	for(int i=1; i<n;i++) os<<","<<v[i];
	os << ")";
	return os;
}

#endif
