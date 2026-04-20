#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include <cstdlib>
#include <cmath>
#include "rand_norm.h"





// Class that represents a mathematical matrix
class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

	void initialise_normal()
	{
		for(int i=0; i<nRows;i++)
		{
			for(int j=0; j<nCols;j++)
			{
				A[j + i * nCols]=rand_normal()/nRows;
			}
		}

	}

	void normalise_columns()
{
    for(int j = 0; j < nCols; j++)
    {
        // compute column norm
        double norm = 0.0;
        for(int i = 0; i < nRows; i++)
        {
            double val = A[j + i*nCols];
            norm += val * val;
        }
        norm = std::sqrt(norm);

        // avoid division by zero
        if(norm > 0.0)
        {
            for(int i = 0; i < nRows; i++)
            {
                A[j + i*nCols] /= norm;
            }
        }
    }
}


private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};

inline MVector operator*(const MMatrix& A, const MVector& x)
{
	if (x.size()!=A.Cols())
	{
		std::cout<<"incompatible operator"<<std::endl;
		return {0};
	}
	MVector v(A.Rows());
	for (int i=0; i<A.Rows(); i++)
	{
		v[i]=0;
		for(int j=0; j<A.Cols();j++)
		{
			v[i]+=A(i,j)*x[j];
		}
	}
	return v;
}

inline MMatrix operator*(const MMatrix& A, const MMatrix& B)
{
	if(A.Cols() !=B.Rows())
	{
		std::cout<<"incompatible matrix operation"<<std::endl;
		MMatrix C(0,0,0);
		return C;
	}
	MMatrix D(A.Rows(),B.Cols(),0);
	for(int i=0;i<A.Rows();i++)
	{
		for(int j=0; j<B.Cols();j++)
		{
			for(int k=0; k<A.Cols();k++)
			{
				D(i,j)+=A(i,k)*B(k,j);
			}
		}
	}
	return D;
}

inline MMatrix Transpose(const MMatrix& A)
{
	MMatrix B(A.Cols(), A.Rows(),0);
	for(int i=0; i<A.Cols(); i++)
	{
		for(int j=0; j<A.Rows(); j++)
		{
			B(i,j)=A(j,i);
		}
	}
	return B;
}

std::ostream& operator<<(std::ostream& os, const MMatrix& A)
{
	os.width(1);
	os << "(";
	os.width(9);
	os.precision(9);
	os<<A(0,0);
	os<<",";

	for(int i=0;i<A.Rows();i++)
	{
		for(int j=0;j<A.Cols();j++)
		{
			if (i==0 && j==0) continue;
			os.width(10);
			os.precision(10);
			os<<A(i,j);
			if(j<A.Cols()-1) os<<",";
		}
		if(i<A.Rows()-1) os<<std::endl;
	}
	os << ")";
	return os;
}

#endif
