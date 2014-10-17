#include <stdio.h>
#include "tensors.h"
#include "assert.h"

///////////////////////////////////////////////////////////////////////////////
// Copy the contents of from basetensor to target tensor.  There is absolutely
// NO error checking being performed here
void basetensor::CopyFrom( basetensor* from )
{

    double* from_vec = from->vec;

    #pragma omp parallel for
    for( int i=0; i < this->length; i++ )
    { vec[i] = from_vec[i]; }

}

// function that sets each entry = x
void basetensor::setall( double x )
{

    #pragma omp parallel for
    for( int i=0; i < this->length; i++ )
    { vec[i] = x; }

}

///////////////////////////////////////////////////////////////////////////////
// 1D tensor. This is an array with base 1 indexing
///////////////////////////////////////////////////////////////////////////////
dTensor1::dTensor1(int length )
{

    // constructor
    this->vec = new double[length];
    this->length = length;

}

dTensor1::~dTensor1()
{
    // destructor
    delete[] vec;
}

#ifdef TENSOR_DEBUG
void dTensor1::set( int i, double value )
{
    int index = i-1;
    assert( index > -1 && index < this->length );
    vec[index] = value;
}
double dTensor1::get( int i ) const
{
    int index = i-1;
    return vec[index];
}
#endif

int dTensor1::get_size( ) const
{
    return this->length;
}


///////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// 1D vector with boundary cells
//////////////////////////////////////////////////////////////////////////////
dTensorBC1::dTensorBC1(int length, int mbc )
{
    l1        = length + 2*mbc;
    this->mbc = mbc;
    this->length = l1;
    vec       = new double[ l1 ];
}
            
dTensorBC1::~dTensorBC1()
{
    delete[] vec;
}

#ifdef TENSOR_DEBUG
double dTensorBC1::get( int i ) const
{
    int index = i + mbc - 1;
    assert( index > -1 && index < this->length );
    return vec[index];
}

void dTensorBC1::set( int i, double value )
{


    int index  = i + mbc-1;

    assert( index > -1 && index < this->length );
    if( index < 0 || index > l1-1 )
    {
        perror("invalid index  ");
        printf("index = %d\n", index );
    }

    vec[index] = value;
}
#endif

int dTensorBC1::get_size( ) const
{
    return l1;
}

//////////////////////////////////////////////////////////////////////////////
// 2D vector with no boundary cells
//////////////////////////////////////////////////////////////////////////////
dTensor2::dTensor2(int mx, int meqn )
{

    this->l1  = mx;
    this->l2  = meqn;
    this->length = l1*l2;
    assert( this->length > 0 );
    vec = new double[ l1*l2 ];

}

dTensor2::~dTensor2()
{
    delete[] vec;
}

#ifdef TENSOR_DEBUG
double dTensor2::get( int i, int j ) const
{

    int index = (i-1)*l2 + (j-1);
    assert( index > -1 && index < this->length );
    return vec[index];

}
void dTensor2::set( int i, int j, double value )
{
    int index  = (i-1)*l2 + (j-1);
    vec[index] = value;
}
#endif

int dTensor2::get_size( int col ) const
{

    switch( col )
    {
        case 1:
        return l1;

        case 2:
        return l2;

        default:
        perror("invalid column number\n");
        return -1;
    }


}

//////////////////////////////////////////////////////////////////////////////
// 3D vector 
//////////////////////////////////////////////////////////////////////////////
dTensor3::dTensor3(int mx, int my, int meqn )
{

    this->l1  = mx;
    this->l2  = my;
    this->l3  = meqn;

    this->length = l1*l2*l3;
    vec = new double[ l1*l2*l3 ];

}

dTensor3::~dTensor3()
{
    delete[] vec;
}

#ifdef TENSOR_DEBUG
double dTensor3::get( int i, int j, int m ) const
{

    int index = ( (i-1)*l2 + (j-1) )*l3 + (m-1);
    assert( index > -1 && index < this->length );
    return vec[index];

}

void dTensor3::set( int i, int j, int m, double value )
{
    int index = ( (i-1)*l2 + (j-1) )*l3 + (m-1);
    assert( index > -1 && index < this->length );
    vec[index] = value;
}
#endif

int dTensor3::get_size( int col ) const 
{

    switch( col )
    {
        case 1:
        return l1;

        case 2:
        return l2;

        case 3:
        return l3;

        default:
        perror("invalid column number\n");
        return -1;
    }


}

//////////////////////////////////////////////////////////////////////////////
// 2D vector with boundary cells on the first index
//////////////////////////////////////////////////////////////////////////////
dTensorBC2::dTensorBC2(int mx, int meqn, int mbc )
{

    this->l1  = mx + 2*mbc;
    this->l2  = meqn;
    this->mbc = mbc;

    this->length = l1*l2;
    vec = new double[ l1*l2 ];

}

dTensorBC2::~dTensorBC2()
{
    delete[] vec;
}

#ifdef TENSOR_DEBUG
double dTensorBC2::get( int i, int j ) const
{

    int index = (i+mbc-1)*l2 + (j-1);
    assert( index > -1 && index < this->length );
    return vec[index];

}

void dTensorBC2::set( int i, int j, double value )
{
    int index  = (i+mbc-1)*l2 + (j-1);
    assert( index > -1 && index < this->length );
    vec[index] = value;
}
#endif

int dTensorBC2::get_size( int col ) const 
{

    switch( col )
    {
        case 1:
        return l1-2*mbc;

        case 2:
        return l2;

        default:
        perror("invalid column number\n");
        return -1;
    }


}

//////////////////////////////////////////////////////////////////////////////
// 3D vector with boundary cells on the first two indices
//////////////////////////////////////////////////////////////////////////////
dTensorBC3::dTensorBC3(int mx, int my, int meqn, int mbc )
{

    this->l1  = mx + 2*mbc;
    this->l2  = my + 2*mbc;
    this->l3  = meqn;
    this->mbc = mbc;

    this->length = l1*l2*l3;
    vec = new double[ l1*l2*l3 ];

}

dTensorBC3::~dTensorBC3()
{
    delete[] vec;
}

#ifdef TENSOR_DEBUG
double dTensorBC3::get( int i, int j, int m ) const
{

    int index = ( (i+mbc-1)*l2 + (j+mbc-1) )*l3 + (m-1);
    assert( index > -1 && index < this->length );
    return vec[index];

}

void dTensorBC3::set( int i, int j, int m, double value )
{
    int index = ( (i+mbc-1)*l2 + (j+mbc-1) )*l3 + (m-1);
    assert( index > -1 && index < this->length );
    vec[index] = value;
}
#endif

int dTensorBC3::get_size( int col ) const 
{

    switch( col )
    {
        case 1:
        return l1-2*mbc;

        case 2:
        return l2-2*mbc;

        case 3:
        return l3;

        default:
        perror("invalid column number\n");
        return -1;
    }


}
