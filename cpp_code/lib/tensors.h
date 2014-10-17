#ifndef _TENSORS_H_
#define _TENSORS_H_

#define TENSOR_DEBUG

class basetensor
{
    public:
        void CopyFrom( basetensor* from );
        void setall( double x );

    protected:
        double* vec;
        int length; // length of the base vector

};

class basetensorBC : public basetensor
{
    public:
        int get_mbc() const {return mbc;}

    protected:
        int mbc;
};

//class dTensor1 : public basetensor
class dTensor1 : public basetensor
{

    // 1d vectors without ghost cells.  These use indices starting at 1 rather
    // than 0.
    public:

        dTensor1(int length );
        ~dTensor1();

        #ifdef TENSOR_DEBUG
        void   set( int index, double value );
        double get( int index ) const;
        #else
        inline void   set( int ind, double val ){ vec[ind-1] = val; }
        inline double get( int index ) const    { return vec[index-1]; }
        #endif
        int    get_size( ) const;

};

class dTensor2 : public basetensor
{

    private:
        
        // lengths of the two indices.  Note: boundary cells are only applied to
        // the first index!
        int l1, l2;

    public:
    
        dTensor2(int i, int m );
        ~dTensor2();

        #ifdef TENSOR_DEBUG
        double get( int i, int j ) const;
        void   set( int i, int m, double value );
        #else
        inline double get( int i, int j ) const
        { return vec[(i-1)*l2 + (j-1)];}
        inline void   set( int i, int j, double value )
        { vec[(i-1)*l2 + (j-1)] = value ;}
        #endif

        int    get_size(int col) const;
        int    get_mbc() const;

};

class dTensor3 : public basetensor
{

    public:
    
        dTensor3(int i, int j, int m );
        ~dTensor3();

        #ifdef TENSOR_DEBUG
        double get( int i, int j, int m ) const;
        void   set( int i, int j, int m, double value );
        #else
        inline double get( int i, int j, int m ) const
        { return vec[ ( (i+m-1)*l2 + (j-1) )*l3 + (m-1) ]; }

        inline void set( int i, int j, int m, double value )
        { vec[ ( (i+m-1)*l2 + (j-1) )*l3 + (m-1) ] = value; }
        #endif

        int    get_size(int col) const;

    private:
        
        // lengths of the two indices.  Note: boundary cells are only applied to
        // the first index!
        int l1, l2, l3;


};


class dTensorBC1 : public basetensorBC
{

    private:

        int l1;     // length of the first index

    public:

        dTensorBC1(int length, int mbc );
        ~dTensorBC1();

        #ifdef TENSOR_DEBUG
        double get( int index ) const;
        void   set( int index, double value );
        #else
        inline double get( int i ) const { return vec[ i+mbc-1 ]; }
        inline void set( int i, double val ) const { vec[ i+mbc-1 ] = val; }
        #endif
        int    get_size( ) const;

};

class dTensorBC2 : public basetensorBC
{

    public:
    
        dTensorBC2(int i, int m, int mbc );
        ~dTensorBC2();

        #ifdef TENSOR_DEBUG
        double get( int i, int j ) const;
        void   set( int i, int m, double value );
        #else
        inline double get( int i, int j ) const
        { return vec[(i+mbc-1)*l2 + (j-1)]; }
        inline void set( int i, int j, double value )
        { vec[(i+mbc-1)*l2 + (j-1)] = value; }
        #endif

        int    get_size(int col) const;

    private:
        
        // lengths of the two indices.  Note: boundary cells are only applied to
        // the first index!
        int l1, l2;


};

class dTensorBC3 : public basetensorBC
{

    public:
    
        dTensorBC3(int i, int j, int m, int mbc );
        ~dTensorBC3();

        #ifdef TENSOR_DEBUG
        double get( int i, int j, int m ) const;
        void   set( int i, int j, int m, double value );
        #else
        inline double get( int i, int j, int m ) const
        { return vec[ ( (i+mbc-1)*l2 + (j+mbc-1) )*l3 + (m-1) ]; }

        inline void set( int i, int j, int m, double value )
        { vec[ ( (i+mbc-1)*l2 + (j+mbc-1) )*l3 + (m-1) ] = value; }
        #endif

        int    get_size(int col) const;

    private:
        
        // lengths of the two indices.  Note: boundary cells are only applied to
        // the first index!
        int l1, l2, l3;


};


#endif
