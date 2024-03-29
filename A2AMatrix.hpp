/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AMatrix.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef A2A_Matrix_hpp_
#define A2A_Matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>
//#include <Hadrons/DiskVector.hpp>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

#ifndef HADRONS_A2AM_NAME 
#define HADRONS_A2AM_NAME "a2aMatrix"
#endif

#ifndef HADRONS_A2AM_IO_TYPE
#define HADRONS_A2AM_IO_TYPE ComplexF
#endif

#define HADRONS_A2AM_PARALLEL_IO

BEGIN_HADRONS_NAMESPACE

// general A2A matrix set based on Eigen tensors and Grid-allocated memory
// Dimensions:
//   0 - ext - external field (momentum, EM field, ...)
//   1 - str - spin-color structure
//   2 - t   - timeslice
//   3 - i   - left  A2A mode index
//   4 - j   - right A2A mode index
template <typename T>
using A2AMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 5, Eigen::RowMajor>>;

template <typename T>
using A2AMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

template <typename T>
using A2AMatrixTr = Eigen::Matrix<T, -1, -1, Eigen::ColMajor>;

// supposing 0 - index1 (igamma)
//           1 - index2 (imode)
//           2 - 4d location
template <typename SiteT>
using FieldMatrix = Eigen::Tensor<SiteT, 3>;

/******************************************************************************
 *                      Abstract class for A2A kernels                        *
 ******************************************************************************/
template <typename T, typename FImpl, typename Field>
class A2AKernel
{
public:
  A2AKernel(void) = default;
  virtual ~A2AKernel(void) = default;
  virtual void operator()(A2AMatrixSet<T> &m,
			  const Field *left, const Field *right,
			  const unsigned int orthogDim, double &time) = 0;
  typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSinglet<vector_type> Singlet_v;

  // New Extended Meson Field for type1
  virtual void operator()(A2AMatrixSet<T> &m,
			  const FieldMatrix<SpinColourVector_v> &left_loop,
			  const Field *right, const int offset,
                          const unsigned int orthogDim, double &time) = 0;

  // Extended Meson Field in type0
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<Singlet_v> &loop_cfs,
			  const Field *left, const Field *right,
                          const unsigned int orthogDim, double &time) = 0;

  // Extended Meson Field in type1
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinColourMatrix_v > &loop_m,
			  const Field *left, const Field *right,
                          const unsigned int orthogDim, double &time) = 0;

  // Extended Meson Field in type2
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<ColourMatrix_v > &loop_m,
			  const Field *left, const Field *right,
                          const unsigned int orthogDim, double &time) = 0;

  // Extended Meson Field in type3
  virtual void operator()(A2AMatrixSet<T> &m,
			  const std::vector<SpinMatrix_v > &loop_m,
			  const Field *left, const Field *right,
                          const unsigned int orthogDim, double &time) = 0;

  // QuarkLoop in type0, ExtendedMesonField
  virtual void operator()(std::vector<Singlet_v > &m,
			  const Field *v, const Field *w,
			  const int Nblock) = 0;

  // QuarkLoop in type1, ExtendedMesonField
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const Field *v, const Field *w,
			  const int Nblock) = 0;

  // QuarkLoop in type2, ExtendedMesonField
  virtual void operator()(std::vector<ColourMatrix_v > &m,
			  const Field *v, const Field *w,
			  const int Nblock) = 0;

  // QuarkLoop in type3, ExtendedMesonField
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const Field *v, const Field *w,
			  const int Nblock) = 0;

  // VpiW in type0, Monster Loop
  virtual void operator()(std::vector<Singlet_v > &m,
			  const Field *v, const Field *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj ) = 0;

  // VpiW in type1, Monster Loop
  virtual void operator()(std::vector<SpinColourMatrix_v > &m,
			  const Field *v, const Field *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj ) = 0;

  // VpiW in type2, Monster Loop
  virtual void operator()(std::vector<ColourMatrix_v > &m,
			  const Field *v, const Field *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj ) = 0;

  // VpiW in type3, Monster Loop
  virtual void operator()(std::vector<SpinMatrix_v > &m,
			  const Field *v, const Field *w,
			  const std::vector<A2AMatrix<Complex> > &mes,
			  const int i,  const int j,
			  const int Ni, const int Nj ) = 0;

  // Left * Loop for type 1
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<SpinColourMatrix_v > &m,
			  const Field *left, const int Nblock,
			  const int offset) = 0;

  // Left * Loop for type 2
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<ColourMatrix_v > &m,
			  const Field *left, const int Nblock,
			  const int offset) = 0;

  // Left * Loop for type 3
  virtual void operator()(FieldMatrix<SpinColourVector_v > &left_loop,
			  const std::vector<SpinMatrix_v > &m,
			  const Field *left, const int Nblock,
			  const int offset) = 0;

  virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
  virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
};

/******************************************************************************
 *                  Class to handle A2A matrix block HDF5 I/O                 *
 ******************************************************************************/
template <typename T>
class A2AMatrixIo
{
public:
    // constructors
    A2AMatrixIo(void) = default;
    A2AMatrixIo(std::string filename, std::string dataname, 
                const unsigned int nt, const unsigned int ni = 0,
                const unsigned int nj = 0);
    // destructor
    ~A2AMatrixIo(void) = default;
    // access
    unsigned int getNi(void) const;
    unsigned int getNj(void) const;
    unsigned int getNt(void) const;
    size_t       getSize(void) const;
    // file allocation
    template <typename MetadataType>
    void initFile(const MetadataType &d, const unsigned int chunkSize);
    // block I/O
    void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                   const unsigned int blockSizei, const unsigned int blockSizej);
    void saveBlock(const A2AMatrixSet<T> &m, const unsigned int ext, const unsigned int str,
                   const unsigned int i, const unsigned int j);
    template <template <class> class Vec, typename VecT>
    void load(Vec<VecT> &v, double *tRead = nullptr, GridBase *grid = nullptr);
private:
    std::string  filename_{""}, dataname_{""};
    unsigned int nt_{0}, ni_{0}, nj_{0};
};

/******************************************************************************
 *                  Wrapper for A2A matrix block computation                  *
 ******************************************************************************/
  template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo = T>
class A2AMatrixBlockComputation
{
private:
    struct IoHelper
    {
        A2AMatrixIo<TIo> io;
        MetadataType     md;
        unsigned int     e, s, i, j;
    };
    typedef std::function<std::string(const unsigned int, const unsigned int)>  FilenameFn;
    typedef std::function<std::string(const unsigned int)>  FilenameFnM;
    typedef std::function<MetadataType(const unsigned int, const unsigned int)> MetadataFn;
    typedef std::function<MetadataType(const unsigned int)> MetadataFnM;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
public:
    // constructor
    A2AMatrixBlockComputation(GridBase *grid,
                              const unsigned int orthogDim,
                              const unsigned int next,
                              const unsigned int nstr,
                              const unsigned int blockSize,
                              const unsigned int cacheBlockSize,
                              TimerArray *tArray = nullptr);
    // execution
    void execute(const std::vector<Field> &left, 
                 const std::vector<Field> &right,
                 A2AKernel<T, FImpl, Field> &kernel,
                 const FilenameFn &ionameFn,
                 const FilenameFn &filenameFn,
                 const MetadataFn &metadataFn);

  //template <typename FImpl>
    void execute0(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn,
		  const int ngamma);
    void execute1(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn);
    void execute2(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn,
		  const int ngamma);
    void execute3(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn);
    void execute0(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
		  const std::vector<A2AMatrix<Complex> > &meson,
		  const int &it_mes,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn,
		  const int ngamma);
    void execute1(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<SpinColourMatrix_v> &mat,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn);
    void execute2(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
		  const std::vector<A2AMatrix<Complex> > &meson,
		  const int &it_mes,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn,
		  const int ngamma);
    void execute3(const std::vector<Field> &left, 
                  const std::vector<Field> &right,
		  const std::vector<Field> &loop1,
		  const std::vector<Field> &loop2,
		  const std::vector<A2AMatrix<Complex> > &meson,
		  const int &it_mes,
                  A2AKernel<T, FImpl, Field> &kernel,
                  const FilenameFnM &ionameFn,
                  const FilenameFnM &filenameFn,
                  const MetadataFnM &metadataFn);
private:
    // I/O handler
    void saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h);
private:
    TimerArray            *tArray_;
    GridBase              *grid_;
    unsigned int          orthogDim_, nt_, next_, nstr_, blockSize_, cacheBlockSize_;
    Vector<T>             mCache_;
    Vector<TIo>           mBuf_;
    std::vector<IoHelper> nodeIo_;
};

/******************************************************************************
 *                       A2A matrix contraction kernels                       *
 ******************************************************************************/
class A2AContraction
{
public:
    // accTrMul(acc, a, b): acc += tr(a*b)
    template <typename C, typename MatLeft, typename MatRight>
    static inline void accTrMul(C &acc, const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;
        if ((MatLeft::Options  == RowMajor) and
            (MatRight::Options == ColMajor))
        {
  	  thread_for(r,a.rows(),
            {
                C tmp;
#ifdef USE_MKL
                dotuRow(tmp, r, a, b);
#else
                tmp = a.row(r).conjugate().dot(b.col(r));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
        else
	  {
            thread_for(c,a.cols(),
            {
                C tmp;
#ifdef USE_MKL 
                dotuCol(tmp, c, a, b);
#else
                tmp = a.col(c).conjugate().dot(b.row(c));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
    }

    template <typename MatLeft, typename MatRight>
    static inline double accTrMulFlops(const MatLeft &a, const MatRight &b)
    {
        double n = a.rows()*a.cols();

        return 8.*n;
    }

    // mul(res, a, b): res = a*b
#ifdef USE_MKL
    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexD, Opts...> &res, 
                           const Mat<ComplexD, Opts...> &a, 
                           const Mat<ComplexD, Opts...> &b)
    {
        static const ComplexD one(1., 0.), zero(0., 0.);
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexD, Opts...>::Options == RowMajor)
        {
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexD, Opts...>::Options == ColMajor)
        {
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }

    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexF, Opts...> &res, 
                           const Mat<ComplexF, Opts...> &a, 
                           const Mat<ComplexF, Opts...> &b)
    {
        static const ComplexF one(1., 0.), zero(0., 0.);
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexF, Opts...>::Options == RowMajor)
        {
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexF, Opts...>::Options == ColMajor)
        {
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }
#else
    template <typename Mat>
    static inline void mul(Mat &res, const Mat &a, const Mat &b)
    {
        res = a*b;
    }
#endif
    template <typename Mat>
    static inline double mulFlops(const Mat &a, const Mat &b)
    {
        double nr = a.rows(), nc = a.cols();

        return nr*nr*(6.*nc + 2.*(nc - 1.));
    }
private:
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotRowPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aRow, 
                                    const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if (MatLeft::Options == RowMajor)
        {
            aPt  = a.data() + aRow*a.cols();
            aInc = 1;
        }
        else if (MatLeft::Options == ColMajor)
        {
            aPt  = a.data() + aRow;
            aInc = a.rows();
        }
        if (MatRight::Options == RowMajor)
        {
            bPt  = b.data() + aRow;
            bInc = b.cols();
        }
        else if (MatRight::Options == ColMajor)
        {
            bPt  = b.data() + aRow*b.rows();
            bInc = 1;
        }
    }

#ifdef USE_MKL
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotColPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aCol, 
                                    const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;
        if (MatLeft::Options == RowMajor)
        {
            aPt  = a.data() + aCol;
            aInc = a.cols();
        }
        else if (MatLeft::Options == ColMajor)
        {
            aPt  = a.data() + aCol*a.rows();
            aInc = 1;
        }
        if (MatRight::Options == RowMajor)
        {
            bPt  = b.data() + aCol*b.cols();
            bInc = 1;
        }
        else if (MatRight::Options == ColMajor)
        {
            bPt  = b.data() + aCol;
            bInc = b.rows();
        }
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexF &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_cdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexF &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_cdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexD &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_zdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexD &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_zdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }
#endif
};

/******************************************************************************
 *                     A2AMatrixIo template implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
A2AMatrixIo<T>::A2AMatrixIo(std::string filename, std::string dataname, 
                            const unsigned int nt, const unsigned int ni,
                            const unsigned int nj)
: filename_(filename), dataname_(dataname)
, nt_(nt), ni_(ni), nj_(nj)
{}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
unsigned int A2AMatrixIo<T>::getNt(void) const
{
    return nt_;
}

template <typename T>
unsigned int A2AMatrixIo<T>::getNi(void) const
{
    return ni_;
}

template <typename T>
unsigned int A2AMatrixIo<T>::getNj(void) const
{
    return nj_;
}

template <typename T>
size_t A2AMatrixIo<T>::getSize(void) const
{
    return nt_*ni_*nj_*sizeof(T);
}

// file allocation /////////////////////////////////////////////////////////////
template <typename T>
template <typename MetadataType>
void A2AMatrixIo<T>::initFile(const MetadataType &d, const unsigned int chunkSize)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t>    dim = {static_cast<hsize_t>(nt_), 
                                   static_cast<hsize_t>(ni_), 
                                   static_cast<hsize_t>(nj_)},
      //chunk = {static_cast<hsize_t>(nt_), 
                            chunk = {static_cast<hsize_t>(1), 
                                     static_cast<hsize_t>(chunkSize), 
                                     static_cast<hsize_t>(chunkSize)};
    if ( ni_ < chunkSize ) chunk[1] = ni_;
    if ( nj_ < chunkSize ) chunk[2] = nj_;
    H5NS::DataSpace         dataspace(dim.size(), dim.data());
    H5NS::DataSet           dataset;
    H5NS::DSetCreatPropList plist;
    
    // create empty file just with metadata
    {
        Hdf5Writer writer(filename_);
        write(writer, dataname_, d);
    }

    // create the dataset
    Hdf5Reader reader(filename_, false);

    push(reader, dataname_);
    auto &group = reader.getGroup();
    std::cout << "setChunk" << std::endl;
    plist.setChunk(chunk.size(), chunk.data());
    std::cout << "setFletcher32" << std::endl;
    plist.setFletcher32();
    std::cout << "createDateSet" << std::endl;
    dataset = group.createDataSet(HADRONS_A2AM_NAME, Hdf5Type<T>::type(), dataspace, plist);
    std::cout << "Done" << std::endl;
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

// block I/O ///////////////////////////////////////////////////////////////////
template <typename T>
void A2AMatrixIo<T>::saveBlock(const T *data, 
                               const unsigned int i, 
                               const unsigned int j,
                               const unsigned int blockSizei,
                               const unsigned int blockSizej)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_, false);
    std::vector<hsize_t> count = {nt_, blockSizei, blockSizej},
                         offset = {0, static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1, 1},
                         block  = {1, 1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    H5NS::DataSet        dataset;
    //    size_t               shift;

    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AM_NAME);
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void A2AMatrixIo<T>::saveBlock(const A2AMatrixSet<T> &m,
                               const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j)
{
    unsigned int blockSizei = m.dimension(3);
    unsigned int blockSizej = m.dimension(4);
    unsigned int nstr       = m.dimension(1);
    size_t       offset     = (ext*nstr + str)*nt_*blockSizei*blockSizej;

    saveBlock(m.data() + offset, i, j, blockSizei, blockSizej);
}

template <typename T>
template <template <class> class Vec, typename VecT>
void A2AMatrixIo<T>::load(Vec<VecT> &v, double *tRead, GridBase *grid)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;

    if (!(grid) || grid->IsBoss())
    {
        Hdf5Reader reader(filename_);
        push(reader, dataname_);
        auto &group = reader.getGroup();
        dataset = group.openDataSet(HADRONS_A2AM_NAME);
        datatype = dataset.getCompType();
        dataspace = dataset.getSpace();
        hdim.resize(dataspace.getSimpleExtentNdims());
        dataspace.getSimpleExtentDims(hdim.data());
        if ((nt_ * ni_ * nj_ != 0) and
            ((hdim[0] != nt_) or (hdim[1] != ni_) or (hdim[2] != nj_)))
        {
            HADRONS_ERROR(Size, "all-to-all matrix size mismatch (got "
                + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + "x"
                + std::to_string(hdim[2]) + ", expected "
                + std::to_string(nt_) + "x" + std::to_string(ni_) + "x"
                + std::to_string(nj_));
        }
        else if (ni_*nj_ == 0)
        {
            if (hdim[0] != nt_)
            {
                HADRONS_ERROR(Size, "all-to-all time size mismatch (got "
                    + std::to_string(hdim[0]) + ", expected "
                    + std::to_string(nt_) + ")");
            }
            ni_ = hdim[1];
            nj_ = hdim[2];
        }
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), &ni_, sizeof(unsigned int));
        grid->Broadcast(grid->BossRank(), &nj_, sizeof(unsigned int));
    }

    A2AMatrix<T>         buf(ni_, nj_);
    int broadcastSize =  sizeof(T) * buf.size();
    std::vector<hsize_t> count    = {1, static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)},
                         stride   = {1, 1, 1},
                         block    = {1, 1, 1},
                         memCount = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)};
    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    std::cout << "Loading timeslice";
    std::cout.flush();
    *tRead = 0.;
    for (unsigned int tp1 = nt_; tp1 > 0; --tp1)
    {
        unsigned int         t      = tp1 - 1;
        std::vector<hsize_t> offset = {static_cast<hsize_t>(t), 0, 0};
        
        if (t % 10 == 0)
        {
            std::cout << " " << t;
            std::cout.flush();
        }
        if (!(grid) || grid->IsBoss())
        {
            dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                      stride.data(), block.data());
        }
        if (tRead) *tRead -= usecond();
        if (!(grid) || grid->IsBoss())
        {
            dataset.read(buf.data(), datatype, memspace, dataspace);
        }
        if (grid)
        {
            grid->Broadcast(grid->BossRank(), buf.data(), broadcastSize);
        }
        if (tRead) *tRead += usecond();
        v[t] = buf.template cast<VecT>();
    }
    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

/******************************************************************************
 *               A2AMatrixBlockComputation template implementation            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::A2AMatrixBlockComputation(GridBase *grid,
                            const unsigned int orthogDim,
                            const unsigned int next, 
                            const unsigned int nstr,
                            const unsigned int blockSize, 
                            const unsigned int cacheBlockSize,
                            TimerArray *tArray)
: grid_(grid), nt_(grid->GlobalDimensions()[orthogDim]), orthogDim_(orthogDim)
, next_(next), nstr_(nstr), blockSize_(blockSize), cacheBlockSize_(cacheBlockSize)
, tArray_(tArray)
{
    mCache_.resize(nt_*next_*nstr_*cacheBlockSize_*cacheBlockSize_);
    mBuf_.resize(nt_*next_*nstr_*blockSize_*blockSize_);
}

#define START_TIMER(name) if (tArray_) tArray_->startTimer(name)
#define STOP_TIMER(name)  if (tArray_) tArray_->stopTimer(name)
#define GET_TIMER(name)   ((tArray_ != nullptr) ? tArray_->getDTimer(name) : 0.)

// execution ///////////////////////////////////////////////////////////////////
template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute(const std::vector<Field> &left, const std::vector<Field> &right,
          A2AKernel<T, FImpl, Field> &kernel, const FilenameFn &ionameFn,
          const FilenameFn &filenameFn, const MetadataFn &metadataFn)
{
    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over blockSize_ factors
    // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
  //std::cout << left << std::endl;
    int    N_i = left.size();
    int    N_j = right.size();
    double flops, bytes, t_kernel;
    double nodes = grid_->NodeCount();

    std::cout << "GLBDIM: " << grid_->GlobalDimensions() << std::endl;
    
    int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
    int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

    for(int i=0;i<N_i;i+=blockSize_)
    for(int j=0;j<N_j;j+=blockSize_)
    {
        // Get the W and V vectors for this block^2 set of terms
        int N_ii = MIN(N_i-i,blockSize_);
        int N_jj = MIN(N_j-j,blockSize_);
        A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

        LOG(Message) << "All-to-all matrix block " 
                     << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
                     << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
                     << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
                     << std::endl;
        // Series of cache blocked chunks of the contractions within this block
        flops    = 0.0;
        bytes    = 0.0;
        t_kernel = 0.0;
        for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
        for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
        {
            double t;
            int N_iii = MIN(N_ii-ii,cacheBlockSize_);
            int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
            A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);
	    /*
	  LOG(Message) << "All-to-all matrix cashe block " 
		       << jj/cacheBlockSize_ + NBlock_j*ii/cacheBlockSize_ + 1 
		       << "/" << pow(blockSize_/cacheBlockSize_,2) << " [" << ii <<" .. " 
		       << ii+N_iii-1 << ", " << jj <<" .. " << jj+N_jjj-1 << "]" 
		       << std::endl;
	    */
            START_TIMER("kernel");
            kernel(mCacheBlock, &left[i+ii], &right[j+jj], orthogDim_, t);
            STOP_TIMER("kernel");
            t_kernel += t;
            flops    += kernel.flops(N_iii, N_jjj);
            bytes    += kernel.bytes(N_iii, N_jjj);

            START_TIMER("cache copy");
            thread_for_collapse( 5,e,next_,{
              for(int s =0;s< nstr_;s++)
              for(int t =0;t< nt_;t++)
              for(int iii=0;iii< N_iii;iii++)
              for(int jjj=0;jjj< N_jjj;jjj++)
              {
                mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
              }
            });
            STOP_TIMER("cache copy");
        }

        // perf
        LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
                     << " Gflop/s/node " << std::endl;
        LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
                     << " GB/s/node "  << std::endl;

        // IO
        double       blockSize, ioTime;
        unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();
    
        LOG(Message) << "Writing block to disk" << std::endl;
        ioTime = -GET_TIMER("IO: write block");
        START_TIMER("IO: total");
        makeFileDir(filenameFn(0, 0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid_->Barrier();
        // make task list for current node
        nodeIo_.clear();
        for(int f = myRank; f < next_*nstr_; f += nRank)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.e  = f/nstr_;
            h.s  = f % nstr_;
            h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                    ionameFn(h.e, h.s), nt_, N_i, N_j);
            h.md = metadataFn(h.e, h.s);
            nodeIo_.push_back(h);
        }
        // parallel IO
        for (auto &h: nodeIo_)
        {
            saveBlock(mBlock, h);
        }
        grid_->Barrier();
#else
        // serial IO, for testing purposes only
        for(int e = 0; e < next_; e++)
        for(int s = 0; s < nstr_; s++)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.e  = e;
            h.s  = s;
            h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                    ionameFn(h.e, h.s), nt_, N_i, N_j);
            h.md = metadataFn(h.e, h.s);
            saveBlock(mfBlock, h);
        }
#endif
        STOP_TIMER("IO: total");
        blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
        ioTime    += GET_TIMER("IO: write block");
        LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                     << ioTime  << " us (" 
                     << blockSize/ioTime*1.0e6/1024/1024
                     << " MB/s)" << std::endl;
    }
}

// Partially contracted 4-quark field used for C13, C23
template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute0(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn,
	   const int ngamma)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSinglet<vector_type> Singlet_v;

  Singlet_v SGL0 = Zero();
  std::vector<Singlet_v> simdlp(MFrvol*ngamma,SGL0);
  std::cout << "Allocated " << MFrvol << "loop scalar" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(simdlp, &loop1[i+ii], &loop2[i+ii], N_iii);
      /*             (a,α)x(β,a)
	 simdlp = Γαβ    / \     for each spacetime x
	                 \_/
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_)
  {
    // Get the W and V vectors for this block^2 set of terms
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

    LOG(Message) << "All-to-all matrix block " 
		 << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		 << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		 << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		 << std::endl;
    // Series of cache blocked chunks of the contractions within this block
    flops    = 0.0;
    bytes    = 0.0;
    t_kernel = 0.0;
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
    {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

      START_TIMER("kernel");
      kernel(mCacheBlock, simdlp, &left[i+ii], &right[j+jj], orthogDim_, t);
      //kernel(mCacheBlock, &left[i+ii], &right[j+jj], orthogDim_, t);
      STOP_TIMER("kernel");
      t_kernel += t;
      flops    += kernel.flops(N_iii, N_jjj);
      bytes    += kernel.bytes(N_iii, N_jjj);

      START_TIMER("cache copy");
      thread_for_collapse( 5,e,next_,{
	  for(int s =0;s< nstr_;s++)
	  for(int t =0;t< nt_;t++)
	  for(int iii=0;iii< N_iii;iii++)
	  for(int jjj=0;jjj< N_jjj;jjj++)
	  {
	    mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	  }
	});
      STOP_TIMER("cache copy");
    }

    // perf
    LOG(Message) << "Kernel time " << t_kernel << " us" << std::endl;
    LOG(Message) << "Kernel flop " << flops << " flops" << std::endl;
    LOG(Message) << "       node " << nodes << std::endl;
    LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		 << " Gflop/s/node " << std::endl;
    LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		 << " GB/s/node "  << std::endl;

    // IO
    double       blockSize, ioTime;
    unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

    LOG(Message) << "Writing block to disk" << std::endl;
    ioTime = -GET_TIMER("IO: write block");
    START_TIMER("IO: total");
    //makeFileDir(filenameFn(0, 0), grid_);
    makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
    grid_->Barrier();
    // make task list for current node
    nodeIo_.clear();
    for(int f = myRank; f < next_*nstr_; f += nRank)
    {
      IoHelper h;

      h.i  = i;
      h.j  = j;
      h.e  = f/nstr_;
      h.s  = f % nstr_;
      //h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
      //		      ionameFn(h.e, h.s), nt_, N_i, N_j);
      h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
			      ionameFn(h.s), nt_, N_i, N_j);
      //h.md = metadataFn(h.e, h.s);
      h.md = metadataFn(h.s);
      nodeIo_.push_back(h);
    }
    // parallel IO
    for (auto &h: nodeIo_)
    {
      saveBlock(mBlock, h);
    }
    grid_->Barrier();
#else
    // serial IO, for testing purposes only
    for(int e = 0; e < next_; e++)
    for(int s = 0; s < nstr_; s++)
    {
      IoHelper h;

      h.i  = i;
      h.j  = j;
      h.e  = e;
      h.s  = s;
      //h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
      //		      ionameFn(h.e, h.s), nt_, N_i, N_j);
      h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
			      ionameFn(h.s), nt_, N_i, N_j);
      //h.md = metadataFn(h.e, h.s);
      h.md = metadataFn(h.s);
      saveBlock(mfBlock, h);
    }
#endif
    STOP_TIMER("IO: total");
    blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
    ioTime    += GET_TIMER("IO: write block");
    LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		 << ioTime  << " us (" 
		 << blockSize/ioTime*1.0e6/1024/1024
		 << " MB/s)" << std::endl;
  }// loop for i, j (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute1(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
           A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
           const FilenameFnM &filenameFn, const MetadataFnM &metadataFn)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  SpinColourMatrix_v SCM0 = Zero();
  std::vector<SpinColourMatrix_v> mat(MFrvol,SCM0);
  std::cout << "Allocated " << MFrvol << "loop matrices" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(mat, &loop1[i+ii], &loop2[i+ii], N_iii);
      /*           (a,α)x(β,b)
	 cfields =     / \     for each spacetime x
	               \_/
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  for(int i=0;i<N_i;i+=blockSize_){
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, mat, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");
    
    for(int j=0;j<N_j;j+=blockSize_){
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_){
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	      mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      //makeFileDir(filenameFn(0, 0), grid_);
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//			ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
	saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
      // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//			ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }// loop for j (Blocks)
  }// loop for i (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute2(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn,
	   const int ngamma)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  ColourMatrix_v CM0 = Zero();
  std::vector<ColourMatrix_v> simdlp(MFrvol*ngamma,CM0);
  std::cout << "Allocated " << MFrvol << "loop scalar" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(simdlp, &loop1[i+ii], &loop2[i+ii], N_iii);
      /*             (a,α)x(β,b)
	 simdlp = Γαβ    / \     for each spacetime x
	                 \_/
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, simdlp, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");
    for(int j=0;j<N_j;j+=blockSize_)
    {
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
      {
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	      mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel time " << t_kernel << " us" << std::endl;
      LOG(Message) << "Kernel flop " << flops << " flops" << std::endl;
      LOG(Message) << "       node " << nodes << std::endl;
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      //makeFileDir(filenameFn(0, 0), grid_);
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//		      ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
			      ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
	saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
    // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//			ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }// loop for j
  }// loop for i (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute3(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
           A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
           const FilenameFnM &filenameFn, const MetadataFnM &metadataFn)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  SpinMatrix_v SM0 = Zero();
  std::vector<SpinMatrix_v> smat(MFrvol,SM0);
  std::cout << "Allocated " << MFrvol << "loop spin matrices" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(smat, &loop1[i+ii], &loop2[i+ii], N_iii);
      /*           (a,α)x(β,a)
	 cfields =     / \     for each spacetime x
	               \_/
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_)
  {
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, smat, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");
    for(int j=0;j<N_j;j+=blockSize_)
    {
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
      {
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	      mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      //makeFileDir(filenameFn(0, 0), grid_);
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//			ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
	saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
      // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	//h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
	//			ionameFn(h.e, h.s), nt_, N_i, N_j);
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	//h.md = metadataFn(h.e, h.s);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }// loop for j (Blocks)
  }// loop for i
}


// Monstor Meson Fields
template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute0(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
	   const std::vector<A2AMatrix<Complex> > &mes, const int &it_mes,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn,
	   const int ngamma)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  //assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSinglet<vector_type> Scalar_v;

  Scalar_v SC0 = Zero();
  std::vector<Scalar_v> cf(MFrvol*ngamma,SC0);
  std::cout << "Allocated " << MFrvol << "loop scalars" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      kernel(cf, &loop1[i+ii], &loop2[j+jj], mes, i+ii, j+jj, N_iii, N_jjj);
      /*               .
	  cf    =     / \     for each spacetime x
	              \ /
		       *
		       Γ
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_)
  {
    // Get the W and V vectors for this block^2 set of terms
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

    LOG(Message) << "All-to-all matrix block " 
		 << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		 << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		 << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		 << std::endl;
    // Series of cache blocked chunks of the contractions within this block
    flops    = 0.0;
    bytes    = 0.0;
    t_kernel = 0.0;
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
    {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

      START_TIMER("kernel");
      kernel(mCacheBlock, cf, &left[i+ii], &right[j+jj], orthogDim_, t);
      STOP_TIMER("kernel");
      t_kernel += t;
      flops    += kernel.flops(N_iii, N_jjj);
      bytes    += kernel.bytes(N_iii, N_jjj);

      START_TIMER("cache copy");
      thread_for_collapse( 5,e,next_,{
	  for(int s =0;s< nstr_;s++)
	  for(int t =0;t< nt_;t++)
	  for(int iii=0;iii< N_iii;iii++)
	  for(int jjj=0;jjj< N_jjj;jjj++)
	  {
	    mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	  }
	});
      STOP_TIMER("cache copy");
    }

    // perf
    LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		 << " Gflop/s/node " << std::endl;
    LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		 << " GB/s/node "  << std::endl;

    // IO
    double       blockSize, ioTime;
    unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

    LOG(Message) << "Writing block to disk" << std::endl;
    ioTime = -GET_TIMER("IO: write block");
    START_TIMER("IO: total");
    makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
    grid_->Barrier();
    // make task list for current node
    nodeIo_.clear();
    for(int f = myRank; f < next_*nstr_; f += nRank)
    {
      IoHelper h;

      h.i  = i;
      h.j  = j;
      h.e  = f/nstr_;
      h.s  = f % nstr_;
      h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
			      ionameFn(h.s), nt_, N_i, N_j);
      h.md = metadataFn(h.s);
      nodeIo_.push_back(h);
    }
    // parallel IO
    for (auto &h: nodeIo_)
    {
      saveBlock(mBlock, h);
    }
    grid_->Barrier();
#else
    // serial IO, for testing purposes only
    for(int e = 0; e < next_; e++)
    for(int s = 0; s < nstr_; s++)
    {
      IoHelper h;

      h.i  = i;
      h.j  = j;
      h.e  = e;
      h.s  = s;
      h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
			      ionameFn(h.s), nt_, N_i, N_j);
      h.md = metadataFn(h.s);
      saveBlock(mfBlock, h);
    }
#endif
    STOP_TIMER("IO: total");
    blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
    ioTime    += GET_TIMER("IO: write block");
    LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		 << ioTime  << " us (" 
		 << blockSize/ioTime*1.0e6/1024/1024
		 << " MB/s)" << std::endl;
  }// loop for i, j (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute1(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<SpinColourMatrix_v> &mat,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  //int    N_i = loop1.size();
  //int    N_j = loop2.size();
  //assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  //int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  //int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  //typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

#if 0
  SpinColourMatrix_v SCM0 = Zero();
  std::vector<SpinColourMatrix_v> mat(MFrvol,SCM0);
  std::cout << "Allocated " << MFrvol << "loop matrices" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      //std::cout << "AAA " << mat[0]()(0,0)(0,0) << std::endl;
      ////kernel(mat, &loop1[i+ii], &loop2[i+ii], N_iii);
      kernel(mat, &loop1[i+ii], &loop2[j+jj], mes, i+ii, j+jj, N_iii, N_jjj);
      /*               .
	 mat    =     / \     for each spacetime x
	              \ /
		  (a,α)x(β,b)
      */
    }
  }
  STOP_TIMER("kernel loop");
  //std::cout << mat << std::endl;
#endif

  int N_i = left.size();
  int N_j = right.size();
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, mat, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");

    for(int j=0;j<N_j;j+=blockSize_) {
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_) {
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	//kernel(mCacheBlock, mat, &left[i+ii], &right[j+jj], orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	      mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
      saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
    // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }
  }// loop for i, j (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute2(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
	   const std::vector<A2AMatrix<Complex> > &mes, const int &it_mes,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn,
	   const int ngamma)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  //assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iColourMatrix<vector_type> ColourMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  ColourMatrix_v CM0 = Zero();
  std::vector<ColourMatrix_v> cmat(MFrvol*ngamma,CM0);
  std::cout << "Allocated " << MFrvol << "loop scalars" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      kernel(cmat, &loop1[i+ii], &loop2[j+jj], mes, i+ii, j+jj, N_iii, N_jjj);
      /*               .
	  cf    =     / \     for each spacetime x
	              \ /
		       *
		       Γ
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, cmat, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");
    for(int j=0;j<N_j;j+=blockSize_)
    {
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
      {
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	//kernel(mCacheBlock, cmat, &left[i+ii], &right[j+jj], orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	    mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
	saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
      // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.e), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }
  }// loop for i, j (Blocks)
}

template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::execute3(const std::vector<Field> &left,  const std::vector<Field> &right,
	   const std::vector<Field> &loop1, const std::vector<Field> &loop2,
	   const std::vector<A2AMatrix<Complex> > &mes, const int &it_mes,
	   A2AKernel<T, FImpl, Field> &kernel, const FilenameFnM &ionameFn,
	   const FilenameFnM &filenameFn, const MetadataFnM &metadataFn)
{
  //////////////////////////////////////////////////////////////////////////
  // i,j   is first  loop over blockSize_ factors
  // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
  // iii,jjj are loops within cacheBlock
  // Total index is sum of these  i+ii+iii etc...
  //////////////////////////////////////////////////////////////////////////
  int    N_i = loop1.size();
  int    N_j = loop2.size();
  //assert(N_i==N_j);
  double flops, bytes, t_kernel;
  double nodes = grid_->NodeCount();
    
  int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

  // MFrvol = local volume / Nsimd
  int MFrvol = grid_->_rdimensions[0]
    *          grid_->_rdimensions[1]
    *          grid_->_rdimensions[2]
    *          grid_->_rdimensions[3];
  //typedef typename FImpl::ComplexField ComplexField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::vector_type vector_type;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinColourVector<vector_type> SpinColourVector_v;

  SpinMatrix_v SM0 = Zero();
  std::vector<SpinMatrix_v> smat(MFrvol,SM0);
  std::cout << "Allocated " << MFrvol << "loop matrices" << std::endl;
  START_TIMER("kernel loop");
  for(int i=0;i<N_i;i+=blockSize_)
  for(int j=0;j<N_j;j+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    int N_jj = MIN(N_j-j,blockSize_);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
    for(int jj=0;jj<N_jj;jj+=cacheBlockSize_) {
      double t;
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
      kernel(smat, &loop1[i+ii], &loop2[j+jj], mes, i+ii, j+jj, N_iii, N_jjj);
      /*               .
	 mat    =     / \     for each spacetime x
	              \ /
		  (a,α)x(β,a)
      */
    }
  }
  STOP_TIMER("kernel loop");

  N_i = left.size();
  N_j = right.size();
  NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
  NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
  for(int i=0;i<N_i;i+=blockSize_) {
    int N_ii = MIN(N_i-i,blockSize_);
    START_TIMER("kernel left * loop");
    FieldMatrix<SpinColourVector_v> left_loop(nstr_,N_ii,MFrvol);
    for(int ii=0;ii<N_ii;ii+=cacheBlockSize_) {
      int N_iii = MIN(N_ii-ii,cacheBlockSize_);
      kernel(left_loop, smat, &left[i+ii], N_iii, ii);
    }
    STOP_TIMER("kernel left * loop");
    for(int j=0;j<N_j;j+=blockSize_) {
      // Get the W and V vectors for this block^2 set of terms
      int N_ii = MIN(N_i-i,blockSize_);
      int N_jj = MIN(N_j-j,blockSize_);
      A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

      LOG(Message) << "All-to-all matrix block " 
		   << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
		   << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
		   << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
		   << std::endl;
      // Series of cache blocked chunks of the contractions within this block
      flops    = 0.0;
      bytes    = 0.0;
      t_kernel = 0.0;
      for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
      for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
      {
	double t;
	int N_iii = MIN(N_ii-ii,cacheBlockSize_);
	int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
	A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

	START_TIMER("kernel");
	kernel(mCacheBlock, left_loop, &right[j+jj], ii, orthogDim_, t);
	//kernel(mCacheBlock, smat, &left[i+ii], &right[j+jj], orthogDim_, t);
	STOP_TIMER("kernel");
	t_kernel += t;
	flops    += kernel.flops(N_iii, N_jjj);
	bytes    += kernel.bytes(N_iii, N_jjj);

	START_TIMER("cache copy");
	thread_for_collapse( 5,e,next_,{
	    for(int s =0;s< nstr_;s++)
	    for(int t =0;t< nt_;t++)
	    for(int iii=0;iii< N_iii;iii++)
	    for(int jjj=0;jjj< N_jjj;jjj++)
	    {
	      mBlock(e,s,t,ii+iii,jj+jjj) = mCacheBlock(e,s,t,iii,jjj);
	    }
	  });
	STOP_TIMER("cache copy");
      }

      // perf
      LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
		   << " Gflop/s/node " << std::endl;
      LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
		   << " GB/s/node "  << std::endl;

      // IO
      double       blockSize, ioTime;
      unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();

      LOG(Message) << "Writing block to disk" << std::endl;
      ioTime = -GET_TIMER("IO: write block");
      START_TIMER("IO: total");
      makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AM_PARALLEL_IO
      grid_->Barrier();
      // make task list for current node
      nodeIo_.clear();
      for(int f = myRank; f < next_*nstr_; f += nRank)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = f/nstr_;
	h.s  = f % nstr_;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	nodeIo_.push_back(h);
      }
      // parallel IO
      for (auto &h: nodeIo_)
      {
	saveBlock(mBlock, h);
      }
      grid_->Barrier();
#else
    // serial IO, for testing purposes only
      for(int e = 0; e < next_; e++)
      for(int s = 0; s < nstr_; s++)
      {
	IoHelper h;

	h.i  = i;
	h.j  = j;
	h.e  = e;
	h.s  = s;
	h.io = A2AMatrixIo<TIo>(filenameFn(h.s), 
				ionameFn(h.s), nt_, N_i, N_j);
	h.md = metadataFn(h.s);
	saveBlock(mfBlock, h);
      }
#endif
      STOP_TIMER("IO: total");
      blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
      ioTime    += GET_TIMER("IO: write block");
      LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
		   << ioTime  << " us (" 
		   << blockSize/ioTime*1.0e6/1024/1024
		   << " MB/s)" << std::endl;
    }
  }// loop for i, j (Blocks)
}

// I/O handler /////////////////////////////////////////////////////////////////
template <typename T, typename FImpl, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, FImpl, Field, MetadataType, TIo>
::saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h)
{
    if ((h.i == 0) and (h.j == 0))
    {
        START_TIMER("IO: file creation");
        h.io.initFile(h.md, blockSize_);
        STOP_TIMER("IO: file creation");
    }
    START_TIMER("IO: write block");
    h.io.saveBlock(m, h.e, h.s, h.i, h.j);
    STOP_TIMER("IO: write block");
}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

END_HADRONS_NAMESPACE

#endif // A2A_Matrix_hpp_
