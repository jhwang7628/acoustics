#ifndef _MLS_MODE_INTERPOLATOR_HPP
#define _MLS_MODE_INTERPOLATOR_HPP

#include <boost/math/special_functions/factorials.hpp>
#include <Eigen/Dense>
#include "MLSConstantsAndTypes.h"

template<int DIM, typename T, int ROWS, int COLS>
class PolyGen
{
};

template<typename T, int ROWS, int COLS>
class PolyGen<1, T, ROWS, COLS>
{
    public:
        inline void gen(const Eigen::Matrix<T, 1, 1> &p,
                        const Eigen::Matrix<T, 1, 1> &pp,
                        int lookupOrder,
                        int row,
                        Eigen::Matrix<T, ROWS, COLS> &A) const
        {
            A(row, 0) = T(1);

            if(lookupOrder == 0)
            {
                return;
            }

            T diffX = pp(0) - p(0);

            A(row, 1) = diffX;

            if(lookupOrder == 1)
            {
                return;
            }

            A(row, 2) = diffX * diffX;

            if(lookupOrder == 2)
            {
                return;
            }

            T diffXSq = diffX * diffX;

            A(row, 3) = diffX * diffXSq;
        }
};

template<typename T, int ROWS, int COLS>
class PolyGen<2, T, ROWS, COLS>
{
    public:
        inline void gen(const Eigen::Matrix<T, 2, 1> &p,
                        const Eigen::Matrix<T, 2, 1> &pp,
                        int lookupOrder,
                        int row,
                        Eigen::Matrix<T, ROWS, COLS> &A) const
        {
            A(row, 0) = T(1);

            if(lookupOrder == 0)
            {
                return;
            }

            T diffX = pp(0) - p(0);
            T diffY = pp(1) - p(1);

            A(row, 1) = diffX;
            A(row, 2) = diffY;

            if(lookupOrder == 1)
            {
                return;
            }

            T diffXSq = diffX * diffX;
            T diffYSq = diffY * diffY;

            A(row, 3) = diffX * diffY;
            A(row, 4) = diffXSq;
            A(row, 5) = diffYSq;

            if(lookupOrder == 2)
            {
                return;
            }

            A(row, 6) = diffY * diffXSq;
            A(row, 7) = diffX * diffYSq;
            A(row, 8) = diffX * diffXSq;
            A(row, 9) = diffY * diffYSq;
        }
};

template<typename T, int ROWS, int COLS>
class PolyGen<3, T, ROWS, COLS>
{
    public:
        inline void gen(const Eigen::Matrix<T, 3, 1> &p,
                        const Eigen::Matrix<T, 3, 1> &pp,
                        int lookupOrder,
                        int row,
                        Eigen::Matrix<T, ROWS, COLS> &A) const
        {
            A(row, 0) = T(1);

            if(lookupOrder == 0)
            {
                return;
            }

            T diffX = pp(0) - p(0);
            T diffY = pp(1) - p(1);
            T diffZ = pp(2) - p(2);

            A(row, 1) = diffX;
            A(row, 2) = diffY;
            A(row, 3) = diffZ;

            if(lookupOrder == 1)
            {
                return;
            }

            T diffXSq = diffX * diffX;
            T diffYSq = diffY * diffY;
            T diffZSq = diffZ * diffZ;

            A(row, 4) = diffX * diffY;
            A(row, 5) = diffX * diffZ;
            A(row, 6) = diffY * diffZ;
            A(row, 7) = diffXSq;
            A(row, 8) = diffYSq;
            A(row, 9) = diffZSq;

            if(lookupOrder == 2)
            {
                return;
            }

            A(row, 10) = diffY * diffXSq;
            A(row, 11) = diffZ * diffXSq;
            A(row, 12) = diffX * diffYSq;
            A(row, 13) = diffX * diffZSq;
            A(row, 14) = diffZ * diffYSq;
            A(row, 15) = diffY * diffZSq;
            A(row, 16) = diffX * diffY * diffZ;
            A(row, 17) = diffXSq * diffX;
            A(row, 18) = diffYSq * diffY;
            A(row, 19) = diffZSq * diffZ;
        }
};

/**
 * Adaptation of the original MLSMode class, templated
 * and updated to use Eigen library. Needed to make it usable with
 * automatic differentiation of ceres library.
 */

/**
 * Class to store a downsampled version of a mode shape.
 * It interpolates the samples to upsample to the full resolution
 * using Moving Least Squares. A good reference is Nealen 2004.
 * Details of the equation we are solving and our implementation are
 * given in the lookup function.
 */
template<typename T, int POINT_DIM = 3, int VAL_DIM = 3, int ROWS = Eigen::Dynamic, int COLS = Eigen::Dynamic>
class MLSModeInterpolator
{
    public:
        typedef Eigen::Matrix<T, POINT_DIM, 1> MLSPoint;
        typedef Eigen::Matrix<T, VAL_DIM, 1> MLSVal;
        typedef Eigen::Matrix<T, ROWS, COLS> MLSMatrix;
        //typedef Eigen::Matrix<T, 1, COLS> MLSMatrixRow;
        typedef Eigen::Matrix<T, ROWS, 1> MLSVector;
        typedef Eigen::Matrix<T, VAL_DIM, COLS> MLSPolyResult;

        MLSModeInterpolator()
            : m_fleishmanRad(MLS_FLEISHMAN_RAD),
              m_compactRad(MLS_COMPACT_RAD),
              m_lookupType(MLS_LOOKUP_TYPE),
              m_lookupOrder(MLS_LOOKUP_ORDER)
        {
        }

        // Copy constructor
        MLSModeInterpolator(const MLSModeInterpolator &rhs)
            : m_fleishmanRad(rhs.m_fleishmanRad),
              m_compactRad(rhs.m_compactRad),
              m_lookupType(rhs.m_lookupType),
              m_lookupOrder(rhs.m_lookupOrder)
        {
        }

        MLSModeInterpolator & operator= (const MLSModeInterpolator &rhs)
        {
            m_fleishmanRad = rhs.m_fleishmanRad;
            m_compactRad = rhs.m_compactRad;
            m_lookupType = rhs.m_lookupType;
            m_lookupOrder = rhs.m_lookupOrder;

            return *this;
        }

        template<typename Y, int RP_DIM, int RV_DIM, int RROWS, int RCOLS>
        void copyParams(const MLSModeInterpolator<Y, RP_DIM, RV_DIM, RROWS, RCOLS> &rhs)
        {
            setFleishmanRadius(T(rhs.getFleishmanRadius()));
            setCompactRadius(T(rhs.getCompactRadius()));
            setLookupType(rhs.getLookupType());
            setLookupOrder(rhs.getLookupOrder());
        }

        /// Getters
        T getFleishmanRadius() const
        {
            return m_fleishmanRad;
        }

        T getCompactRadius() const
        {
            return m_compactRad;
        }

        MLS_lookup_t getLookupType() const
        {
            return m_lookupType;
        }

        int getLookupOrder() const
        {
            return m_lookupOrder;
        }

        /// Setters
        void setFleishmanRadius(T radius)
        {
            m_fleishmanRad = radius;
        }

        void setCompactRadius(T radius)
        {
            m_compactRad = radius;
        }

        void setLookupType(MLS_lookup_t lt)
        {
            m_lookupType = lt;
        }

        void setLookupOrder(int lo)
        {
            m_lookupOrder = lo;
        }

        /// converts the lookup order to the number of terms
        int orderToSize() const
        {
            using boost::math::unchecked_factorial;

            if (m_lookupOrder > 3)
            {
                std::ostringstream os;
                os << "unknown lookup order " << m_lookupOrder;
                throw std::runtime_error(os.str());
            }

            // (d + m)! / (d! m!)
            return static_cast<int>(unchecked_factorial<double>(POINT_DIM + m_lookupOrder) /
                    (unchecked_factorial<double>(POINT_DIM) * unchecked_factorial<double>(m_lookupOrder)));
        }

        /**
         * Solves for modal displacement at surface point p.
         * Comments refer to ASAP Intro to LS/WLS/MLS by Andrew Nealen.
         * Generally, this is the operation in equation 5.
         */
        template<class P_ALLOCATOR, class V_ALLOCATOR>
        MLSVal
        lookup(const MLSPoint &p,
               const std::vector<MLSPoint, P_ALLOCATOR> &points,
               const std::vector<MLSVal, V_ALLOCATOR> &attributes,
               T knnRadius = T(-1),
               T *condNum = NULL,
               const std::vector<int> *indices = NULL) const
        {
            // Polynomial degree (m) -> num coefficients (k).
            // k = (d+m)!/m!/d! (d = num degrees)
            int size = orderToSize();

            T radius;
            if (FLEISHMAN == m_lookupType)
            {
                radius = m_fleishmanRad;
            }
            else if (COMPACT == m_lookupType)
            {
                radius = m_compactRad;
            }
            else
            {
                throw std::runtime_error("unknown lookup type");
            }

            if (knnRadius > 0)
            {
                radius = knnRadius / T(3);
            }

            int numPts = indices ? indices->size() : points.size();

            // We are solving A c = b. A is the least squares matrix of all the samples,
            // c is unknown, and b is the sample values.
            // This is overconstrained (A is tall and skinny). We currently solve it with
            // a QR decomp of A.
            // This needs to be done 3 times (mode values are vectors), but the only part
            // which changes each time is b, so we precompute the factorization of A.

            // Initialize the weights
            MLSVector weights(numPts);
            initializeWeights(p,
                              points,
                              radius * radius,
                              weights,
                              indices);

            // Form the A matrix. The ith 1 x k row vector contains evaluated
            // polynomial terms (such as 1, x_i, y_i, z_i...) weighted by the
            // ith point's distance from p. i=0..n-1 where n is the number of
            // sampled points.
            MLSMatrix A(numPts,
                        size);

            assignLsqrMatrix(A,
                             p,
                             points,
                             indices);

            // Multiply each row of A by the corresponding weight
            A.array().colwise() *= weights.array();

            if (condNum)
            {
                // Find the condition number of A
                Eigen::JacobiSVD<MLSMatrix> svd(A,
                                                Eigen::ComputeThinU | Eigen::ComputeThinV);

                if (svd.nonzeroSingularValues() < A.cols())
                {
                    // Singular, infinite condition number
                    throw std::runtime_error("infinite condition number");
                }
                else
                {
                    *condNum = svd.singularValues()(0) / svd.singularValues()(A.cols() - 1);
                }
            }

            // Init b, c
            //MLSVector b(numPts);
            Eigen::Matrix<T, COLS, 1> c(size);

            Eigen::Matrix<T, ROWS, VAL_DIM> b(numPts, VAL_DIM);

            assignAttributes(b,
                             points,
                             attributes,
                             p,
                             indices);

            b.array().colwise() *= weights.array();

            //return A.fullPivLu().solve((b.array().colwise() * weights.array()).matrix()).row(0);

            MLSVal output;

            // Full pivot LU always works
            // A^T A is always positive semidefinite, so the LDLT (robust Cholesky) works
            // Doesn't seem to be much speed difference between them
            // According to the Eigen documentation, the full pivot LU is
            // more accurate
            //
            // QR'ing the A matrix for now, avoids squaring the condition number.
            //
            // Could also do a truncated SVD solve on A, but I was worried about the
            // robustness when using auto diff.
            // SVD needs to use sqrt I think, which is a problem if we ever take sqrt(0).
            //Eigen::FullPivLU<MLSMatrix> decomp((A.transpose() * A));
            //Eigen::LDLT<MLSMatrix> decomp((A.transpose() * A));
            //Eigen::ColPivHouseholderQR<MLSMatrix> decomp(A);
            //Eigen::FullPivHouseholderQR<MLSMatrix> decomp(A);
            //Eigen::HouseholderQR<MLSMatrix> decomp(A);
            //Eigen::FullPivLU<MLSMatrix> decomp(A);
            //Eigen::JacobiSVD<MLSMatrix> decomp(A,
                                               //Eigen::ComputeThinU | Eigen::ComputeThinV);
            //Eigen::JacobiSVD<MLSMatrix> decomp(A,
                                               //Eigen::ComputeFullU | Eigen::ComputeFullV);
           Eigen::JacobiSVD<MLSMatrix, Eigen::FullPivHouseholderQRPreconditioner> decomp(A,
                                               Eigen::ComputeFullU | Eigen::ComputeFullV);

            for(int i = 0; i < VAL_DIM; ++i)
            {
                // Assign b, the n x 1 RHS vector.
                // These are the function values. There's three of them but
                // they get evaluated one at a time.
                //assignAttributes(b,
                                 //i,
                                 //points,
                                 //attributes,
                                 //p,
                                 //indices);

                //// Multiply each entry of b by the corresponding weight
                //b.array() *= weights.array();

                // Solve the linear system to get the coefficients.
                //c = decomp.solve(A.transpose() * b);
                c = decomp.solve(b.col(i));

                // Get the modal displacement coordinate. Since we assigned A
                // in a coordinate system translated so that point p is at the origin,
                // x=y=z=0 - in other words only c[0], the constant, is non-zero.
                output(i) = c(0);
            }

            return output;
        }

        template<class P_ALLOCATOR, class V_ALLOCATOR>
        MLSPolyResult
        lookupFullPolynomial(const MLSPoint &p,
                             const std::vector<MLSPoint, P_ALLOCATOR> &points,
                             const std::vector<MLSVal, V_ALLOCATOR> &attributes,
                             T knnRadius = T(-1),
                             T *condNum = NULL,
                             const std::vector<int> *indices = NULL) const
        {
            // Polynomial degree (m) -> num coefficients (k).
            // k = (d+m)!/m!/d! (d = num degrees, 3)
            int size = orderToSize();

            T radius;
            if (FLEISHMAN == m_lookupType)
            {
                radius = m_fleishmanRad;
            }
            else if (COMPACT == m_lookupType)
            {
                radius = m_compactRad;
            }
            else
            {
                throw std::runtime_error("unknown lookup type");
            }

            if (knnRadius > 0)
            {
                radius = knnRadius / T(3);
            }

            int numPts = indices ? indices->size() : points.size();

            // We are solving A c = b. A is the least squares matrix of all the samples,
            // c is unknown, and b is the sample values.
            // This is overconstrained (A is tall and skinny). We currently solve it with
            // a QR decomp of A.
            // This needs to be done 3 times (mode values are vectors), but the only part
            // which changes each time is b, so we precompute the factorization of A.

            // Initialize the weights
            MLSVector weights(numPts);
            initializeWeights(p,
                              points,
                              radius * radius,
                              weights,
                              indices);

            // Form the A matrix. The ith 1 x k row vector contains evaluated
            // polynomial terms (such as 1, x_i, y_i, z_i...) weighted by the
            // ith point's distance from p. i=0..n-1 where n is the number of
            // sampled points.
            MLSMatrix A(numPts,
                        size);

            assignLsqrMatrix(A,
                             p,
                             points,
                             indices);

            // Multiply each row of A by the corresponding weight
            A.array().colwise() *= weights.array();

            if (condNum)
            {
                // Find the condition number of A
                Eigen::JacobiSVD<MLSMatrix> svd(A,
                                                Eigen::ComputeThinU | Eigen::ComputeThinV);

                if (svd.nonzeroSingularValues() < A.cols())
                {
                    // Singular, infinite condition number
                    throw std::runtime_error("infinite condition number");
                }
                else
                {
                    *condNum = svd.singularValues()(0) / svd.singularValues()(A.cols() - 1);
                }
            }

            // Init b, c
            MLSVector b(numPts);
            Eigen::Matrix<T, COLS, 1> c(size);

            MLSPolyResult output(VAL_DIM, size);

            // Full pivot LU always works
            // A^T A is always positive semidefinite, so the LDLT (robust Cholesky) works
            // Doesn't seem to be much speed difference between them
            // According to the Eigen documentation, the full pivot LU is
            // more accurate
            //
            // QR'ing the A matrix for now, avoids squaring the condition number.
            //
            // Could also do a truncated SVD solve on A, but I was worried about the
            // robustness when using auto diff.
            // SVD needs to use sqrt I think, which is a problem if we ever take sqrt(0).
            //Eigen::FullPivLU<MLSMatrix> decomp((A.transpose() * A));
            //Eigen::LDLT<MLSMatrix> decomp((A.transpose() * A));
            //Eigen::ColPivHouseholderQR<MLSMatrix> decomp(A);
            Eigen::HouseholderQR<MLSMatrix> decomp(A);
            //Eigen::FullPivLU<MLSMatrix> decomp(A);
            //Eigen::JacobiSVD<MLSMatrix> decomp(A,
                                               //Eigen::ComputeThinU | Eigen::ComputeThinV);

            for(int i = 0; i < VAL_DIM; ++i)
            {
                // Assign b, the n x 1 RHS vector.
                // These are the function values. There's three of them but
                // they get evaluated one at a time.
                assignAttributes(b,
                                 i,
                                 points,
                                 attributes,
                                 p,
                                 indices);

                // Multiply each entry of b by the corresponding weight
                b.array() *= weights.array();

                // Solve the linear system to get the coefficients.
                //c = decomp.solve(A.transpose() * b);
                c = decomp.solve(b);

                // Get the modal displacement coordinate. Since we assigned A
                // in a coordinate system translated so that point p is at the origin,
                // x=y=z=0 - in other words only c[0], the constant, is non-zero.
                output.row(i) = c;
            }

            return output;
        }

        //template <typename Y> friend class MLSModeInterpolator;

    protected:
        /// lookup radius, if relevant
        T m_fleishmanRad;
        T m_compactRad;

        /// lookup type (changes the distance function)
        MLS_lookup_t m_lookupType;

        /// order of the polynomial fit
        int m_lookupOrder;

        /// Generates polynomial basis functions of points
        PolyGen<POINT_DIM, T, ROWS, COLS> m_polyGenerator;

        /**
         * Computes the distance weight, based on the lookup type.
         * Returns the sqrt of the distance weight.
         */
        T
        distanceWeight(T dSquared,
                       T radiusSquared) const
        {
            switch(m_lookupType)
            {
                case FLEISHMAN:
                {
                    // The division by 2 returns the square root of the distance weight
                    //T radius = m_fleishmanRad;
                    return exp(-(dSquared) / (2. * radiusSquared));
                }

                case COMPACT:
                {
                    //T radius = m_compactRad;
                    throw std::runtime_error("unstable due to sqrt");
                    T distance = sqrt(dSquared);
                    T radius = sqrt(radiusSquared);
                    if(distance <= radius)
                    {
                        return sqrt(pow((T(1) - distance / radius), T(4))
                               * (T(4) * distance / radius + T(1)));
                    }
                    else
                    {
                        return T(0.0);
                    }
                }

                // I don't think this is reachable, but it stops compiler warnings
                default:
                {
                    throw std::runtime_error("Unknown lookup type");
                }
            }
        }

        /// fills in the A matrix
        template<class ALLOCATOR>
        void assignLsqrMatrix(MLSMatrix &A,
                              const MLSPoint &p,
                              const std::vector<MLSPoint, ALLOCATOR> &points,
                              const std::vector<int> *indices = NULL) const
        {
            assert(static_cast<size_t>(A.rows()) == (indices ? indices->size() : points.size())
                   && A.cols() == orderToSize());

            if (indices)
            {
                for(size_t i = 0; i < indices->size(); i++)
                {
                    m_polyGenerator.gen(p,
                                        points[(*indices)[i]],
                                        m_lookupOrder,
                                        i,
                                        A);
                }
            }
            else
            {
                for(size_t i = 0; i < points.size(); i++)
                {
                    m_polyGenerator.gen(p,
                                        points[i],
                                        m_lookupOrder,
                                        i,
                                        A);
                }
            }
        }

        /// assigns the b vector
        template<class P_ALLOCATOR, class V_ALLOCATOR>
        void assignAttributes(MLSVector &b,
                              int k,
                              const std::vector<MLSPoint, P_ALLOCATOR> &points,
                              const std::vector<MLSVal, V_ALLOCATOR> &attributes,
                              MLSPoint p,
                              const std::vector<int> *indices = NULL) const
        {
            assert(static_cast<size_t>(b.rows()) == (indices ? indices->size() : points.size()));

            if (indices)
            {
                for(size_t i = 0; i < indices->size(); i++)
                {
                    b(i) = T(attributes[(*indices)[i]](k));
                }
            }
            else
            {
                for(size_t i = 0; i < points.size(); i++)
                {
                    b(i) = T(attributes[i](k));
                }
            }
        }

        template<class P_ALLOCATOR, class V_ALLOCATOR>
        void assignAttributes(Eigen::Matrix<T, ROWS, VAL_DIM> &b,
                              const std::vector<MLSPoint, P_ALLOCATOR> &points,
                              const std::vector<MLSVal, V_ALLOCATOR> &attributes,
                              MLSPoint /*p*/,
                              const std::vector<int> *indices = NULL) const
        {
            assert(static_cast<size_t>(b.rows()) == (indices ? indices->size() : points.size()));

            if (indices)
            {
                for(size_t i = 0; i < indices->size(); i++)
                {
                    b.row(i) = (attributes[(*indices)[i]]);
                }
            }
            else
            {
                for(size_t i = 0; i < points.size(); i++)
                {
                    b.row(i) = (attributes[i]);
                }
            }
        }

        template<class ALLOCATOR>
        void initializeWeights(const MLSPoint &p,
                               const std::vector<MLSPoint, ALLOCATOR> &points,
                               T radiusSquared,
                               MLSVector &weights,
                               const std::vector<int> *indices = NULL) const
        {
            int sz = indices ? indices->size() : points.size();

            assert(weights.rows() == sz);

            if (indices)
            {
                for (int i = 0; i < sz; ++i)
                {
                    weights(i) = distanceWeight((p - points[(*indices)[i]]).squaredNorm(),
                                                radiusSquared);
                }
            }
            else
            {
                for (int i = 0; i < sz; ++i)
                {
                    weights(i) = distanceWeight((p - points[i]).squaredNorm(),
                                                radiusSquared);
                }
            }
        }
};

#endif // _MLS_MODE_INTERPOLATOR_HPP

