#ifndef MATHBITS_SPARSELU_HPP
#define MATHBITS_SPARSELU_HPP

#include <mathbits/sparsematrix.hpp>
#include <cmath>

namespace mathbits {

    struct SparseLU
    {
        SparseMatrix lu;
        std::vector<double> inv_diag;
        std::vector<unsigned int> src_row, dst_row;

        void compute(const SparseMatrix& m)
        {
            size_t n = m.num_cols();
            lu = m;
            inv_diag.resize(n);
            src_row.resize(n);
            dst_row.resize(n);
    
            for (size_t i=0; i<n; ++i)
                dst_row[i] = src_row[i] = i;

            int swaps = 0;
            unsigned int row = 0;
            for (size_t i=0; i<n; ++i) {
                SparseMatrix::col_iterator argmax;
                double maxval = 0;
                for (SparseMatrix::col_iterator it = lu.col_begin(i); it != lu.col_end(i); ++it)
                {
                    unsigned int r = dst_row[it.row()];
                    if (r < row)
                        continue;
            
                    double ji = fabs(*it);
                    if (ji > maxval) {
                        maxval = ji;
                        argmax = it;
                    }
                }
        
                if (maxval < 1e-15) {
                    inv_diag[i] = 0;
                    continue;
                }

                unsigned int max_row = dst_row[argmax.row()];
                if (max_row != row) {
                    unsigned int src_lo = src_row[row];
                    
                    dst_row[src_lo] = max_row;
                    src_row[max_row] = src_lo;

                    dst_row[argmax.row()] = row;                    
                    src_row[row] = argmax.row();
                    ++swaps;
                }
                        
                double inv_ii = 1 / *argmax;
                inv_diag[i] = inv_ii;

                unsigned sr = argmax.row();
                for (SparseMatrix::col_iterator it = lu.col_begin(i); it != lu.col_end(i); ++it)
                {
                    unsigned int rj = it.row();
                    unsigned int r = dst_row[rj];
                    if (r <= row)
                        continue;
                    double factor = *it * inv_ii;            
                    for (SparseMatrix::row_iterator rit = lu.row_begin(sr); rit != lu.row_end(sr); ++rit) {
                        unsigned int ck = rit.col();
                        if (ck <= i)
                            continue;
                        lu(rj, ck) -= factor * *rit;
                    }
                    *it = factor;
                }
                ++row;                                               
            }
        }

        void inverse_times(const std::vector<double>& v,
                           std::vector<double>& x) const
        {
            const size_t n = inv_diag.size();
            assert(v.size() == n);
            x.resize(n);
    
            for (size_t i=0; i<n; ++i) {
                unsigned int r = src_row[i];
                double yi = v[src_row[i]];
                for (SparseMatrix::const_row_iterator it = lu.row_begin(r);
                     it != lu.row_end(r);
                     ++it)
                {
                    if (it.col() >= i)
                        continue;
                    yi -= *it * x[it.col()];
                }
                x[i] = yi;
            }
            for (int i=n-1; i>=0; --i) {
                double yi = x[i];
                unsigned int r = src_row[i];
                for (SparseMatrix::const_row_iterator it = lu.row_begin(r);
                     it != lu.row_end(r);
                     ++it)
                {
                    if ((int)it.col() <= i)
                        continue;
                    yi -= *it * x[it.col()];
                }
                x[i] = inv_diag[i] * yi;
            }
        }    
    };
    
    
}

#endif
