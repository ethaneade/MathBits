#ifndef MATHBITS_SPARSEMATRIX_HPP
#define MATHBITS_SPARSEMATRIX_HPP

#include <vector>
#include <algorithm>
#include <iostream>

namespace mathbits {

    class SparseMatrix
    {
    private:
        struct ColEntry
        {
            unsigned int row;
            double value;
            ColEntry() {}
            ColEntry(unsigned int r, double v) : row(r), value(v) {}
            bool operator<(const ColEntry& e) const { return e.row < row; }
        };

        struct RowEntry
        {
            unsigned int col;
            unsigned int index;
            RowEntry() {}
            RowEntry(unsigned int c, unsigned int i) : col(c), index(i) {}
            bool operator<(const RowEntry& e) const { return e.col < col; }
        };

        typedef std::vector<ColEntry> Column;
        typedef std::vector<RowEntry> Row;
    
        std::vector<Column> cols;
        std::vector<Row> rows;
        size_t num_entries;

    public:
        SparseMatrix() : num_entries(0) {}

        void swap(SparseMatrix& m) {
            cols.swap(m.cols);
            rows.swap(m.rows);
            std::swap(num_entries, m.num_entries);
        }
    
    
        size_t num_rows() const { return rows.size(); }
        size_t num_cols() const { return cols.size(); }
        size_t fill_count() const { return num_entries; }
    
        class row_iterator;
        class const_row_iterator;
        class col_iterator;
        class const_col_iterator;
    
        class row_iterator
        {
        private:
            SparseMatrix *m;
            Row::iterator it;
            friend class SparseMatrix;
            friend class const_row_iterator;
        public:        
            row_iterator& operator++() { ++it; return *this; }
            row_iterator operator++(int) { row_iterator rit = *this; ++rit.it; return rit; }
            row_iterator& operator--() { --it; return *this; }
            row_iterator operator--(int) { row_iterator rit = *this; --rit.it; return rit; }

            bool operator==(const row_iterator& other) const { return it == other.it; }
            bool operator!=(const row_iterator& other) const { return !(*this == other); }

            unsigned int col() const { return it->col; }
            double& operator*() { return m->cols[it->col][it->index].value; }
        };

        class const_row_iterator
        {
        private:
            const SparseMatrix *m;
            Row::const_iterator it;        
            friend class SparseMatrix;
        public:
            const_row_iterator() {}
            const_row_iterator(const row_iterator& rit) : m(rit.m), it(rit.it) {}

            const_row_iterator& operator++() { ++it; return *this; }
            const_row_iterator operator++(int) { const_row_iterator rit = *this; ++rit.it; return rit; }
            const_row_iterator& operator--() { --it; return *this; }
            const_row_iterator operator--(int) { const_row_iterator rit = *this; --rit.it; return rit; }

            bool operator==(const const_row_iterator& other) const { return it == other.it; }
            bool operator!=(const const_row_iterator& other) const { return !(*this == other); }
        
            unsigned int col() const { return it->col; }
            double operator*() const { return m->cols[it->col][it->index].value; }
        };

        class col_iterator
        {
        private:
            Column::iterator it;
            friend class SparseMatrix;
            friend class const_col_iterator;
        public:
            col_iterator& operator++() { ++it; return *this; }
            col_iterator operator++(int) { col_iterator cit = *this; ++cit.it; return cit; }
            col_iterator& operator--() { --it; return *this; }
            col_iterator operator--(int) { col_iterator cit = *this; --cit.it; return cit; }

            bool operator==(const col_iterator& other) const { return it == other.it; }
            bool operator!=(const col_iterator& other) const { return !(*this == other); }
        
            unsigned int row() const { return it->row; }
            double& operator*() { return it->value; }
        };

        class const_col_iterator
        {
        private:
            SparseMatrix *m;
            Column::const_iterator it;
            friend class SparseMatrix;
        public:
            const_col_iterator() {}
            const_col_iterator(const col_iterator& cit) : it(cit.it) {}
        
            const_col_iterator& operator++() { ++it; return *this; }
            const_col_iterator operator++(int) { const_col_iterator cit = *this; ++cit.it; return cit; }
            const_col_iterator& operator--() { --it; return *this; }
            const_col_iterator operator--(int) { const_col_iterator cit = *this; --cit.it; return cit; }

            bool operator==(const const_col_iterator& other) const { return it == other.it; }
            bool operator!=(const const_col_iterator& other) const { return !(*this == other); }
        
            unsigned int row() const { return it->row; }
            double operator*() const { return it->value; }
        };

        row_iterator row_begin(unsigned int i) {
            row_iterator rit;
            rit.m = this;
            rit.it = rows[i].begin();
            return rit;
        }
    
        const_row_iterator row_begin(unsigned int i) const {
            const_row_iterator rit;
            rit.m = this;
            rit.it = rows[i].begin();
            return rit;
        }
    
        row_iterator row_end(unsigned int i) {
            row_iterator rit;
            rit.m = this;
            rit.it = rows[i].end();
            return rit;
        }
    
        const_row_iterator row_end(unsigned int i) const {
            const_row_iterator rit;
            rit.m = this;
            rit.it = rows[i].end();
            return rit;
        }

        col_iterator col_begin(unsigned int i) {
            col_iterator cit;
            cit.it = cols[i].begin();
            return cit;
        }
    
        const_col_iterator col_begin(unsigned int i) const {
            const_col_iterator cit;
            cit.it = cols[i].begin();
            return cit;
        }
    
        col_iterator col_end(unsigned int i) {
            col_iterator cit;
            cit.it = cols[i].end();
            return cit;
        }
    
        const_col_iterator col_end(unsigned int i) const {
            const_col_iterator cit;
            cit.it = cols[i].end();
            return cit;
        }

        void reserve(unsigned int m, unsigned int n)
        {
            if (rows.size() < m)
                rows.resize(m);
            if (cols.size() < n)
                cols.resize(n);
        }
    
        double get(unsigned int i, unsigned int j) const
        {
            if (j >= cols.size() || i >= rows.size())
                return 0;
            if (rows[i].size() < cols[j].size()) {
                for (size_t k=0; k<rows[i].size(); ++k)
                    if (rows[i][k].col == j)
                        return cols[j][rows[i][k].index].value;
            } else {
                for (size_t k=0; k<cols[j].size(); ++k)
                    if (cols[j][k].row == i)
                        return cols[j][k].value;
            }
            return 0;
        }
    
        double& operator()(unsigned int i, unsigned int j)
        {
            reserve(i+1, j+1);
        
            if (rows[i].size() < cols[j].size()) {            
                for (Row::iterator it = rows[i].begin(); it != rows[i].end(); ++it)
                    if (it->col == j)
                        return cols[j][it->index].value;
            } else {
                for (Column::iterator it = cols[j].begin(); it != cols[j].end(); ++it)
                    if (it->row == i)
                        return it->value;
            }
            rows[i].push_back(RowEntry(j,cols[j].size()));
            cols[j].push_back(ColEntry(i,0));
            ++num_entries;
            return cols[j].back().value;
        }

        void insert(unsigned int i, unsigned int j, double value)
        {
            reserve(i+1,j+1);
            rows[i].push_back(RowEntry(j,cols[j].size()));
            cols[j].push_back(ColEntry(i,value));
            ++num_entries;
        }

        void sort_entries()
        {
            for (size_t i=0; i<rows.size(); ++i)
                rows[i].clear();
        
            for (size_t i=0; i<cols.size(); ++i)
            {
                std::sort(cols[i].begin(), cols[i].end());
                for (size_t j=0; j<cols[i].size(); ++j)
                {
                    unsigned int r = cols[i][j].row;
                    assert(r < rows.size());
                    rows[r].push_back(RowEntry(i, j));
                }
            }
        }

        void print_sparse(std::ostream& out) const
        {
            for (size_t i=0; i<rows.size(); ++i)
            {
                if (rows[i].empty())
                    continue;
                out.width(6);
                out << i << ":: ";
                for (size_t j=0; j<rows[i].size(); ++j)
                {
                    unsigned int c = rows[i][j].col;
                    unsigned int ix = rows[i][j].index;
                    out << "(" << c << ": ";
                    assert(cols[c][ix].row == i);
                    double v = cols[c][ix].value;
                    out << v << ") ";
                }
                out << std::endl;
            }
        }
    
        void print_dense(std::ostream& out) const
        {
            int w = (int)out.precision() + 7;
            for (size_t i=0; i<rows.size(); ++i) {
                for (size_t j=0; j<cols.size(); ++j) {
                    out.width(w);
                    out << get(i,j);
                }
                out << std::endl;
            }
        }
    };

    void multiply(const SparseMatrix& m, const std::vector<double>& x,
                  std::vector<double>& y)
    {
        const size_t n = m.num_cols();
        assert(x.size() == n);
        y.resize(m.num_rows());

        for (size_t i=0; i<y.size(); ++i)
        {
            double sum = 0;
            for (SparseMatrix::const_row_iterator it = m.row_begin(i);
                 it != m.row_end(i);
                 ++it)
            {
                sum += *it * x[it.col()];
            }
            y[i] = sum;
        }
    }
    
}

#endif
