#ifndef HJ_ZJUAD_MATRIX_MN_H_
#define HJ_ZJUAD_MATRIX_MN_H_

#include <zjucad/matrix/matrix_expression.h>

namespace zjucad { namespace matrix {

template <typename T, size_t M, size_t N>
class matrix_mn : public matrix_expression<matrix_mn<T, M, N> >
{
public:
  typedef matrix_mn<T, M, N> expression_type;
	typedef zjucad::matrix::size_type size_type;
	typedef T value_type;
	typedef const T& const_reference;
	typedef T& reference;
	typedef const T* const_pointer;
	typedef T* pointer;
  typedef const_pointer const_iterator;
  typedef pointer iterator;

    struct stack_memory {
        T dat_[M*N];

        iterator begin() { return dat_; }
        iterator end() { return dat_ + M*N; }
        const_iterator begin() const { return dat_; }
        const_iterator end() const { return dat_ + M*N; }
    };

    typedef stack_memory raw_data_type;
    typedef const stack_memory const_raw_data_type;

    matrix_mn(){}
	template <typename E>
    matrix_mn(const matrix_expression<E> &e){
		*this = e();
	}

	// size
    size_type size(void) const {return M*N;}
    size_type size(int dim) const {return (dim==1)?M:N;}

    void resize(size_type size) {
    assert(M == size && N == 1);
	}
	void resize(size_type nrows, size_type ncols) {
    assert(M == nrows && N == ncols);
	}

	// element access
    const_reference operator[](idx_type i) const {return dat_.dat_[i];}
    const_reference operator()(idx_type i) const {return dat_.dat_[i];}
    const_reference operator()(idx_type row, idx_type col) const {
        return dat_.dat_[row+col*M];
	}
    reference operator[](idx_type i) {return dat_.dat_[i];}
    reference operator()(idx_type i) {return dat_.dat_[i];}
    reference operator()(idx_type row, idx_type col) {
        return dat_.dat_[row+col*M];
	}

	// iterator access
    const_iterator begin(void) const {return dat_.begin();}
    const_iterator end(void) const {return dat_.end();}
    iterator begin(void) {return dat_.begin();}
    iterator end(void) {return dat_.end();}

	template <typename E>
    const matrix_mn<T, M, N>& operator=(const matrix_expression<E> &e){
		assert(size(1) == e().size(1) || size(2) == e().size(2));
#ifdef __CUDACC__
        gpu_copy(e().begin(), e().end(), begin());
#else
        std::copy(e().begin(), e().end(), begin());
#endif
		return *this;
	}

    const_raw_data_type &data(void) const {return dat_;}
    raw_data_type &data(void) {return dat_;}

	PROXY_ACCESS;
	MATRIX_SELF_OP;

private:
  raw_data_type dat_;
};

}}

#endif
