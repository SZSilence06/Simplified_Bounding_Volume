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
    typedef pointer raw_data_type;
    typedef const_pointer const_raw_data_type;

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
	const_reference operator[](idx_type i) const {return dat_[i];}
	const_reference operator()(idx_type i) const {return dat_[i];}
	const_reference operator()(idx_type row, idx_type col) const {
		return dat_[row+col*M];
	}
	reference operator[](idx_type i) {return dat_[i];}
	reference operator()(idx_type i) {return dat_[i];}
	reference operator()(idx_type row, idx_type col) {
		return dat_[row+col*M];
	}

	// iterator access
	const_iterator begin(void) const {return dat_;}
	const_iterator end(void) const {return dat_+M*N;}
	iterator begin(void) {return dat_;}
	iterator end(void) {return dat_+M*N;}

	template <typename E>
	const matrix_mn<T, M, N>& operator=(const matrix_expression<E> &e){
		assert(size(1) == e().size(1) || size(2) == e().size(2));
		std::copy(e().begin(), e().end(), begin());
		return *this;
	}

    const_raw_data_type data(void) const {return dat_;}
    raw_data_type data(void) {return dat_;}

	PROXY_ACCESS;
	MATRIX_SELF_OP;

private:
  T dat_[M*N];
};

}}

#endif
