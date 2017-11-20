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

        ZJUCAD_CPU_GPU_FUNC iterator begin() { return dat_; }
        ZJUCAD_CPU_GPU_FUNC iterator end() { return dat_ + M*N; }
        ZJUCAD_CPU_GPU_FUNC const_iterator begin() const { return dat_; }
        ZJUCAD_CPU_GPU_FUNC const_iterator end() const { return dat_ + M*N; }
    };

    typedef stack_memory raw_data_type;
    typedef const stack_memory const_raw_data_type;

    ZJUCAD_CPU_GPU_FUNC matrix_mn(){}
	template <typename E>
    ZJUCAD_CPU_GPU_FUNC matrix_mn(const matrix_expression<E> &e){
		*this = e();
	}

	// size
    ZJUCAD_CPU_GPU_FUNC size_type size(void) const {return M*N;}
    ZJUCAD_CPU_GPU_FUNC size_type size(int dim) const {return (dim==1)?M:N;}

    void resize(size_type size) {
    assert(M == size && N == 1);
	}
	void resize(size_type nrows, size_type ncols) {
    assert(M == nrows && N == ncols);
	}

	// element access
    ZJUCAD_CPU_GPU_FUNC const_reference operator[](idx_type i) const {return dat_.dat_[i];}
    ZJUCAD_CPU_GPU_FUNC const_reference operator()(idx_type i) const {return dat_.dat_[i];}
    ZJUCAD_CPU_GPU_FUNC const_reference operator()(idx_type row, idx_type col) const {
        return dat_.dat_[row+col*M];
	}
    ZJUCAD_CPU_GPU_FUNC reference operator[](idx_type i) {return dat_.dat_[i];}
    ZJUCAD_CPU_GPU_FUNC reference operator()(idx_type i) {return dat_.dat_[i];}
    ZJUCAD_CPU_GPU_FUNC reference operator()(idx_type row, idx_type col) {
        return dat_.dat_[row+col*M];
	}

	// iterator access
    ZJUCAD_CPU_GPU_FUNC const_iterator begin(void) const {return dat_.begin();}
    ZJUCAD_CPU_GPU_FUNC const_iterator end(void) const {return dat_.end();}
    ZJUCAD_CPU_GPU_FUNC iterator begin(void) {return dat_.begin();}
    ZJUCAD_CPU_GPU_FUNC iterator end(void) {return dat_.end();}

	template <typename E>
    ZJUCAD_CPU_GPU_FUNC const matrix_mn<T, M, N>& operator=(const matrix_expression<E> &e){
		assert(size(1) == e().size(1) || size(2) == e().size(2));
#ifdef __CUDACC__
        gpu_copy(e().begin(), e().end(), begin());
#else
        std::copy(e().begin(), e().end(), begin());
#endif
		return *this;
	}

    ZJUCAD_CPU_GPU_FUNC const_raw_data_type &data(void) const {return dat_;}
    ZJUCAD_CPU_GPU_FUNC raw_data_type &data(void) {return dat_;}

	PROXY_ACCESS;
	MATRIX_SELF_OP;

private:
  raw_data_type dat_;
};

}}

#endif
