#include <stdlib.h>
#include <vector>
using namespace std;

//(done)убрать &
//(fixed)обработка ошибок
//(fixed)разные типы
//(fixed)return types

class temp_matrix;

template <typename T = double>
class matrix {
	size_t _nrow;
	size_t _ncol;
	typedef vector<T> _Trow;
	vector<_Trow> _matr;

	public:
		matrix();
		matrix(size_t);
		matrix(size_t, size_t);		
		template <typename U> matrix(const matrix<U> &);
		matrix(temp_matrix &);
//		~matrix();

		template <typename U> matrix & operator = (const matrix<U> &);
		void scan ();
		void print () const;

		size_t get_row_num () const;
		size_t get_col_num () const;
		const _Trow operator [] (const size_t id) const { return _matr[id]; }

		matrix operator * (T) const;

		matrix operator + (const matrix &) const;
		matrix operator - (const matrix &) const;
		matrix operator * (const matrix &) const;

		T trace() const;
		T det() const;

		matrix transp() const;
		matrix<double> inv() const;

	private:


};


class temp_matrix {
	public:
		size_t _nrow;
		size_t _ncol;
		typedef vector<double> _Trow;
		vector<_Trow> _matr;

		temp_matrix(size_t nrow, size_t ncol) {
			_nrow = nrow;
			_ncol = ncol;
			_matr.resize(_nrow);
			for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);
		}

		template <typename U> temp_matrix(const matrix<U> & m) {
			_nrow = m.get_row_num();
			_ncol = m.get_col_num();

			_matr.resize(_nrow);
			for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);

			for (size_t i = 0; i < _nrow; i++) {
				for (size_t j = 0; j < _ncol; j++) {
					_matr[i][j] = m[i][j];
				}
			}
		}

		void strswap(size_t x, size_t y)
		{
			_Trow temp = _matr[x];
			_matr[x] = _matr[y];
			_matr[y] = temp;
		}

		void strmult(size_t x, double p)
		{
			for (size_t i = 0; i < _ncol; i++) _matr[x][i] *= p;
		}

		void strsub(size_t x, size_t y, double p)
		{
			for (size_t i = 0; i < _ncol; i++) _matr[x][i] -= _matr[y][i] * p;
		}
};

//-----constructors---------------------------------------------------------------

template <typename T>
matrix<T>::matrix()
{
	_nrow = 0;
	_ncol = 0;
}


template <typename T>
matrix<T>::matrix(size_t n)
{
	_nrow = n;
	_ncol = n;
	_matr.resize(_nrow);
	for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);
}


template <typename T>
matrix<T>::matrix(size_t nrow, size_t ncol)
{
	_nrow = nrow;
	_ncol = ncol;
	_matr.resize(_nrow);
	for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);
}


template <typename T>
template <typename U>
matrix<T>::matrix(const matrix<U> & m)
{
	_nrow = m.get_row_num();
	_ncol = m.get_col_num();

	_matr.resize(_nrow);
	for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);

	for (size_t i = 0; i < _nrow; i++) {
		for (size_t j = 0; j < _ncol; j++) {
			_matr[i][j] = m[i][j];
		}
	}
}

template <typename T>
matrix<T>::matrix(temp_matrix & m)
{
	_nrow = m._nrow;
	_ncol = m._ncol;

	_matr.resize(_nrow);
	for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);

	for (size_t i = 0; i < _nrow; i++) {
		for (size_t j = 0; j < _ncol; j++) {
			_matr[i][j] = m._matr[i][j];
		}
	}
}


//-----assign_input_output--------------------------------------------------------

template <typename T>
template <typename U> 
matrix<T> & matrix<T>::operator = (const matrix<U> & m)
{
	_nrow = m.get_row_num();
	_ncol = m.get_col_num();

	_matr.resize(_nrow);
	for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);

	for (size_t i = 0; i < _nrow; i++) {
		for (size_t j = 0; j < _ncol; j++) {
			_matr[i][j] = m[i][j];
		}
	}
}


template <typename T>
std::istream & operator >> (std::istream & is, matrix<T> & m)
{
	m.scan();
	return is;
}


template <typename T>
std::ostream & operator << (std::ostream & os, const matrix<T> & m)
{
	m.print();
	return os;
}


template <typename T>
void matrix<T>::scan ()
{
	if (_nrow == 0) {
		cout << "Enter count of rows and columns: ";
		cin >> _nrow >> _ncol;
		_matr.resize(_nrow);
		for (size_t i = 0; i < _nrow; i++) _matr[i].resize(_ncol);
	}

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < _ncol; j++) cin >> _matr[i][j];
}


template <typename T>
void matrix<T>::print () const
{
	cout << "size: " << _nrow << "x" << _ncol << endl;

	for (size_t i = 0; i < _nrow; i++) {
		for (size_t j = 0; j < _ncol; j++) cout << _matr[i][j] << ' ';
		cout << endl;
	}
}



//-----get_size_fuctions----------------------------------------------------------

template <typename T>
size_t matrix<T>::get_row_num () const
{
	return _nrow;
}


template <typename T>
size_t matrix<T>::get_col_num () const
{
	return _ncol;
}



//-----arifmetic_fuctions---------------------------------------------------------

template <typename T>
matrix<T> matrix<T>::operator * (T c) const
{
	matrix * ans = new matrix(*this);

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < _ncol; j++) ans->_matr[i][j] *= c;

	return *ans;
}


template <typename T>
matrix<T> matrix<T>::operator + (const matrix & m) const
{
	matrix * ans = new matrix(_nrow, _ncol);

	if ((_nrow != m._nrow) || (_ncol != m._ncol)) {
		cerr << "Incorrect size" << endl;
		return *ans;
	}

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < _ncol; j++) ans->_matr[i][j] = _matr[i][j] + m._matr[i][j];

	return *ans;
}


template <typename T>
matrix<T> matrix<T>::operator - (const matrix & m) const
{
	matrix * ans = new matrix(_nrow, _ncol); //<double> ли как??????

	if ((_nrow != m._nrow) || (_ncol != m._ncol)) {
		cerr << "Incorrect size" << endl;
		return *ans;
	}

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < _ncol; j++) ans->_matr[i][j] = _matr[i][j] - m._matr[i][j];

	return *ans;
}


template <typename T>
matrix<T> matrix<T>::operator * (const matrix & m) const
{
	matrix * ans = new matrix(_nrow, m._ncol);

	if (_ncol != m._nrow) {
		cerr << "Incorrect size" << endl;
		return *ans;
	}

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < m._ncol; j++) {
			T sum = 0;
			for (size_t k = 0; k < _ncol; k++) sum += _matr[i][k] * m._matr[k][j];
			ans->_matr[i][j] = sum;
		}

	return *ans;
}


//-----matrix_numbers_and_representations-----------------------------------------

template <typename T>
T matrix<T>::trace() const
{
	if (_nrow != _ncol) {
		cerr << "Incorrect size" << endl;
		return 0;
	}

	T tr = 0;
	for (size_t i = 0; i < _nrow; i++) tr += _matr[i][i];

	return tr;
}


template <typename T>
T matrix<T>::det() const
{
 	if (_nrow != _ncol) {
		cerr << "Incorrect size" << endl;
		return 0;
	}

	temp_matrix temp(*this);

	double d = 1;

	for (size_t i = 0; i < _nrow; i++) {
		size_t j;
		for (j = i; j < _nrow; j++)
			if (temp._matr[j][i] != 0) break;

		if (j == _nrow) return 0;

		if (i != j) {
			temp.strswap(i, j);
			d *= -1;
		}

		for (j = i + 1; j < _nrow; j++) temp.strsub(j, i, temp._matr[j][i] / temp._matr[i][i]);

		d *= temp._matr[i][i];
	}

	return (T (d));
}

template <typename T>
matrix<T> matrix<T>::transp() const
{
	matrix ans(_ncol, _nrow);

	for (size_t i = 0; i < _nrow; i++)
		for (size_t j = 0; j < _ncol; j++) ans._matr[j][i] = _matr[i][j];

	return ans;
}

template <typename T>
matrix<double> matrix<T>::inv() const
{
	temp_matrix ans(_nrow, _nrow);

 	if (_nrow != _ncol) {
		cerr << "Incorrect size" << endl;
		return matrix<double>(ans);
	}

	temp_matrix temp = *this;

	for (size_t i = 0; i < _nrow; i++) ans._matr[i][i] = 1;

	for (size_t i = 0; i < _nrow; i++) {
		size_t j;
		for (j = i; j < _nrow; j++)
			if (temp._matr[j][i] != 0) break;

		if (j == _nrow) {
			cerr << "Inverse matrix doesn't exist" << endl;
			exit(EXIT_FAILURE);
		}

		if (i != j) {
			ans.strswap(i, j);
			temp.strswap(i, j);
		}

		ans.strmult(i, ((double) 1) / temp._matr[i][i]);
		temp.strmult(i, ((double) 1) / temp._matr[i][i]);

		for (j = i + 1; j < _nrow; j++) {
			ans.strsub(j, i, temp._matr[j][i]);
			temp.strsub(j, i, temp._matr[j][i]);
		}
	}

	for (size_t i = _nrow - 1; i > 0; i--)
		for (size_t j = i; j > 0; j--) {
			ans.strsub(j - 1, i, temp._matr[j-1][i]);
			temp.strsub(j - 1, i, temp._matr[j-1][i]);
		}

	return matrix<double>(ans);
}
