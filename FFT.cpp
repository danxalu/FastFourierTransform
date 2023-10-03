#define PI 3.141592653589793238462643383279502884197169399375105820974944
#include <complex>
#include <vector>
#include <iostream>
#include <random>

typedef std::complex<double> complex;
typedef std::vector<complex> complex_vector;


complex_vector operator - (const complex_vector& c1, const complex_vector& c2)
{
	int n = (int)c1.size();
	if (n != (int)c2.size()) return c1;
	complex_vector answer(n);
	for (int i = 0; i < n; i++) {
		answer[i] = c1[i] - c2[i];
	}
	return answer;
}


complex_vector copy_vector(const complex_vector& c1) {
	complex_vector c2;
	for (auto i : c1) {
		c2.push_back(i);
	}
	return c2;
}


void print_vector(const complex_vector& vec) {
	for (auto i : vec) {
		std::cout << i << " ";
	}
	std::cout << "\n";
}


double taxi_norm_vector(const complex_vector& vec) {
	double norm = 0;
	for (auto i : vec) {
		norm += abs(i);
	}
	return norm;
}


void fft_5(complex_vector& a, bool inv) {
	int n = (int)a.size();
	if (n % 5 != 0) return;

	complex_vector a0(n / 5), a1(n / 5), a2(n / 5), a3(n / 5), a4(n / 5);
	for (int i = 0, j = 0; i < n; i += 5, ++j) {
		a0[j] = a[i];
		a1[j] = a[i + 1];
		a2[j] = a[i + 2];
		a3[j] = a[i + 3];
		a4[j] = a[i + 4];
	}
	fft_5(a0, inv);
	fft_5(a1, inv);
	fft_5(a2, inv);
	fft_5(a3, inv);
	fft_5(a4, inv);

	double ang = 2 * PI / n * (inv ? -1 : 1);
	double ang_5 = 2 * PI / 5 * (inv ? -1 : 1);
	complex w(1, 0), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			a[i + j*n / 5] = a0[i] + complex(cos(ang_5*j), sin(ang_5*j)) * w * a1[i] + pow(complex(cos(ang_5*j), sin(ang_5*j)) * w, 2) * a2[i] + pow(complex(cos(ang_5*j), sin(ang_5*j)) * w, 3) * a3[i] + pow(complex(cos(ang_5*j), sin(ang_5*j)) * w, 4) * a4[i];
			if (inv) {
				a[i + j*n / 5] /= 5;
			}
		}
		w *= wn;
	}
}


void fft_3(complex_vector& a, bool inv) {
	int n = (int)a.size();
	if (n % 3 != 0) {
		if (n % 5 == 0) {
			fft_5(a, inv);
		}
		return;
	}
	
	complex_vector a0(n / 3), a1(n / 3), a2(n / 3);
	for (int i = 0, j = 0; i < n; i += 3, ++j) {
		a0[j] = a[i];
		a1[j] = a[i + 1];
		a2[j] = a[i + 2];
	}
	fft_3(a0, inv);
	fft_3(a1, inv);
	fft_3(a2, inv);

	double ang = 2 * PI / n * (inv ? -1 : 1);
	double ang_3 = 2 * PI / 3 * (inv ? -1 : 1);
	complex w(1, 0), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 3; ++i) {
		a[i] = a0[i] + w * a1[i] + pow(w, 2) * a2[i];
		a[i + n / 3] = a0[i] + complex(cos(ang_3), sin(ang_3)) * w * a1[i] + complex(cos(ang_3 * 2), sin(ang_3 * 2)) * pow(w, 2) * a2[i];
		a[i + 2 * n / 3] = a0[i] + complex(cos(ang_3 * 2), sin(ang_3 * 2)) * w * a1[i] + complex(cos(ang_3), sin(ang_3)) * pow(w, 2) * a2[i];
		if (inv) {
			a[i] /= 3;
			a[i + n / 3] /= 3;
			a[i + 2*n / 3] /= 3;
		}
		w *= wn;
	}
}


void fft_2(complex_vector& a, bool inv) {
	int n = (int)a.size();
	if (n % 2 != 0) {
		if (n % 3 == 0) {
			fft_3(a, inv);
		}
		else {
			if (n % 5 == 0) {
				fft_5(a, inv);
			}
		}
		return;
	}

	complex_vector a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i < n; i += 2, ++j) {
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}
	fft_2(a0, inv);
	fft_2(a1, inv);

	double ang = 2 * PI / n * (inv ? -1 : 1);
	complex w(1, 0), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 2; ++i) {
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (inv) {
			a[i] /= 2;
			a[i + n / 2] /= 2;
		}
		w *= wn;
	}
}


void fft(complex_vector& a, bool inv) {
	int n = (int)a.size();
	if (n == 1) return;
	if (n % 2 == 0) {
		fft_2(a, inv);
	}
	else {
		if (n % 3 == 0) {
			fft_3(a, inv);
		}
		else {
			if (n % 5 == 0) {
				fft_5(a, inv);
			}
		}
	}
}

int main() {
	complex_vector vec;
	srand(time(NULL));
	int len = pow(2, (rand() % 3)) * pow(3, (rand() % 3)) * pow(5, (rand() % 3));
	for (int i = 0; i < len; i++) {
		bool sign1 = rand() % 2;
		bool sign2 = rand() % 2;
		vec.push_back(complex(pow(-1, sign1) * (rand() % 100), pow(-1, sign2) * (rand() % 100)));
	}

	std::cout << "Random vector (length: " << len << "): ";
	complex_vector vec_begin = copy_vector(vec);
	print_vector(vec_begin);

	std::cout << "\nDirect FFT vector: ";
	fft(vec, false);
	complex_vector vec_direct = copy_vector(vec);
	print_vector(vec_direct);

	std::cout << "\nInverse FFT vector: ";
	fft(vec, true);
	complex_vector vec_inverse = copy_vector(vec);
	print_vector(vec_inverse);

	std::cout << "\nError: ";
	complex_vector vec_error = vec_begin - vec_inverse;
	std::cout << taxi_norm_vector(vec_error);
}
