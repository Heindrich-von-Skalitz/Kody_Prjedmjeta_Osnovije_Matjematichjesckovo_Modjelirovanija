#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include <cmath>
#define pi 3.1415926535897932

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> meshgrid(std::vector<double> inp1, std::vector<double> inp2) {
	std::vector<std::vector<double>> out1(inp1.size(), std::vector<double>(inp2.size(), 0)), out2(inp1.size(), std::vector<double>(inp2.size(), 0));
	for (size_t i = 0; i < inp1.size(); i++) {
		for (size_t j = 0; j < inp2.size(); j++) {
			out1[i][j] = inp1[i];
			out2[i][j] = inp2[j];
		}
	}
	return { out1, out2 };
}

std::vector<double> linspace(double minim, double maxim, size_t numb) {
	std::vector<double> out(numb);
	double d = (maxim - minim) / (numb - 1);
	for (size_t i = 0; i < numb; i++) {
		out[i] = minim + d * i;
	}
	return out;
}

std::vector<double> freqspace(double minim, double maxim, size_t numb, bool shifted = false) {
	std::vector<double> out(numb);
	if (numb < 2) {
		if (numb == 1) out[0] = 0.0;
		return out;
	}
	double dt = (maxim - minim) / static_cast<double>(numb - 1), fs = 1.0 / dt, df = fs / static_cast<double>(numb);

	if (shifted) {
		// 何移位：负频率在前，零频在中心 (若 numb 为偶数：-fs/2 到 fs/2-df)
		int half = static_cast<int>(numb) / 2;
		for (size_t i = 0; i < numb; ++i) {
			out[i] = (static_cast<int>(i) - half) * df;
		}
	}
	else {
		// 未移位：0, df, 2*df, ... , (numb-1)*df
		for (size_t i = 0; i < numb; ++i) {
			out[i] = static_cast<double>(i) * df;
		}
	}
	return out;
}

void save(std::vector<double> inp, const std::string filename, const std::string interval = " ") {
	std::ofstream file(filename);
	if (!file.is_open()) {
		return;
	}
	for (size_t i = 0; i < inp.size(); ++i) {
		file << inp[i];
		if (i != inp.size() - 1) {
			file << interval;
		}
	}
}

void save(std::vector<std::vector<double>> inp, const std::string filename, const std::string interval = " ") {
	std::ofstream file(filename);
	if (!file.is_open()) {
		return;
	}
	for (const auto& row : inp) {
		for (size_t j = 0; j < row.size(); ++j) {
			file << row[j];
			if (j != row.size() - 1) {
				file << interval;
			}
		}
		file << '\n';
	}
}

double xmin = 0, xmax = pi / 2, ymin = 0, ymax = pi / 2, tmin = 0, tmax = 0.1;
size_t Nx = 256, Ny = 256, Nt = 2000;
int plotN = 10;
auto X = linspace(xmin, xmax, Nx), Y = linspace(ymin, ymax, Ny), T = linspace(tmin, tmax, 2*Nt);
auto [XM, YM] = meshgrid(X, Y);
double dx = X[1] - X[0], dy = Y[1] - Y[0], dt = T[1] - T[0];

std::vector<double>chasingMethod(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double>d) {
	size_t N = b.size();
	if (a.size() != N - 1 || c.size() != N - 1) {
		throw(matherr, "维度错误！");
	}
	std::vector<double>result(N), y(N);
	std::vector<double>l(N - 1), U(N), C(N - 1);
	U[0] = b[0];
	y[0] = d[0];
	for (size_t i = 1; i < N; i++) {
		l[i - 1] = a[i - 1] / U[i - 1];
		U[i] = b[i] - l[i - 1] * c[i - 1];
		y[i] = d[i] - l[i - 1] * y[i - 1];
	}
	result[N - 1] = y[N - 1] / U[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		result[i] = (y[i] - c[i] * result[i + 1]) / U[i];
	}
	return result;
}

std::vector<std::vector<double>> alternatingDirectionMethod(std::vector<std::vector<double>> u_prev) {
	const double f = 0, ax = -1 / (dx * dx), bx = 2 / (dx * dx) + 1 / dt, cx = -1 / (dx * dx), ay = -1 / (dy * dy), by = 2 / (dy * dy) + 1 / dt, cy = -1 / (dy * dy);

	std::vector<std::vector<double>> u_k = u_prev;

	// X-Direction
	std::vector<std::vector<double>> u_mid(Nx, std::vector<double>(Ny, 0.0));

	for (size_t j = 0; j < Ny; ++j) {
		std::vector<double> a_sub(Nx - 1), b_vec(Nx), c_sub(Nx - 1), d_vec(Nx);

		for (size_t i = 0; i < Nx; ++i) {
			double d2u_dy2 = 0.0;
			if (j == 0) {
				d2u_dy2 = 2.0 * (u_k[i][1] - u_k[i][0]) / (dy * dy);
			}
			else if (j == Ny - 1) {
				d2u_dy2 = 2.0 * (u_k[i][Ny - 2] - u_k[i][Ny - 1]) / (dy * dy);
			}
			else {
				d2u_dy2 = (u_k[i][j + 1] - 2.0 * u_k[i][j] + u_k[i][j - 1]) / (dy * dy);
			}
			double g_xi = d2u_dy2 + f + u_k[i][j] / dt;

			// Construction of Matrix
			if (i == 0) {
				b_vec[i] = bx + cx;
				c_sub[i] = ax;
				d_vec[i] = g_xi;
			}
			else if (i == Nx - 1) {
				a_sub[i - 1] = cx;
				b_vec[i] = bx + ax;
				d_vec[i] = g_xi;
			}
			else {
				a_sub[i - 1] = cx;
				b_vec[i] = bx;
				c_sub[i] = ax;
				d_vec[i] = g_xi;
			}
		}
		std::vector<double> col = chasingMethod(a_sub, b_vec, c_sub, d_vec);
		for (size_t i = 0; i < Nx; ++i)
			u_mid[i][j] = col[i];
	}

	// y-Direction
	std::vector<std::vector<double>> u_next(Nx, std::vector<double>(Ny, 0.0));

	for (size_t i = 0; i < Nx; ++i) {
		std::vector<double> a_sub(Ny - 1), b_vec(Ny), c_sub(Ny - 1), d_vec(Ny);

		for (size_t j = 0; j < Ny; ++j) {
			double d2u_dx2 = 0.0;
			if (i == 0) {
				d2u_dx2 = 2.0 * (u_mid[1][j] - u_mid[0][j]) / (dx * dx);
			}
			else if (i == Nx - 1) {
				d2u_dx2 = 2.0 * (u_mid[Nx - 2][j] - u_mid[Nx - 1][j]) / (dx * dx);
			}
			else {
				d2u_dx2 = (u_mid[i + 1][j] - 2.0 * u_mid[i][j] + u_mid[i - 1][j]) / (dx * dx);
			}
			double g_yj = d2u_dx2 + f + u_mid[i][j] / dt;

			// Construction of Matrix
			if (j == 0) {
				b_vec[j] = by + cy;
				c_sub[j] = ay;
				d_vec[j] = g_yj;
			}
			else if (j == Ny - 1) {
				a_sub[j - 1] = cy;
				b_vec[j] = by + ay;
				d_vec[j] = g_yj;
			}
			else {
				a_sub[j - 1] = cy;
				b_vec[j] = by;
				c_sub[j] = ay;
				d_vec[j] = g_yj;
			}
		}
		std::vector<double> row = chasingMethod(a_sub, b_vec, c_sub, d_vec);
		for (size_t j = 0; j < Ny; ++j)
			u_next[i][j] = row[j];
	}
	return u_next;
}

int main() {
	std::vector<std::vector<double>> Analytics(Nx, std::vector<double>(Ny, 0.0)), Result(Nx, std::vector<double>(Ny, 0.0));
	std::string analy = "analytics", solu = "solution";

	save(XM, "XM.txt");
	save(YM, "YM.txt");

	for (size_t i = 0; i < Nx; i++) {
		for (size_t j = 0; j < Ny; j++) {
			Result[i][j] = cos(4 * XM[i][j]) * cos(2 * YM[i][j]);
		}
	}
	save(Result, solu + "0.txt");
	save(Result, analy + "0.txt");

	int saveInterval = Nt / plotN;
	for (size_t k = 1; k < Nt; k++) {
		Result = alternatingDirectionMethod(Result);   // 从 u^{k} -> u^{k+2}

		if (k % saveInterval == 0) {
			save(Result, solu + std::to_string(T[2*k]) + ".txt");
			std::cout << "已计算: " << float(k) / float(Nt) * 100 << "%" << std::endl;
		}
	}
	return 0;
}