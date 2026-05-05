#include<iostream>
#include<vector>
#include<utility>
#include<fstream>
#include<string>
#include<cmath>
#include<chrono>

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

void save(std::vector<double> inp, const std::string filename, const std::string interval = " ") {
	std::ofstream file(filename);
	if (!file.is_open()) {
		return;  // 可选：抛出异常或输出错误信息
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

int main() {
	//初始化 Инициаляция
	double minx, maxx, mint, maxt, tau, tCentre, A, epsilon, dt, dx;
	int maxTimes = 1E5, stopped = 0, Nx, Nt;
	minx = 0;
	maxx = 1;
	mint = 0;
	maxt = 10;
	Nx = 2048;
	Nt = 32768;
	std::vector<double> X, T, utemp(3);
	std::vector<std::vector<double>>u(Nx, std::vector<double>(Nt, 0));//Результат вычисления, автоматически совпадаемый ГУ u(x=0)=0 u函数，自动满足u(x=0)=0边界条件

	X = linspace(minx, maxx, Nx);
	T = linspace(mint, maxt, Nt);
	auto[Xm, Tm] = meshgrid(X, T);
	save(Xm, "Xm.txt");
	save(Tm, "Tm.txt");
	dx = X[2] - X[1];
	dt = T[2] - T[1];

	//边界条件 Граничное условие для t=0
	for (int i = 0; i < Nx; i++) {
		u[i][0] = X[i];
	}

	//计算 Вычисление
	double tanval, F, dF, v;
	epsilon = 1E-10;
	auto total_start = std::chrono::steady_clock::now();
	for (size_t i = 1; i < Nx; i++) {
		stopped = 0;
		auto row_start = std::chrono::steady_clock::now();
		for (size_t j = 0; j < Nt-1; j++) {
			//第一次计算 Первое вычисление
			v = u[i][j];
			tanval = atan(u[i - 1][j + 1]); 
			F = v - u[i][j] + dt * (atan(v) - tanval) / dx;//F = 0
			dF = 1 + dt / (dx * (1 + v *v));
			v -= F / dF;
			for (size_t k = 1; k < maxTimes; k++) {
				utemp[0] = v;
				F = v - u[i][j] + dt * (atan(v) - tanval) / dx;
				dF = 1 + dt / (dx * (1 + v * v));
				v -= F / dF;
				utemp[1] = v;
				F = v - u[i][j] + dt * (atan(v) - tanval) / dx;
				dF = 1 + dt / (dx * (1 + v * v));
				utemp[2] = v;
				if (abs((utemp[2] - utemp[1]) / (1 - (utemp[2] - utemp[1]) / (utemp[1] - utemp[0] + 1E15))) < epsilon) {
					u[i][j + 1] = utemp[2];
					stopped++;
					break;
				}
			}
			u[i][j + 1] = utemp[2];
		}

		auto row_end = std::chrono::steady_clock::now();
		std::chrono::duration<double> row_elapsed = row_end - row_start;
		std::chrono::duration<double> total_elapsed = row_end - total_start;

		// 计算进度百分比
		double progress = 100.0 * (i + 1) / Nx;   // i+1 因为 i 从 1 开始，但第一行是 i=1，共 Nx-1 行
		// 预估剩余时间（基于已用总时间）
		double remaining = (total_elapsed.count() / (i)) * (Nx - 1 - i);  // i 是当前处理的行数（已完成 i 行）
		std::cout << "Row " << i << "/" << Nx - 1
			<< " (" << progress << "%) | "
			<< "Row time: " << row_elapsed.count() << "s | "
			<< "Total: " << total_elapsed.count() << "s | "
			<< "ETA: " << remaining << "s | "
			<< "Early Stopped: " << stopped << std::endl;
	}
	save(u, "u.txt");
	double balan = 0;
	for (size_t m = 0; m < Nx; ++m) {
		balan += dx * (u[m][Nt - 1] - u[m][0]);
	}
	for (size_t n = 0; n < Nt - 1; ++n) {
		balan += dt * std::atan(u[Nx - 1][n]);
	}
	std::cout << balan;
	return 0;
}