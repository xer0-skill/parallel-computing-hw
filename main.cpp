#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <thread>
#include <time.h>
#include <stdlib.h>
#include <chrono>

using namespace std;
using namespace galois;

#define PRINT_ELAPSED_TIME_WITH_PREFIX(prefix, code)														\
start = std::chrono::high_resolution_clock::now();															\
code																										\
finish = std::chrono::high_resolution_clock::now();															\
std::cout << prefix"Elapsed time: " << std::chrono::duration<double>(finish - start).count() << " s\n";	\

#define PRINT_ELAPSED_TIME(code) PRINT_ELAPSED_TIME_WITH_PREFIX("", code)


shared_ptr<GaloisField> _pGf;

void initField(vector<unsigned int> primPoly) {
	_pGf = make_shared<GaloisField>((int)(primPoly.size() - 1), primPoly.data());
}

bool powOf2(unsigned int val) {
	return ((val - 1) & val) == 0;
}

GaloisFieldPolynomial trace(GaloisField* _gf, unsigned char i) {
	unsigned int m = _gf->pwr();
	auto betaI = GaloisFieldElement(_gf, _gf->alpha(i));
	auto nullElem = betaI + betaI;
	size_t deg = (size_t)0x01 << (m - 1);
	vector<GaloisFieldElement> cPoly(deg + 1);
	for (unsigned int i = 0; i <= deg; ++i) {
		if (powOf2(i)) {
			cPoly[i] = betaI ^ i;
		}
		else {
			cPoly[i] = nullElem;
		}
	}
	auto tracePoly = GaloisFieldPolynomial(_gf, (unsigned int)deg, cPoly.data());
	return tracePoly;
}

// 0 deg		- 1
// 1 deg		- alpha
// 2^m - 2 deg	- alpha^(2^m - 2) = alpha^(-1)
// 2^m - 1 deg	- 0
GaloisFieldPolynomial createPoly(GaloisField* _gf, vector<GFSymbol> rootDegs) {
	auto poly = GaloisFieldPolynomial(_gf, 0, &GaloisFieldElement(_gf, _gf->alpha(0))); // poly = 1
	for (auto i : rootDegs) {
		GaloisFieldElement cFirstDegPoly[] = {				// x + alpha^i == x - alpha^i
			GaloisFieldElement(_gf, _gf->alpha(i)),
			GaloisFieldElement(_gf, _gf->alpha(0)) };
		auto firstDegPoly = GaloisFieldPolynomial(_gf, 1, cFirstDegPoly);
		poly *= firstDegPoly;
	}
	return poly;
}

void findPolyRoots(const GaloisFieldPolynomial& poly, vector<GFSymbol>& roots, unsigned int i = 0) {
	if (poly.deg() == 1) {
		//roots.push_back((poly / poly[1])[0].index());
		roots.push_back((poly[0] * poly[1].inverse()).index());
		return;
	}
	if (poly.deg() == 2) {
		int a = 0;
	}
	GaloisField* pGf = poly.field();
	unsigned int power = pGf->pwr();

	for (; i < power; ++i) {
		auto poly1 = trace(pGf, i);
		//auto poly2 = poly1;
		//poly2[0] += GaloisFieldElement(pGf, pGf->alpha(0));
		auto gcd1 = gcd(poly, poly1);
		if (gcd1.deg() == 0 || gcd1.deg() == poly.deg()) {
			continue;
		}
		auto gcd2 = poly / gcd1;

		findPolyRoots(gcd1, roots, i + 1);
		findPolyRoots(gcd2, roots, i + 1);		

		return;
	}

	abort();
}

tuple<GaloisFieldPolynomial, GaloisFieldPolynomial, unsigned int> splitPolyBy2(const GaloisFieldPolynomial& poly, unsigned int i) {
	assert(poly.deg() >= 2);
	GaloisField* pGf = poly.field();
	unsigned int power = pGf->pwr();

	for (; i < power; ++i) {
		auto tracePoly = trace(pGf, i);
		auto gcd1 = gcd(poly, tracePoly);
		if (gcd1.deg() == 0 || gcd1.deg() == poly.deg()) {
			continue;
		}
		auto gcd2 = poly / gcd1;

		return { gcd1, gcd2, i};
	}

	abort();
}

vector<pair<GaloisFieldPolynomial, unsigned int>> splitPoly(const GaloisFieldPolynomial& poly, unsigned int parts) {
	assert(poly.deg() >= parts);
	typedef GaloisFieldPolynomial gfp;
	typedef pair<gfp, unsigned int> gfp_i;
	vector<gfp_i> res;
	res.push_back({ poly, 0 });
	auto cmp = [](const gfp_i p1, const gfp_i p2) {return p1.first.deg() < p2.first.deg(); };

	while (res.size() < parts) {
		sort(res.begin(), res.end(), cmp);
		gfp_i p = res.back();
		res.pop_back();
		auto [p1 ,p2, i] = splitPolyBy2(p.first, p.second);
		res.push_back({ p1, i });
		res.push_back({ p2, i });
	}
	return res;
}

vector<GFSymbol> findPolyRootsSync(const GaloisFieldPolynomial& poly) {
	vector<GFSymbol> roots;
	findPolyRoots(poly, roots);
	return roots;
}

vector<GFSymbol> findPolyRootsAsync(const GaloisFieldPolynomial& poly) {
	decltype(std::chrono::high_resolution_clock::now()) start, finish;

	unsigned int cpuCount = 2;
	cout << "\tN of threads: " << cpuCount << endl;

	assert(poly.deg() >= cpuCount);

	vector<vector<GFSymbol>> vRoots(cpuCount);
	cout << "\tBut need to split to " << cpuCount << " factors first..." << endl;
	PRINT_ELAPSED_TIME_WITH_PREFIX("\t", auto polys = splitPoly(poly, cpuCount););
	cout << "\tTheir degrees:";
	for (const auto& i : polys) {
		cout << " " << i.first.deg();
	}
	cout << endl;

	cout << "\tFinally running threads..." << endl;
	PRINT_ELAPSED_TIME_WITH_PREFIX("\t",
		vector<thread> ths(cpuCount);
		for (unsigned int i = 0; i < cpuCount; ++i) {
			ths[i] = thread(&findPolyRoots, std::ref(polys[i].first), std::ref(vRoots[i]), polys[i].second);
		}
		for (unsigned int i = 0; i < cpuCount; ++i) {
			ths[i].join();
		}
	);

	vector<GFSymbol> roots;
	size_t nRoots = 0;
	for (const auto& i : vRoots) {
		nRoots += i.size();
	}
	roots.reserve(nRoots);
	for (const auto& i : vRoots) {
		roots.insert(roots.end(), i.begin(), i.end());
	}

	return roots;
}

//vector<unsigned int> primPoly = { 1,1,0,0,1};
//vector<unsigned int> primPoly = { 1,0,0,0,0,0,1,1 }; //deg 7
//vector<unsigned int> primPoly = { 1,0,0,1,0,0,0,0,0,0,1 }; //deg 10
//vector<unsigned int> primPoly = { 1,0,0,0,0,0,0,0,0,1,0,1 }; //deg 11
//vector<unsigned int> primPoly = { 1,0,1,1,0,0,0,0,0,1,0,0,1 }; //deg 12
//vector<unsigned int> primPoly = { 1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1 }; //deg 15

int main() {
	decltype(std::chrono::high_resolution_clock::now()) start, finish;
	vector<unsigned int> primPoly = { 1,0,0,1,0,0,0,0,0,0,1 }; //deg 10
	cout << "m = " << primPoly.size() - 1 << endl;

	cout << "Initializing field..." << endl;
	initField(primPoly);


	vector<GFSymbol> rootsSet((size_t)1 << (primPoly.size() - 1));
	for (unsigned int i = 0; i < rootsSet.size(); ++i) {
		rootsSet[i] = i;
	}
	size_t nRoots = rootsSet.size() / 2;
	vector<GFSymbol> roots(nRoots);
	srand((unsigned int)time(NULL));
	for (unsigned int i = 0; i < nRoots; ++i) {
		unsigned int index = rand() % rootsSet.size();
		roots[i] = rootsSet[index];
		rootsSet.erase(rootsSet.begin() + index);
	}

	cout << "Creating polynomial with " << roots.size() << " roots..." << endl;
	auto poly = createPoly(_pGf.get(), roots);

	//cout << poly;

	cout << "Running multithreaded root finding..." << endl;
	auto roots2 = findPolyRootsAsync(poly);

	cout << "Running one-threaded root finding..." << endl;
	PRINT_ELAPSED_TIME(auto roots1 = findPolyRootsSync(poly););

	return 0;
}
