#include <iostream>
#include <random>
using namespace std;

std::mt19937 rng(10);

int main(){
	cout << rng() << endl;
	cout << rng() << endl;
	return 0;
}
