#include <cstdint>
#include <iostream>
using namespace std;

uint16_t linear_congruence(uint16_t a, uint16_t c, uint16_t prev) {
    return a * prev + c;
}

int main() {
    cout << "Hello World!" << endl;
}