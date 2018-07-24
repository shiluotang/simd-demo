#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>

#include "../src/vector.hpp"

int main(int argc, char* argv[]) {
    using namespace std;
    using namespace org::sqg;

    size_t const N = 0x7ffffff;
    vector<float, 0xff> a, b;
    for (size_t i = 0, n = a.size(); i != n; ++i)
        a[i] = b[i] = i;

    double dot = 0.0;
    time_t t1, t2;
    t1 = time(NULL);
    for (size_t i = 0, n = N; i != n; ++i)
        dot = a * b;
    t2 = time(NULL);
    cout
        << "N = " << N
        << ", size(a) = " << a.size()
        << ", total = " << (t2 - t1) << " seconds"
        << ", average = " << static_cast<double>(t2 - t1) / N << " seconds"
        << endl;
    cout << "dot = " << dot << endl;
    return EXIT_SUCCESS;
}
