#include "gps_gausscoord_transform.h"

#include <iostream>
#include <iomanip>

using std::cin;
using std::cout;

int main(int argc, char* argv[]) {
    auto latitude{0.0f};
    auto longitude{0.0f};
    while(cin >> latitude >> longitude) {
        auto [x, y] = GisConverter::gisToGauss(latitude, longitude, WGS84);
        cout << std::setprecision(15) << x << ' ' << std::setprecision(15) << y << '\n';
    }
    return 0;
}
