#ifndef GISCONVERT_HHH
#define GISCONVERT_HHH

#include <math.h> // for pow and sin
// #include <tuple> // for return value
// #include <utility> 
#include <gsl/gsl_integration.h>

#define EARTH_LONG 6378137.0f
#define FALSE_EASTING 500000
#define BANDS 6

enum ModelType {WGS84, CGCS2000};

// template<typename T> 
// constexpr T PI = T{3.1415926535897932385L};
template<typename T> 
constexpr T PI{3.1415926535897932385L};

typedef struct EarthParams {
    double long_radius;
    double flatteness;
    double minor_radius;
    double ecc1;
    double ecc2;
} EarthParams;

class EarthModel {
    // member
    auto _calMinorRadius() {
        _minor_radius = _long_radius * (1.0f - _flatteness);
    }

    auto _calEccentricity1() {
        _eccentricity1 = 2.0f * _flatteness - _flatteness * _flatteness;
    }

    auto _calEccentricity2() {
        _eccentricity2 = (_long_radius * _long_radius - _minor_radius * _minor_radius) / (_minor_radius * _minor_radius);
    }
public:
    // constructor by model_type
    EarthModel(ModelType model_type = WGS84) : 
       _long_radius{EARTH_LONG}, _flatteness{0.0f}, _minor_radius(0.0f),
       _eccentricity1{0.0f}, _eccentricity2{0.0f} {
        calEarthParams(model_type);
    }

    // get ell model EarthParams;
    auto getEarthParams(ModelType model_type = WGS84) -> EarthParams {
        return {_long_radius, _flatteness, _minor_radius, _eccentricity1, _eccentricity2};
    }

    auto calEarthParams(ModelType model_type = WGS84) -> void {
        // cal 
        switch(model_type) {
            //'wgs84'
            case WGS84:
                _flatteness = 1.0f/ 298.257223563;
                break;
            // 'cgcs2000'
            case CGCS2000:
                _flatteness = 1.0f / 298.257222101;
                break;
            default:
                break;
        }
        _calMinorRadius();
        _calEccentricity1();
        _calEccentricity2();
    }
    

private:
    // member
    double _long_radius;
    double _flatteness;
    double _minor_radius;
    double _eccentricity1;
    double _eccentricity2;
};


typedef struct Coor2D {
    double x;
    double y;
} Coor2D;
class GisConverter {
private:
    static auto area(double x, void* p) -> double {
        EarthParams params = *((EarthParams *)p);
        return params.long_radius * (1 - params.ecc2) / std::pow(1 - params.ecc2 * std::pow(std::sin(x), 2), 3.0f / 2);
    }
public:
    static auto gisToGauss(double latitude, double longitude, ModelType type) -> Coor2D {
        auto lat_rad = latitude / 180 * PI<double>;
        // auto lon_rad = longitude / 180 * PI<double>;

        // six degree band model assumption
        auto band_order = int(longitude / BANDS) + 1;
        auto l = (longitude - (band_order * BANDS - 3)) / 180 * PI<double>;
        
        auto params = EarthModel(type).getEarthParams();
        auto long_radius = params.long_radius;
        auto flatteness = params.flatteness;
        auto minor_raduis = params.minor_radius;
        auto ecc1 = params.ecc1;
        auto ecc2 = params.ecc2;
        // gsl_integration
        gsl_function f;
        f.function = &area; 
        f.params = &params;

        // gsl EarthParams
        auto ret{0.0}; 
        auto err{0.0};
        auto neval = size_t{0};

        gsl_integration_qng(&f, 0.0, lat_rad, 1e-11, 1e-11, &ret, &err, &neval);

        // gauss plane coordinates
        auto N = long_radius / std::sqrt(1 - flatteness * (2 - flatteness) * std::sin(lat_rad) * std::sin(lat_rad));
        auto latitude_tan = std::tan(lat_rad);
        auto niu2 = ecc2 * std::pow(std::cos(lat_rad), 2);
        
        auto x_coor = ret + N / 2.0f * std::sin(lat_rad) * cos(lat_rad) * l * l + N / 24.0f * std::sin(lat_rad) * std::pow(cos(lat_rad), 3) * (5 - std::pow(latitude_tan, 2) + 9 * niu2 + 4 * std::pow(niu2, 2)) * std::pow(l, 4) + N / 720.0f * std::sin(lat_rad) * std::pow(std::cos(lat_rad), 5) * (61 - 58 * std::pow(latitude_tan, 2) + std::pow(latitude_tan, 4)) * std::pow(l, 6);

        auto y_coor = N * std::cos(lat_rad) * l + double(N) / BANDS * std::pow(std::cos(lat_rad), 3) * (1 - std::pow(latitude_tan, 2) + niu2) * std::pow(l, 3)  + N / 120.0f * std::pow(std::cos(lat_rad), 5) * (5 - 18 * std::pow(latitude_tan, 2) + std::pow(latitude_tan, 4) + 14 * niu2 - 58 * std::pow(latitude_tan, 2) * niu2) * std::pow(l, 5);
        y_coor += FALSE_EASTING;
        y_coor += band_order * std::pow(10, std::floor(std::log10(y_coor) + 1));

        return {x_coor, y_coor};
    }
};

#endif 
