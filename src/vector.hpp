#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED

#include <cstddef>

#ifdef _MSC_VER
#   include <intrin.h>
#elif defined(__MINGW32__) || defined(__CYGYWIN__)
#   include <w32api.h>
#   include <intrin.h>
#elif __GNUC__
#   include <x86intrin.h>
#endif

namespace org {
    namespace sqg {

        struct float4 {
            __m128 xmm;
            float4() {}
            float4(__m128 const &x) { xmm = x; }
            float4 &operator=(__m128 const &x) {
                xmm = x;
                return *this;
            }

            float4& load(float const *p) {
                xmm = _mm_loadu_ps(p);
                return *this;
            }

            float4 const& store(float *p) const {
                _mm_storeu_ps(p, xmm);
                return *this;
            }

            float sum() const {
                float values[4];
                store(&values[0]);
                return values[0] + values[1] + values[2] + values[3];
            }

            operator __m128() const { return xmm; }
        };

        static inline float4 operator+(float4 const &a, float4 const &b) {
            return _mm_add_ps(a, b);
        }
        static inline float4 operator*(float4 const &a, float4 const &b) {
            return _mm_mul_ps(a, b);
        }

        struct block3 {
            float4 x, y, z;
        };

        struct block4 {
            float4 x, y, z, w;
        };

        static inline float4 dot(block3 const &a, block3 const &b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }

        static inline float4 dot(block4 const &a, block4 const &b) {
            return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
        }

        template <typename T, size_t N> class vector {
        public:
            T const &operator[](size_t idx) const { return _M_data[idx]; }
            T &operator[](size_t idx) { return _M_data[idx]; }

            size_t size() const { return N; }

            T operator*(vector<T, N> const &other) const { return dot(other); }

        private:
            T dot(vector<T, N> const &other) const {
                T sum = T(0);
                for (size_t i = 0; i != N; ++i)
                    sum += _M_data[i] * other._M_data[i];
                return sum;
            }

        private:
            T _M_data[N];
        };

        template <size_t N> class vector<double, N> {
        public:
            double const &operator[](size_t idx) const { return _M_data[idx]; }
            double &operator[](size_t idx) { return _M_data[idx]; }

            size_t size() const { return N; }

            double operator*(vector<double, N> const &other) const {
                return dot(other);
            }

        private:
            double dot(vector<double, N> const &other) const {
                double sum = double(0);
#ifdef __SSE3__
                size_t const UNROLL = 2;
                for (size_t i = 0, n = size() / UNROLL * UNROLL; i != n;
                     i += UNROLL) {
                    sum += _M_data[i + 0] * other._M_data[i + 0] +
                           _M_data[i + 1] * other._M_data[i + 1];
                }
                for (size_t i = size() / UNROLL * UNROLL, n = size(); i != n;
                     ++i)
                    sum += _M_data[i] * other._M_data[i];
#else
                for (size_t i = 0, n = size(); i != n; ++i)
                    sum += _M_data[i] * other._M_data[i];
#endif
                return sum;
            }

        private:
            double _M_data[N];
        };

        template <size_t N> class vector<float, N> {
        public:
            float const &operator[](size_t idx) const { return _M_data[idx]; }
            float &operator[](size_t idx) { return _M_data[idx]; }

            size_t size() const { return N; }

            float operator*(vector<float, N> const &other) const {
                return dot(other);
            }

        private:
            float dot(vector<float, N> const &other) const {
                float sum = float(0);
#ifdef __SSE3__
                size_t const UNROLL = 4;
                float4 a, b, c;
                for (size_t i = 0, n = size() / UNROLL * UNROLL; i != n;
                     i += UNROLL) {
                    a.load(&_M_data[i]);
                    b.load(&other._M_data[i]);
                    c = c + a * b;
                    //sum += _M_data[i + 0] * other._M_data[i + 0] +
                    //       _M_data[i + 1] * other._M_data[i + 1] +
                    //       _M_data[i + 2] * other._M_data[i + 2] +
                    //       _M_data[i + 3] * other._M_data[i + 3];
                }
                sum = c.sum();
                for (size_t i = size() / UNROLL * UNROLL, n = size(); i != n;
                     ++i)
                    sum += _M_data[i] * other._M_data[i];
#else
                for (size_t i = 0, n = size(); i != n; ++i)
                    sum += _M_data[i] * other._M_data[i];
#endif
                return sum;
            }

        private:
            float _M_data[N];
        };
    }
}

#endif // VECTOR_HPP_INCLUDED
