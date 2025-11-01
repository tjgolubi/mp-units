// The MIT License (MIT)
//
// Copyright (c) 2018 Mateusz Pusz
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "Vector3D.hpp"

#include <catch2/catch_test_macros.hpp>

#include <mp-units/compat_macros.h>
#include <mp-units/ext/format.h>
#ifdef MP_UNITS_IMPORT_STD
import std;
#else
#include <matrix>
#include <ostream>
#endif
#ifdef MP_UNITS_MODULES
import mp_units;
#else
#include <mp-units/math.h>
#include <mp-units/systems/isq/mechanics.h>
#include <mp-units/systems/isq/space_and_time.h>
#include <mp-units/systems/isq/electromagnetism.h>
#include <mp-units/systems/si.h>
#endif

#include <array>
#include <cmath>

template<typename Rep = double>
using vector = STD_LA::fixed_size_column_vector<Rep, 3>;

template<typename Rep = double>
using matrix = STD_LA::fixed_size_matrix<Rep, 3, 3>;

namespace STD_LA {

template<typename Rep>
[[nodiscard]] inline Rep magnitude(const fixed_size_column_vector<Rep, 3>& v) noexcept
{
  using std::hypot;
  return static_cast<Rep>(hypot(v(0), v(1), v(2)));
} // magnitude

template<typename Rep>
[[nodiscard]] inline auto direction(const fixed_size_column_vector<Rep, 3>& v) noexcept
{
  using std::atan2;
  return atan2(v(1), v(0));
} // direction

template<typename Rep>
[[nodiscard]] inline auto elevation(const fixed_size_column_vector<Rep, 3>& v)
{
  using std::hypot;
  using std::atan2;
  return atan2(v(2), hypot(v(1), v(0)));
} // elevation

template<typename Rep>
[[nodiscard]] inline Rep magnitude(const fixed_size_matrix<Rep, 3, 3>& m) noexcept
{
  using std::sqrt;
  auto sum = Rep{};
  for (std::size_t i=0; i!=3; ++i) {
    for (std::size_t j=0; j!=3; ++j) {
      const auto& e = m(i, j);
      sum += e * e;
    }
  }
  return sqrt(sum);
} // magnitude

} // STD_LA

template<typename Rep>
std::ostream& operator<<(std::ostream& os, const vector<Rep>& v)
{
  os << "|";
  for (auto i = 0U; i < v.size(); ++i) {
    os << MP_UNITS_STD_FMT::format(" {:>9}", v(i));
  }
  os << " |";
  return os;
}

template<typename Rep>
std::ostream& operator<<(std::ostream& os, const matrix<Rep>& m)
{
  for (auto i = 0U; i != m.rows(); ++i) {
    os << '|';
    for (auto j = 0U; j != m.columns(); ++j)
      os << MP_UNITS_STD_FMT::format(" {:>9}", m(i, j));
    os << " |\n";
  }
  return os;
}

namespace {

using namespace mp_units;
using namespace mp_units::si::unit_symbols;

constexpr struct resistance_vec final
  : quantity_spec<isq::resistance, quantity_character::vector>
  { } resistance_vec;

constexpr struct current_vec final
  : quantity_spec<isq::electric_current, quantity_character::vector>
  { } current_vec;

template<QuantitySpec auto QS, QuantityOf<QS> Q>
  requires(Q::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar) &&
          requires (const Q& q) { q.numerical_value_ref_in(q.unit)(0); }
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const Q& q)
{
  const auto& v = q.numerical_value_ref_in(q.unit);
  return hypot(v(0) * QS[Q::unit], v(1) * QS[Q::unit], v(2) * QS[Q::unit]);
}

template<QuantitySpec auto QS, QuantityOf<QS> T>
  requires(T::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar)
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const vector<T>& v)
{
  return hypot(QS(v(0)), QS(v(1)), QS(v(2)));
}

template<QuantitySpec auto QS, QuantityOf<QS> Q>
  requires(Q::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar) &&
          requires (const Q& q) { q.numerical_value_ref_in(q.unit)(0)(0); }
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const Q& q)
{
  const auto& m = q.numerical_value_ref_in(q.unit);
  return magnitude(m) * QS[Q::unit];
}

template<QuantitySpec auto QS, QuantityOf<QS> T>
  requires(T::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar)
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const matrix<T>& m)
{
  return QS(magnitude(m));
}

template<QuantitySpec auto QS, QuantityOf<QS> Q>
  requires(Q::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar) &&
          requires (const Q& q) { q.numerical_value_ref_in(q.unit)[0]; }
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const Q& q)
{
  const auto& v = q.numerical_value_ref_in(q.unit);
  return hypot(v[0] * QS[Q::unit], v[1] * QS[Q::unit], v[2] * QS[Q::unit]);
}

template<QuantitySpec auto QS, QuantityOf<QS> T>
  requires(T::quantity_spec.character == quantity_character::vector) &&
          (QS.character == quantity_character::real_scalar)
[[nodiscard]] constexpr QuantityOf<QS> auto get_magnitude(const tjg::Vector3D<T>& v)
{
  return hypot(QS(v[0]), QS(v[1]), QS(v[2]));
}

template<typename T, typename U>
[[nodiscard]] constexpr vector<decltype(T{} * U{})> cross_product(const vector<T>& a, const vector<U>& b)
{
  return {a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0)};
}

template<Quantity Q1, Quantity Q2>
[[nodiscard]] constexpr QuantityOf<Q1::quantity_spec * Q2::quantity_spec> auto cross_product(const Q1& q1, const Q2& q2)
{
  return cross_product(q1.numerical_value_ref_in(q1.unit), q2.numerical_value_ref_in(q2.unit)) *
         (Q1::reference * Q2::reference);
}

}  // namespace

TEST_CASE("vector quantity", "[la]")
{
  SECTION("cast of unit")
  {
    SECTION("non-truncating")
    {
      const auto v = vector<int>{3, 2, 1} * isq::displacement[km];
      CHECK(v.numerical_value_in(m) == vector<int>{3000, 2000, 1000});
    }

    SECTION("truncating")
    {
      const auto v = vector<int>{1001, 1002, 1003} * isq::displacement[m];
      CHECK(v.force_numerical_value_in(km) == vector<int>{1, 1, 1});
    }
  }

  SECTION("to scalar magnitude")
  {
    const auto v = vector<int>{2, 3, 6} * isq::velocity[km / h];
    const auto speed = get_magnitude<isq::speed>(v);
    CHECK(speed.numerical_value_ref_in(km / h) == 7);
  }

  SECTION("multiply by scalar value")
  {
    const auto v = vector<int>{1, 2, 3} * isq::displacement[m];

    SECTION("integral")
    {
      SECTION("scalar on LHS") { CHECK((2 * v).numerical_value_in(m) == vector<int>{2, 4, 6}); }
      SECTION("scalar on RHS") { CHECK((v * 2).numerical_value_in(m) == vector<int>{2, 4, 6}); }
    }

    SECTION("floating-point")
    {
      SECTION("scalar on LHS") { CHECK((0.5 * v).numerical_value_in(m) == vector<double>{0.5, 1., 1.5}); }
      SECTION("scalar on RHS") { CHECK((v * 0.5).numerical_value_in(m) == vector<double>{0.5, 1., 1.5}); }
    }
  }

  SECTION("divide by scalar value")
  {
    const auto v = vector<int>{2, 4, 6} * isq::displacement[m];

    SECTION("integral") { CHECK((v / 2).numerical_value_in(m) == vector<int>{1, 2, 3}); }
    SECTION("floating-point") { CHECK((v / 0.5).numerical_value_in(m) == vector<double>{4., 8., 12.}); }
  }

  SECTION("add")
  {
    const auto v = vector<int>{1, 2, 3} * isq::displacement[m];

    SECTION("same unit")
    {
      const auto u = vector<int>{3, 2, 1} * isq::displacement[m];
      CHECK((v + u).numerical_value_in(m) == vector<int>{4, 4, 4});
    }
    SECTION("different units")
    {
      const auto u = vector<int>{3, 2, 1} * isq::displacement[km];
      CHECK((v + u).numerical_value_in(m) == vector<int>{3001, 2002, 1003});
    }
  }

  SECTION("subtract")
  {
    const auto v = vector<int>{1, 2, 3} * isq::displacement[m];

    SECTION("same unit")
    {
      const auto u = vector<int>{3, 2, 1} * isq::displacement[m];
      CHECK((v - u).numerical_value_in(m) == vector<int>{-2, 0, 2});
    }
    SECTION("different units")
    {
      const auto u = vector<int>{3, 2, 1} * isq::displacement[km];
      CHECK((v - u).numerical_value_in(m) == vector<int>{-2999, -1998, -997});
    }
  }

  SECTION("multiply by scalar quantity")
  {
    const auto v = vector<int>{1, 2, 3} * isq::velocity[m / s];

    SECTION("integral")
    {
      const auto mass = 2 * isq::mass[kg];

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK((mass * v).numerical_value_in(kg * m / s) == vector<int>{2, 4, 6}); }
        SECTION("scalar on RHS") { CHECK((v * mass).numerical_value_in(kg * m / s) == vector<int>{2, 4, 6}); }
      }
      SECTION("quantity_cast to momentum")
      {
        SECTION("scalar on LHS")
        {
          CHECK(quantity_cast<isq::momentum>(mass * v).numerical_value_in(N * s) == vector<int>{2, 4, 6});
        }
        SECTION("scalar on RHS")
        {
          CHECK(quantity_cast<isq::momentum>(v * mass).numerical_value_in(N * s) == vector<int>{2, 4, 6});
        }
      }
      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const quantity<isq::momentum[N * s], vector<int>> momentum = mass * v;
          CHECK(momentum.numerical_value_ref_in(N * s) == vector<int>{2, 4, 6});
        }
        SECTION("scalar on RHS")
        {
          const quantity<isq::momentum[N * s], vector<int>> momentum = v * mass;
          CHECK(momentum.numerical_value_ref_in(N * s) == vector<int>{2, 4, 6});
        }
      }
    }

    SECTION("floating-point")
    {
      const auto mass = 0.5 * isq::mass[kg];

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK((mass * v).numerical_value_in(kg * m / s) == vector<double>{0.5, 1., 1.5}); }
        SECTION("scalar on RHS") { CHECK((v * mass).numerical_value_in(kg * m / s) == vector<double>{0.5, 1., 1.5}); }
      }
      SECTION("quantity_cast to momentum")
      {
        SECTION("scalar on LHS")
        {
          CHECK(quantity_cast<isq::momentum>(mass * v).numerical_value_in(N * s) == vector<double>{0.5, 1., 1.5});
        }
        SECTION("scalar on RHS")
        {
          CHECK(quantity_cast<isq::momentum>(v * mass).numerical_value_in(N * s) == vector<double>{0.5, 1., 1.5});
        }
      }
      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const quantity<isq::momentum[N * s], vector<double>> momentum = mass * v;
          CHECK(momentum.numerical_value_ref_in(N * s) == vector<double>{0.5, 1., 1.5});
        }
        SECTION("scalar on RHS")
        {
          const quantity<isq::momentum[N * s], vector<double>> momentum = v * mass;
          CHECK(momentum.numerical_value_ref_in(N * s) == vector<double>{0.5, 1., 1.5});
        }
      }
    }
  }

  SECTION("divide by scalar quantity")
  {
    const auto pos = vector<int>{30, 20, 10} * isq::displacement[km];

    SECTION("integral")
    {
      const auto dur = 2 * isq::duration[h];

      SECTION("derived_quantity_spec") { CHECK((pos / dur).numerical_value_in(km / h) == vector<int>{15, 10, 5}); }
      SECTION("quantity_cast to velocity")
      {
        CHECK(quantity_cast<isq::velocity>(pos / dur).numerical_value_in(km / h) == vector<int>{15, 10, 5});
      }
      SECTION("quantity of velocity")
      {
        const quantity<isq::velocity[km / h], vector<int>> v = pos / dur;
        CHECK(v.numerical_value_ref_in(km / h) == vector<int>{15, 10, 5});
      }
    }

    SECTION("floating-point")
    {
      const auto dur = 0.5 * isq::duration[h];

      SECTION("derived_quantity_spec") { CHECK((pos / dur).numerical_value_in(km / h) == vector<double>{60, 40, 20}); }
      SECTION("quantity_cast to velocity")
      {
        CHECK(quantity_cast<isq::velocity>(pos / dur).numerical_value_in(km / h) == vector<double>{60, 40, 20});
      }
      SECTION("quantity of velocity")
      {
        const quantity<isq::velocity[km / h], vector<double>> v = pos / dur;
        CHECK(v.numerical_value_ref_in(km / h) == vector<double>{60, 40, 20});
      }
    }
  }

  SECTION("cross product with a vector quantity")
  {
    const auto r = vector<int>{3, 0, 0} * isq::displacement[m];
    const auto f = vector<int>{0, 10, 0} * isq::force[N];

    CHECK(cross_product(r, f) == vector<int>{0, 0, 30} * isq::moment_of_force[N * m]);
  }

#if 0
  SECTION("operator()")
  {
    SECTION("get")
    {
      const auto v = vector<int>{3, 2, 1} * isq::displacement[km];
      auto v0 = v(0);
      // auto v1 = v(1);
      // auto v2 = v(2);
      CHECK(v0.numerical_value_in(m) == 3000);
      // CHECK(v(1).numerical_value_in(m) == 2000);
      // CHECK(v(2).numerical_value_in(m) == 1000);
    }

    SECTION("set")
    {
      auto v = vector<int>{} * isq::displacement[m];
      v(0) =    7 * km;
      v(1) = 8000 * m;
      v(2) =  900 * cm;
      CHECK(v(0).numerical_value_in(m) == 7000);
      CHECK(v(1).numerical_value_in(m) == 8000);
      CHECK(v(2).numerical_value_in(m) ==    9);
    }
  }
#endif
}

TEST_CASE("vector of quantities", "[la]")
{
  SECTION("cast of unit")
  {
    SECTION("non-truncating")
    {
      const vector<quantity<isq::displacement[km], int>> v = {3 * km, 2 * km, 1 * km};

      CHECK(vector<quantity<isq::displacement[m], int>>(v) ==
            vector<quantity<isq::displacement[m], int>>{3000 * m, 2000 * m, 1000 * m});
    }

    // truncating not possible (no way to apply quantity_cast to sub-components of a vector)
  }

  SECTION("to scalar magnitude")
  {
    const vector<quantity<isq::velocity[km / h], int>> v = {2 * km / h, 3 * km / h, 6 * km / h};
    const auto speed = get_magnitude<isq::speed>(v);
    CHECK(speed.numerical_value_ref_in(km / h) == 7);
  }

  SECTION("element access")
  {
    SECTION("vector")
    {
      auto v = vector<int>{3, 1, 4} * isq::displacement[km];
      CHECK(v.unit == km);
      CHECK(v(0) == 3 * km);
      CHECK(v(1) == 1 * km);
      CHECK(v(2) == 4 * km);
      v(1) = 2 * km;
      CHECK(v(1) == 2 * km);
    }

    SECTION("3-port impedance matrix produces correct voltages")
    {
      // Make vectors from scalars

      // Z-matrix in ohms
      const auto raw_z = matrix<double>{ { 50.0,  2.0,  1.0 },
                                         {  2.0, 55.0,  3.0 },
                                         {  1.0,  3.0, 45.0 } };
      CHECK(raw_z(0, 0) == 50.0);
      CHECK(raw_z(1, 1) == 55.0);
      CHECK(raw_z(2, 2) == 45.0);

      auto z = quantity{raw_z, resistance_vec[ohm]};

      CHECK(z.unit == ohm);
      CHECK(z(0, 0) == 50.0 * ohm);
      CHECK(z(ohm)(1, 1) == 55.0);

      z(1, 1) = .045 * kohm;
      CHECK(z(1, 1) == 45.0 * ohm);

      CHECK(z(2, 2) == 45.0 * ohm);
      z(ohm)(2, 2) = 55.0;
      // z(kohm)(2, 2) == 55.0/1000 // won't compile
      CHECK(z(2, 2) == 55 * ohm);

      // Current vector in amperes
      constexpr auto raw_i = vector<double>{ {  2.0 },
                                             { -1.0 },
                                             {  0.5 } };

      constexpr auto i = quantity{raw_i, current_vec[si::ampere]};

      auto v = z * i; // Compute voltage vector

      constexpr auto expect = vector<double>{ {  98.5 },
                                              { -39.5 },
                                              {  26.5 } };
      // Verify units and numerical results
      CHECK(v(V) == expect);
      CHECK(v(2)(mV) == expect(2) * 1000);
      CHECK(v(1) == expect(1) * si::volt);
    }
  }

  SECTION("multiply by scalar value")
  {
    const vector<quantity<isq::displacement[m], int>> v = {1 * m, 2 * m, 3 * m};

    SECTION("integral")
    {
      const vector<quantity<isq::displacement[m], int>> result = {2 * m, 4 * m, 6 * m};

      SECTION("scalar on LHS") { CHECK(2 * v == result); }
      SECTION("scalar on RHS") { CHECK(v * 2 == result); }
    }

    SECTION("floating-point")
    {
      const vector<quantity<isq::displacement[m], double>> result = {0.5 * m, 1. * m, 1.5 * m};

      SECTION("scalar on LHS") { CHECK(0.5 * v == result); }
      SECTION("scalar on RHS") { CHECK(v * 0.5 == result); }
    }
  }

  SECTION("divide by scalar value")
  {
    const vector<quantity<isq::displacement[m], int>> v = {2 * m, 4 * m, 6 * m};

    SECTION("integral") { CHECK(v / 2 == vector<quantity<isq::displacement[m], int>>{1 * m, 2 * m, 3 * m}); }
    SECTION("floating-point")
    {
      CHECK(v / 0.5 == vector<quantity<isq::displacement[m], double>>{4. * m, 8. * m, 12. * m});
    }
  }

  SECTION("add")
  {
    const vector<quantity<isq::displacement[m], int>> v = {1 * m, 2 * m, 3 * m};

    SECTION("same unit")
    {
      const vector<quantity<isq::displacement[m], int>> u = {3 * m, 2 * m, 1 * m};

      CHECK(v + u == vector<quantity<isq::displacement[m], int>>{4 * m, 4 * m, 4 * m});
    }
    SECTION("different units")
    {
      const vector<quantity<isq::displacement[km], int>> u = {3 * km, 2 * km, 1 * km};

      CHECK(v + u == vector<quantity<isq::displacement[m], int>>{3001 * m, 2002 * m, 1003 * m});
    }
  }

  SECTION("subtract")
  {
    const vector<quantity<isq::displacement[m], int>> v = {1 * m, 2 * m, 3 * m};

    SECTION("same unit")
    {
      const vector<quantity<isq::displacement[m], int>> u = {3 * m, 2 * m, 1 * m};
      CHECK(v - u == vector<quantity<isq::displacement[m], int>>{-2 * m, 0 * m, 2 * m});
    }
    SECTION("different units")
    {
      const vector<quantity<isq::displacement[km], int>> u = {3 * km, 2 * km, 1 * km};
      CHECK(v - u == vector<quantity<isq::displacement[m], int>>{-2999 * m, -1998 * m, -997 * m});
    }
  }

#if 0
  SECTION("multiply by scalar quantity")
  {
    const vector<quantity<isq::velocity[m / s], int>> v = {1 * m / s, 2 * m / s, 3 * m / s};

    SECTION("integral")
    {
      const auto mass   = 2 * isq::mass[kg];
      const auto result = vector<quantity<isq::momentum[N * s], int>>{2 * N * s, 4 * N * s, 6 * N * s};

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK(mass * v == result); }
        SECTION("scalar on RHS") { CHECK(v * mass == result); }
      }

      // no way to apply quantity_cast to sub-components
      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const vector<quantity<isq::momentum[N * s], int>> momentum = mass * v;
          CHECK(momentum == result);
        }
        SECTION("scalar on RHS")
        {
          const vector<quantity<isq::momentum[N * s], int>> momentum = v * mass;
          CHECK(momentum == result);
        }
      }
    }

    SECTION("floating-point")
    {
      const auto mass = 0.5 * isq::mass[kg];
      const auto result = vector<quantity<isq::momentum[N * s], double>>{0.5 * N * s, 1. * N * s, 1.5 * N * s};

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK(mass * v == result); }
        SECTION("scalar on RHS") { CHECK(v * mass == result); }
      }

      // no way to apply quantity_cast to sub-components

      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const vector<quantity<isq::momentum[N * s], double>> momentum = mass * v;
          CHECK(momentum == result);
        }
        SECTION("scalar on RHS")
        {
          const vector<quantity<isq::momentum[N * s], double>> momentum = v * mass;
          CHECK(momentum == result);
        }
      }
    }
  }

  SECTION("divide by scalar quantity")
  {
    const vector<quantity<isq::displacement[km], int>> pos = {30 * km, 20 * km, 10 * km};

    SECTION("integral")
    {
      const auto dur = 2 * isq::duration[h];

      SECTION("derived_quantity_spec")
      {
        CHECK(pos / dur == vector<quantity<isq::velocity[km / h], int>>{15 * km / h, 10 * km / h, 5 * km / h});
      }

      // no way to apply quantity_cast to sub-components

      SECTION("quantity of velocity")
      {
        const vector<quantity<isq::velocity[km / h], int>> v = pos / dur;
        CHECK(v == vector<quantity<isq::velocity[km / h], int>>{15 * km / h, 10 * km / h, 5 * km / h});
      }
    }

    SECTION("floating-point")
    {
      const auto dur = 0.5 * isq::duration[h];

      SECTION("derived_quantity_spec")
      {
        CHECK(pos / dur == vector<quantity<isq::velocity[km / h], double>>{60. * km / h, 40. * km / h, 20. * km / h});
      }

      // no way to apply quantity_cast to sub-components

      SECTION("quantity of velocity")
      {
        const vector<quantity<isq::velocity[km / h], double>> v = pos / dur;
        CHECK(v == vector<quantity<isq::velocity[km / h], double>>{60. * km / h, 40. * km / h, 20. * km / h});
      }
    }
  }
#endif

  SECTION("cross product with a vector of quantities")
  {
    const vector<quantity<isq::displacement[m], int>> r = {3 * m, 0 * m, 0 * m};
    const vector<quantity<isq::force[N], int>> f = {0 * N, 10 * N, 0 * N};

    CHECK(cross_product(r, f) == vector<quantity<isq::moment_of_force[N * m], int>>{0 * N * m, 0 * N * m, 30 * N * m});
  }
}

TEST_CASE("vector quantity", "[tjg]")
{
  SECTION("cast of unit")
  {
    SECTION("non-truncating")
    {
      const auto v = tjg::Vector3D<int>{3, 2, 1} * isq::displacement[km];
      CHECK(v.numerical_value_in(m) == tjg::Vector3D<int>{3000, 2000, 1000});
    }

    SECTION("truncating")
    {
      const auto v = tjg::Vector3D<int>{1001, 1002, 1003} * isq::displacement[m];
      CHECK(v.force_numerical_value_in(km) == tjg::Vector3D<int>{1, 1, 1});
    }
  }
  SECTION("to scalar magnitude")
  {
    const auto v = tjg::Vector3D<int>{2, 3, 6} * isq::velocity[km / h];
    const auto speed = get_magnitude<isq::speed>(v);
    CHECK(speed.numerical_value_ref_in(km / h) == 7);
  }
  SECTION("element access")
  {
    auto v = tjg::Vector3D<int>{3, 1, 4} * isq::velocity[km / h];
    CHECK(v[0] == 3 * km / h);
    CHECK(v[1] == 1 * km / h);
    CHECK(v[2] == 4 * km / h);
    v[1] = 2 * km / h;
    CHECK(v[1] == 2 * km / h);
  }
  SECTION("multiply by scalar value")
  {
    const auto v = tjg::Vector3D<int>{1, 2, 3} * isq::displacement[m];

    SECTION("integral")
    {
      SECTION("scalar on LHS") { CHECK((2 * v).numerical_value_in(m) == tjg::Vector3D<int>{2, 4, 6}); }
      SECTION("scalar on RHS") { CHECK((v * 2).numerical_value_in(m) == tjg::Vector3D<int>{2, 4, 6}); }
    }

    SECTION("floating-point")
    {
      SECTION("scalar on LHS") { CHECK((0.5 * v).numerical_value_in(m) == tjg::Vector3D<double>{0.5, 1., 1.5}); }
      SECTION("scalar on RHS") { CHECK((v * 0.5).numerical_value_in(m) == tjg::Vector3D<double>{0.5, 1., 1.5}); }
    }
  }

  SECTION("divide by scalar value")
  {
    const auto v = tjg::Vector3D<int>{2, 4, 6} * isq::displacement[m];

    SECTION("integral") { CHECK((v / 2).numerical_value_in(m) == tjg::Vector3D<int>{1, 2, 3}); }
    SECTION("floating-point") { CHECK((v / 0.5).numerical_value_in(m) == tjg::Vector3D<double>{4., 8., 12.}); }
  }

  SECTION("add")
  {
    const auto v = tjg::Vector3D<int>{1, 2, 3} * isq::displacement[m];

    SECTION("same unit")
    {
      const auto u = tjg::Vector3D<int>{3, 2, 1} * isq::displacement[m];
      CHECK((v + u).numerical_value_in(m) == tjg::Vector3D<int>{4, 4, 4});
    }
    SECTION("different units")
    {
      const auto u = tjg::Vector3D<int>{3, 2, 1} * isq::displacement[km];
      CHECK((v + u).numerical_value_in(m) == tjg::Vector3D<int>{3001, 2002, 1003});
    }
  }

  SECTION("subtract")
  {
    const auto v = tjg::Vector3D<int>{1, 2, 3} * isq::displacement[m];

    SECTION("same unit")
    {
      const auto u = tjg::Vector3D<int>{3, 2, 1} * isq::displacement[m];
      CHECK((v - u).numerical_value_in(m) == tjg::Vector3D<int>{-2, 0, 2});
    }
    SECTION("different units")
    {
      const auto u = tjg::Vector3D<int>{3, 2, 1} * isq::displacement[km];
      CHECK((v - u).numerical_value_in(m) == tjg::Vector3D<int>{-2999, -1998, -997});
    }
  }

  SECTION("multiply by scalar quantity")
  {
    const auto v = tjg::Vector3D<int>{1, 2, 3} * isq::velocity[m / s];

    SECTION("integral")
    {
      const auto mass = 2 * isq::mass[kg];

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK((mass * v).numerical_value_in(kg * m / s) == tjg::Vector3D<int>{2, 4, 6}); }
        SECTION("scalar on RHS") { CHECK((v * mass).numerical_value_in(kg * m / s) == tjg::Vector3D<int>{2, 4, 6}); }
      }
      SECTION("quantity_cast to momentum")
      {
        SECTION("scalar on LHS")
        {
          CHECK(quantity_cast<isq::momentum>(mass * v).numerical_value_in(N * s) == tjg::Vector3D<int>{2, 4, 6});
        }
        SECTION("scalar on RHS")
        {
          CHECK(quantity_cast<isq::momentum>(v * mass).numerical_value_in(N * s) == tjg::Vector3D<int>{2, 4, 6});
        }
      }
      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const quantity<isq::momentum[N * s], tjg::Vector3D<int>> momentum = mass * v;
          CHECK(momentum.numerical_value_ref_in(N * s) == tjg::Vector3D<int>{2, 4, 6});
        }
        SECTION("scalar on RHS")
        {
          const quantity<isq::momentum[N * s], tjg::Vector3D<int>> momentum = v * mass;
          CHECK(momentum.numerical_value_ref_in(N * s) == tjg::Vector3D<int>{2, 4, 6});
        }
      }
    }

    SECTION("floating-point")
    {
      const auto mass = 0.5 * isq::mass[kg];

      SECTION("derived_quantity_spec")
      {
        SECTION("scalar on LHS") { CHECK((mass * v).numerical_value_in(kg * m / s) == tjg::Vector3D<double>{0.5, 1., 1.5}); }
        SECTION("scalar on RHS") { CHECK((v * mass).numerical_value_in(kg * m / s) == tjg::Vector3D<double>{0.5, 1., 1.5}); }
      }
      SECTION("quantity_cast to momentum")
      {
        SECTION("scalar on LHS")
        {
          CHECK(quantity_cast<isq::momentum>(mass * v).numerical_value_in(N * s) == tjg::Vector3D<double>{0.5, 1., 1.5});
        }
        SECTION("scalar on RHS")
        {
          CHECK(quantity_cast<isq::momentum>(v * mass).numerical_value_in(N * s) == tjg::Vector3D<double>{0.5, 1., 1.5});
        }
      }
      SECTION("quantity of momentum")
      {
        SECTION("scalar on LHS")
        {
          const quantity<isq::momentum[N * s], tjg::Vector3D<double>> momentum = mass * v;
          CHECK(momentum.numerical_value_ref_in(N * s) == tjg::Vector3D<double>{0.5, 1., 1.5});
        }
        SECTION("scalar on RHS")
        {
          const quantity<isq::momentum[N * s], tjg::Vector3D<double>> momentum = v * mass;
          CHECK(momentum.numerical_value_ref_in(N * s) == tjg::Vector3D<double>{0.5, 1., 1.5});
        }
      }
    }
  }

  SECTION("divide by scalar quantity")
  {
    const auto pos = tjg::Vector3D<int>{30, 20, 10} * isq::displacement[km];

    SECTION("integral")
    {
      const auto dur = 2 * isq::duration[h];

      SECTION("derived_quantity_spec") { CHECK((pos / dur).numerical_value_in(km / h) == tjg::Vector3D<int>{15, 10, 5}); }
      SECTION("quantity_cast to velocity")
      {
        CHECK(quantity_cast<isq::velocity>(pos / dur).numerical_value_in(km / h) == tjg::Vector3D<int>{15, 10, 5});
      }
      SECTION("quantity of velocity")
      {
        const quantity<isq::velocity[km / h], tjg::Vector3D<int>> v = pos / dur;
        CHECK(v.numerical_value_ref_in(km / h) == tjg::Vector3D<int>{15, 10, 5});
      }
    }

    SECTION("floating-point")
    {
      const auto dur = 0.5 * isq::duration[h];

      SECTION("derived_quantity_spec") { CHECK((pos / dur).numerical_value_in(km / h) == tjg::Vector3D<double>{60, 40, 20}); }
      SECTION("quantity_cast to velocity")
      {
        CHECK(quantity_cast<isq::velocity>(pos / dur).numerical_value_in(km / h) == tjg::Vector3D<double>{60, 40, 20});
      }
      SECTION("quantity of velocity")
      {
        const quantity<isq::velocity[km / h], tjg::Vector3D<double>> v = pos / dur;
        CHECK(v.numerical_value_ref_in(km / h) == tjg::Vector3D<double>{60, 40, 20});
      }
    }
  }

  SECTION("cross product with a vector quantity")
  {
    const auto r = tjg::Vector3D<int>{3, 0, 0} * isq::displacement[m];
    const auto f = tjg::Vector3D<int>{0, 10, 0} * isq::force[N];

    CHECK(cross_product(r, f) == tjg::Vector3D<int>{0, 0, 30} * isq::moment_of_force[N * m]);
  }

}
