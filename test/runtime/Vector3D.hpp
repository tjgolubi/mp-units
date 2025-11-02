#pragma once
#include <iostream>
#include <format>
#include <array>
#include <type_traits>
#include <cmath>

namespace tjg {

template<typename T1, typename T2>
using RetT = std::common_type_t<T1, T2>;

template<typename T = double>
struct Vector3D: std::array<T, 3>
{
  using Base = std::array<T, 3>;
  using Base::Base;  // keep aggregate-like construction, e.g. {x,y,z}

  Vector3D() = default;
  Vector3D(const T& x, const T& y, const T& z)
    : Base{x, y, z} { }

  Vector3D(const Vector3D&) = default;
  Vector3D& operator=(const Vector3D&) = default;

  template<typename U>
  constexpr explicit(!std::convertible_to<U, T>)
  Vector3D(const Vector3D<U>& v) noexcept
    : Vector3D{static_cast<T>(v[0]),
               static_cast<T>(v[1]),
               static_cast<T>(v[2])}
  { }

  // Converting assignment: Vector3D<U> -> Vector3D<T>
  template<class U> requires std::convertible_to<U, T>
  constexpr Vector3D& operator=(const Vector3D<U>& v) noexcept {
    return *this = Vector3D{static_cast<T>(v[0]),
                            static_cast<T>(v[1]),
                            static_cast<T>(v[2])};
  }
}; // Vector3D

template<typename T, typename U>
constexpr bool operator==(const Vector3D<T>& lhs, const Vector3D<U>& rhs) noexcept
{
  return lhs[0] == rhs[0]
      && lhs[1] == rhs[1]
      && lhs[2] == rhs[2];
}

template<typename T> requires std::is_arithmetic_v<T>
[[nodiscard]] constexpr Vector3D<T> operator-(const Vector3D<T>& v) noexcept
  { return {-v[0], -v[1], -v[2]}; }

template<typename T, typename S> requires std::is_arithmetic_v<S>
[[nodiscard]] constexpr Vector3D<RetT<T, S>> operator*(const Vector3D<T>& v, const S& s) noexcept
{
  using R = RetT<T, S>;
  return { static_cast<R>(v[0] * s),
           static_cast<R>(v[1] * s),
           static_cast<R>(v[2] * s)};
}

template<typename S, typename T> requires std::is_arithmetic_v<S>
[[nodiscard]] constexpr Vector3D<RetT<T, S>> operator*(const S& s, const Vector3D<T>& v) noexcept
  { return v * s; }

template<typename T, typename S> requires std::is_arithmetic_v<S>
[[nodiscard]] constexpr Vector3D<RetT<T, S>> operator/(const Vector3D<T>& v, const S& s) noexcept
{
  using R = RetT<T, S>;
  return {static_cast<R>(v[0] / s),
          static_cast<R>(v[1] / s),
          static_cast<R>(v[2] / s)};
}

template<typename T, typename U>
[[nodiscard]] constexpr Vector3D<RetT<T, U>> operator+(const Vector3D<T>& a, const Vector3D<U>& b) noexcept
{
  using R = RetT<T, U>;
  return {static_cast<R>(a[0] + b[0]),
          static_cast<R>(a[1] + b[1]),
          static_cast<R>(a[2] + b[2])};
}

template<typename T, typename U>
[[nodiscard]] constexpr Vector3D<RetT<T, U>> operator-(const Vector3D<T>& a, const Vector3D<U>& b) noexcept
{
  using R = RetT<T, U>;
  return {static_cast<R>(a[0] - b[0]),
          static_cast<R>(a[1] - b[1]),
          static_cast<R>(a[2] - b[2])};
}

template<typename T, typename U>
[[nodiscard]] constexpr RetT<T, U> dot_product(const Vector3D<T>& a, const Vector3D<U>& b)
{
  return static_cast<RetT<T, U>>(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

template<typename T, typename U>
[[nodiscard]] constexpr Vector3D<RetT<T, U>> cross_product(const Vector3D<T>& a, const Vector3D<U>& b)
{
  using R = RetT<T, U>;
  return {static_cast<R>(a[1] * b[2] - a[2] * b[1]),
          static_cast<R>(a[2] * b[0] - a[0] * b[2]),
          static_cast<R>(a[0] * b[1] - a[1] * b[0])};
}

template<typename T>
[[nodiscard]] constexpr T magnitude(const Vector3D<T>& v) noexcept
{
  using std::hypot;
  return static_cast<T>(hypot(v[0], v[1], v[2]));
} // magnitude

template<typename T>
[[nodiscard]] constexpr auto direction(const Vector3D<T>& v) noexcept
{
  using std::atan2;
  return atan2(v[1], v[0]);
} // direction

template<typename T>
[[nodiscard]] constexpr auto elevation(const Vector3D<T>& v) noexcept
{
  using std::hypot;
  using std::atan2;
  return atan2(v[2], hypot(v[1], v[0]));
} // elevation

template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector3D<T>& v)
{
  os << "|";
  for (auto i = 0U; i < v.size(); ++i)
    os << std::format(" {:>9}", v[i]);
  os << " |";
  return os;
}

} // tjg
