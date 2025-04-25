#include <algorithm>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <vector>

namespace mp = boost::multiprecision;

/**
 * @brief Split a double number into two parts: high and low
 * @param a - number to split
 * @return a pair of high and low parts
 *
 * The algorithm is based on the idea of multiplying a number by a specific
 * value (splitter) and then subtracting the integer part of the result from
 * the original number. The integer part of the result is the high part of the
 * number, and the low part is the result of the subtraction.
 */
std::pair<double, double> Split(double a) {
  constexpr double splitter =
      (1 << (std::numeric_limits<double>::digits / 2 + 1)) + 1;

  if (std::isinf(a)) {
    return {a, 0.0};
  } else if (std::isnan(a)) {
    return {a, a};
  }

  double temp = a * splitter;
  double a_high = temp - (temp - a);
  double a_low = a - a_high;

  return {a_high, a_low};
}

/**
 * @brief Split a double number into two parts: high and low and then
 * multiply each part of the first number with each part of the second number
 * @param a - first number
 * @param b - second number
 * @return a pair of the product and the error
 *
 * The algorithm is based on the idea of multiplying each part of the first
 * number with each part of the second number and then summing up all the
 * products. The result is the product of the two numbers, and the error is
 * the difference between the exact product and the result.
 */
std::pair<double, double> Product(double a, double b) {
  auto [a_high, a_low] = Split(a);
  auto [b_high, b_low] = Split(b);

  double production = a * b;
  if (std::isinf(production)) {
    return {production, 0.0};
  } else if (std::isnan(production)) {
    return {production, production};
  }

  double error =
      ((a_high * b_high - production) + a_high * b_low + a_low * b_high) +
      a_low * b_low;

  return {production, error};
}

/**
 * @brief Compute the sum of two numbers with error correction
 * @param a - first number
 * @param b - second number
 * @return a pair of the sum and the error
 *
 * The algorithm is based on the idea of computing the sum and then subtracting
 * each number from the result in order to compute the round-off error. The
 * error is then added to the result to get the corrected sum.
 */
std::pair<double, double> Sum(double a, double b) {
  double sum = a + b;

  if (std::isinf(sum)) {
    return {sum, 0.0};
  } else if (std::isnan(sum)) {
    return {sum, sum};
  }

  double b_virtual = sum - a;
  double a_virtual = sum - b_virtual;
  double b_roundoff = b - b_virtual;
  double a_roundoff = a - a_virtual;
  double error = a_roundoff + b_roundoff;

  return {sum, error};
}

/**
 * @brief Compute the sum of all numbers in a vector using precise floating
 * point arithmetic.
 *
 * The algorithm is based on the idea of splitting each number into two parts:
 * high and low. Then the addition is performed for each part of the first
 * number with each part of the second number. The results are then summed up.
 * The error is the difference between the exact sum and the result.
 *
 * @param numbers The vector of numbers to sum.
 * @return A pair of the sum and the error.
 */
std::pair<double, double> Sum(const std::vector<double>& numbers) {
  bool has_infinity = false;
  bool has_negative_infinity = false;
  std::vector<double> partials;

  for (double number : numbers) {
    if (std::isnan(number))
      return {std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
    else if (std::isinf(number)) {
      if (number > 0)
        has_infinity = true;
      else
        has_negative_infinity = true;

      if (has_infinity && has_negative_infinity)
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    for (double& partial : partials) {
      if (std::abs(number) > std::abs(partial)) {
        std::swap(number, partial);
      }

      auto [sum, error] = Sum(partial, number);
      partial = sum;
      number = error;

      if (error == 0.0) break;
    }

    if (number != 0.0) {
      partials.push_back(number);
    }
  }

  if (has_infinity)
    return {std::numeric_limits<double>::infinity(), 0.0};
  else if (has_negative_infinity)
    return {-std::numeric_limits<double>::infinity(), 0.0};

  double total = 0.0;
  double compensation = 0.0;

  for (double p : partials) {
    double y = p - compensation;
    double t = total + y;
    compensation = (t - total) - y;
    total = t;
  }

  return {total, -compensation};
}

/**
 * @brief Compute the dot product of two vectors using precise floating point
 * arithmetic
 *
 * @param a First vector.
 * @param b Second vector.
 * @return The dot product of @a a and @a b.
 *
 * The algorithm is based on the idea of splitting each number into two parts:
 * high and low. Then the multiplication is performed for each part of the
 * a number with each part of the b number. The results are then
 * summed up. The error is the difference between the exact product and the
 * result.
 */
double DotProduct(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  std::vector<double> partials;
  for (size_t i = 0; i < a.size(); ++i) {
    auto [product, error] = Product(a[i], b[i]);
    partials.push_back(product);
    partials.push_back(error);
  }

  auto [dot_product, _] = Sum(partials);
  return dot_product;
}

/**
 * @brief Compute the dot product of two vectors using Boost.Multiprecision.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return The dot product of @a a and @a b.
 */
double BoostDotProduct(const std::vector<double>& a,
                       const std::vector<double>& b) {
  mp::cpp_dec_float_100 sum = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    auto product = mp::cpp_dec_float_100(a[i]) * mp::cpp_dec_float_100(b[i]);
    auto sum_double = sum.convert_to<double>();
    auto product_double = product.convert_to<double>();
    if (std::isnan(sum_double + product_double)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    sum += product;
  }
  return sum.convert_to<double>();
}

/**
 * @brief Check if two floating point numbers are equal up to a certain number
 * of ULPs
 * @param a - first number
 * @param b - second number
 * @param max_ulp_diff - maximum allowed difference in ULPs (default is 1)
 * @return true if numbers are equal up to specified number of ULPs, false
 * otherwise
 */
bool BitwiseAlmostEqual(double a, double b, int max_ulp_diff = 0) {
  if (a == b) return true;

  if (std::isnan(a) && std::isnan(b)) {
    return true;
  }

  if (a == 0.0 && b == 0.0) {
    return std::signbit(a) == std::signbit(b);
  }

  int64_t a_bits, b_bits;
  memcpy(&a_bits, &a, sizeof(double));
  memcpy(&b_bits, &b, sizeof(double));

  int64_t diff = std::abs(a_bits - b_bits);
  return diff <= max_ulp_diff;
}

/**
 * @brief Generate a set of test cases for vector dot product computations.
 *
 * This function generates a variety of test cases that include vectors with
 * different characteristics, such as medium-sized numbers, very small numbers,
 * alternating signs, denormalized numbers, combinations of very large and very
 * small numbers, self-destructive giants with small numbers, and vectors
 * containing many small numbers with one large number. These test cases are
 * used to evaluate the precision and robustness of dot product calculations.
 *
 * @return A vector of test cases, where each test case is a pair of vectors.
 */
std::vector<std::pair<std::vector<double>, std::vector<double>>>
GenerateTestCases() {
  std::vector<std::pair<std::vector<double>, std::vector<double>>> test_cases;
  std::mt19937_64 rng(std::random_device{}());

  std::uniform_real_distribution<double> big_dist(-1e300, 1e300);
  std::uniform_real_distribution<double> small_dist(-1e-300, 1e-300);

  auto precise_dist = [&rng]() {
    std::uniform_int_distribution<int64_t> int_dist(
        -9999999999999999ll,  // -1e16
        9999999999999999ll    // +1e16
    );
    // Scaling to [-1.0, 1.0]
    return static_cast<double>(int_dist(rng)) * 1e-16;
  };

  // Medium-sized numbers
  for (int i = 0; i < 100; ++i) {
    std::vector<double> a(1000), b(1000);
    for (int j = 0; j < 1000; ++j) {
      // 1e-15 - 1e15
      double scale = std::pow(10.0, static_cast<int>(rng() % 30) - 15);
      a[j] = precise_dist() * scale;
      b[j] = precise_dist() * scale;
    }
    test_cases.emplace_back(a, b);
  }

  /* Very small numbers */
  for (int i = 0; i < 100; ++i) {
    std::vector<double> a(1000), b(1000);
    for (int j = 0; j < 1000; ++j) {
      // 1e-30 - 1e-15
      a[j] = precise_dist() * std::pow(10.0, static_cast<int>(rng() % 30) - 30);
      b[j] = precise_dist() * std::pow(10.0, static_cast<int>(rng() % 30) - 30);
    }
    test_cases.emplace_back(a, b);
  }

  /* Alternating signs */
  {
    std::vector<double> a(1000), b(1000);
    for (int j = 0; j < 1000; ++j) {
      a[j] = (j % 2) ? 1.0 : -1.0;
      b[j] = 1.0 / static_cast<double>(j + 1);
    }
    test_cases.emplace_back(a, b);
  }

  /* Denormalized numbers */
  {
    std::vector<double> a{std::numeric_limits<double>::denorm_min()};
    std::vector<double> b{std::numeric_limits<double>::denorm_min()};
    test_cases.emplace_back(a, b);
  }

  /* Very large and very small */
  {
    std::vector<double> a{big_dist(rng), small_dist(rng)};
    std::vector<double> b{small_dist(rng), big_dist(rng)};
    test_cases.emplace_back(a, b);
  }

  /* Self-destructive giants with small numbers */
  {
    std::vector<double> a{1e100, 1.0, small_dist(rng), -1e100, 1.0};
    std::vector<double> b{1.0, 1e100, big_dist(rng), 1.0, -1e100};
    test_cases.emplace_back(a, b);
  }

  /* Many small numbers and one big */
  {
    std::vector<double> a(1000, 1e-100);
    std::vector<double> b(1000, 1e-100);
    a.back() = 1e100;
    b.back() = 1e100;
    test_cases.emplace_back(a, b);
  }

  return test_cases;
}

/**
 * @brief Runs tests for the DotProduct function
 *
 * The function generates test cases using GenerateTestCases and then runs
 * the DotProduct function on each test case and checks the result against
 * the result of BoostDotProduct. If the results are equal up to a certain
 * number of ULPs (default is 1), the test is considered passed. The results
 * of the tests are then printed to the console.
 *
 * @see GenerateTestCases
 * @see DotProduct
 * @see BoostDotProduct
 * @see BitwiseAlmostEqual
 */
void RunTests() {
  auto test_cases = GenerateTestCases();

  int passed = 0;
  int total = 0;

  std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  for (const auto& [a, b] : test_cases) {
    double expected = BoostDotProduct(a, b);
    double result = DotProduct(a, b);

    std::cout << "======================================================"
              << std::endl;
    std::cout << "Expected: " << expected << std::endl;
    std::cout << "Got:      " << result << std::endl;

    if (BitwiseAlmostEqual(result, expected)) {
      std::cout << "Test passed!" << std::endl;
      ++passed;
    } else {
      std::cout << "Test failed!" << std::endl;
    }

    std::cout << "First 5 elements:" << std::endl;
    for (int i = 0; i < std::min(5ul, a.size()); ++i) {
      std::cout << a[i] << " * " << b[i] << " = "
                << mp::cpp_dec_float_100(a[i]) * mp::cpp_dec_float_100(b[i])
                << std::endl;
    }

    ++total;
  }

  std::cout << "Passed: " << passed << "/" << total << " ("
            << (100.0 * passed / total) << "%)" << std::endl;
}

int main() {
  RunTests();
  return 0;
}