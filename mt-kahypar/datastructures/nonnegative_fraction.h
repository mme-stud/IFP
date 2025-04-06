
#pragma once

namespace mt_kahypar {
namespace ds {



/**
 * @brief Non-negative Fraction for the conductance priority queue.
 * Numerator and denominator are nonnegative numbers.
 * 0 / 0 = -inf, x / 0 = +inf
 */
template <typename Numerator>
class NonnegativeFraction {
  using Denominator = Numerator;
private:
  Numerator numerator;
  Denominator denominator;
  
  static const uintmax_t MAX_QUICK = std::numeric_limits<uint32_t>::max();
public:
  // infinity per default
  NonnegativeFraction() :
    numerator(1),
    denominator(0) { }

  NonnegativeFraction(const Numerator& n, const Denominator& d) :
    numerator(n),
    denominator(d) { 
      ASSERT(d >= 0);
      ASSERT(n >= 0);
  }
  
  // ############## Redefined Operators (<, ==, >, <<)  #####################

  bool operator< (const NonnegativeFraction& other) const {
    uint8_t null_result = isLess_nullCases(other);
    if (null_result == 2) {
      if (small() && other.small())
        return isLessQuick(other);
      return isLessSlow(*this, other);
    }
    return static_cast<bool>(null_result);
  }

  bool operator== (const NonnegativeFraction& other) const {
    uint8_t null_result = isEqual_nullCases(other);
    if (null_result == 2) {
      if (small() && other.small())
        return isEqualQuick(other);
      return isEqualSlow(*this, other);
    }
    return static_cast<bool>(null_result);
  }

  bool operator> (const NonnegativeFraction& other) const {
    uint8_t null_result = other.isLess_nullCases(*this);
    if (null_result == 2) {
      if (small() && other.small())
        return isGreaterQuick(other);
      return isLessSlow(other, *this);
    }
    return static_cast<bool>(null_result);
  }

  friend std::ostream& operator<< (std::ostream& s, const NonnegativeFraction& f) {
    return s << f.numerator << " / " << f.denominator;
  }

  // ####################### Getters and Setters ########################

  // Numerator must be non-negative 
  void setNumerator(const Numerator& n) {
    ASSERT(n >= 0);
    numerator = n;
  }
  // ! Denominator must be greater than 0
  void setDenominator(const Denominator& d) {
    ASSERT(d >= 0);
    denominator = d;
  }
  Numerator getNumerator() const {
    return numerator;
  }
  Denominator getDenominator() const {
    return denominator;
  }

  double_t value() const {
    if (denominator == 0) {
      return std::numeric_limits<double_t>::max();
    }
    return static_cast<double_t>(numerator) / static_cast<double_t>(denominator);
  }
private:  
// ############ Support of operators < = > for big / strange numbers #############

  // ! True if both numerator and denominator are 0
  bool null() const {
    return numerator == 0 && denominator == 0;
  }

  // ! 0 = false, 1 = true, 2 = not clear, but no zeroes
  int8_t isLess_nullCases(const NonnegativeFraction& other) const {
    bool this_null = null();
    bool other_null = other.null();
    if (this_null) {      // this = -inf
      if (other_null) {   // this = other
        return 0; // false
      } else {            // this = -inf < other
        return 1; // true
      }
    } // ##########################################
    if (other_null) {     // this > -inf = other
      return 0; // false
    }
    return 2; // not clear, but no nulls
  }
  
  // ! 0 = false, 1 = true, 2 = not clear, but no zeroes
  int8_t isEqual_nullCases(const NonnegativeFraction& other) const {
    bool this_null = null();
    bool other_null = other.null();
    if (this_null) {      // this = -inf
      if (other_null) {   // this = other
        return 1; // true
      } else {            // this = -inf < other
        return 0; // false
      }
    } // #########################################
    if (other_null) {     // this > -inf = other
      return 0; // false
    }
    return 2; // not clear, but no nulls
  }

  // ! True if both numerator and denominator are not greater than MAX_QUICK
  bool small() const {
    ASSERT(std::numeric_limits<uintmax_t>::max() >= MAX_QUICK * MAX_QUICK);
    return numerator <= MAX_QUICK && denominator <= MAX_QUICK;
  }

  Numerator gcd(Numerator a, Numerator b) const {
    while (b != 0) {
      a %= b;
      std::swap(a, b);
    }
    return a;
  }

  // ! Only for operator< : The fractions are never reduced automaticly
  void reduce() {
    Numerator d = gcd(numerator, denominator);
    ASSERT(d != 0);
    numerator /= d;
    denominator /= d;
  }

  // ! 0 = false, 1 = true, 2 = not clear, but no zeroes
  int8_t isLess_zeroCases(NonnegativeFraction lhs, NonnegativeFraction rhs) const {
    if (rhs.denominator == 0) {   // rhs = +inf
      if (lhs.denominator == 0) { // lhs = +inf = rhs
        return 0; // false
      } else {                    // lhs < +inf = rhs
        return 1; // true
      }
    } // ###################################################################
    if (lhs.denominator == 0) { // lhs = +inf , rhs != +inf 
      return 0; // false        // lhs = +inf !< rhs
    } // ###################################################################
    if (lhs.numerator == 0) {   // lhs = 0
      if (rhs.numerator == 0) { // rhs = 0
        return 0; // false
      } else {                  // rhs != 0
        return 1; // true
      }
    } // ###################################################################
    if (rhs.numerator == 0) {   // rhs = 0 , lhs != 0
      return 0; // false        // lhs !< 0 = rhs
    }
    return 2; // not clear, no zeroes
  }

  // ! 0 = false, 1 = true, 2 = not clear, but no zeroes
  int8_t isEqual_zeroCases(NonnegativeFraction lhs, NonnegativeFraction rhs) const {
    if (rhs.denominator == 0) { // rhs = +inf
      if (lhs.denominator == 0){// lhs = +inf = rhs
        return 1; // true
      } else {                  // lhs != +inf = rhs
        return 0; // false
      }
    } // ###################################################################
    if (lhs.denominator == 0) { // lhs = +inf , rhs != +inf 
      return 0; // false        // lhs = +inf != rhs
    } // ###################################################################
    if (lhs.numerator == 0) {   // lhs = 0
      if (rhs.numerator == 0) { // rhs = 0
        return 1; // true
      } else {                  // rhs != 0 = lhs
        return 0; // false
      }
    } // ###################################################################
    if (rhs.numerator == 0) {   // rhs = 0 , lhs != 0
      return 0; // false        // lhs != 0 = rhs
    }
    return 2; // not clear, no zeroes
  }  

  // ! Quick version of operator< for fractions with small numerator and denominator
  bool isLessQuick(const NonnegativeFraction& other) const {
    ASSERT(small() && other.small());
    uintmax_t lhs = static_cast<uintmax_t>(numerator) * static_cast<uintmax_t>(other.denominator);
    uintmax_t rhs = static_cast<uintmax_t>(other.numerator) * static_cast<uintmax_t>(denominator);
    return lhs < rhs;
  }

  // ! Quick version of operator< for fractions with small numerator and denominator
  bool isGreaterQuick(const NonnegativeFraction& other) const {
    ASSERT(small() && other.small());
    uintmax_t lhs = static_cast<uintmax_t>(numerator) * static_cast<uintmax_t>(other.denominator);
    uintmax_t rhs = static_cast<uintmax_t>(other.numerator) * static_cast<uintmax_t>(denominator);
    return lhs > rhs;
  }
  
  // ! Quick version of operator< for fractions with small numerator and denominator
  bool isEqualQuick(const NonnegativeFraction& other) const {
    ASSERT(small() && other.small());
    uintmax_t lhs = static_cast<uintmax_t>(numerator) * static_cast<uintmax_t>(other.denominator);
    uintmax_t rhs = static_cast<uintmax_t>(other.numerator) * static_cast<uintmax_t>(denominator);
    return lhs == rhs;
  }

  // ! Slower version of operator< for fractions with big numerator or denominator
  bool isLessSlow(NonnegativeFraction lhs, NonnegativeFraction rhs) const {
    int8_t result_zero_cases = isLess_zeroCases(lhs, rhs);
    while (result_zero_cases == 2) { // 2: not clear but no zeroes or +inf
      if (lhs.small() && rhs.small()) {
        return lhs.isLessQuick(rhs);
      }
      Numerator int_lhs = lhs.numerator / lhs.denominator;
      Numerator int_rhs = rhs.numerator / rhs.denominator;
      if (int_lhs != int_rhs) {
        return int_lhs < int_rhs;
      }
      lhs.reduce();   rhs.reduce();
      lhs.numerator %= lhs.denominator;
      rhs.numerator %= rhs.denominator;
      std::swap(lhs.numerator, rhs.denominator);
      std::swap(lhs.denominator, rhs.numerator);
      result_zero_cases = isLess_zeroCases(lhs, rhs);
    }
    return static_cast<bool>(result_zero_cases); // 0 = false , 1 = true
  } 

  // ! Slower version of operator= for fractions with big numerator or denominator
  bool isEqualSlow(NonnegativeFraction lhs, NonnegativeFraction rhs) const {
    int8_t result_zero_cases = isEqual_zeroCases(lhs, rhs);
    while (result_zero_cases == 2) { // 2: not clear but no zeroes or +inf
      if (lhs.small() && rhs.small()) {
        return lhs.isEqualQuick(rhs);
      }
      Numerator int_lhs = lhs.numerator / lhs.denominator;
      Numerator int_rhs = rhs.numerator / rhs.denominator;
      if (int_lhs != int_rhs) {
        return false;
      }
      lhs.reduce();   rhs.reduce();
      lhs.numerator %= lhs.denominator;
      rhs.numerator %= rhs.denominator;
      std::swap(lhs.numerator, rhs.denominator);
      std::swap(lhs.denominator, rhs.numerator);
      result_zero_cases = isEqual_zeroCases(lhs, rhs);
    }
    return static_cast<bool>(result_zero_cases); // 0 = false , 1 = true
  }
};

}  // namespace ds
}  // namespace mt_kahypar