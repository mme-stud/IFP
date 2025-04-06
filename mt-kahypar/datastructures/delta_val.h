
#pragma once

#include "mt-kahypar/macros.h" // for assertions

namespace mt_kahypar {
namespace ds {


/**
 * @brief Delta Value for adjusting keys in the conductance pririty queue 
 * during parallel calls of changeNodePart(u, from, to, ..)
 * 
 * Value is always non-negative. Boolean flag tells the sign.
 */
template <typename NonnegativeNumT>
class DeltaValue {
private:
  bool _negative = false;
  NonnegativeNumT _val;
public:
  DeltaValue() :
    _negative(false),
    _val(0)  { }

  DeltaValue(const NonnegativeNumT& val, const bool& negative = false) :
    _negative(negative),
    _val(val)  { }

  // ############################# Getters ################################
  bool isNegative() const {
    return _negative;
  }
  NonnegativeNumT abs() const {
    return _val;
  }
  friend std::ostream& operator<< (std::ostream& s, const DeltaValue& f) {
    if (f._negative) s << "-";
    return s << f._val;
  }

  // ########################### Arithmerics ##############################
  void operator+=(const NonnegativeNumT& add_val) {
    if (_negative) {
        if (add_val < _val) { // still negative
            _val -= add_val;
        } else {              // now non-negative
            _val = add_val - _val;
            _negative = false;
        }
    } else {
        ASSERT(static_cast<uintmax_t>(std::numeric_limits<NonnegativeNumT>::max()) 
                    - static_cast<uintmax_t>(_val)
                > static_cast<uintmax_t>(add_val), "Overflow");
        _val += add_val;      // still non-negative
    }
  }

  void operator-=(const NonnegativeNumT& sub_val) {
    if (_negative) {
        ASSERT(static_cast<uintmax_t>(std::numeric_limits<NonnegativeNumT>::max()) 
                    - static_cast<uintmax_t>(_val)
                > static_cast<uintmax_t>(sub_val), "Overflow");
        _val += sub_val;      // still negative
    } else {
        if (sub_val <= _val) {// still non-negative
            _val -= sub_val;
        } else {              // now negative
            _val = sub_val - _val;
            _negative = true;
        }
    }
  }
  
  void operator+=(const DeltaValue& other) {
    if (other._negative) {
      operator-= (other._val);
    } else { 
      operator+= (other._val);
    }
  }

  void operator-=(const DeltaValue& other) {
    if (other._negative) {
      operator+= (other._val);
    } else { 
      operator-= (other._val);
    }
  }

  DeltaValue operator-(const DeltaValue& other) const {
    DeltaValue result = *this;
    result -= other;
    return result;
  }

  DeltaValue operator+(const DeltaValue& other) const {
    DeltaValue result = *this;
    result += other;
    return result;
  }

  DeltaValue operator+(const NonnegativeNumT& add_val) const {
    DeltaValue result = *this;
    result += add_val;
    return result;
  }

  DeltaValue operator-(const NonnegativeNumT& sub_val) const {
    DeltaValue result = *this;
    result -= sub_val;
    return result;
  }

  // ########################## Comarators ################################
  bool operator<(const NonnegativeNumT& other_val) const {
    return (_negative || _val < other_val);
  }
  bool operator>(const NonnegativeNumT& other_val) const {
    return (!_negative && _val > other_val);
  }
  bool operator==(const NonnegativeNumT& other_val) const {
    return (!_negative && _val == other_val);
  }
  bool operator!=(const NonnegativeNumT& other_val) const {
    return (_negative || _val != other_val);
  }

  bool operator<(const DeltaValue& other) const {
    if (_negative != other._negative) return _negative;
    if (_negative) return _val > other._val;
    return _val < other._val;
  }
  bool operator>(const DeltaValue& other) const {
    if (_negative != other._negative) return !_negative;
    if (_negative) return _val < other._val;
    return _val > other._val;
  }
  bool operator==(const DeltaValue& other) const {
    return (_val == other._val && _negative == other._negative);
  }
  bool operator!=(const DeltaValue& other) const {
    return (_val != other._val || _negative != other._negative);
  }
};

}  // namespace ds
}  // namespace mt_kahypar