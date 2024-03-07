#include <float.h> 
#include <limits>

template<typename T>
bool is_infinite( const T &value )
{

    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;
 
    return ! ( min_value <= value && value <= max_value );
}
 
template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}
 
template<typename T>
bool is_valid( const T &value )
{
    return ! is_infinite(value) && ! is_nan(value);
}