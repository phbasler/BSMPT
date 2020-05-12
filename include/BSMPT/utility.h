#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <iostream>
#include <vector>

/**
 * @file
 */
namespace BSMPT {
/**
 * @brief seperator used to write into output files
 */
const std::string sep = "\t";

/**
 * Overload to print out vectors with the << operator
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    bool first = true;
    for (const auto& el : vec)
    {
        if(not first) {
            os << sep;
        }
        else{
            first = false;
        }
        os << el;
    }
    return os;
}


/**
 * @brief operator << overload for the model parameter
 */
namespace ModelID {
  enum class ModelIDs;
}
std::ostream& operator<<(std::ostream& os, const ModelID::ModelIDs& Model);


}

#endif // UTILITY_H
