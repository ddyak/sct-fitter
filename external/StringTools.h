/**************************************************************************
 * AURORA (SCT Framework)                                                 *
 * Copyright(C) 2019 - SCT Collaboration                                  *
 *                                                                        *
 * Author: The SCT Collaboration                                          *
 * Contributors: Vitaly Vorobyev                                          *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/
#pragma once

#include <string_view>
#include <optional>
#include <utility>
#include <vector>

namespace sct::comm {

class ST {
 public:
    static std::string_view Strip(std::string_view sw);
    static std::pair<std::optional<std::string_view>, std::optional<std::string_view>> SplitAndStripTwo(std::string_view s, std::string_view delimiter);
    static std::vector<std::string_view> SplitAndStrip(std::string_view s, std::string_view delimiter);
    /** Returns position of the matched closing parenthesis if the first character in the given
      * string contains an opening parenthesis. Otherwise return 0. */
    static size_t findMatchedParenthesis(std::string_view str, char open='[', char close=']');
    /** Returns the position of a pattern in a string ignoring everything that is in parenthesis. */
    static size_t findIgnoringParenthesis(std::string_view str, std::string_view pattern, size_t begin=0);
    /** Split into std::vector on delimiter ignoring delimiters between parenthesis */
    static std::vector<std::string_view> splitOnDelimiterAndConserveParenthesis(std::string_view str, char delimiter, char open, char close);
};

}  // namespace sct::comm
