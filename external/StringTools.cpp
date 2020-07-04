#include "StringTools.h"

#include <stdexcept>

using namespace std;

namespace sct::comm {

string_view ST::Strip(string_view sw) {
    while (sw.back() == ' ') sw.remove_suffix(1);
    while (sw.front() == ' ') sw.remove_prefix(1);
    return sw;
}

pair<optional<string_view>, optional<string_view>> ST::SplitAndStripTwo(string_view s, string_view delimiter) {
    const size_t pos = s.find(delimiter);
    if (pos == string::npos) {
        return {Strip(s), nullopt};
    } else if (pos == 0) {
        return {nullopt, Strip(s.substr(delimiter.length()))};
    } else {
        return {Strip(s.substr(0, pos)), Strip(s.substr(pos + delimiter.length()))};
    }
}

vector<string_view> ST::SplitAndStrip(string_view s, string_view delimiter) {
    vector<string_view> parts;
    if (s.empty()) {
        return parts;
    }
    while (true) {
        const auto [lhs_opt, rhs_opt] = SplitAndStripTwo(s, delimiter);
        if (lhs_opt.has_value()) {
            parts.push_back(Strip(lhs_opt.value()));
        }
        if (!rhs_opt) {
            break;
        }
        s = *rhs_opt;
    }
    return parts;
}

size_t ST::findMatchedParenthesis(string_view str, char open, char close) {
    size_t end = 1;
    if (str[0] == open) {
        size_t count = 1;
        for (end = 1; end < str.size() and count > 0; ++end) {
            if (str[end] == open) ++count;
            else if (str[end] == close) --count;
        }

        if (count > 0)
            throw runtime_error("Variable string has an invalid format: " + string(str));
    }
    return end - 1;
}

size_t ST::findIgnoringParenthesis(string_view str, string_view pattern, size_t begin) {
    if (str.size() < pattern.size())
        return string::npos;

    for (size_t i = begin; i < str.size() - pattern.size(); ++i) {
        if (str[i] == '[') {
            i += findMatchedParenthesis(str.substr(i), '[', ']');
            continue;
        }
        if (str[i] == '(') {
            i += findMatchedParenthesis(str.substr(i), '(', ')');
            continue;
        }
        if (str[i] == '{') {
            i += findMatchedParenthesis(str.substr(i), '{', '}');
            continue;
        }

        for (size_t j = 0; j < pattern.size(); ++j) {
            if (str[i + j] != pattern[j]) {
                break;
            }
            if (j == pattern.size() - 1) {
            return i;
            }
        }
    }
    return string::npos;
}

vector<string_view> ST::splitOnDelimiterAndConserveParenthesis(string_view str, char delimiter, char open, char close) {
    vector<string_view> result;
    size_t lastdelimiter = 0;
    for (size_t i = 0; i < str.size(); ++i) {
        if (str[i] == open) {
            i += findMatchedParenthesis(str.substr(i), open, close);
            continue;
        }
        if (str[i] == delimiter) {
            result.push_back(str.substr(lastdelimiter, i - lastdelimiter));
            lastdelimiter = i + 1;
        }
    }
    string_view last = str.substr(lastdelimiter);
    if (last.size() != 0) {
        result.push_back(last);
    }
    return result;
}

}  // namespace sct::comm
