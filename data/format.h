#ifndef __FORMAT_H__

#define __FORMAT_H__

#include <iostream>
#include <algorithm>
#include <string>

enum class Format {
    GREG,
    JT, // Jackson's data
    RABBOTT, // My data
    RAW, // Output from measurment program
    RAW_V2, // Output from the newer, more modular measurment program
    RAW_WITH_MIXED, // Output from measurment program with additional measurments on color-mixed diagrams
    ZIYUAN,
};

inline std::ostream& operator<< (std::ostream& os, const Format& f) {
    switch (f) {
        case Format::GREG: os << "Greg"; break;
        case Format::JT: os << "JT"; break;
        case Format::RABBOTT: os << "Rabbott"; break;
        case Format::RAW: os << "Raw"; break;
        case Format::RAW_V2: os << "Raw (Version 2)"; break;
        case Format::RAW_WITH_MIXED: os << "Raw (with mixed)"; break;
        case Format::ZIYUAN: os << "Ziyuan"; break;
    }
    return os;
}

static inline Format parse_format(std::string str) {
    // Make the string all lowercase
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

    if (str == "rabbott" || str == "rab") {
        return Format::RABBOTT;
    }
    if (str == "greg") {
        return Format::GREG;
    }
    if (str == "jt") {
        return Format::JT;
    }
    if (str == "raw") {
        return Format::RAW;
    }
    if (str == "raw_v2") {
        return Format::RAW_V2;
    }
    if (str == "ziyuan" || str == "daiqian") {
        return Format::ZIYUAN;
    }

    std::cerr << "Error: unkown format '" << str << "'" << std::endl;
    exit(1);
}


#endif
