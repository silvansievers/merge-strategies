/** @file */
#ifndef OPTION_STRUCT_H
#define OPTION_STRUCT_H

#include <vector>
#include <string>

/** The generic declaration of a struct wrapping some basic type, e.g., int,
 * that is used as options in fast-downward */
#define DECLARE_OPT(T, OPT_STRUCT) \
    struct OPT_STRUCT { \
        static std::string opt_key(); \
        static std::string def_val(); \
        OPT_STRUCT(const T value_ = 0); \
        operator T() const; \
        const T value; \
    }

/** The generic definition of a struct wrapping some basic type, e.g., int,
 * that is used as options in fast-downward */
#define DEFINE_OPT(T, OPT_STRUCT, KEY, DEF) \
    OPT_STRUCT::OPT_STRUCT(const T value_) : value(value_) {} \
    OPT_STRUCT::operator T() const {return value; } \
    std::string OPT_STRUCT::opt_key() {return KEY; } \
    std::string OPT_STRUCT::def_val() {return DEF; }

/** @brief The generic declaration of a struct wrapping an enum
 * that is used as options in fast-downward */
#define DECLARE_ENUM_OPT(ENUM_OPT_STRUCT) \
    struct ENUM_OPT_STRUCT { \
        enum E { \
            ENUM_OPT_STRUCT ## Table \
        }; \
        static const char *C[]; \
        static std::vector<std::string> S(); \
        static std::string option_key(); \
        static std::string default_value(); \
        ENUM_OPT_STRUCT(const int e_int = 0); \
        ENUM_OPT_STRUCT(const E e_ = E(0)); \
        operator E() const; \
        const E e; \
    }

/** @brief The generic definition of a struct wrapping an enum
 * that is used as options in fast-downward */
#define DEFINE_ENUM_OPT(ENUM_OPT_STRUCT, KEY, DEF) \
    const char *ENUM_OPT_STRUCT::C[] = {ENUM_OPT_STRUCT ## Table}; \
    ENUM_OPT_STRUCT::ENUM_OPT_STRUCT(const int e_int) : e(E(e_int)) {} \
    ENUM_OPT_STRUCT::ENUM_OPT_STRUCT(const E e_) : e(e_) {} \
    ENUM_OPT_STRUCT::operator E() const {return e; } \
    std::string ENUM_OPT_STRUCT::option_key() {return KEY; } \
    std::string ENUM_OPT_STRUCT::default_value() {return C[DEF]; } \
    std::vector<std::string> ENUM_OPT_STRUCT::S() { \
        return std::vector<std::string>(C, C + sizeof(C) / sizeof(char *)); \
    }

#endif // OPTION_STRUCT_H
