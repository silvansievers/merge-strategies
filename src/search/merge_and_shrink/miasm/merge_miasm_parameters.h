/** @file */
#ifndef MERGE_MIASM_PARAMETERS_H
#define MERGE_MIASM_PARAMETERS_H

#include "option_struct.h"

#define X(a) a

#define MiasmInternalTable \
    /** by variable level */ \
    X(LEVEL), \
    /** by reverse variable level */ \
    X(REVERSE_LEVEL)

/** @brief \link #DECLARE_ENUM_OPT Enum-Option-Struct \endlink
 * defines MIASM's internal merging strategies */
DECLARE_ENUM_OPT(MiasmInternal);


#define MiasmExternalTable \
    /** non-singleton set before singleton set; <br> */ \
    /** causal graph predecssor before goal; <br> */ \
    /** incresing number of variables; <br> */ \
    /** level of the smallest variable */ \
    X(NUM_VAR_CGL), \
    /** TODO */ \
    X(RNR_SIZE_CGL), \
    /** causal graph predecssor before goal; <br> */ \
    /** reverse level of the smallest variable */ \
    X(CGRL)

/** @brief \link #DECLARE_ENUM_OPT Enum-Option-Struct \endlink
 * defines MIASM's external merging strategies */
DECLARE_ENUM_OPT(MiasmExternal);

#undef X

#endif // MERGE_MIASM_PARAMETERS_H
