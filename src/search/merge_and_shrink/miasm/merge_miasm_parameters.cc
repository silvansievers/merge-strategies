#include "merge_miasm_parameters.h"

#define X(a) # a

DEFINE_ENUM_OPT(MiasmInternal, "miasm_internal", LEVEL)

DEFINE_ENUM_OPT(MiasmExternal, "miasm_external", NUM_VAR_CGL)

#undef X
