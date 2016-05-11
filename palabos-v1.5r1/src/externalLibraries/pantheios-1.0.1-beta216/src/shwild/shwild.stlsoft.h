/* /////////////////////////////////////////////////////////////////////////
 * File:        shwild.stlsoft.h
 *
 * Purpose:     Includes the STLSoft library header(s), and verifies
 *              required version.
 *
 * Created:     11th February 2009
 * Updated:     10th August 2010
 *
 * Author:      Matthew Wilson
 *
 * Home:        http://www.shwild.org/
 *
 * Copyright (c) 2009-2010, Matthew Wilson
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted in accordance with the license and warranty
 * information described in shwild.h (included in this distribution, or
 * available from http://shwild.org/)
 *
 * ////////////////////////////////////////////////////////////////////// */


#ifndef SHWILD_INCL_H_SHWILD_STLSOFT
#define SHWILD_INCL_H_SHWILD_STLSOFT

/** \file shwild.stlsoft.h Includes STLSoft
 */

/* /////////////////////////////////////////////////////////////////////////
 * Includes
 */

#include <stlsoft/stlsoft.h>

/* /////////////////////////////////////////////////////////////////////////
 * Compatibility
 */

#include <stlsoft/stlsoft.h>    /* If the compiler cannot find this, you are not using STLSoft 1.9 or later, as required. */
#if !defined(_STLSOFT_VER) || \
   _STLSOFT_VER < 0x010964ff
# error Requires STLSoft 1.9.100, or later. (www.stlsoft.org/downloads.html)
#endif /* STLSoft version */

/* ////////////////////////////////////////////////////////////////////// */

#endif /* SHWILD_INCL_H_SHWILD_STLSOFT */

/* ///////////////////////////// end of file //////////////////////////// */
