
    /***************************************************************************/
    /* THIS IS PROPRIETARY SOURCE CODE OF THE UNIVERSITY OF WALES, SWANSEA     */
    /*                      All Rights Reserved                                */
    /*     Duplication or use of any form of this file (source, object,        */
    /*     or executable) is prohibited without the express written            */
    /*     consent of Prof. N. P. Weatherill                                   */
    /*     Modifications to this software or use without the stated permission */
    /*     is an infringement of copyright and appropriate legal action        */
    /*     may be taken.                                                       */
    /* THIS MODULE IS MADE AVAILABLE TO THE CAESAR PROJECT,FUNDED BY THE CEC   */
    /*        TEL. (+U.K.) 1792-295598                                         */
    /*                                                                         */
    /*        Date of initial work : 1st January 1995                          */
    /*        Reference : CAESAR/H/UWS/y1/0.6/1.7.96                           */
    /*        Module Name : gb_localmath                                       */
    /*        Description : header -- gb_localmath for GB                      */
    /***************************************************************************/




/*
 *      Local Math Header
 */

#define SQ(f)  ((f)*(f))
#define MAX(a,b)  ( (a) < (b) ? (b) : (a) )
#define MIN(a,b)  ( (a) > (b) ? (b) : (a) )
#define ABS(a) (((a) < 0) ? -(a):(a))


#define ONE(a)		(MIN(MAX((a), 0.0), 1.0))

#define MAX3(a, b, c) \
((a) > (b) ? \
 ((a) > (c) ? (a) : (c)) : \
 ((b) > (c) ? (b) : (c)))

#define MIN3(a, b, c) \
((a) < (b) ? \
 ((a) < (c) ? (a) : (c)) : \
 ((b) < (c) ? (b) : (c)))
