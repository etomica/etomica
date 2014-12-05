/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos.bigdecimal ;

import java.lang.* ;
import java.util.* ;
import java.math.* ;

/** Harmonic numbers.
* H(n) is the sum of the inverses of the integers from 1 to n.
* @since 2008-10-19
* @author Richard J. Mathar
*/
public class Harmonic
{
        /** ctor()
        * Does nothing.
        */
        public Harmonic()
        {
        }

        /** The Harmonic number at the index specified
        * @param n the index, non-negative.
        * @return the H_1=1 for n=1, H_2=3/2 for n=2 etc.
        *   For values of n less than 1, zero is returned.
        */
        public Rational at(int n)
        {
                if ( n < 1)
                        return(new Rational(0,1)) ;
                else
                {
                        /* start with 1 as the result
                        */
                        Rational a = new Rational(1,1) ;

                        /* add 1/i for i=2..n
                        */
                        for( int i=2 ; i <=n ; i++)
                                a = a.add(new Rational(1,i)) ;
                        return a ;
                }
        }
} /* Harmonic */
