/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.random;

/*
    A C-program for MT19937, with initialization improved 2002/1/26.
    Coded by Takuji Nishimura and Makoto Matsumoto.
    
    Before using, initialize the state by using init_genrand(seed)  
    or init_by_array(init_key, key_length).
    
    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
    All rights reserved.                          
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    
      1. Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
    
      2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
         documentation and/or other materials provided with the distribution.
    
      3. The names of its contributors may not be used to endorse or promote 
         products derived from this software without specific prior written 
         permission.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Any feedback is very welcome.
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
    email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

import etomica.meta.annotations.IgnoreProperty;

/**
 * Mersenne Twister RNG.  seed initialization and nextInt() methods written in
 * C by Takuji Nishimura and Makoto Matsumoto as described above, and
 * translated into Java.
 * <p>
 * Other methods added for the Java version.
 *
 * @author Takuji Nishimura
 * @author Makoto Matsumoto
 * @author Andrew Schultz
 */
public class RandomMersenneTwister implements IRandom {

    /* Period parameters */
    static final protected int N = 624;
    static final protected int M = 397;
    static final protected int MATRIX_A = 0x9908b0df;   /* constant vector a */
    static final protected int UPPER_MASK = 0x80000000; /* most significant w-r bits */
    static final protected int LOWER_MASK = 0x7fffffff; /* least significant r bits */
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    static final int[] mag01 = new int[]{0x0, MATRIX_A};
    protected final int[] mt = new int[N]; /* the array for the state vector  */
    protected int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */

    protected boolean hasNextGaussian = false;
    protected double nextGaussian;
    protected int[] savedSeedArray;
    protected int savedSeed;

    /**
     * Creates a Mersenne Twister with the given seed.
     */
    public RandomMersenneTwister(int s) {
        setSeed(s);
    }

    /**
     * Creates a Mersenne Twister with the given array of seeds.
     */
    public RandomMersenneTwister(int[] s) {
        this(5);
        setSeedArray(s);
    }

    @IgnoreProperty
    public int getSeed() {
        if (savedSeedArray != null) throw new RuntimeException("it's a seed array");
        return savedSeed;
    }

    /**
     * Configures the RNG with the given seed.
     */
    public void setSeed(int s) {
        savedSeed = s;
        savedSeedArray = null;
        /* initializes mt[N] with a seed */
        mt[0] = s;
        for (mti = 1; mti < N; mti++) {
            mt[mti] =
                    (1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >>> 30)) + mti);
            /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
            /* In the previous versions, MSBs of the seed affect   */
            /* only MSBs of the array mt[].                        */
            /* 2002/01/09 modified by Makoto Matsumoto             */
        }
        hasNextGaussian = false;
    }

    public int[] getSeedArray() {
        return savedSeedArray;
    }

    /**
     * Configures the RNG with an array of seeds
     */
    public void setSeedArray(int init_key[]) {
        /* initialize by an array with array-length */
        /* init_key is the array for initializing keys */
        /* slight change for C++, 2004/2/26 */
        int key_length = init_key.length;
        int i, j, k;
        setSeed(19650218);
        savedSeedArray = init_key.clone();
        i = 1;
        j = 0;
        k = (N > key_length ? N : key_length);
        for (; k > 0; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >>> 30)) * 1664525))
                    + init_key[j] + j; /* non linear */
            i++;
            j++;
            if (i >= N) {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if (j >= key_length) j = 0;
        }
        for (k = N - 1; k > 0; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >>> 30)) * 1566083941))
                    - i; /* non linear */
            i++;
            if (i >= N) {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }

        mt[0] = 0x80000000; /* MSB is 1; assuring non-zero initial array */
    }

    /* generates a random number on [0,0xffffffff]-interval */
    public int nextInt() {
        int y;

        if (mti >= N) { /* generate N words at one time */
            int kk;

            for (kk = 0; kk < N - M; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >>> 1) ^ mag01[y & 0x1];
            }
            for (; kk < N - 1; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + (M - N)] ^ (y >>> 1) ^ mag01[y & 0x1];
            }
            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >>> 1) ^ mag01[y & 0x1];

            mti = 0;
        }

        y = mt[mti++];

        /* Tempering */
        y ^= (y >>> 11);
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= (y >>> 18);

        return y;
    }

    public long nextLong() {
        return (((long) nextInt()) << 32) | (nextInt() & ((1L << 32) - 1));
    }

    /**
     * Returns a random int ranging from 0 to max-1.
     */
    public int nextInt(int max) {
        if (max < 1) {
            throw new RuntimeException("max must be positive");
        }

        // we take a random int s and return s%max.  But our s must be less
        // than the largest integer multiple of max
        // maxRand+1 is the largest integer value that is a multiple of max
        int maxRand = Integer.MAX_VALUE - (int) ((Integer.MAX_VALUE + 1L) % max);
        int s;
        do {
            s = nextInt() & LOWER_MASK;
        } while (s > maxRand);
        return s % max;
    }

    public double nextFixedDouble() {
        // yes, yes, only 53 bits matter.  who cares.
        // if you need that, call nextDouble
        return (nextLong() >>> 1) / ((double) (-1L >>> 1));
    }

    public double nextDouble() {
        // nominally, we'd generate a 53-bit long (x) and then divide by 2^53
        // but that's sub-optimal for small values of x -- we can generate
        // 0, 1/2^53, 2/2^53, 3/2^53, etc.
        // but there are many other numbers within that range that we would not
        // be able to return.  Instead we effectively generate random bits
        // until we have leading 0's followed by a 1 and then the 52 random
        // bits came after the 1.  We then divide that by 2^a where a is the
        // total number of bits we generated (leading zero's + 53).

        double shiftFac = 9007199254740992.0; //9007199254740992 = 1<<53
        // we want y to be an int with 10 leading zeros (so that it has 22 bits)
        // then we'll grab the lowest 31 bits from z  (22+31 = 53)
        int y = 0;
        while (y == 0) {
            y = nextInt();
            if (y == 0) {
                shiftFac *= (1L << 32);
                continue;
            }
            int p = Integer.numberOfLeadingZeros(y);
            // shift the leading 1 to the 10th bit
            if (p < 11) {
                // dump unneeded low bits
                y = y >>> (10 - p);
                shiftFac *= 1 << p;
            } else {
                y = y << (p - 10);
                shiftFac *= 1L << p;
                // we left a hole in the low bits.  fill that in.
                y |= nextInt() >>> (42 - p);
            }
        }

        int z = nextInt();

        return ((((long) y) << 31) | (z & LOWER_MASK)) / shiftFac;
    }

    public double nextGaussian() {
        if (hasNextGaussian) {
            hasNextGaussian = false;
            return nextGaussian;
        }

        double x1, x2, w;
        do {
            x1 = nextDouble();
            x2 = nextDouble();
            w = x1 * x1 + x2 * x2;
        } while (w >= 1);
        w = Math.sqrt(-2 * Math.log(w) / w);
        int signs = nextInt();
        if ((signs & 0x00000001) == 1) x1 = -x1;
        if ((signs & 0x00000002) == 2) x2 = -x2;
        nextGaussian = x2 * w;
        hasNextGaussian = true;
        return x1 * w;
    }

    public float nextFloat() {
        // see nextDouble

        float shiftFac = 16777216.0f;  // 16777216 = 1<<24
        // we want y to be an int with 8 leading zeros (so that it has 24 bits)
        int y = 0;
        while (y == 0) {
            y = nextInt();
            if (y == 0) {
                shiftFac *= (1L << 32);
                continue;
            }
            int p = Integer.numberOfLeadingZeros(y);
            // shift the leading 1 to the 8th bit
            if (p < 9) {
                // dump unneeded low bits
                y = y >>> (8 - p);
                shiftFac *= 1 << p;
            } else {
                y = y << (p - 8);
                shiftFac *= 1L << p;
                // we left a hole in the low bits.  fill that in.
                y |= nextInt() >>> (40 - p);
            }
        }

        return y / shiftFac;
    }
}
