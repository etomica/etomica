/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.random;

import java.io.FileInputStream;
import java.io.IOException;

/**
 * Random number generator that takes random bits from files.  This can be used
 * with the unix special files /dev/random or /dev/urandom.
 * <p>
 * Due to file access, this is likely to be slower than other RNG, but could be
 * used to provide an initial seed.
 *
 * @author Andrew Schultz
 */
public class RandomNumberGeneratorUnix implements IRandom {

    protected final FileInputStream fReader;
    protected final byte[] buf = new byte[8];
    // temporary storage for nextGaussian, which generates 2 Gaussians at a time
    protected boolean hasNextGaussian = false;
    protected double nextGaussian = 0;

    /**
     * Create an RNG using /dev/urandom as a source of random bits.
     */
    public RandomNumberGeneratorUnix() {
        this("/dev/urandom");
    }

    /**
     * Create an RNG using the given file as a source of random bits.
     * <p>
     * /dev/urandom as a randFile will not block
     * /dev/random as a randFile may block if the entropy pool is depleted
     */
    public RandomNumberGeneratorUnix(String randFile) {
        try {
            fReader = new FileInputStream(randFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Convenience method that returns an array of integers, useful as seeds to
     * another RNG.  If an instance of this class can be instantiated, it is
     * used to return 4 integers.  Otherwise, System.nanotime is split into 2
     * integers.
     */
    public static int[] getRandSeedArray() {
        int[] rv;
        RandomNumberGeneratorUnix rng = null;
        try {
            rng = new RandomNumberGeneratorUnix();
        } catch (RuntimeException e) {
        }
        if (rng == null) {
            // no real problem, probably just not on a unix box.
            rv = new int[2];
            long t = System.nanoTime();
            rv[0] = (int) t; // lower 32 bits
            rv[1] = (int) (t >> 32); // upper 32 bits
        } else {
            rv = new int[4];
            for (int i = 0; i < rv.length; i++) {
                rv[i] = rng.nextInt();
            }
            rng.dispose();
        }
        return rv;
    }

    public static void main(String[] args) {
        long l = (1L << 51) + 3;
        System.out.println(l / 2.0 + " " + (l + 1) / 2.0 + " " + ((l / 2.0) - (l + 1) / 2.0));
        l = l << 1;
        System.out.println(l / 2.0 + " " + (l + 1) / 2.0 + " " + ((l / 2.0) - (l + 1) / 2.0));
        l = l << 1;
        System.out.println(l / 2.0 + " " + (l + 1) / 2.0 + " " + ((l / 2.0) - (l + 1) / 2.0));
        RandomNumberGeneratorUnix random = new RandomNumberGeneratorUnix();
        float min = 1;
        float max = 0;
        float t = 0;
        // min should be ~1e-7, never 0
        // max should be ~0.9999999, never 1
        // avg should be ~0.5
        for (int i = 0; i < 10000000; i++) {
            float f = random.nextFloat();
            if (f < min) min = f;
            if (f > max) max = f;
            t += f;
            if ((i + 1) % 1000000 == 0) System.out.println(min + " " + max + " " + t + " " + t / (i + 1));
        }

        System.out.println(((1L << 50) - 1) + " " + random.nextLong((1L << 50) - 1));
        int[] bins = new int[3];
        for (int i = 0; i < 10000000; i++) {
            bins[random.nextInt(bins.length)]++;
        }

        // bins should have (more or less) equal hits
        for (int i = 0; i < bins.length; i++) {
            System.out.print(bins[i] + " ");
        }
    }

    /**
     * Cleanup, close the file.
     */
    public void dispose() {
        try {
            fReader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Fills the byte buffer with random bytes up to numBytes, starting at
     * offset.
     */
    protected void fill(int offset, int numBytes) {
        try {
            int n = offset;
            while (n < numBytes) {
                n += fReader.read(buf, n, numBytes - n);
            }
            if (n != numBytes) {
                throw new RuntimeException("Could not read from file (" + n + "/" + numBytes + ")");
            }
        } catch (IOException e) {
            throw new RuntimeException("problem reading file");
        }
    }

    /**
     * Returns a random long (64-bit integer).
     */
    public long nextLong() {
        fill(0, 8);
        long s = 0;
        for (int i = 0; i < 8; i++) {
            s += (buf[i] + 128L) << (i * 8);
        }
        return s;
    }

    /**
     * Returns a random int (32-bit integer).
     */
    public int nextInt() {
        fill(0, 4);
        int s = 0;
        for (int i = 0; i < 4; i++) {
            s += (buf[i] + 128) << (i * 8);
        }
        return s;
    }

    /**
     * Returns a random int ranging from 0 to max-1.
     */
    public int nextInt(int max) {
        if (max < 1) {
            throw new RuntimeException("max must be positive");
        } else if (max == 1) return 0;
        // number of random bits we need
        int maxBit = 32 - Integer.numberOfLeadingZeros(max - 1);
        // number of random bytes we need
        int maxByte = (maxBit + 7) / 8;
        // number of extra bits in the last byte
        int extra = maxByte * 8 - maxBit;
        while (true) {
            fill(0, maxByte);
            int s = ((buf[0] + 128) >> extra) << ((maxByte - 1) * 8);
            for (int i = 1; i < maxByte; i++) {
                s += (buf[i] + 128) << ((maxByte - i - 1) * 8);
            }
            if (s < max) {
                return s;
            }
        }
    }

    /**
     * Returns a random long ranging from 0 to max-1.
     */
    public long nextLong(long max) {
        if (max < 1) {
            throw new RuntimeException("max must be positive");
        } else if (max == 1) return 0;
        // number of random bits we need
        int maxBit = 64 - Long.numberOfLeadingZeros(max - 1);
        // number of random bytes we need
        int maxByte = (maxBit + 7) / 8;
        // number of extra bits in the last byte
        int extra = maxByte * 8 - maxBit;
        while (true) {
            fill(0, maxByte);
            long s = ((buf[0] + 128L) >> extra) << ((maxByte - 1) * 8);
            for (int i = 1; i < maxByte; i++) {
                s += (buf[i] + 128L) << ((maxByte - i - 1) * 8);
            }
            if (s < max) {
                return s;
            }
        }
    }

    /**
     * Returns a random float, uniformly distributed between 0 and 1.
     * Neither 0 nor 1 will ever be returned.
     */
    public float nextFloat() {
        fill(0, 3);
        float byteOffsetFac = 1;
        while (buf[0] == -128) {
            // we have 0s for the highest 8 bits.
            // Generate another 8 low bits
            // we'll divide our final result by 256 to compensate
            byteOffsetFac *= 256;
            buf[0] = buf[1];
            buf[1] = buf[2];
            fill(2, 3);
        }
        int s = 0;
        for (int i = 0; i < 3; i++) {
            s += (buf[2 - i] + 128) << (i * 8);
        }
        return s / ((1 << 24) * byteOffsetFac);
    }

    /**
     * Returns a random double, uniformly distributed between 0 and 1.
     * Neither 0 nor 1 will ever be returned.
     */
    public double nextDouble() {
        fill(0, 7);
        double byteOffsetFac = 1;
        while ((buf[0] + 128) >> 3 == 0) {
            // we have 0s for the highest 5 bits
            // Generate (at least) another 8 low bits
            int shift = 1;
            for (; shift < 8; shift++) {
                if ((buf[shift] + 128) >> 3 != 0) {
                    // we found a non-zero byte
                    break;
                }
                byteOffsetFac *= 256;
            }
            // shift the non-zero bits up
            for (int i = 0; i < 7 - shift; i++) {
                buf[i] = buf[i + shift];
            }
            // grab some new low bits
            fill(7 - shift, 7);
        }
        long s = ((buf[0] + 128L) >> 3) << 48;
        for (int i = 0; i < 6; i++) {
            s += (buf[6 - i] + 128L) << (i * 8);
        }
        return s / ((1L << 53) * byteOffsetFac);
    }

    public double nextFixedDouble() {
        // yes, yes, only 53 bits matter.  who cares.
        // if you need that, call nextDouble
        return (nextLong() >>> 1) / ((double) (-1L >>> 1));
    }

    /**
     * Returns a double taken from a Gaussian distribution with 0 mean and
     * variance of 1.
     */
    public double nextGaussian() {
        if (hasNextGaussian) {
            hasNextGaussian = false;
            return nextGaussian;
        }
        double r1 = nextDouble();
        double r2 = nextDouble();
        double arg = 2.0 * Math.PI * r2;
        double sqrt = Math.sqrt(-2.0 * Math.log(r1));
        nextGaussian = sqrt * Math.cos(arg);
        return sqrt * Math.sin(arg);
    }
}
