/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

/**
 * Base class for all unit prefixes, such as kilo, micro, nano, etc.
 * A prefix class is specified (with a Unit) to form a PrefixedUnit object.
 * All standard prefixes are available as static fields of this class.  These go from
 * YOCTO (10^-24) to YOTTA (10^+24).
 *
 * @see PrefixedUnit
 */
public enum Prefix {
    YOCTO   (1.0e-24,   "yocto",  "y"),
    ZEPTO   (1.0e-21,   "zepto",  "z"),
    ATTO    (1.0e-18,   "atto",   "a"),
    FEMTO   (1.0e-15,   "femto",  "f"),
    PICO    (1.0e-12,   "pico",   "p"),
    NANO    (1.0e-09,   "nano",   "n"),
    MICRO   (1.0e-06,   "micro",  "\u00B5"), //unicode for micro sign
    MILLI   (1.0e-03,   "milli",  "m"),
    CENTI   (1.0e-02,   "centi",  "c"),
    DECI    (1.0e-01,   "deci",   "d"),
    NULL    (1.0e+00,   "",       ""),
    DEKA    (1.0e+01,   "deka",   "da"),
    HECTO   (1.0e+02,   "hecto",  "h"),
    KILO    (1.0e+03,   "kilo",   "k"),
    MEGA    (1.0e+06,   "mega",   "M"),
    GIGA    (1.0e+09,   "giga",   "G"),
    TERA    (1.0e+12,   "tera",   "T"),
    PETA    (1.0e+15,   "peta",   "P"),
    EXA     (1.0e+18,   "exa",    "E"),
    ZETTA   (1.0e+21,   "zetta",  "Z"),
    YOTTA   (1.0e+24,   "yotta",  "Y");

    /**
     * Constant used when invoking constructor of SimpleUnit class to indicate
     * that the unit permits the use of a prefix.
     */
    public static final boolean ALLOWED = true;
    /**
     * Constant used when invoking constructor of SimpleUnit class to indicate
     * that the unit does not permit the use of a prefix.
     */
    public static final boolean NOT_ALLOWED = false;

    private final double value;
    private final String name, symbol;

    Prefix(double value, String name, String symbol) {
        this.value = value;
        this.name = name;
        this.symbol = symbol;
    }

    /**
     * Returns a prefix corresponding to the given key.
     * For example, 'm' gives MILLI, 'M' give MEGA, 'f' gives FEMTO, etc.
     * Note that " " (space) returns the NULL prefix, and 'D' returns DEKA
     * (which has a two-letter symbol "da").
     * Returns null if the char does not correspond to any Prefix.
     */
    public static Prefix keySelect(char aKey) {
        if(aKey == ' ') return NULL;
        if(aKey == 'D') return DEKA;//DEKA.symbol() is "da", so use 'D' for it here
        String sKey = Character.toString(aKey);
        for(Prefix prefix : values()) {
            if(prefix.symbol().equals(sKey)) return prefix;
        }
        return null;
    }


    /**
     * @return the value associated with the prefix. For example, KILO.value() returns 1000.
     */
    public double value() {
        return value;
    }

    /**
     * @return the name of the prefix. For example, KILO.toString() returns "kilo".
     * All letters in returned string are lower case.
     */
    public String toString() {
        return name;
    }

    /**
     * Symbol used to show that prefix is part of the unit. For example the "k" in kg.
     * @return the symbol representing the prefix.
     */
    public String symbol() {
        return symbol;
    }

/*
    public static void main(String[] args) {
        System.out.println(keySelect('k'));
        System.out.println(keySelect('M'));
        System.out.println(keySelect(' '));
        System.out.println(keySelect('D'));
        System.out.println(keySelect('Y'));
        System.out.println(keySelect('y'));
    }
*/

}