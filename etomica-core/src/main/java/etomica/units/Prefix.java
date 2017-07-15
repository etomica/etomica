/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import java.io.ObjectStreamException;

/**
 * Base class for all unit prefixes, such as kilo, micro, nano, etc.
 * A prefix class is specified (with a Unit) to form a PrefixedUnit object.
 * All standard prefixes are available as static fields of this class.  These go from
 * YOCTO (10^-24) to YOTTA (10^+24).
 * 
 * @see PrefixedUnit
 */
public abstract class Prefix implements java.io.Serializable {
 
    /**
     * Private constructor to prevent instantiation externally.
     */
    private Prefix() {}
    
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
    
    /**
     * Prefix indicating 10<sup>-24</sup>.
     */
    public static final Prefix YOCTO = new Yocto();
    
    /**
     * Prefix indicating 10<sup>-21</sup>.
     */
    public static final Prefix ZEPTO = new Zepto();
    
    /**
     * Prefix indicating 10<sup>-18</sup>.
     */
    public static final Prefix ATTO = new Atto();
    
    /**
     * Prefix indicating 10<sup>-15</sup>.
     */
    public static final Prefix FEMTO = new Femto();
    
    /**
     * Prefix indicating 10<sup>-12</sup>.
     */
    public static final Prefix PICO = new Pico();
    
    /**
     * Prefix indicating 10<sup>-9</sup>.
     */
    public static final Prefix NANO = new Nano();
    
    /**
     * Prefix indicating 10<sup>-6</sup>.
     */
    public static final Prefix MICRO = new Micro();
    
    /**
     * Prefix indicating 10<sup>-3</sup> = 0.001.
     */
    public static final Prefix MILLI = new Milli();
    
    /**
     * Prefix indicating 10<sup>-2</sup> = 0.01.
     */
    public static final Prefix CENTI = new Centi();
    
    /**
     * Prefix indicating 10<sup>-1</sup> = 0.1.
     */
    public static final Prefix DECI = new Deci();
    
    /**
     * Prefix indicating 10<sup>0</sup> = 1.
     */
    public static final Prefix NULL = new Null();
    
    /**
     * Prefix indicating 10<sup>+1</sup> = 10.
     */
    public static final Prefix DEKA = new Deka();
    
    /**
     * Prefix indicating 10<sup>+2</sup> = 100.
     */
    public static final Prefix HECTO = new Hecto();
    
    /**
     * Prefix indicating 10<sup>+3</sup> = 1000.
     */
    public static final Prefix KILO = new Kilo();
    
    /**
     * Prefix indicating 10<sup>+6</sup>.
     */
    public static final Prefix MEGA = new Mega();
    
    /**
     * Prefix indicating 10<sup>+9</sup>.
     */
    public static final Prefix GIGA = new Giga();
    
    /**
     * Prefix indicating 10<sup>+12</sup>.
     */
    public static final Prefix TERA = new Tera();
    
    /**
     * Prefix indicating 10<sup>+15</sup>.
     */
    public static final Prefix PETA = new Peta();
    
    /**
     * Prefix indicating 10<sup>+18</sup>.
     */
    public static final Prefix EXA = new Exa();
    
    /**
     * Prefix indicating 10<sup>+21</sup>.
     */
    public static final Prefix ZETTA = new Zetta();
    
    /**
     * Prefix indicating 10<sup>+24</sup>.
     */
    public static final Prefix YOTTA = new Yotta();
    
    public static final Prefix[] ALL = new Prefix[] {
        YOCTO, ZEPTO, ATTO, FEMTO, PICO, NANO, MICRO, MILLI, CENTI, DECI, NULL,
        DEKA, HECTO, KILO, MEGA, GIGA, TERA, PETA, EXA, ZETTA, YOTTA};

    public abstract double value();
    public abstract String toString();
    public abstract String symbol();
    
    private final static class Yocto extends Prefix {
        public double value() {return 1.0e-24;}
        public String toString() {return "yocto";}
        public String symbol() {return "y";}
    }
    private final static class Zepto extends Prefix {
        public double value() {return 1.0e-21;}
        public String toString() {return "zepto";}
        public String symbol() {return "z";}
   }
    private final static class Atto extends Prefix {
        public double value() {return 1.0e-18;}
        public String toString() {return "atto";}
        public String symbol() {return "a";}
    }
    private final static class Femto extends Prefix {
        public double value() {return 1.0e-15;}
        public String toString() {return "femto";}
        public String symbol() {return "f";}
    }
    private final static class Pico extends Prefix {
        public double value() {return 1.0e-12;}
        public String toString() {return "pico";}
        public String symbol() {return "p";}
    }
    private final static class Nano extends Prefix {
        public double value() {return 1.0e-9;}
        public String toString() {return "nano";}
        public String symbol() {return "n";}
   }
    private final static class Micro extends Prefix {
        public double value() {return 1.0e-6;}
        public String toString() {return "micro";}
        public String symbol() {return "\u00B5";} //unicode for micro sign
   }
    private final static class Milli extends Prefix {
        public double value() {return 1.0e-3;}
        public String toString() {return "milli";}
        public String symbol() {return "m";}
    }
    private final static class Centi extends Prefix {
        public double value() {return 1.0e-2;}
        public String toString() {return "centi";}
        public String symbol() {return "c";}
    }
    private final static class Deci extends Prefix {
        public double value() {return 1.0e-1;}
        public String toString() {return "deci";}
        public String symbol() {return "d";}
    }
    private final static class Null extends Prefix {
        public double value() {return 1.0;}
        public String toString() {return "";}
        public String symbol() {return "";}
    }
    private final static class Deka extends Prefix {
        public double value() {return 1.0e+1;}
        public String toString() {return "deka";}
        public String symbol() {return "da";}
    }
    private final static class Hecto extends Prefix {
        public double value() {return 1.0e+2;}
        public String toString() {return "hecto";}
        public String symbol() {return "h";}
    }
    private final static class Kilo extends Prefix {
        public double value() {return 1.0e+3;}
        public String toString() {return "kilo";}
        public String symbol() {return "k";}
    }
    private final static class Mega extends Prefix {
        public double value() {return 1.0e+6;}
        public String toString() {return "mega";}
        public String symbol() {return "M";}
    }
    private final static class Giga extends Prefix {
        public double value() {return 1.0e+9;}
        public String toString() {return "giga";}
        public String symbol() {return "G";}
    }
    private final static class Tera extends Prefix {
        public double value() {return 1.0e+12;}
        public String toString() {return "tera";}
        public String symbol() {return "T";}
    }
    private final static class Peta extends Prefix {
        public double value() {return 1.0e+15;}
        public String toString() {return "peta";}
        public String symbol() {return "P";}
    }
    private final static class Exa extends Prefix {
        public double value() {return 1.0e+18;}
        public String toString() {return "exa";}
        public String symbol() {return "E";}
    }
    private final static class Zetta extends Prefix {
        public double value() {return 1.0e+21;}
        public String toString() {return "zetta";}
        public String symbol() {return "Z";}
    }
    private final static class Yotta extends Prefix {
        public double value() {return 1.0e+24;}
        public String toString() {return "yotta";}
        public String symbol() {return "Y";}
    }

    /**
     * Returns a prefix corresponding to the given key.
     * For example, 'm' gives MILLI, 'M' give MEGA, 'f' gives FEMTO, etc.
     * Returns null if the char does not correspond to any prefix.
     */
    public static Prefix keySelect(char aKey) {
        int i = intKeySelect(aKey);
        return i!=-1 ? ALL[i] : null;
    }
    /**
     * Returns an integer code corresponding to the prefix for the given key.
     * The integers correspond to the index of the prefix in the ALL array.
     * Returns -1 if the char does not correspond to any prefix.
     */
	public static int intKeySelect(char aKey) {
	    switch(aKey) {
	        case 'y': return 0;
	        case 'z': return 1;
	        case 'a': return 2;
	        case 'f': return 3;
	        case 'p': return 4;
	        case 'n': return 5;
	        case 'u': return 6;
	        case 'm': return 7;
	        case 'c': return 8;
	        case 'd': return 9;
	        case ' ': return 10;
	        case 'D': return 11;
	        case 'h': return 12;
	        case 'H': return 12;
	        case 'k': return 13;
	        case 'K': return 13;
	        case 'M': return 14;
	        case 'G': return 15;
	        case 'T': return 16;
	        case 'P': return 17;
	        case 'E': return 18;
	        case 'Z': return 19;
	        case 'Y': return 20;
	        default: return -1;
	    }
	}

}