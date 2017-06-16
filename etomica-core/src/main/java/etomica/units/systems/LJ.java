/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.systems;

import etomica.units.Count;
import etomica.units.Decimal;
import etomica.units.dimensions.Dimension;
import etomica.units.Electron;
import etomica.units.Radian;
import etomica.units.Unit;
import etomica.util.Strings;

/**
 * Lennard-Jones system of units, such that all quantities are made
 * dimensionless with respect to a characteristic size (sigma), energy (epsilon)
 * and mass. Values of sigma, epsilon and mass reside with the instance of
 * the LJ unit system, to which each unit refers when making its conversion.
 * Thus any changes in sigma, epsilon, or mass are immediately recognized
 * by the unit conversions.
 * 
 * Values of sigma, epsilon, and mass are specified in simulation units.
 * Default values are each 1.0.
 */

// consider whether LJ units for charge, dipole, etc should be reduced by 4 pi eps0 rather than by electron charge 

public class LJ extends UnitSystem  {

    /**
     * Constructs LJ unit system with default values of sigma, epsilon and mass each equal to unity.
     * Default uses unicode for symbols
     */
    public LJ() {
        this(1.0, 1.0, 1.0, true);
    }

    /**
     * Constructs LJ unit system with given values of sigma, epsilon, mass, given in simulation units.
     * Unit classes can be used to convert to simulation units if needed. For example, if epsilon is
     * given as epsilon/kB in units of Kelvin, the argument would take Kelvin.UNIT.toSim(epsilon).
     *
     * @param sigma LJ length parameter
     * @param epsilon LJ energy parameter
     * @param mass characteristic mass, nominally the mass of a LJ atom
     * @param unicodeSymbols
     *            if true, Unit symbols will be constructed using unicode
     *            expressions for sigma and epsilon and other symbols; otherwise
     *            substitutes are used.
     */
    public LJ(double sigma, double epsilon, double mass, boolean unicodeSymbols) {
        this.sigma = sigma;
        this.epsilon = epsilon;
        this.mass = mass;
        if(unicodeSymbols) {
            EPS = "\u03B5";
            SIG = "\u03C3";
            HALF = "\u00BD";
            DOT = "\u22C5";
        } else {
            EPS = "epsilon";
            SIG = "sigma";
            HALF = "1/2";
            DOT = "-";
        }
    }

    //accessor and mutator methods
    public double getSigma() {
        return sigma;
    }

    public void setSigma(double s) {
        sigma = s;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public void setEpsilon(double e) {
        epsilon = e;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double m) {
        mass = m;
    }

    //implementation of UnitSystem interface
    public Unit quantity() {
        return Count.UNIT;
    }

    public Unit fraction() {
        return Decimal.UNIT;
    }

    public Unit mass() {
        return massUnit;
    }

    public Unit length() {
        return lengthUnit;
    }

    public Unit time() {
        return timeUnit;
    }

    public Unit angle() {
        return Radian.UNIT;
    }

    public Unit charge() {
        return Electron.UNIT;
    }
    
    public Unit current() {
        return currentUnit;
    }

    public Unit dipole() {
        return dipoleUnit;
    }
    
    public Unit force() {
        return forceUnit;
    }

    public Unit energy() { return energyUnit; }

    public Unit power() { return powerUnit; }

    public Unit temperature() {
        return temperatureUnit;
    }

    public Unit pressure() {
        return pressureUnit;
    }

    public Unit volume() {
        return volumeUnit;
    }

    public Unit area() {
        return areaUnit;
    }
    
    public Unit viscosity() {
        return viscosityUnit;
    }
    
    protected double sigma = 1.0;
    protected double epsilon = 1.0;
    protected double mass = 1.0;
    private static final long serialVersionUID = 1;
    
    // \u03B5 is unicode for epsilon
    // \u03C3 is unicode for sigma
    // \u00BD is unicode for "1/2"
    private final String EPS;
    private final String SIG;
    private final String HALF;
    private final String DOT;

    private final Unit massUnit = new Mass(this);
    private final Unit lengthUnit = new Length(this);
    private final Unit timeUnit = new Time(this);
    private final Unit dipoleUnit = new Dipole(this);
    private final Unit forceUnit = new Force(this);
    private final Unit energyUnit = new Energy(this);
    private final Unit powerUnit = new Power(this);
    private final Unit temperatureUnit = new Temperature(this);
    private final Unit pressureUnit = new Pressure(this);
    private final Unit volumeUnit = new Volume(this);
    private final Unit areaUnit = new Area(this);
    private final Unit viscosityUnit = new Viscosity(this);
    private final Unit currentUnit = new Current(this);

    //definitions of the LJ units
    

    //parent class for the LJ units
    public static abstract class LJUnit implements Unit {

        /**
         * Unit must provide instance of LJ unit system at construction.
         * Conversions should be defined for each unit in reference to 
         * lj.sigma, lj.epsilon, and lj.mass.
         */
        protected LJUnit(LJ ljSystem) {
            lj = ljSystem;
        }

        public abstract Dimension dimension();

        public abstract String symbol();

        public abstract double fromSim(double x);

        public double toSim(double x) {
            if (x == 0.0) {
                return 0.0;
            }
            return 1.0 / fromSim(1.0 / x);
        }

        public boolean prefixAllowed() {
            return false;
        }

        protected final LJ lj;
    }

    private static final class Mass extends LJUnit {
        private Mass(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Mass.DIMENSION;
        }

        public String symbol() {
            return "m";
        }

        public double fromSim(double x) {
            return x / lj.mass;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Length extends LJUnit {
        private Length(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Length.DIMENSION;
        }

        public String symbol() {
            return lj.SIG;
        }

        public double fromSim(double x) {
            return x / lj.sigma;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Time extends LJUnit {
        private Time(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Time.DIMENSION;
        }

        public String symbol() {
            return "("+lj.EPS+Strings.exponent(-1)+lj.DOT+"m"+lj.DOT+lj.SIG+Strings.exponent(2)+")"+Strings.exponent(lj.HALF);
        }

        public double fromSim(double x) {
            return x * Math.sqrt(lj.epsilon / lj.mass) / lj.sigma;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Current extends LJUnit {
        private Current(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Current.DIMENSION;
        }

        public String symbol() {
            return "e"+lj.DOT+"("+lj.EPS+lj.DOT+"m"+lj.DOT+lj.SIG+Strings.exponent(2)+")"+Strings.exponent("-"+lj.HALF);
        }

        public double fromSim(double x) {
            return x * Electron.UNIT.fromSim(1.0)/(Math.sqrt(lj.epsilon / lj.mass) / lj.sigma);
        }

        private static final long serialVersionUID = 1;
    }


    private static final class Dipole extends LJUnit {
        private Dipole(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Dipole.DIMENSION;
        }

        public String symbol() {
            return "e"+lj.DOT+lj.SIG;
        }

        public double fromSim(double x) {
            return Electron.UNIT.fromSim(x) / lj.sigma;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Force extends LJUnit {
        private Force(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Force.DIMENSION;
        }

        public String symbol() { // \u03B5/\u03C3";// epsilon/sigma
            return lj.EPS+"/"+lj.SIG;
        }

        public double fromSim(double x) {
            return x * lj.sigma / lj.epsilon;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Energy extends LJUnit {
        private Energy(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Energy.DIMENSION;
        }

        public String symbol() {
            return lj.EPS;//"\u03B5";
        }

        public double fromSim(double x) {
            return x / lj.epsilon;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Power extends LJUnit {
        private Power(LJ ljSystem) { super(ljSystem); }

        public Dimension dimension() {
            return etomica.units.dimensions.Power.DIMENSION;
        }

        public String symbol() {
            return  "("+lj.EPS+Strings.exponent(1.5)+lj.DOT+"m"+Strings.exponent(-0.5)+lj.DOT+lj.SIG+Strings.exponent(-1)+")";
        }

        public double fromSim(double x) {
            return x/(lj.epsilon * Math.sqrt(lj.epsilon / lj.mass) / lj.sigma );
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Temperature extends LJUnit {
        private Temperature(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Temperature.DIMENSION;
        }

        public String symbol() {
            return lj.EPS+"/k";//"\u03B5/k";
        }

        public double fromSim(double x) {
            return x / lj.epsilon;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Pressure extends LJUnit {
        private Pressure(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Pressure.DIMENSION;
        }

        public String symbol() {
            return lj.EPS+"/"+lj.SIG+Strings.exponent(3);//"\u03B5/\u03C3^3";
        }

        public double fromSim(double x) {
            return x * lj.sigma * lj.sigma * lj.sigma / lj.epsilon;
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Volume extends LJUnit {
        private Volume(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Volume.DIMENSION;
        }

        public String symbol() {
            return lj.SIG+Strings.exponent(3);//"\u03C3^3";
        }

        public double fromSim(double x) {
            return x / (lj.sigma * lj.sigma * lj.sigma);
        }

        private static final long serialVersionUID = 1;
    }

    private static final class Area extends LJUnit {
        private Area(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Area.DIMENSION;
        }

        public String symbol() {
            return lj.SIG+Strings.exponent(2);
        }

        public double fromSim(double x) {
            return x / (lj.sigma * lj.sigma);
        }

        private static final long serialVersionUID = 1;
    }
    
    private static final class Viscosity extends LJUnit {
        private Viscosity(LJ ljSystem) {
            super(ljSystem);
        }

        public Dimension dimension() {
            return etomica.units.dimensions.Viscosity.DIMENSION;
        }

        public String symbol() {
            return lj.SIG+Strings.exponent(2)+"/(m"+lj.DOT+lj.EPS+")"+Strings.exponent(lj.HALF);//\u03C3^2/(m\u03B5)^\u00BD";
        }

        public double fromSim(double x) {
            return x * lj.sigma * lj.sigma / Math.sqrt(lj.epsilon * lj.mass);
        }

        private static final long serialVersionUID = 1;
    }


}
