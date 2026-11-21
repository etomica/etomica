/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import etomica.units.Joule;
import etomica.units.systems.LJ;

/**
 * Collection of assorted physical constants.  All values
 * are in simulations units in which time is in picoseconds, length is in Angstroms,
 * and mass is in Daltons (amu).  Also defines several enumerated-type constants.
 * 
 * The values for most of these constants can be found at
 *   http://physics.nist.gov/cuu/Constants/index.html
 */
public final class Constants {

    public static final double TWO_PI = 2.0*Math.PI;
    /**
     * Avogadro's number, 6.022140857(74)e23 molecules/mole.
     */
    // reference:  http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=abbr_in!
    public static final double AVOGADRO = 6.02214076e23;
    
    /*             
       Units for internal calculations (simulation units)
            time: ps
          length: Angstroms
            mass: Daltons, or atomic mass units (amu)
    */
    /**
     * The standard acceleration of gravity on Earth.
     */
    //  acceleration of gravity (on Earth), in A/ps^2
    public static final double G_EARTH = 9.8 * 1e10 / 1e24;
    /**
     * Boltzmann's constant, in (sim units)/kelvin.  Specifically,
     * (daltons)(angstroms^2)(ps^-2)(kelvin^-1).
     */
    //  Boltzmann's constant, converted from J/K to D-A^2/ps^2/K (where it equals ~0.8314)
    // (~1.38e-23 kg-m^2/s^2/K/molecule)(1000 g/kg)(N_avo D/g)(10^10 A/m)^2 (10^-12 s/ps)^2
    // reference:  http://physics.nist.gov/cgi-bin/cuu/Value?k
    public static final double BOLTZMANN_K = 1.380649e-23 * 1000 * AVOGADRO * 1e20 * 1e-24;
    /**
     * Planck's constant, in simulation units.  Specifically, equal to approximately 39.903127 D-A^2/ps.
     */
    //convert from J-s to simulation units
    public static final double PLANCK_H = Joule.UNIT.toSim(6.62607015e-34) * 1e12;
    /**
     * Rydberg's constant, in Angstroms.
     */
    public static final double RYDBERG_R = 1.0973731568160E-3;
    /**
     * The permittivity of a vacuum, in (sim units) (electron charge)^2. Specifically,
     * (ps^2)(daltons^-1)(angstroms^-3)(electrons^2).
     */
    //epsilon0, converted from C^2/(N-m^2) to e^2 ps^2/(D-A^3)
    // (~8.854e-12 C^2 s^2/(kg-m^3)) (1/1.60217653e-12 e/C)^2 (10^12 ps/s)^2 (10^-3 kg/g) (1/Avo g/D) (10^-10 m/A)^3\
    // https://pml.nist.gov/cgi-bin/cuu/Value?ep0
    public static final double EPSILON_0 = 8.8541878188e-12 * Math.pow(1.602176634e-19, -2)
            * 1e24 * 1e-3 / AVOGADRO * 1e-30;
    /**
     * The speed of light, in simulation units.  Equal to 2997924.58 Angstroms/picosecond
     */
    public static final double LIGHT_SPEED = 299792458 * (1e10 * 1e-12);//convert m/s to A/ps
    /**
     * The gravitational constant, in simulation units.
     * Equal to approximately 1.1e-31 A^3/ps^2/D.
     */
    // ~6.674e-11 m^3 s^-2 kg^-1 (1e10 A/m)^3 (10-12 s/ps)^2 (1kg/1000g) (1g/Avo D)
    // reference http://physics.nist.gov/cgi-bin/cuu/Value?bg
    public static final double G = 6.67430e-11 * (1e30 * 1e-24 * 1e-3) / AVOGADRO;
    
    //private constructor to prevent instantiation
    private Constants() {
    }
    
    public static void main(String arg[]) {
        System.out.println("Avogadro's number: "+AVOGADRO);
        System.out.println("Boltzmann's constant: "+BOLTZMANN_K);
        System.out.println("C: "+LIGHT_SPEED);
        System.out.println("Planck's constant: "+PLANCK_H);
        System.out.println("Epsilon0: "+EPSILON_0);
        System.out.println("G: "+G);
        System.out.println("1.0/sqrt(4 Pi Epsilon0): "+1.0/Math.sqrt(4.*Math.PI*EPSILON_0));
        System.out.println("unit toSim: "+new etomica.units.systems.MKS().viscosity().toSim(1.0));
        System.out.println("unit toSim: "+etomica.units.Volt.UNIT.toSim(1.0));
        System.out.println("symbol: "+new LJ(1,1,1,false).viscosity().symbol());
    }

    /**
     * Enumerated type for the compass directions NORTH, SOUTH, EAST, WEST.
     * Used to express the orientation of an object.
     */
   public enum CompassDirection {
        NORTH("North"), SOUTH("South"), WEST("West"), EAST("East");

        private final String label;

        CompassDirection(String label) {
            this.label = label;
        }

        @Override
        public String toString() {
            return label;
        }
    }

}
    
    
