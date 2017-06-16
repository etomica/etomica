/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.statmech;
import etomica.math.function.Function;
import etomica.space.Space;
import etomica.units.Mass;

/**
 * Non-instantiable class with static methods for implementing features of the 
 * Maxwell-Boltzmann velocity distribution.
 */
 
public final class MaxwellBoltzmann {
    
    //private to prevent instantiation
    private MaxwellBoltzmann() {}
    
    /**
     * Function giving the Maxwell-Boltzmann distribution of speeds.
     */
    public static class Distribution implements Function, java.io.Serializable {
        private double temperature;
        private double mass;
        private Space space;
        private double a, c;

        public Distribution(Space space, double temperature, double mass) {
            this.temperature = temperature;
            this.mass = mass;
            this.space = space;
            update();
        }
        public double f(double v) {
            if(v == 0.0) return (space.D() == 1) ? c : 0.0;
            double v2 = v*v;
            return c*space.powerD(v)*Math.exp(-a*v2)/v;
        }
        //not well checked.
        public double dfdx(double v) {
            throw new RuntimeException("MaxwellBoltzmann.Distribution.dfdx not implemented");
          //  double v2 = v*v;
          //  return c*Math.exp(-a*v2)*2*v*(1 - a*v2);
        }
        public double inverse(double x) {
            throw new RuntimeException("MaxwellBoltxmann.Distribution.inverse not defined");
        }
        public void setTemperature(double t) {
            temperature = t;
            update();
            
        }
        private void update() {
            //c = Math.pow(0.5*mass/Math.PI/temperature,1.5);
            a = 0.5*mass/temperature;
            switch(space.D()) {
                case 1: c = Math.sqrt(2.0/Math.PI*mass/temperature); break;
                case 2: c = mass/temperature; break;
                case 3: c = Math.sqrt(2.0/Math.PI)*Math.pow(mass/temperature, 1.5); break;
            }
         //   c /= space.sphereArea(1.0);
        }
        public double getTemperature() {return temperature;}
        public etomica.units.Dimension getTemperatureDimension() {return etomica.units.Temperature.DIMENSION;}

        public void setMass(double m) {
            mass = m;
            update();
        }
        public double getMass() {return mass;}
        public etomica.units.Dimension getMassDimension() {return Mass.DIMENSION;}
    }//end of Distribution
        
    
    
}
