/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.Constants;
import etomica.Phase;
import etomica.Constants.Alignment;
import etomica.units.Dimension;
import etomica.units.Kelvin;


public final class AtomTypeWall extends AtomType {
        
        int thickness = 4;  //thickness when drawn to screen (if horizontal or vertical)
        int[] drawShift; //specifies simple shift (in pixels) of wall when drawn to screen
        private boolean vertical, horizontal, wide;
        private double cosX, sinX, tanX, cosY, sinY, tanY, cosZ, sinZ, tanZ;
//        double[] f = new double[Space.D];   //force on wall
//        double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
        protected double xAngle, yAngle, zAngle;   //orientation (radians) with respect to positive x-axis, taking r at the origin
        protected double length;  //length of line in simulation units
        protected boolean longWall;  //set to true if wall is infinite in length
        //  these may belong in integratorAgent
        protected double temperature = Kelvin.UNIT.toSim(300.);
        protected boolean adiabatic = true;
        private Constants.Alignment alignment;
        public double pAccumulator, qAccumulator; //net sum of momentum and heat transferred to wall
        
        public AtomTypeWall(AtomFactory creator, double m, double l, double x, double y, double z) {
            super(creator, m);
            setLength(l);
            setXAngle(x);
            setYAngle(y);
            setZAngle(z);
            drawShift = new int[creator.space.D()];
        }
        
        private void checkAlignment() {
          horizontal = (xAngle == 0.0) && (zAngle == 0.0);
          vertical = (xAngle == 0.0) && (zAngle == Math.PI/2);
          wide = (yAngle == Math.PI/2) && (zAngle == Math.PI/2);
          if(horizontal) alignment = Constants.HORIZONTAL;
          if(vertical) alignment = Constants.VERTICAL;
          if(wide) alignment = Constants.WIDTH;
        }
        
        public void setXAngle(double t) {
            double pi = Math.PI;
            xAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosX = Math.cos(xAngle);
            sinX = Math.sin(xAngle);
            tanX = Math.tan(xAngle);
            checkAlignment();
        }
        public void setYAngle(double t) {
            double pi = Math.PI;
            yAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosY = Math.cos(yAngle);
            sinY = Math.sin(yAngle);
            tanY = Math.tan(yAngle);
            checkAlignment();
        }
        public void setZAngle(double t) {
            double pi = Math.PI;
            zAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosZ = Math.cos(zAngle);
            sinZ = Math.sin(zAngle);
            tanZ = Math.tan(zAngle);
            checkAlignment();
        }
        public final double getXAngle() {return(xAngle);}
        public final double getYAngle() {return(yAngle);}
        public final double getZAngle() {return(zAngle);}
        public final double getSinX() {return(sinX);}
        public final double getSinY() {return(sinY);}
        public final double getSinZ() {return(sinZ);}
        public final double getCosX() {return(cosX);}
        public final double getCosY() {return(cosY);}
        public final double getCosZ() {return(cosZ);}
        public final double getTanX() {return(tanX);}
        public final double getTanY() {return(tanY);}
        public final double getTanZ() {return(tanZ);}
        
        public final void setAlignment(Constants.Alignment a) {
            alignment = a;
            if(a == Constants.VERTICAL) {
              setXAngle(0.0);
              setYAngle(0.0);
              setZAngle(Math.PI/2);
            } else if(a == Constants.HORIZONTAL) {
              setXAngle(0.0);
              setYAngle(0.0);
              setZAngle(0.0);
            } else {
              setXAngle(Math.PI/2);
              setYAngle(0.0);
              setZAngle(0.0);
            }
        }
        public final Constants.Alignment getAlignment() {return alignment;}
        
        public final int[] getDrawShift() {return drawShift;}
            
        public final boolean isVertical() {return vertical;}
        public final boolean isHorizontal() {return horizontal;}
        public final boolean isWide() {return wide;}
        
        public final int getThickness() {return thickness;}
        public final void setThickness(int thickness) {this.thickness = thickness;}
        
        public final double getLength() {return length;}
        public final void setLength(double length) {
            this.length = length; 
        }
        
        public final void setLongWall(boolean b) {longWall = b;}
        public final boolean isLongWall() {return longWall;}
        
        public final double getTemperature() {return temperature;}
//        public final void setTemperature(int t) {setTemperature((double)t);}  //for connection to sliders, etc.
        public final void setTemperature(double t) {temperature = t;}
        
        public final boolean isAdiabatic() {return adiabatic;}
        public final void setAdiabatic(boolean a) {adiabatic = a;}
        
        //Meter that evaluates pressure based on accumulated momentum from collisions with wall
        public class MeterPressure extends etomica.data.meter.MeterScalar {
            
            private double timeSum;
            public MeterPressure() {
               super();
            }
            public etomica.units.Dimension getDimension() {return etomica.units.Dimension.PRESSURE;}

            // Phase argument is not used
            public double getDataAsScalar(Phase p) {
                double flux = pAccumulator/timeSum;   //divide by time interval
                flux /= length; //divide by area
                return flux;
            }
        }//end of MeterPressure
    }