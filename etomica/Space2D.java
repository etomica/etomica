package simulate;
import java.awt.Graphics;
import java.awt.Color;
import java.util.Random;

public class Space2D extends Space {
    
    public static final int D = 2;
    public final int D() {return D;}
    
    public Space.Coordinate makeCoordinate(Space.Occupant o) {return new Coordinate(o);}
    public Space.Vector makeVector() {return new Vector();}
    public Space.CoordinatePair makeCoordinatePair(Space.Coordinate c1, Space.Coordinate c2, Space.Boundary boundary) {return new CoordinatePair(c1, c2, (Boundary)boundary);}
    public Space.CoordinatePair makeCoordinatePair(Space.Boundary boundary) {return new CoordinatePair((Boundary)boundary);}
    public Space.Boundary makeBoundary(int b) {
        switch(b) {
            case(Boundary.NONE):            return new BoundaryNone();
            case(Boundary.PERIODIC_SQUARE): return new BoundaryPeriodicSquare();
            default:                        return null;
        }
    }
        
    public final static class Vector extends Space.Vector {  //declared final for efficient method calls
        public static final Random random = new Random();
        public static final Vector ORIGIN = new Vector(0.0,0.0);
        double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double a1, double a2) {x = a1; y = a2;}
        public double component(int i) {return (i==0) ? x : y;}
        public void setComponent(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Space.Vector u) {E((Vector)u);}
        public void PE(Space.Vector u) {PE((Vector)u);}
        public void ME(Space.Vector u) {ME((Vector)u);}
        public void TE(Space.Vector u) {TE((Vector)u);}
        public void DE(Space.Vector u) {DE((Vector)u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double a) {x = a; y = a;}
        public void E(int i, double a) {if(i==0) x = a; else y = a;}  //assumes i = 0 or 1
        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y;}
        public void PE(Vector u) {x += u.x; y += u.y;}
        public void ME(Vector u) {x -= u.x; y -= u.y;}
        public void PE(int i, double a) {if(i==0) x += a; else y += a;}
        public void TE(double a) {x *= a; y *= a;}
        public void TE(int i, double a) {if(i==0) x *= a; else y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
        public void randomStep(double d) {x += (2.*random.nextDouble()-1.0)*d; y+= (2.*random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = random.nextDouble()*d; y = random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = random.nextDouble()*dx; y = random.nextDouble()*dy;}
        public void setRandom(Vector u) {setRandom(u.x,u.y);}
        public void setRandomCube() {
            x = random.nextDouble() - 0.5; 
            y = random.nextDouble() - 0.5;
        }
        public void setRandomSphere() {
            x = Math.cos(2*Math.PI*random.nextDouble()); 
            y = Math.sqrt(1.0 - x*x);
        }
        public void setToOrigin() {x = ORIGIN.x; y = ORIGIN.y;}
        public void randomDirection() {
            x = Math.cos(2*Math.PI*random.nextDouble()); 
            y = Math.sqrt(1.0 - x*x);
        }
    }
    
    protected static final class CoordinatePair extends Space.CoordinatePair {  
        Coordinate c1;
        Coordinate c2;
        final Boundary boundary;
        final Vector dimensions;   //assumes this is not transferred between phases
        private final Vector dr = new Vector();
        private double drx, dry, dvx, dvy;
        public CoordinatePair() {boundary = new BoundaryNone(); dimensions = (Vector)boundary.dimensions();}
        public CoordinatePair(Space.Boundary b) {boundary = (Boundary)b; dimensions = (Vector)boundary.dimensions();}
        public CoordinatePair(Space.Coordinate c1, Space.Coordinate c2, Space.Boundary b) {
            boundary = (Boundary)b;
            dimensions = (Vector)boundary.dimensions();           
            if(c1 != null && c2 != null) {reset(c1,c2);}
        }
        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
//                dr.x -= (dr.x > 0.0) ? Math.floor(dr.x+0.5) : Math.ceil(dr.x-0.5);
//                dr.y -= (dr.y > 0.0) ? Math.floor(dr.y+0.5) : Math.ceil(dr.y-0.5);
            boundary.apply(dr);
            drx = dr.x; 
            dry = dr.y;
            r2 = drx*drx + dry*dry;
            double rm1 = c1.parent().rm();
            double rm2 = c2.parent().rm();
            dvx = (rm2*c2.p.x - rm1*c1.p.x);  
            dvy = (rm2*c2.p.y - rm1*c1.p.y);  
        }
                
//        public double r2() {
//            return r2;
//        }
        public double v2() {
            double rm1 = c1.parent().rm();
            double rm2 = c2.parent().rm();
            return dvx*dvx + dvy*dvy;
        }
        public double vDotr() {
            double rm1 = c1.parent().rm();
            double rm2 = c2.parent().rm();
            return drx*dvx + dry*dvy;
        }
        public void push(double impulse) {  //changes momentum in the direction joining the atoms
            c1.p.x += impulse*drx;
            c1.p.y += impulse*dry;
            c2.p.x -= impulse*drx;
            c2.p.y -= impulse*dry;
        }
        public void setSeparation(double r2New) {
            double ratio = c2.parent().mass()*c1.parent().rm();  // (mass2/mass1)
            double delta = (Math.sqrt(r2New/r2()) - 1.0)/(1+ratio);
            c1.r.x -= ratio*delta*drx;
            c1.r.y -= ratio*delta*dry;
            c2.r.x += delta*drx;
            c2.r.y += delta*dry;
            //need call reset?
        }
    }

    static class Coordinate extends Space.Coordinate {
        public final Vector r = new Vector();  //Cartesian coordinates, scaled by volume
        public final Vector p = new Vector();  //Momentum vector
        public final Vector posn = new Vector();  //Unscaled coordinates
        public Coordinate(Space.Occupant o) {super(o);}
        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.component(i);}
        public double momentum(int i) {return p.component(i);}
        public Space.Vector makeVector() {return new Vector();}
        public final double kineticEnergy(double mass) {return 0.5*p.squared()/mass;}
    } 
    
    public interface Boundary extends Space.Boundary {
        public static final int NONE = 0;
        public static final int PERIODIC_SQUARE = 1;
        public void apply(Vector dr);
        public void centralImage(Vector r);
    }

    /**
     * Class for implementing no periodic boundary conditions
     */
    protected static final class BoundaryNone implements Boundary {
        private final Vector temp = new Vector();
        private final double[][] shift0 = new double[0][D];
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        public static final Random random = new Random();
        public void apply(Space.Vector dr) {}
        public void centralImage(Space.Vector r) {}
        public void apply(Vector dr) {}
        public void centralImage(Vector r) {}
        public double volume() {return Double.MAX_VALUE;}
        public void inflate(double s) {}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
            temp.x = random.nextDouble(); 
            temp.y = random.nextDouble(); 
            return temp;
        }
        public void draw(Graphics g, int[] origin, double scale) {}
    }

    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    protected static final class BoundaryPeriodicSquare implements Boundary {
        private final Vector temp = new Vector();
        public static final Random random = new Random();
        private final double[][] shift0 = new double[0][D];
        private final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
        private final double[][] shift3 = new double[3][D];
        public BoundaryPeriodicSquare() {this(1.0,1.0);}
        public BoundaryPeriodicSquare(double lx, double ly) {dimensions.x = lx; dimensions.y = ly;}
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*random.nextDouble(); 
            temp.y = dimensions.y*random.nextDouble(); 
            return temp;}
        public void apply(Space.Vector dr) {apply((Vector)dr);}
        public void apply(Vector dr) {
            dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
            dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
        }
        public void centralImage(Space.Vector r) {centralImage((Vector)r);}
        public void centralImage(Vector r) {
            r.x -= dimensions.x * ((r.x > 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            r.y -= dimensions.y * ((r.y > 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
        }
        public void inflate(double scale) {dimensions.TE(scale);}
        public double volume() {return dimensions.x * dimensions.y;}
        public void draw(Graphics g, int[] origin, double scale) {
            g.setColor(Color.gray);
            double toPixels = scale*DisplayConfiguration.SIM2PIXELS;
            g.drawRect(origin[0],origin[1],(int)(toPixels*dimensions.component(0))-1,(int)(toPixels*dimensions.component(1))-1);
            }
        /** Computes origins for periodic images
        */
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells+1)*(2*nShells+1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                for(int j=-nShells; j<=nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x;
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }

        /** Returns coordinate shifts needed to draw all images that overflow into central image
         * 0, 1, or 3 shifts may be returned
         */
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {
            Vector r = (Vector)rr;
            int shiftX = 0;
            int shiftY = 0;
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if(shiftX == 0) {
                if(shiftY == 0) {
                    return shift0;
                }
                else {
                    shift1[0][0] = 0.0;
                    shift1[0][1] = shiftY*dimensions.y;
                    return shift1;
                }
            }
            else { //shiftX != 0
                if(shiftY == 0) {
                    shift1[0][0] = shiftX*dimensions.x;
                    shift1[0][1] = 0.0;
                    return shift1;
                }
                else {
                    shift3[0][0] = shiftX*dimensions.x;
                    shift3[0][1] = 0.0;
                    shift3[1][0] = 0.0;
                    shift3[1][1] = shiftY*dimensions.y;
                    shift3[2][0] = shift3[0][0];
                    shift3[2][1] = shift3[1][1];
                    return shift3;
                }
            }
        } //end of getOverflowShifts
    }  //end of BoundarySquarePeriodic
            
}