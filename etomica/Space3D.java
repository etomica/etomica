package etomica;
import etomica.units.*;

/**
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
 
 /* History of changes
  * 7/24/02 started recording change history
  * 7/24/02 (DW) modified RotationTensor
  *         (DW) change to CoordinatePair to use vector operators, needed to ensure proper working with CoordinateGroup
  *         (?) unknown changes to Orientation
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  */

public class Space3D extends Space implements EtomicaElement {

    public static String version() {return "Space3D:01.07.11/"+Space.VERSION;}
    public static final int D = 3;
    public final int D() {return D;}
    public final int powerD(int n) {return n*n*n;}
    public final double powerD(double a) {return a*a*a;}
    public int[] makeArrayD(int i) {return new int[] {i, i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d, d};}
    
    public static final Vector ORIGIN = new Vector();
    public final Space.Vector origin() {return ORIGIN;}
    
    public Space3D() {super(3);}
    
    public double sphereVolume(double r) {return (Math.PI*4.0*r*r*r/3.0);}
    public double sphereArea(double r)  {return (Math.PI*4*r*r);}
    public Space.Vector makeVector() {return new Vector();}
    public Space.Orientation makeOrientation() {return new Orientation();}
    public Space.Tensor makeTensor() {return new Tensor();}
    public Space.Tensor makeRotationTensor() {return new RotationTensor();}
    public Space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        else return new Coordinate(a);
    }
    public Space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}

    public Space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
   //     else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Three-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        Vector.WORK.y = u1.y - u2.y;
        Vector.WORK.z = u1.z - u2.z;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x + Vector.WORK.y*Vector.WORK.y + Vector.WORK.z*Vector.WORK.z;
    }
   
    
    public static final class Vector extends Space.Vector{
        public static final Vector ORIGIN = new Vector(0.0, 0.0, 0.0);
        public static final Vector WORK = new Vector();
        private double x, y, z;
        public int length() {return D;}
        public int D() {return D;}
        public Vector () {x = 0.0; y = 0.0; z = 0.0;}
        public Vector (double a1, double a2, double a3) {x = a1; y = a2; z = a3;}
        public Vector (double[] a) {x = a[0]; y = a[1]; z = a[2];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public String toString() {return "("+x+", "+y+", "+z+")";}
        public double x(int i) {return((i==0) ? x : (i==1) ? y : z);}
        public double[] toArray() {return new double[] {x, y, z};}
        public boolean equals(Space.Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x) && (y == v.y) && (z == v.z);}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y + z*z);
            result[1] = Math.acos(z/result[0]); //theta
            result[2] = Math.atan2(x, y);  //phi
        }
        public void E(Vector u) {x = u.x; y = u.y; z = u.z;}
        public void E(double a) {x = a; y = a; z = a;}
        public void E(double a, double b, double c) {x = a; y = b; z = c;}
        public void E(int i, double a) {if (i==0) x=a; else if (i==1) y=a; else z=a;}
        public void E(double[] u) {x = u[0]; y = u[1]; z = u[2];}  //should check length of array for exception
        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y; z = a1*u1.z;}
        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y; z += a1*u1.z;}
        public void PE(Vector u) {x += u.x; y += u.y; z += u.z;}
        public void PE(double a) {x += a; y += a; z += a;}
        public void ME(Vector u) {x -= u.x; y -= u.y; z -= u.z;}
        public void PE(int i, double a) {if (i==0) x+=a; else if (i==1) y+=a; else z+=a;}
        public void TE(double a) {x *= a; y *= a; z *= a;}
        public void TE(int i, double a) {if (i==0) x*=a; else if (i==1) y*=a; else z*=a;}
        public void TE(Vector u) {x *= u.x; y *= u.y; z *= u.z;}
        public void DE(double a) {x /= a; y /= a; z /= a;}
        public void DE(Vector u) {x /= u.x; y /= u.y ; z /= u.z;}
        public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
            y = v1.y + v2.y;
            z = v1.z + v2.z;
        }
        public void Ev1Mv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
            y = v1.y - v2.y;
            z = v1.z - v2.z;
        }
        public void mod(Space.Vector u) {
            mod((Vector)u);
        }
        public void mod(Vector u) {
            while(x > u.x) x -= u.x;
            while(x < 0.0) x += u.x;
            while(y > u.y) y -= u.y;
            while(y < 0.0) y += u.y;
            while(z > u.z) z -= u.z;
            while(z < 0.0) z += u.z;
        }
        public void mod(double a) {
            while(x > a)   x -= a;
            while(x < 0.0) x += a;
            while(y > a)   y -= a;
            while(y < 0.0) y += a;
            while(z > a)   z -= a;
            while(z < 0.0) z += a;
        }

        
        public Space.Vector P(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; WORK.y = y+u1.y; WORK.z = z+u1.z; return WORK;}
        public Space.Vector M(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; WORK.y = y-u1.y; WORK.z = z-u1.z; return WORK;}
        public Space.Vector T(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; WORK.y = y*u1.y; WORK.z = z*u1.z; return WORK;}
        public Space.Vector D(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; WORK.y = y/u1.y; WORK.z = z/u1.z; return WORK;}
        public Space.Vector abs() {WORK.x = (x>0)?x:-x; WORK.y = (y>0)?y:-y; WORK.z = (z>0)?z:-z; return WORK;}
        public double min() {return (x < y) ? (x<z)?x:z : (y<z)?y:z;}
        public double max() {return (x > y) ? (x>z)?x:z : (y>z)?y:z;}
        public double squared() {return x*x + y*y + z*z;}
        public double dot(Vector u) {return x*u.x + y*u.y + z*u.z;}
        public void transform(Space.Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {
            x = A.xx*x + A.xy*y + A.xz*z; 
            y = A.yx*x + A.yy*y + A.yz*z;
            z = A.zx*x + A.zy*y + A.zz*z;
        }
        public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
        public void transform(Boundary b, Vector r0, Tensor A) {
            WORK.x = x - r0.x; WORK.y = y - r0.y; WORK.z = z - r0.z;
            b.nearestImage(WORK);
            x = r0.x + A.xx*WORK.x + A.xy*WORK.y + A.xz*WORK.z;
            y = r0.y + A.yx*WORK.x + A.yy*WORK.y + A.yz*WORK.z;
            z = r0.z + A.zx*WORK.x + A.zy*WORK.y + A.zz*WORK.z;
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1)*d; y +=(2.*Simulation.random.nextDouble()-1)*d; z +=(2.*Simulation.random.nextDouble()-1)*d;}
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d; y = Simulation.random.nextDouble()*d; z = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy, double dz) {x = Simulation.random.nextDouble()*dx; y = Simulation.random.nextDouble()*dy; z = Simulation.random.nextDouble()*dz;}
        public void setRandom(Vector u) {setRandom(u.x, u.y, u.z);}
        public void setRandomCube() {
            x = Simulation.random.nextDouble() - 0.5;
            y = Simulation.random.nextDouble() - 0.5;
            z = Simulation.random.nextDouble() - 0.5;
        }
        public void setX(int a, double d) { if(a==0) x=d; else if(a==1) y=d; else z=d;}
        public void setRandomSphere() {//check before using
            double z1 = 0.0;
            double z2 = 0.0;
            double z3 = 0.0;
            double rsq = Double.MAX_VALUE;
            while(rsq > 1.0) {
                
                z1 = 1.0 - 2.0*Simulation.random.nextDouble();
                z2 = 1.0 - 2.0*Simulation.random.nextDouble();
                z3 = 1.0 - 2.0*Simulation.random.nextDouble();
        
                rsq = z1*z1+z2*z2+z3*z3;
            }
            double r = Math.sqrt(rsq);
            x = z1/r;
            y = z2/r;
            z = z3/r;
        }
        public void randomRotate(double thetaStep){//check before using
            //could be made more efficient by merging with setRandomSphere
            if(thetaStep == 0.0) return;
            if(Math.abs(thetaStep) > Math.PI) {setRandomSphere(); return;}
            double xOld = x; double yOld = y; double zOld = z;
            double r = Math.sqrt(x*x + y*y + z*z);
            double dotMin = r*Math.cos(thetaStep);
            do {
                setRandomSphere();
            } while(xOld*x + yOld*y + zOld*z < dotMin);
            x *= r; y *= r; z*=r;
        }
        
        public void setToOrigin() {x = ORIGIN.x; y = ORIGIN.y; z = ORIGIN.z;}

        public void E(Space.Vector u) {E((Vector) u);}
        public void PE(Space.Vector u) {PE((Vector) u);}
        public void TE(Space.Vector u) {TE((Vector) u);}
        public void ME(Space.Vector u) {ME((Vector) u);}
        public void DE(Space.Vector u) {DE((Vector) u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
        public Space3D.Vector cross(Space2D.Vector u) {
            Space3D.Vector.WORK.x = -z*u.x(1);
            Space3D.Vector.WORK.y = z*u.x(0);
            Space3D.Vector.WORK.z = x*u.x(1) - y*u.x(0);
            return Space3D.Vector.WORK;
        }
        public Space3D.Vector cross(Space3D.Vector u) {//not thread safe
            Space3D.Vector.WORK.x = y*u.z - z*u.y;
            Space3D.Vector.WORK.y = z*u.x - x*u.z;
            Space3D.Vector.WORK.z = x*u.y - y*u.x;
            return Space3D.Vector.WORK;
        }
        public void XE(Vector u) {//cross product
            double xNew = y*u.z - z*u.y;
            double yNew = z*u.x - x*u.z;
            z = x*u.y - y*u.x;
            y = yNew;
            x = xNew;
        }
        public void normalize() {
            double norm = 1./Math.sqrt(x*x + y*y + z*z);
            x *= norm;
            y *= norm;
            z *= norm;
        }
    }
    
    public static class Tensor implements Space.Tensor {
        double xx, xy, xz, yx, yy, yz, zx, zy, zz;
        public static final Tensor ORIGIN = new Tensor();
        public static final Tensor WORK = new Tensor();
        public Tensor () {xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;}
        public Tensor (double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
            this.xx=xx; this.xy=xy; this.xz=xz; this.yx=yx; this.yy=yy; this.yz=yz; this.zx=zx; this.zy=zy; this.zz=zz;
        }
        public double component(int i, int j) {
            return ( i==0 ) ? ( (j==0) ? xx : ( j==1 ? xy : xz ) ) : ( (i==1) ? ( (j==0) ? yx : ( (j==1) ? yy : yz ) ) : ( (j==0) ? zx : ((j==1) ? zy : zz)));
        }
        public int length() {return D;}
        public void setComponent(int i, int j, double d) {
            if (i==0) {if (j==0) {xx = d;} else if (j==1) {xy = d;} else xz = d;}
            else if (i==1) {if (j==0) {yx = d;} else if (j==1) {yy=d;} else yz = d;}
            else {if (j==0) {zx = d;} else if (j==1) {zy = d;} else zz = d;}
        }
        public void E(Tensor t) {xx=t.xx; xy=t.xy; xz=t.xz; yx=t.yx; yy=t.yy; yz=t.yz; zx=t.zx; zy=t.zy; zz=t.zz;}
        public void E(Vector u1, Vector u2) {xx=u1.x*u2.x; xy=u1.x*u2.y; xz=u1.x*u2.z; yx=u1.y*u2.x; yy=u1.y*u2.y; yz=u1.y*u2.z; zx=u1.z*u2.x; zy=u1.z*u2.y; zz=u1.z*u2.z;}
        public void E(double a) {xx=xy=xz=yx=yy=yz=zx=zy=zz=a;}
        public void PE(Tensor t) {xx+=t.xx; xy+=t.xy; xz+=t.xz; yx+=t.yx; yy+=t.yy; yz+=t.yz; zx+=t.zx; zy+=t.zy; zz+=t.zz;}
        public void PE(int i, int j, double d) {
            if (i==0) {if (j==0) {xx += d;} else if (j==1) {xy += d;} else xz += d;}
            else if (i==1) {if (j==0) {yx += d;} else if (j==1) {yy += d;} else yz += d;}
            else {if (j==0) {zx += d;} else if (j==1) {zy += d;} else zz += d;}
        }
        public double trace() {return xx+yy+zz;}
        public void E(Space.Tensor t) {E((Tensor)t);}
        public void E(Space.Vector u1, Space.Vector u2) {E((Vector)u1, (Vector)u2);}
        public void PE(Space.Tensor t) {PE((Tensor) t);}
        public void PE(Space.Vector u1, Space.Vector u2) {PE((Vector)u1, (Vector)u2);}
        public void TE(double a) {xx*=a; xy*=a; xz*=a; yx*=a; yy*=a; yz*=a; zx*=a; zy*=a; zz*=a;}
    }
    
    public static class RotationTensor extends Tensor implements Space.RotationTensor {
        public RotationTensor() {super(); reset();}
        public void reset() {
            xx = 1.0; xy = 0.0; xz = 0.0;
            yx = 0.0; yy = 1.0; yz = 0.0;
            zx = 0.0; zy = 0.0; zz = 1.0;
        }
        /**
         * Sets tensor to for rotation by the given angle about a randomly selected axis.
         */
        public void setAxial(double theta) {
            int n = (int)(Simulation.random.nextDouble()*3);
            setAxial(n, theta);
        }
        /**
         * Sets tensor for rotation about the indicated axis (0=x,1=y,2=z) by 
         * the given angle.
         */
        public void setAxial(int i, double theta) {
            double st = Math.sin(theta);
            double ct = Math.cos(theta);
            switch(i) {
                case 0: xx = 1.; xy = 0.; xz = 0.;
                        yx = 0.; yy = ct; yz = -st;
                        zx = 0.; zy = st; zz = ct;
                        return;
                case 1: xx = ct; xy = 0.; xz = -st;
                        yx = 0.; yy = 1.; yz = 0.;
                        zx = st; zy = 0.; zz = ct;
                        return;
                case 2: xx = ct; xy = -st; xz = 0.;
                        yx = st; yy = ct;  yz = 0.;
                        zx = 0.; zy = 0.;  zz = 1.;
                        return;
                default: throw new IllegalArgumentException("Improper axis specified for Space3D.RotationTensor.setAxial");
            }
        }
        /**
         * Not yet implemented.
         */
        public void setAngles(double[] angles) {
            throw new RuntimeException("Space3D.CoordinateGroup.setAngles() not yet implemented");
        }
        public void invert() {
            double det = xx*yy*zz - xx*yz*zy - yx*xy*zz + yx*xz*zy + zx*xy*yz - zx*xz*yy ;
            double xx1 = (yy*zz - yz*zy)/det;
            double xy1 = (-xy*zz + xz*zy)/det;
            double xz1 = (xy*yz - xz*yy)/det;
            
            double yx1 = (-yx*zz + yz*zx)/det;
            double yy1 = (xx*zz - xz*zx)/det;
            double yz1 = (-xx*yz + xz*yx)/det;
            
            double zx1 = (yx*zy - yy*zx)/det;
            double zy1 = (-xx*zy + xy*zx)/det;
            double zz1 = (xx*yy - xy*yx)/det;
            
            this.xx = xx1; this.xy = xy1; this.xz = xz1;
            this.yx = yx1; this.yy = yy1; this.yz = yz1;
            this.zx = zx1; this.zy = zy1; this.zz = zz1;
            
        }
        /**
         * Method to test rotation tensor.
         */
        public static void main (String[] args) {
            
            Vector r1 = new Vector(2,2,3);
            System.out.println("r1_before" + r1.toString());
            Tensor tensor = new Tensor(1,2,0,1,1,2,0,0,1);
            RotationTensor tensor2 = new RotationTensor();
            tensor2.E(tensor);
            System.out.println("tensor2_before " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
            System.out.println();
        
            r1.transform(tensor2);
            System.out.println("r1_transform(tensor2)" + r1.toString());
            tensor2.invert();
            System.out.println("tensor2_invert " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
            //tensor2.setAxial(1, 2*Math.PI);
            //System.out.println("tensor2_rotate_360 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
            //System.out.println();
        
            //r1.transform(tensor2);
            //System.out.println("r1_afterInvert_andRotate360 " + r1.toString());
        }//end of main 
        
    }//end of Space3D.RotationTensor
    
    public static final class CoordinatePair extends Space.CoordinatePair {
        Coordinate c1;
        Coordinate c2;
        private final Vector dr = new Vector();
        private double dvx, dvy, dvz; //drx, dry, drz;
        public CoordinatePair() {super();}

        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.Ev1Mv2(c2.position(),c1.position());
         //   c2.position(); c1.position();
         //   dr.x = c2.r.x - c1.r.x;
         //   dr.y = c2.r.y - c1.r.y;
         //   dr.z = c2.r.z - c1.r.z;
            c1.atom.node.parentPhase().boundary().nearestImage(dr);
   //         drx = dr.x;
   //         dry = dr.y;
   //         drz = dr.z;
            r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
            /*  comment here if not doing hard dynamics
            double rm1 = c1.rm();
            double rm2 = c2.rm();
            dvx = rm2*c2.p.x - rm1*c1.p.x;
            dvy = rm2*c2.p.y - rm1*c1.p.y;
            dvz = rm2*c2.p.z - rm1*c1.p.z;
          //  */  //end of non-dynamics commenting
        }
            
        public void reset(Space3D.Vector M) {
            dr.x = c2.r.x - c1.r.x + M.x;
            dr.y = c2.r.y - c1.r.y + M.y;
            dr.z = c2.r.z - c1.r.z + M.z;
        //    drx = dr.x;
        //    dry = dr.y;
        //    drz = dr.z;
            r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        }
        public double r2() {
            return r2;
         /*   dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
            dr.z = c2.r.z - c1.r.z;
            c1.atom.node.parentPhase().boundary().nearestImage(dr);
            return dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;*/
        ///    dr.Ev1Mv2(c2.r, c1.r);
        ///    c1.atom.node.parentPhase().boundary().nearestImage(dr);
        ///    return dr.squared();
        }
            
        public Space.Vector dr() {return dr;}
        public double dr(int i) {return (i==0) ? dr.x : ((i==1) ? dr.y : dr.z);}
        public double dv(int i) {return (i==0) ? dvx : ((i==1) ? dvy : dvz);}
        public double v2() {return dvx*dvx + dvy*dvy + dvz*dvz;}
        public double vDot(Space.Vector u) {return vDot((Space3D.Vector)u);}
        public double vDot(Space3D.Vector u) {return dvx*u.x + dvy*u.y + dvz*u.z;}
        public double vDotr() {return dr.x*dvx + dr.y*dvy + dr.z*dvz;}
        public void push(double impulse) {
            c1.p.x += impulse*dr.x;
            c1.p.y += impulse*dr.y;
            c1.p.z += impulse*dr.z;
            c2.p.x -= impulse*dr.x;
            c2.p.y -= impulse*dr.y;
            c2.p.z -= impulse*dr.z;
        }
        public void setSeparation(double r2New) {
            double ratio = c2.mass()*c1.rm();
            double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1 + ratio);
            c1.r.x -= ratio*delta*dr.x;
            c1.r.y -= ratio*delta*dr.y;
            c1.r.z -= ratio*delta*dr.z;
            c2.r.x += ratio*delta*dr.x;
            c2.r.y += ratio*delta*dr.y;
            c2.r.z += ratio*delta*dr.z;
        }
    }
        
    public static class Coordinate extends Space.Coordinate {
        public final Vector r = new Vector();
        public final Vector p = new Vector();
        public final Vector rLast = new Vector();  //vector for saving position
        public final Vector work = new Vector();
        public Coordinate(Atom a) {super(a);}
                        
        public void transform(Space.Vector r0, Space.Tensor A) {
            r.transform((Boundary)atom.node.parentPhase().boundary(),(Vector)r0, (Tensor)A);
            atom.seq.moveNotify();
        }
        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.x(i);}
        public double momentum(int i) {return p.x(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {
            double tM = t*rm(); // t/mass
            r.x += p.x*tM;
            r.y += p.y*tM;
            r.z += p.z*tM;
            atom.seq.moveNotify();
        }
        /**
         * Scales positions of atoms by multiplying by given value.  Does not notify sequencers.
         */
        public void inflate(double s) {r.x *= s; r.y *= s; r.z *= s;}
        /**
         * Scales positions of atoms isotropically by multiplying by coordinate in each
         * direction by the corresponding element of the given vector.  Does not notify sequencers.
         */
        public void inflate(Space.Vector s) {Vector u = (Vector)s; r.x *= u.x; r.y *= u.y; r.z *= u.z;}
        
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(Space.Vector u) {
            r.PE((Vector)u); 
            atom.seq.moveNotify();
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(double d, Space.Vector u) {
            r.PEa1Tv1(d,(Vector)u);
            atom.seq.moveNotify();
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateTo(Space.Vector u) {
            r.E((Vector)u);
            atom.seq.moveNotify();
        }      
        public void replace() {r.E(rLast); atom.seq.moveNotify();}
        
        public void displaceBy(Space.Vector u) {rLast.E(r); translateBy((Vector)u);}
        public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,(Vector)u);}
        public void displaceTo(Space.Vector u) {rLast.E(r); translateTo((Vector)u);}  
        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
    //    public final void inflate(double s) {r.TE(s);}

        public void accelerateBy(Space.Vector u) {p.PE(u);}
        public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}
        public void accelerateTo(Space.Vector u) {p.E(u);}

        public void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
            if(isStationary()) {p.E(0.0); return;}
            double magnitude = Math.sqrt(mass()*temperature*(double)D);  //need to divide by sqrt(m) to get velocity
            momentum().setRandomSphere();
            momentum().TE(magnitude);
        }
    }//end of Coordinate
    
public static class CoordinateGroup extends Coordinate {

    private final AtomIteratorListSimple childIterator = new AtomIteratorListSimple();
    
    public CoordinateGroup(Atom a) {
        super(a);
        childIterator.setBasis(((AtomTreeNodeGroup)a.node).childList);
    }

    /**
        * Applies transformation to COM of group, keeping all internal atoms at same relative
        * positions.
        */
    public void transform(Space.Vector r0, Space.Tensor A) {
        work.E(position()); //work = r
        work.transform(atom.node.parentPhase().boundary(), r0, A);
        work.ME(r);//now work vector contains translation vector for COM
        translateBy(work);
    }
    public Space.Vector position() {
        r.E(0.0); double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            r.PEa1Tv1(a.coord.mass(), a.coord.position()); 
            massSum += a.coord.mass();
        }
        r.DE(massSum);
        return r;
    }
    public Space.Vector momentum() {
        p.E(0.0);
        childIterator.reset();
        while(childIterator.hasNext()) {
            p.PE(childIterator.next().coord.momentum());
        }
        return p;
    }
    public double position(int i) {
        double sum = 0.0; double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            sum += a.coord.mass()*a.coord.position(i); 
            massSum += a.coord.mass();
        }
        sum /= massSum;
        return sum;
    }
    public double momentum(int i) {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            sum += a.coord.mass()*a.coord.momentum(i);
        }
        return sum;
    }
    public double kineticEnergy() {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            sum += childIterator.next().coord.kineticEnergy();
        }
        return sum;
    }
    public void freeFlight(double t) {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.freeFlight(t);
        }
        atom.seq.moveNotify();
    }
    public void inflate(double scale) {
        work.E(position());
        work.TE(scale-1.0);
        translateBy(work);
    }
    public void inflate(Space.Vector scale) {
        scale.PE(-1.0);
        work.E(position());
        work.TE(scale);
        translateBy(work);
        scale.PE(1.0);
    }
    
    public void translateBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.translateBy(u);
        }
        atom.seq.moveNotify();
    }
    public void translateBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.translateBy(d, u);
        }
        atom.seq.moveNotify();
    }
    public void translateTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        translateBy(work);
    }
    public void displaceBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.displaceBy(u);
        }
        atom.seq.moveNotify();
    }
    public void displaceBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.displaceBy(d, u);
        }
        atom.seq.moveNotify();
    }
    public void displaceTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        displaceBy(work);
    }
    public void displaceToRandom(etomica.Phase p) {
        displaceTo(p.boundary().randomPosition());
    }
    public void replace() {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.replace();
        }
        atom.seq.moveNotify();
    }
    public void accelerateBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.accelerateBy(u);
        }
    }
    public void accelerateBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.accelerateBy(d, u);
        }
    }
    public void accelerateTo(Space.Vector u) {
        work.Ea1Tv1(-1.0/childIterator.size(),momentum());//probably need this first
        work.PE(u);
        accelerateBy(work);
    }
    public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
    public void randomizeMomentum(double temperature) {
        switch(((AtomTreeNodeGroup)atom.node).childAtomCount()) {
            case 0: return;
            case 1: ((AtomTreeNodeGroup)atom.node).firstChildAtom().coord.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                    return;
            default://multi-atom group
   //             work.E(0.0); double sum=0.0;
                childIterator.reset();
                while(childIterator.hasNext()) {
                    Atom a = childIterator.next();
                    a.coord.randomizeMomentum(temperature);
   //                 work.PE(a.coord.momentum());
   //                 sum++;
                }
   //             work.DE(-sum);
   //             childIterator.reset();
   //             while(childIterator.hasNext()) {
    //                childIterator.next().coord.accelerateBy(work);
   //             }
        }//end switch
    }//end randomizeMomentum
    
}//end of CoordinateGroup

    /*public static class CoordinateGroup extends Coordinate {
        public CoordinateGroup(AtomGroup a) {super(a);}
                
/*        public void updateMass() {
            double mass = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                mass += coord.mass();
            }
            rm = 1.0/mass;
        } 
     //   public double rm() {return 1.0/mass();}
        /**
         * Applies transformation to COM of group, keeping all internal atoms at same relative positions.
         * /
        public void transform(Space.Vector r0, Space.Tensor A) {
            work.E(position()); //work = r
            work.transform((Boundary)atom.node.parentPhase().boundary(),(Vector)r0, (Tensor)A);
            work.ME(r);//now work vector contains translation vector for COM
            translateBy(work);
        }
        public Space.Vector position() {
            r.E(0.0); double massSum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                r.PEa1Tv1(coord.mass(), coord.position()); massSum += coord.mass();
                if(coord == lastChild) break;
            }
            r.DE(massSum);
            return r;
        }
        public Space.Vector momentum() {
            p.E(0.0);
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                p.PE(coord.momentum());
                if(coord == lastChild) break;
            }
            return p;
        }
        public double position(int i) {
            double sum = 0.0; double massSum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.mass()*coord.position(i); massSum += coord.mass();
                if(coord == lastChild) break;
            }
            sum /= massSum;
            return sum;
        }
        public double momentum(int i) {
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.mass()*coord.momentum(i);
                if(coord == lastChild) break;
            }
            return sum;
        }
        public double kineticEnergy() {
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.kineticEnergy();
                if(coord == lastChild) break;
            }
            return sum;
        }
        public void freeFlight(double t) {
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.freeFlight(t);
                if(coord == lastChild) break;
            }
        }
        public void translateBy(Space.Vector u) {translateBy((Vector)u);}
        public void translateBy(Vector u0) {
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.translateBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void translateBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.translateBy(d, u0);
                if(coord == lastChild) break;
            }
        }
        public void translateTo(Space.Vector u) {
            work.Ea1Tv1(-1,position()); //position() uses work, so need this first
            work.PE((Vector)u);
            translateBy(work);
        }
        public void displaceBy(Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.displaceBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void displaceBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.displaceBy(d, u0);
                if(coord == lastChild) break;
            }
        }
        public void displaceTo(Space.Vector u) {
            work.E((Vector)u);
            work.ME(position());
            displaceBy(work);
        }
        public void displaceToRandom(etomica.Phase p) {
            displaceTo((Vector)p.boundary().randomPosition());
        }
        public void replace() {
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.replace();
                if(coord == lastChild) break;
            }
        }
        public void accelerateBy(Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.accelerateBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void accelerateBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.accelerateBy(d, u0);
                if(coord == lastChild) break;
            }
        }
        public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
        public void randomizeMomentum(double temperature) {
            switch(((AtomGroup)atom).node.childAtomCount()) {
                case 0: return;
                case 1: firstChild.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                        return;
                default://multi-atom group
                    work.E(0.0); double sum=0.0;
                    for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                        coord.randomizeMomentum(temperature);
                        work.PE(coord.momentum());
                        sum++;
                        if(coord == lastChild) break;
                    }
                    work.DE(-sum);
                    for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                        coord.accelerateBy(work);
                        if(coord == lastChild) break;
                    }
            }//end switch
        }//end randomizeMomentum
        
    }//end of CoordinateGroup
    */
    public static class OrientedCoordinate extends Coordinate implements Space.Coordinate.Angular {
        private double L = 0.0; //magnitude of angular momentum
        private final Space3D.Vector vector = new Space3D.Vector();//used to return vector quantities (be sure to keep x and y components zero)
        private final double[] I;
        private final Orientation orientation = new Orientation();
        public OrientedCoordinate(Atom a) {
            super(a);
            I = ((AtomType.SphericalTop)a.type).momentOfInertia();
        }
        public Space3D.Vector angularMomentum() {vector.z = L; return vector;}
        public Space3D.Vector angularVelocity() {vector.z = L/I[0]; return vector;}
        public void angularAccelerateBy(Space3D.Vector t) {L += t.z;}
        public Space.Orientation orientation() {return orientation;}
        public double kineticEnergy() {return super.kineticEnergy() + 0.5*L*L/I[0];}
        public void freeFlight(double t) {
            super.freeFlight(t);
            orientation.rotateBy(t*L/I[0]);//all elements of I equal for spherical top
        }
    }
     
     public static class Orientation extends Space.Orientation {
        public Orientation() {
            throw new RuntimeException("Space3D.Orientation should be checked for correctness before using");
        }
         //The rotation matrix A operates on the components of a vector in the space-fixed frame to yield the
         //components in the body-fixed frame
        private final double[][] A = new double[D][D];
        private final Vector[] bodyFrame = new Vector[] {new Vector(1.0,0.0, 0.0), new Vector(0.0,1.0,0.0), new Vector(0.0,0.0,1.0)};
        private final double[] angle = new double[3];
        private final Vector orientVector = new Vector(1.0,0.0,0.0);
        private boolean needToUpdateA = true;
        private transient double x1,y1,z1;//temp variable
        private transient Vector v1 = new Vector();
        public void E(Space.Orientation o) {E((Orientation)o);}
        public Space.Vector getOrientation(){return orientVector;}
        public void setOrientation(Space.Vector vect){orientVector.E(vect);}
        public void E(Orientation o) {
          angle[0] = o.angle[0];
          angle[1] = o.angle[1];
          angle[2] = o.angle[2];
          needToUpdateA = true;
          orientVector.E(o.getOrientation());
        }
        public Space.Vector[] bodyFrame() {return bodyFrame;}
        public double[] angle() {return angle;}
        public final void rotateBy(double dt[]) {
            rotateBy(0, dt[0]);
            rotateBy(1, dt[1]);
            rotateBy(2, dt[2]);
        }
        public final void rotateBy(double dt) {
          rotateBy(0, dt);
          rotateBy(1, dt);
          rotateBy(2, dt);
        }
        public void rotateBy(int i, double dt) {
            angle[i] += dt;
            if(angle[i] > Constants.TWO_PI) angle[i] -= Constants.TWO_PI;
            else if(angle[i] < 0.0) angle[i] += Constants.TWO_PI;
            needToUpdateA = true;
        }
          /**
        * Choose one of the three spaced fixed axes at random and rotate 
        * it by t (in radian).
        */
        public void randomRotation(double stepsize) {
            double t = (2.0*Simulation.random.nextDouble()-1.0)*stepsize;
            int i = (int)(3.0 *(Simulation.random.nextDouble()));
            double ct,st;
            
            ct= Math.cos(t);st =Math.sin(t);
            
            //rotateBy(i,t);
            
            switch(i){
             case 0:
                    x1=orientVector.x(0);
                    y1=ct*orientVector.x(1)+ st*orientVector.x(2);      
                    z1=ct*orientVector.x(2)-st*orientVector.x(1);
                    break;
             case 1:
                    x1=ct*orientVector.x(0)-st*orientVector.x(2);
                    y1=orientVector.x(1);
                    z1=ct*orientVector.x(2)+st*orientVector.x(0);
                    break;
             case 2:
                    x1=ct*orientVector.x(0)+st*orientVector.x(1);
                    y1=ct*orientVector.x(1)-st*orientVector.x(0);
                    z1=orientVector.x(2);
                    break;
            }
                
            orientVector.E(x1,y1,z1);
           // orientVector.normalize(); // Not needed
            
                
                
            
        }
        private final void updateRotationMatrix() {
            double ct0,ct1,ct2;
            double st0,st1,st2;
                        
            ct0=Math.cos(angle[0]);ct1=Math.cos(angle[1]);ct2=Math.cos(angle[2]);
            st0=Math.sin(angle[0]);st1=Math.sin(angle[1]);st2=Math.sin(angle[2]);
                     
            //LEFT HAND SYSTEM  
            A[0][0]= ct1*ct2;
            A[0][1]= -ct1*st2;
            A[0][2]= st1;
            A[1][0]= st0*st1*ct2+ct0*st2;
            A[1][1]= -st0*st1*st2+ct0*ct2;
            A[1][2]= -st0*ct1;
            A[2][0]= -ct0*st1*ct2+st0*st2;
            A[2][1]= ct0*st1*st2+st0*ct2;
            A[2][2]= ct0*ct1;
                      
            bodyFrame[0].E(A[0]);
            bodyFrame[1].E(A[1]);
            bodyFrame[2].E(A[2]);
            needToUpdateA = false;
         }
      //   public double[][] rotationMatrix() {return A;}
          public void convertToBodyFrame(Vector v) {
            if(needToUpdateA) updateRotationMatrix();
            v1.x = A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z;
            v1.y = A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z;
            v1.z = A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z;
            v.E(v1);
         }
         //V_space = A_transpose*V_body
         public void convertToSpaceFrame(Vector v) {
            if(needToUpdateA) updateRotationMatrix();
            v1.x = A[0][0]*v.x + A[1][0]*v.y + A[2][0]*v.z;
            v1.y = A[0][1]*v.x + A[1][1]*v.y + A[2][1]*v.z;
            v1.z = A[0][2]*v.x + A[1][2]*v.y + A[2][2]*v.z;
            
            v.E(v1);
         }
        public void convertToBodyFrame(Space.Vector v) {convertToBodyFrame((Vector)v);}
        public void convertToSpaceFrame(Space.Vector v) {convertToSpaceFrame((Vector)v);}
    }

    /*public static class Orientation extends Space.Orientation {
        public final Vector[] eArray = new Vector[] {new Vector(1.0,0.0,0.0)};
        public final Vector L = new Vector(); //angular momentum
        public final Vector e = eArray[0];
        private final Vector[] bodyFrame = new Vector[] {new Vector(1.0,0.0,0.0), 
                                                         new Vector(0.0,1.0,0.0),
                                                         new Vector(0.0,0.0,1.0)};
        public final double[] angle = new double[1];
        public Space.Vector[] direction() {return eArray;}
        public double[] angle() {angle[0] = Math.atan2(e.y,e.x); return angle;}
        public Space.Vector momentum() {return L;}
        public Space.Vector[] bodyFrame() {return bodyFrame;} 
        public void freeFlight(double I) {//motion of a spherical top
        }
    }
    */
    
    public static abstract class Boundary extends Space.Boundary {
        public static class Type extends Space.Boundary.Type {
            private Type(String label) {super(label);}
            public Constants.TypedConstant[] choices() {return TYPES;}
        }
        public static final String[] TAGS = {"None", "Periodic Square"/*, "Sliding Brick"*/};
        public static final Type NONE = new Type("None");
        public static final Type PERIODIC_SQUARE = new Type("Periodic Square");
  //      public static final Type SLIDING_BRICK = new Type("Sliding Brick");
        public static final Type[] TYPES = {NONE,PERIODIC_SQUARE/*,SLIDING_BRICK*/};
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public abstract void nearestImage(Vector dr);
        public abstract boolean centralImage(Vector r);
        public abstract boolean centralImage(Coordinate c);
    }
    
    
    public static final class BoundaryNone extends Boundary {
        private final Vector temp = new Vector();
        public final Vector dimensions = new Vector(Default.BOX_SIZE, Default.BOX_SIZE, Default.BOX_SIZE);
        public final Space.Vector dimensions() {return dimensions;}
        public BoundaryNone() {super();}
        public BoundaryNone(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.NONE;}
        public void nearestImage(Space.Vector dr) {}
        public boolean centralImage(Space.Vector r) {return false;}
        public void nearestImage(Vector dr) {}
        public boolean centralImage(Vector r) {return false;}
        public boolean centralImage(Coordinate c) {return false;}
        public double volume() {return dimensions.x*dimensions.y*dimensions.z;}
        public void inflate(double s) {dimensions.TE(s);}
        public void inflate(Space.Vector s) {dimensions.TE(s);}
        public void setDimensions(Space.Vector v) {dimensions.E(v);}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble();
            temp.y = dimensions.y*Simulation.random.nextDouble();
            temp.z = dimensions.z*Simulation.random.nextDouble();
            return temp;
        }
    }//end of BoundaryNone
        
    protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic  {
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx, double ly, double lz) {super(p);dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
        public BoundaryPeriodicSquare(double lx, double ly, double lz) {super();dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
        public Space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        private static Space.Tensor zilch = new Tensor();
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector();
        private final Vector dimensionsCopy = new Vector();
        private final Vector dimensionsHalf = new Vector();
        public final Space.Vector dimensions() {return dimensionsCopy;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble();
            temp.y = dimensions.y*Simulation.random.nextDouble();
            temp.z = dimensions.z*Simulation.random.nextDouble();
            return temp;
        }
        private final void updateDimensions() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
        }
        public void nearestImage(Space.Vector dr) {nearestImage((Vector) dr);}
        public void nearestImage(Vector dr) {
      //      dr.x -= dimensions.x*((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x + 0.5) : Math.ceil(dr.x/dimensions.x - 0.5));
      //      dr.y -= dimensions.y*((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y + 0.5) : Math.ceil(dr.y/dimensions.y - 0.5));
      //      dr.z -= dimensions.z*((dr.z > 0.0) ? Math.floor(dr.z/dimensions.z + 0.5) : Math.ceil(dr.z/dimensions.z - 0.5));
      //      final double dimxHalf = 0.5*dimensions.x;
      //      final double dimyHalf = 0.5*dimensions.y;
      //      final double dimzHalf = 0.5*dimensions.z;
        /*    while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
            while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
            while(dr.y > +dimensionsHalf.y) dr.y -= dimensions.y;
            while(dr.y < -dimensionsHalf.y) dr.y += dimensions.y;
            while(dr.z > +dimensionsHalf.z) dr.z -= dimensions.z;
            while(dr.z < -dimensionsHalf.z) dr.z += dimensions.z;*/
            dr.PE(dimensionsHalf);
            dr.mod(dimensions);
            dr.ME(dimensionsHalf);
            //System.out.println("dimesions = "+dimensions);
        //    System.out.print(dr.x+"  ");dr.x %= dimensionsHalf.x; System.out.println(dr.x);
        //    dr.x = ((dr.x + dimensions.x) % dimensions.x) - dimensionsHalf.x;
        //    dr.y = ((dr.y + dimensions.y) % dimensions.y) - dimensionsHalf.y;
        //    dr.z = ((dr.z + dimensions.z) % dimensions.z) - dimensionsHalf.z;
        }
        public boolean centralImage(Coordinate c) {return centralImage(c.r);}
        public boolean centralImage(Space.Vector r) {return centralImage((Vector) r);}
        public boolean centralImage(Vector r) {
            temp.E(r);
            r.mod(dimensions);
            return temp.equals(r);
       /*     while(r.x > dimensions.x) r.x -= dimensions.x;
            while(r.x < 0.0)          r.x += dimensions.x;
            while(r.y > dimensions.y) r.y -= dimensions.y;
            while(r.y < 0.0)          r.y += dimensions.y;
            while(r.z > dimensions.y) r.z -= dimensions.z;
            while(r.z < 0.0)          r.z += dimensions.z;*/
         //   r.x -= dimensions.x* ((r.x>0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x - 1.0));
         //   r.y -= dimensions.y *((r.y>0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y - 1.0));
         //   r.z -= dimensions.z *((r.z>0) ? Math.floor(r.z/dimensions.z) : Math.ceil(r.z/dimensions.z - 1.0));
        }
        public void inflate(double scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void inflate(Space.Vector scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void setDimensions(Space.Vector v) {dimensions.E(v); updateDimensions();}
        public double volume() {return dimensions.x*dimensions.y*dimensions.z;}
                
        /**
         * imageOrigins and getOverFlowShifts are both probably incorrect, if they are
         * even completed.  They should definitely be checked before being implemented.
         */
        
        int shellFormula, nImages, i, j, k, m;
        double[][] origins;
        public double[][] imageOrigins(int nShells) {
            shellFormula = (2 * nShells) + 1;
            nImages = shellFormula*shellFormula*shellFormula-1;
            origins = new double[nImages][D];
            for (k=0,i=-nShells; i<=nShells; i++) {
                for (j=-nShells; j<=nShells; j++) {
                    for (m=-nShells; m<=nShells; m++) {
                        if ((i==0 && j==0) && m==0 ) {continue;}
                        origins[k][0] = i*dimensions.x;
                        origins[k][1] = j*dimensions.y;
                        origins[k][2] = m*dimensions.z;
                        k++;
                    }
                }
            }
            return origins;
        }
        
        
        //getOverflowShifts ends up being called by the display routines quite often
        //so, in the interest of speed, i moved these outside of the function;
        int shiftX, shiftY, shiftZ;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
            shiftX = 0; shiftY = 0; shiftZ = 0;
            r = (Vector)rr;
            
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if(r.z-distance < 0.0) {shiftZ = +1;}
            else if(r.z+distance > dimensions.z) {shiftZ = -1;}
              
            if((shiftX == 0) && (shiftY == 0) && (shiftZ == 0)) {
              shift = shift0;
            } else if((shiftX != 0) && (shiftY == 0) && (shiftZ == 0)) {
              shift = new float[1][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[1][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
            } else if((shiftX == 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[1][D];
              shift[0][2] = (float)(shiftZ*dimensions.z);
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][0] = shift[0][0];
              shift[2][1] = shift[1][1];
            } else if((shiftX != 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][0] = shift[0][0];
              shift[2][2] = shift[1][2];
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[3][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][1] = shift[0][1];
              shift[2][2] = shift[1][2];
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[7][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][2] = (float)(shiftZ*dimensions.z);
              shift[3][0] = shift[0][0];
              shift[3][1] = shift[1][1];
              shift[4][1] = shift[1][1];
              shift[4][2] = shift[2][2];
              shift[5][0] = shift[0][0];
              shift[5][2] = shift[2][2];
              shift[6][0] = shift[0][0];
              shift[6][1] = shift[1][1];
              shift[6][2] = shift[2][2];
            }
            
            return(shift);
        }
    }
    
/*    public static final class BoundarySlidingBrick extends BoundaryPeriodicSquare {
        private double gamma = 0.0;
        private double delvx;
        private IntegratorMD.ChronoMeter timer;
        public BoundarySlidingBrick() {super();}
        public BoundarySlidingBrick(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.SLIDING_BRICK;}
        public void setShearRate(double g) {gamma = g; computeDelvx();}
        public double getShearRate() {return gamma;}
        private void computeDelvx() {delvx = gamma*dimensions.y;}
        
        public void setTimer(IntegratorMD.ChronoMeter t) {timer = t;}
        
        public void nearestImage(Vector dr) {
            double delrx = delvx*timer.currentValue();
            double cory = ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y + 0.5) : Math.ceil(dr.y/dimensions.y - 0.5));
            double corz = ((dr.z > 0.0) ? Math.floor(dr.z/dimensions.z + 0.5) : Math.ceil(dr.z/dimensions.z - 0.5));
            dr.x -= cory*delrx;
            dr.x -= dimensions.x*((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x + 0.5) : Math.ceil(dr.x/dimensions.x - 0.5));
            dr.y -= dimensions.y * cory;
            dr.z -= dimensions.z * corz;
        }
        
        public void centralImage(Vector r) {
            double delrx = delvx*timer.currentValue();
            double cory = ((r.y > 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y - 1.));
            double corz = ((r.z > 0.0) ? Math.floor(r.z/dimensions.z) : Math.ceil(r.z/dimensions.z - 1.));
            r.x -= cory*delrx;
            r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x - 1.0));
            r.y -= dimensions.y * cory;
            r.z -= dimensions.z * corz;
        }
        
        public void centralImage(Coordinate c) {
            Vector r = c.r;
            double cory = ((r.y > 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
            double corx = ((r.x > 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            double corz = ((r.z > 0.0) ? Math.floor(r.z/dimensions.z) : Math.ceil(r.z/dimensions.z-1.0));
            if (corz == 0.0 && (cory == 0.0 && corz == 0.0)) return;
            double delrx = delvx*timer.currentValue();
            Vector p = c.p;
            r.x -= cory*delrx;
            r.x -= dimensions.x * corx;
            r.y -= dimensions.y * cory;
            r.z -= dimensions.z * corz;
            p.x -= corx*delvx;
        }
        
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells + 1)*(2*nShells + 1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for (int i=-nShells; i <= nShells; i++) {
                for (int j=-nShells; j <= nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x + j*delvx*timer.currentValue();
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }
    }//end of BoundarySlidingBrick */

/*    protected final static class BoundaryDeformableCell extends Boundary implements Space.Boundary.Periodic, Space.Boundary.Deformable {

        public BoundaryDeformableCell() {this(Default.BOX_SIZE);}
        public BoundaryDeformableCell(double length) {this(new double[] {length, length, length});}
        public BoundaryDeformableCell(double[] lengths) {
            super();
            dimensions.E(lengths);
            boundaryTensor.diagE(lengths);
            updateCopies();
        }
        public Space.Boundary.Type type() {return Boundary.DEFORMABLE_CELL;}

        private Tensor boundaryTensor = new Tensor();
        private Tensor boundaryTensorCopy = new Tensor();
        
        private final Vector dimensions = new Vector();
        private final Vector dimensionsCopy = new Vector();
        private final Vector dimensionsHalf = new Vector();
        private final Vector workVector = new Vector();
        private final Tensor workTensor = new Tensor();
        
        
        public final Space.Vector dimensions() {return dimensionsCopy;}
        public Space.Tensor boundaryTensor() {return boundaryTensorCopy;}

        public Space.Vector randomPosition() {
            temp.x = Simulation.random.nextDouble();
            temp.y = Simulation.random.nextDouble();
            temp.z = Simulation.random.nextDouble();
            temp.transform(boundaryTensor);
            return temp;
        }
        
        private void updateCopies() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
            boundaryTensorCopy.E(boundaryTensor);
        }

        public void nearestImage(Space.Vector dr) {nearestImage((Vector) dr);}
        public void nearestImage(Vector dr) {
            check(dr);
            if(hey){
            while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
            while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
            while(dr.y > +dimensionsHalf.y) dr.y -= dimensions.y;
            while(dr.y < -dimensionsHalf.y) dr.y += dimensions.y;
            while(dr.z > +dimensionsHalf.z) dr.z -= dimensions.z;
            while(dr.z < -dimensionsHalf.z) dr.z += dimensions.z;}
//            System.out.println(" after " + dr);
        }
        public void centralImage(Coordinate c) {centralImage(c.r);}
        public void centralImage(Space.Vector r) {centralImage((Vector) r);}
        public void centralImage(Vector r) {
            check(r);
            if(hey){
                while(r.x > dimensions.x) r.x -= dimensions.x;
                while(r.x < 0.0)          r.x += dimensions.x;
                while(r.y > dimensions.y) r.y -= dimensions.y;
                while(r.y < 0.0)          r.y += dimensions.y;
                while(r.z > dimensions.y) r.z -= dimensions.z;
                while(r.z < 0.0)          r.z += dimensions.z;
            }
//            System.out.println(" after " + r);
//             System.out.println(" rx, ry, rz "+ (xi*DT.component(0,0)+eta*DT.component(1,0)+zeta*DT.component(2,0))+" "
//                              +(xi*DT.component(0,1)+eta*DT.component(1,1)+zeta*DT.component(2,1))+" "
//                              +(xi*DT.component(0,2)+eta*DT.component(1,2)+zeta*DT.component(2,2))+" ");
                             // +"          "+ "xi eta zeta " + xi+" " +eta+" " +zeta);   
        }

//     public void check(Vector r2){r2.E(r2);hey=false;} //if MCMoveatom is not involed....      
        
// MCMoveatom is involved...
        public void check(Vector r2){ 
//System.out.println("Before x,y,z"+r2);
             A1=r2.x; A2=r2.y; A3=r2.z;
//             System.out.println(" A1 A2 A3 "+ A1+" "+A2+" "+A3);   
//             System.out.println(DT.component(0,1)+" "+DT.component(0,2)+" "+DT.component(1,0)+" "+DT.component(1,2)+" "+DT.component(2,0)+" "+DT.component(2,1));   
             hey = true;
             
         if( DT.component(0,1)!=0.0 || DT.component(0,2)!=0.0 || DT.component(1,0)!=0.0 
               || DT.component(1,2)!=0.0 || DT.component(2,0)!=0.0 || DT.component(2,1)!=0.0 )
          { 
            
         xi = (-DT.component(2,0)*A3*DT.component(1,1)+DT.component(2,0)*DT.component(1,2)*A2 
               -DT.component(2,2)*DT.component(1,0)*A2-DT.component(1,2)*DT.component(2,1)*A1
               +A3*DT.component(2,1)*DT.component(1,0)+ DT.component(2,2)*A1*DT.component(1,1))
             /(-DT.component(2,0)*DT.component(0,2)*DT.component(1,1)+DT.component(2,0)*DT.component(0,1)*DT.component(1,2)
               +DT.component(2,2)*DT.component(0,0)*DT.component(1,1)-DT.component(0,1)*DT.component(2,2)*DT.component(1,0)
               +DT.component(0,2)*DT.component(2,1)*DT.component(1,0)-DT.component(2,1)*DT.component(0,0)*DT.component(1,2)); 

         eta = -(DT.component(0,1)*DT.component(2,2)*A1-DT.component(2,0)*DT.component(0,1)*A3
                +DT.component(2,0)*DT.component(0,2)*A2-DT.component(0,2)*DT.component(2,1)*A1
                -DT.component(2,2)*DT.component(0,0)*A2+DT.component(2,1)*DT.component(0,0)*A3)
             /(-DT.component(2,0)*DT.component(0,2)*DT.component(1,1)+DT.component(2,0)*DT.component(0,1)*DT.component(1,2)
               +DT.component(2,2)*DT.component(0,0)*DT.component(1,1)-DT.component(0,1)*DT.component(2,2)*DT.component(1,0)
               +DT.component(0,2)*DT.component(2,1)*DT.component(1,0)-DT.component(2,1)*DT.component(0,0)*DT.component(1,2)); 
             
         zeta = (DT.component(0,0)*A3*DT.component(1,1)-DT.component(0,0)*DT.component(1,2)*A2
                -DT.component(1,0)*DT.component(0,1)*A3+DT.component(1,0)*DT.component(0,2)*A2
                -A1*DT.component(0,2)*DT.component(1,1)+A1*DT.component(0,1)*DT.component(1,2))
             /(-DT.component(2,0)*DT.component(0,2)*DT.component(1,1)+DT.component(2,0)*DT.component(0,1)*DT.component(1,2)
               +DT.component(2,2)*DT.component(0,0)*DT.component(1,1)-DT.component(0,1)*DT.component(2,2)*DT.component(1,0)
               +DT.component(0,2)*DT.component(2,1)*DT.component(1,0)-DT.component(2,1)*DT.component(0,0)*DT.component(1,2)); 

//             System.out.println("Previous x,y,z"+A1+" "+A2+" "+A3);             
//             System.out.println("Previous xi,eta,zeta"+xi+" "+eta+" "+zeta);

                if(Math.abs(xi)< 1e-10){xi=0;}
                if(Math.abs(eta)< 1e-10){eta=0;}
                if(Math.abs(zeta)< 1e-10){zeta=0;}

//       System.out.println(" here in space3D Previous: " +r2.y+ "  "+xi+" "+DT.component(0,1)+" "+eta+" "+DT.component(1,1)+
//        " "+zeta+" "+DT.component(2,1));

                while(xi<0) {xi+=1;}
                while(xi>1) {xi-=1;}
                while(eta<0) {eta+=1;}
                while(eta>1) {eta-=1;}
                while(zeta<0) {zeta+=1;}
                while(zeta>1) {zeta-=1;}

//       System.out.println(" here in space3D Previous: " +r2.y+ "  "+xi+" "+DT.component(0,1)+" "+eta+" "+DT.component(1,1)+
//        " "+zeta+" "+DT.component(2,1));
         
          r2.x=xi*DT.component(0,0)+eta*DT.component(1,0)+zeta*DT.component(2,0);
          r2.y=xi*DT.component(0,1)+eta*DT.component(1,1)+zeta*DT.component(2,1);
          r2.z=xi*DT.component(0,2)+eta*DT.component(1,2)+zeta*DT.component(2,2);
               hey = false;


          }           
    }//end check
    
    private void updateDimensions() {
        Vector[] columns = boundaryTensor.columns();
        for(i=0; i<D; i++) dimensions.setX(i, Math.sqrt(columns[i].squared()));
    }
    
    public void deform(Tensor deformationTensor) {
        boundaryTensor.TE(deformationTensor);
        updateDimensions();
        updateCopies();
    }

        /**
         * Multiplies all lengths isotropically by the given value.
         * /
        public void inflate(double scale) {
            dimensions.TE(scale);
            boundaryTensor.TE(scale);
            updateCopies();
        }
        
        /**
         * Multiples each dimension (x, y, z) by the corresponding factor in the given vector.
         * Same as deform with a tensor argument formed with the given vector on its diagonal.
         * /
        public void inflate(Space.Vector scale) {
            workTensor.E(0.0);
            workTensor.diagE(scale);
            deform(workTensor);
        }

        /**
         * Sets the length of each boundary edge to the corresponding value in the
         * given vector, keeping the shape of the box unchanged (other than the change in size).
         * /
        public void setDimensions(Space.Vector v) {
            workVector.E(v);
            workVector.DE(dimensions);
            inflate(workVector);
        }
        
        /**
         * Sets the boundary tensor to equal the given tensor.
         * /
        public void setDimensions(Space.Tensor t) {
            boundaryTensor.E(t);
            updateDimensions();
            updateCopies();
        }        
            
        public double volume() {
            Vector[] columns = boundaryTensor.columns();
            return Math.abs(columns[0].dot(columns[1].cross(columns[2])));
        } 
       
                
        /**
         * imageOrigins and getOverFlowShifts are both probably incorrect, if they are
         * even completed.  They should definitely be checked before being implemented.
         * /
        
        int shellFormula, nImages, i, j, k, m;
        double[][] origins;
        public double[][] imageOrigins(int nShells) {
            throw new RuntimeException("imageOrigins not implemented in Space3D.BoundaryDeformableCell");
         /*   shellFormula = (2 * nShells) + 1;
            nImages = shellFormula*shellFormula*shellFormula-1;
            origins = new double[nImages][D];
            for (k=0,i=-nShells; i<=nShells; i++) {
                for (j=-nShells; j<=nShells; j++) {
                    for (m=-nShells; m<=nShells; m++) {
                        if ((i==0 && j==0) && m==0 ) {continue;}
                        origins[k][0] = i*dimensions.x;
                        origins[k][1] = j*dimensions.y;
                        origins[k][2] = m*dimensions.z;
                        k++;
                    }
                }
            }
            return origins;
            * /
        }//end of imageOrigins
        
        
        //getOverflowShifts ends up being called by the display routines quite often
        //so, in the interest of speed, i moved these outside of the function;
        int shiftX, shiftY, shiftZ;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
            throw new RuntimeException("Space3D.BoundaryDeformableCell.getOverflowShifts not implmented");
         /*   shiftX = 0; shiftY = 0; shiftZ = 0;
            r = (Vector)rr;
            
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if(r.z-distance < 0.0) {shiftZ = +1;}
            else if(r.z+distance > dimensions.z) {shiftZ = -1;}
              
            if((shiftX == 0) && (shiftY == 0) && (shiftZ == 0)) {
              shift = shift0;
            } else if((shiftX != 0) && (shiftY == 0) && (shiftZ == 0)) {
              shift = new float[1][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[1][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
            } else if((shiftX == 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[1][D];
              shift[0][2] = (float)(shiftZ*dimensions.z);
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][0] = shift[0][0];
              shift[2][1] = shift[1][1];
            } else if((shiftX != 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][0] = shift[0][0];
              shift[2][2] = shift[1][2];
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[3][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][1] = shift[0][1];
              shift[2][2] = shift[1][2];
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[7][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][2] = (float)(shiftZ*dimensions.z);
              shift[3][0] = shift[0][0];
              shift[3][1] = shift[1][1];
              shift[4][1] = shift[1][1];
              shift[4][2] = shift[2][2];
              shift[5][0] = shift[0][0];
              shift[5][2] = shift[2][2];
              shift[6][0] = shift[0][0];
              shift[6][1] = shift[1][1];
              shift[6][2] = shift[2][2];
            }
            
            return(shift);
            * /
        }//end of getOverflowShifts
    }//end of BoundaryDeformableCell 
    */

    public static void main (String[] args) {
        Vector r1 = new Vector(2,2,3);
        System.out.println("r1_before" + r1.toString());
        Tensor tensor = new Tensor(1,2,0,1,1,2,0,0,1);
        RotationTensor tensor2 = new RotationTensor();
        //r1.transform(tensor2);
        tensor2.E(tensor);
        System.out.println("tensor2_before " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        System.out.println();
        
        r1.transform(tensor2);
        System.out.println("r1_transform(tensor2)" + r1.toString());
        tensor2.invert();
        System.out.println("tensor2_invert " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        tensor2.setAxial(1, 2*Math.PI);
        System.out.println("tensor2_rotate_360 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        System.out.println();
        
        r1.transform(tensor2);
        System.out.println("r1_afterInvert_andRotate360 " + r1.toString());
        //System.out.println("tensor2 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
    }
    
}//end of Space3D