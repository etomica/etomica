package etomica.performance;

import etomica.*;



public class SpaceP extends Space{

    public static String version() {return "SpaceP:01.12.14/"+Space.VERSION;}

    public static int D=3;   

    public SpaceP(int dimension){
        D=dimension;
    }

    public final int D() {return D;}

    public static final Vector ORIGIN = new Vector();

    public final Space.Vector origin() {return ORIGIN;}

    

    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r

    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)

    public Space.Vector makeVector() {return new Vector();}

    public Space.Orientation makeOrientation() {return null;}

    public Space.Tensor makeTensor() {return null;}

    public Space.Tensor makeRotationTensor() {return null;}

    public Space.Coordinate makeCoordinate(Atom a) {
        if(a instanceof AtomGroup) return new CoordinateGroup((AtomGroup)a);

        //else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);

        else {return new Coordinate(a);}
    }

    public Space.CoordinatePair makeCoordinatePair() {return null;}

    public Space.Boundary.Type[] boundaryTypes() {return null;}

    public Space.Boundary makeBoundary() {return makeBoundary(Space3D.Boundary.PERIODIC_SQUARE);}  //default

    public Space.Boundary makeBoundary(Space.Boundary.Type t) {

    /*

        if(t == Boundary.NONE) {return new BoundaryNone();}

        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}

        else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();

        else return null;

        */

        return new BoundaryPeriodicSquare();

    }

    

    public static EtomicaInfo getEtomicaInfo() {

        EtomicaInfo info = new EtomicaInfo("Two-dimensional space");

        return info;

    }



    public static final double r2(Vector u1, Vector u2, Boundary b) {

    System.out.println("checking in r2 ");
       /*

        Vector.WORK.x = u1.x - u2.x;

        Vector.WORK.y = u1.y - u2.y;

        b.nearestImage(Vector.WORK);

        */

        //return Vector.WORK.x*Vector.WORK.x + Vector.WORK.y*Vector.WORK.y;

        return 0.0;

    }

    public final static class Vector extends Space.Vector{  //declared final for efficient method calls
        public double X[];

        public static final Vector WORK = new Vector();

        public Vector () {X = new double[D];
            for(int i=0;i<D;i++) X[i]=0.0;
        }

        public Vector (double x, double y, double z) {X[0] = x; X[1]= y; X[2]=z;}

        public Vector (double[] a) {for(int i=0;i<D;i++)X[i]=a[i];}//should check length of a for exception

        public Vector (Vector u) {this.E(u);}

        public double[] toArray() {return X;}

    	public void setArray(double[] xNew) {X = xNew;}//should check length

        //NEED TO CHANGE

        public void sphericalCoordinates(double[] result) {

        }

        public int length() {return D;}//bad name for this

        public int D() {return D;}

        public double component(int i) {return X[i];}

        public void setComponent(int i, double d) {X[i]=d;}

        public void E(Vector u) {for(int i=0;i<D;i++)X[i]=u.X[i];}
       
        public void E(double[] u) {for(int i=D;i>0;i--)X[i]=u[i];}  //should check length of array for exception

//        public void E(double[] u) {X=u;}  //should check length of array for exception

        public void E(double a) {for(int i=0;i<D;i++)X[i]=a;}


        public void E(Space.Vector u) {
             if (u instanceof Space3D.Vector){
         X=u.toArray();
         }
        else 
                E((Vector)u);
        }

        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u;for(int i=0;i<D;i++)X[i]=a1*u1.X[i];}

        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; for(int i=0;i<D;i++)X[i]+=a1*u1.X[i];}

        public void PE(Vector u) {for(int i=0;i<D;i++)X[i]+=u.X[i];}

        public void PE(double a) {for(int i=0;i<D;i++)X[i]+=a;}

        public void ME(Vector u) {for(int i=0;i<D;i++)X[i]-=u.X[i];}

        public void PE(int i, double a) {X[i]+=a;}

        public void TE(double a) {for(int i=0;i<D;i++)X[i]*=a;}

        public void TE(Vector u) {for(int i=0;i<D;i++)X[i]*=u.X[i];}

        public void TE(int i, double a) {X[i]*=a;}

        public void DE(double a) {for(int i=0;i<D;i++)X[i]/=a;}

        public void DE(Vector u) {for(int i=0;i<D;i++)X[i]/=u.X[i];}

        public double min() {double m=Double.MAX_VALUE;for(int i=0;i<D;i++)if(X[i]<m) m=X[i]; return m;}

        public double max() {double m=Double.MAX_VALUE;for(int i=0;i<D;i++)if(X[i]>m) m=X[i]; return m;}

        public double squared() {double sum=0.0; for(int i=0;i<D;i++)sum+=X[i]*X[i]; return sum;}

        public double dot(Vector u) {double sum=0.0;for(int i=0;i<D;i++)sum+=X[i]*u.X[i]; return sum;}

        public Space.Vector P(Space.Vector u) {return WORK;}//NEED TO CHANGE

        public Space.Vector M(Space.Vector u) { return WORK;}//NEED TO CHANGE

        public Space.Vector T(Space.Vector u) { return WORK;}//NEED TO CHANGE

        public Space.Vector D(Space.Vector u) {return WORK;} //NEED TO CHANGE

        public Space.Vector abs() { return WORK;} // NEED TO CHANGE

       public void setRandomCube(){for(int i=0;i<D;i++)X[i]=Simulation.random.nextDouble()-0.5;}
        /**
         * Replaces this vector with its cross-product with the given 3D vector, with result projected
         * onto the 2D plane.  This vector becomes the result of (this vector) X u.
         */
        public void XE(Space3D.Vector u) {
        }

        //NEED TO CHANGE

         public Space3D.Vector cross(Space3D.Vector u) {//not thread safe

            /*

            Space3D.Vector.WORK.x = y*u.z;

            Space3D.Vector.WORK.y = -x*u.z;

            Space3D.Vector.WORK.z = x*u.y - y*u.x;

            */

            return Space3D.Vector.WORK;

        }

        //NEED TO CHANGE

        public Space3D.Vector cross(Space2D.Vector u) {//not thread safe

            /*

            Space3D.Vector.WORK.x = 0.0;

            Space3D.Vector.WORK.y = 0.0;

            Space3D.Vector.WORK.z = x*u.y - y*u.x;

            */

            return Space3D.Vector.WORK;

        }

         

         

        public void normalize() {

            

        }

        public void transform(Space.Tensor A) {transform((Tensor)A);}

        public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {}

        public void transform(Boundary b, Vector r0, Tensor A) {

            

        }

        public void randomStep(double d) {} //uniformly distributed random step in x and y, within +/- d

        public void setRandom(double d){ }

        public void setRandom(double dx, double dy) {}

        public void setRandom(Vector u) {}

        public void setRandomSphere() {

         }

        public void randomRotate(double thetaStep){

         }

        

       // public void setArray(double [] a){for(int i=0;i<D;i++)X[i]=a[i];}

        public double[] getArray(){return X;}

        

        
        public void PE(Space.Vector u) {PE((Vector)u);}

        public void ME(Space.Vector u) {ME((Vector)u);}

        public void TE(Space.Vector u) {TE((Vector)u);}

        public void DE(Space.Vector u) {DE((Vector)u);}

        public double dot(Space.Vector u) {return dot((Vector)u);}
        
      
        

    }

        

   public static class Coordinate extends Space.Coordinate {

        public Coordinate nextCoordinate, previousCoordinate;

        public final Vector r = new Vector();  //Cartesian coordinates

        public final Vector p = new Vector();  //Momentum vector

        public final Vector rLast = new Vector();  //vector for saving position

        protected final Vector work = new Vector(); 



        public Coordinate(Atom a) {super(a);}

        public Space.Vector position() {return r;}

        public Space.Vector momentum() {return p;}

        public double position(int i) {return r.component(i);}

        public double momentum(int i) {return p.component(i);}

        

        public void setNextAtom(Atom a) {

            if(a == null) nextCoordinate = null;

            else {

                nextCoordinate = (Coordinate)a.coord;

                ((Coordinate)a.coord).previousCoordinate = this;

            }

        }

        

        public Atom nextAtom() {return nextCoordinate!=null ? nextCoordinate.atom : null;}

        public Atom previousAtom() {return previousCoordinate!=null ? previousCoordinate.atom : null;}

        public void clearPreviousAtom() {previousCoordinate = null;}

        

        

        

        public void randomizeMomentum(double temperature) { 

        }

        public double kineticEnergy() {return 0.5*p.squared()*rm();}

       //NEED TO CHANGE

        public void freeFlight(double t) {

            /*double tM = t*rm(); // t/mass

            for(int i=r.length;i>0;i--){

            r.X[i] += p.X[i]*tM;

            }

            */

       }

      

        public void transform(Space.Vector r0, Space.Tensor A) {

            r.transform((Boundary)atom.parentPhase().boundary(), (Vector)r0, (Tensor)A);}

        /**

        * Moves the atom by some vector distance

        * 

        * @param u

        */

        public void translateBy(Space.Vector u) {r.PE(u);}

        /**

        * Moves the atom by some vector distance

        * 

        * @param u

        */

        public void translateBy(double d, Space.Vector u) {r.PEa1Tv1(d,(Vector)u);}

        /**

        * Moves the atom by some vector distance

        * 

        * @param u

        */

        public void translateTo(Space.Vector u) {r.E(u);}      

        public void displaceBy(Space.Vector u) {rLast.E(r); translateBy((Vector)u);}

        public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,(Vector)u);}

        public void displaceTo(Space.Vector u) {rLast.E(r); translateTo((Vector)u);}  

        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}

        public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}

        public void replace() {r.E(rLast);}

    //    public final void inflate(double s) {r.TE(s);}



        public void accelerateBy(Space.Vector u) {p.PE(u);}

        public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}



       

        

    }

    

    

     public static class CoordinateGroup extends Coordinate implements Space.CoordinateGroup {
        public Coordinate firstChild, lastChild;
        public CoordinateGroup(AtomGroup a) {super(a);}

        public final Atom firstAtom() {return (firstChild != null) ? firstChild.atom : null;}
        public final void setFirstAtom(Atom a) {firstChild = (a != null) ? (Coordinate)a.coord : null;}
        public final Atom lastAtom() {return (lastChild != null) ? lastChild.atom : null;}
        public final void setLastAtom(Atom a) {lastChild = (a != null) ? (Coordinate)a.coord : null;}
                
        public double mass() {
            double massSum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                massSum += coord.mass();
            }
            return massSum;
        } 
        public double rm() {return 1.0/mass();}
        /**
         * Applies transformation to COM of group, keeping all internal atoms at same relative positions.
         */
        public void transform(Space.Vector r0, Space.Tensor A) {
            work.E(position()); //work = r
            work.transform((Boundary)atom.parentPhase().boundary(),(Vector)r0, (Tensor)A);
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
            switch(((AtomGroup)atom).childAtomCount()) {
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

    
     //protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic  {
    public static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic  {
        private final Vector temp = new Vector();
//        private static Space.Tensor zilch = new Tensor();
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx, double ly, double lz) {super(p);dimensions.X[0]=lx; dimensions.X[1]=ly; dimensions.X[2]=lz;}
        public BoundaryPeriodicSquare(double lx, double ly, double lz) {super();dimensions.X[0]=lx; dimensions.X[1]=ly; dimensions.X[2]=lz;}
        public Space.Boundary.Type type() {return Space3D.Boundary.PERIODIC_SQUARE;}
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        
        public Space.Vector randomPosition() {
            temp.X[0] = dimensions.X[0]*Simulation.random.nextDouble();
            temp.X[1] = dimensions.X[1]*Simulation.random.nextDouble();
            temp.X[2] = dimensions.X[2]*Simulation.random.nextDouble();
            return temp;
        }
        
        
        public void nearestImage(Space.Vector dr) {nearestImage((Vector) dr);}
        public void nearestImage(Vector dr) {
            System.out.println(" Inside SpaceP.nearestImage");
            /*
            dr.x -= dimensions.x*((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x + 0.5) : Math.ceil(dr.x/dimensions.x - 0.5));
            dr.y -= dimensions.y*((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y + 0.5) : Math.ceil(dr.y/dimensions.y - 0.5));
            dr.z -= dimensions.z*((dr.z > 0.0) ? Math.floor(dr.z/dimensions.z + 0.5) : Math.ceil(dr.z/dimensions.z - 0.5));
            
            */
        }
        
        public void centralImage(Coordinate c) {centralImage(c.r);}
       /* 
        public void centralImage(Coordinate c) {            
        
            Vector r = (Vector)c.position();
            temp.x = (r.x > dimensions.x) ? -dimensions.x : (r.x < 0.0) ? dimensions.x : 0.0;
            temp.y = (r.y > dimensions.y) ? -dimensions.y : (r.y < 0.0) ? dimensions.y : 0.0;
            temp.z =  (r.z >dimensions.z) ? -dimensions.z : (r.z < 0.0) ? dimensions.z : 0.0;
            if(temp.x != 0.0 || temp.y != 0.0 || temp.z!=0.0) c.translateBy(temp);
       
        }

        */
        
        public void centralImage(Space.Vector r) {centralImage((Vector) r);}
        public void centralImage(Vector r) {
            r.X[0] -= dimensions.X[0]*((r.X[0]>0) ? Math.floor(r.X[0]/dimensions.X[0]) : Math.ceil(r.X[0]/dimensions.X[0] - 1.0));
            r.X[1] -= dimensions.X[1]*((r.X[1]>0) ? Math.floor(r.X[1]/dimensions.X[1]) : Math.ceil(r.X[1]/dimensions.X[1] - 1.0));
            r.X[2] -= dimensions.X[2]*((r.X[2]>0) ? Math.floor(r.X[2]/dimensions.X[2]) : Math.ceil(r.X[2]/dimensions.X[2] - 1.0));
        }
        public void inflate(double scale) {dimensions.TE(scale);}
        public double volume() {return dimensions.X[0]*dimensions.X[1]*dimensions.X[2];}
                
        /**
         * imageOrigins and getOverFlowShifts are both probably incorrect, if they are
         * even completed.  They should definitely be checked before being implemented.
         */
        
        int shellFormula, nImages, i, j, k, m;
        double[][] origins;
        public double[][] imageOrigins(int nShells) {
            /*
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
            */
           
           System.out.println("Inside SpaceP.imageOrigins");
           return null;
            //return origins;
        }
        
        //getOverflowShifts ends up being called by the display routines quite often
        //so, in the interest of speed, i moved these outside of the function;
        int shiftX, shiftY, shiftZ;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
         
            System.out.println(" inside SpaceP.getOverflo..");
            return null;
        }
        
      
    }
    
    

    

    

    

    

    

    

    

    

    

}

