package etomica.space3d;

import etomica.Atom;
import etomica.Phase;
import etomica.statmech.MaxwellBoltzmann;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Coordinate extends etomica.space.Coordinate {
        public final Vector r = new Vector();
        public final Vector p = new Vector();
        public final Vector rLast = new Vector();  //vector for saving position
        public final Vector work = new Vector();
        public Coordinate(Atom a) {super(a);}
                        
        public void transform(etomica.space.Vector r0, etomica.space.Tensor A) {
            r.transform((Boundary)atom.node.parentPhase().boundary(),(Vector)r0, (Tensor)A);
        }
        public etomica.space.Vector position() {return r;}
        public etomica.space.Vector truePosition(double falseTime) {
            work.E(r);
            work.PEa1Tv1(falseTime*rm(),p);
            return work;
        }
        public etomica.space.Vector momentum() {return p;}
        public double position(int i) {return r.x(i);}
        public double truePosition(int i, double falseTime) {return r.x(i)+falseTime*rm()*p.x(i);}
        public double momentum(int i) {return p.x(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {
            double tM = t*rm(); // t/mass
            r.x += p.x*tM;
            r.y += p.y*tM;
            r.z += p.z*tM;
        }
        /**
         * Scales positions of atoms by multiplying by given value.  Does not notify sequencers.
         */
        public void inflate(double s) {r.x *= s; r.y *= s; r.z *= s;}
        /**
         * Scales positions of atoms isotropically by multiplying by coordinate in each
         * direction by the corresponding element of the given vector.  Does not notify sequencers.
         */
        public void inflate(etomica.space.Vector s) {Vector u = (Vector)s; r.x *= u.x; r.y *= u.y; r.z *= u.z;}
        
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(etomica.space.Vector u) {translateBy((Vector)u);}
//            r.PE((Vector)u);
//        }
        //for use by other methods/classes in Space3D
        //(not sure if this is used, or the Space.Vector version)
        protected void translateBy(Vector u) {
            r.PE(u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(double d, etomica.space.Vector u) {
            r.PEa1Tv1(d,u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateTo(etomica.space.Vector u) {
            r.E((Vector)u);
        }      
        public void replace() {r.E(rLast);}
        
        public void displaceBy(etomica.space.Vector u) {rLast.E(r); translateBy((Vector)u);}
        public void displaceBy(double d, etomica.space.Vector u) {rLast.E(r); translateBy(d,u);}
        public void displaceTo(etomica.space.Vector u) {rLast.E(r); translateTo(u);}  
        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        public void displaceToRandom(etomica.Phase phase) {rLast.E(r); translateToRandom(phase);}
    //    public final void inflate(double s) {r.TE(s);}

        public void accelerateBy(etomica.space.Vector u) {p.PE(u);}
        public void accelerateBy(double d, etomica.space.Vector u) {p.PEa1Tv1(d,u);}
        public void accelerateTo(etomica.space.Vector u) {p.E(u);}
        public void trueAccelerateTo(etomica.space.Vector u, double falseTime) {
            double tm = falseTime * rm();
            r.x -= tm * (((Vector)u).x - p.x);
            r.y -= tm * (((Vector)u).y - p.y);
            r.z -= tm * (((Vector)u).z - p.z);
            p.x = ((Vector)u).x;
            p.y = ((Vector)u).x;
            p.z = ((Vector)u).z;
        }

        public void randomizeMomentum(double temperature) {
            if(isStationary()) {p.E(0.0); return;}
            p.setX(0,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
            p.setX(1,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
            p.setX(2,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
        }
    }
