package etomica.space3d;

import etomica.NearestImageTransformer;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class CoordinatePair extends etomica.space.CoordinatePair {
        Coordinate c1;
        Coordinate c2;
        private final Vector dr = new Vector();
        private double dvx, dvy, dvz; //drx, dry, drz;
        private NearestImageTransformer nearestImageTransformer = Boundary;

		public void setNearestImageTransformer(NearestImageTransformer b) {this.nearestImageTransformer = b;}
		public NearestImageTransformer getNearestImageTransformer() {return nearestImageTransformer;}
        public void trueReset(Coordinate coord1, Coordinate coord2, double falseTime) {
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            trueReset(falseTime);
        }
        public void trueReset(double falseTime) {
            resetV();
            dr.Ev1Mv2(c2.r,c1.r);

            dr.x += falseTime * dvx;
            dr.y += falseTime * dvy;
            dr.z += falseTime * dvz;
            nearestImageTransformer.nearestImage(dr);
        }

        public void reset(Coordinate coord1, Coordinate coord2) {
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.Ev1Mv2(c2.r,c1.r);
         //   c2.position(); c1.position();
         //   dr.x = c2.r.x - c1.r.x;
         //   dr.y = c2.r.y - c1.r.y;
         //   dr.z = c2.r.z - c1.r.z;
         //   c1.atom.node.parentPhase().boundary().nearestImage(dr);
            nearestImageTransformer.nearestImage(dr);
   //         drx = dr.x;
   //         dry = dr.y;
   //         drz = dr.z;
        }
        
        public void resetV() {
  //       /*    comment here if not doing hard dynamics
            double rm1 = c1.rm();
            double rm2 = c2.rm();
            dvx = rm2*c2.p.x - rm1*c1.p.x;
            dvy = rm2*c2.p.y - rm1*c1.p.y;
            dvz = rm2*c2.p.z - rm1*c1.p.z;
          //  /* */  //end of non-dynamics commenting
        }
            
        public double r2() {
            return dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
         /*   dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
            dr.z = c2.r.z - c1.r.z;
            c1.atom.node.parentPhase().boundary().nearestImage(dr);
            return dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;*/
        ///    dr.Ev1Mv2(c2.r, c1.r);
        ///    c1.atom.node.parentPhase().boundary().nearestImage(dr);
        ///    return dr.squared();
        }
            
        public Vector dr() {return dr;}
        public double dr(int i) {return (i==0) ? dr.x : ((i==1) ? dr.y : dr.z);}
        public double dv(int i) {return (i==0) ? dvx : ((i==1) ? dvy : dvz);}
        public double v2() {return dvx*dvx + dvy*dvy + dvz*dvz;}
        public double vDot(etomica.space.Vector u) {return vDot((Vector)u);}
        public double vDot(Vector u) {return dvx*u.x + dvy*u.y + dvz*u.z;}
        public double vDotr() {return dr.x*dvx + dr.y*dvy + dr.z*dvz;}
        public void push(double impulse) {
            c1.p.x += impulse*dr.x;
            c1.p.y += impulse*dr.y;
            c1.p.z += impulse*dr.z;
            c2.p.x -= impulse*dr.x;
            c2.p.y -= impulse*dr.y;
            c2.p.z -= impulse*dr.z;
        }
        public void truePush(Vector u, double falseTime) {
            c1.p.PE(u);
            c2.p.ME(u);
            
            c1.r.PEa1Tv1(-falseTime*c1.rm(),u);
            c2.r.PEa1Tv1(falseTime*c2.rm(),u);
        }
        public void nudge(double rDelta) {
            double ratio = c2.mass()*c1.rm()*rDelta;
            c1.r.x -= ratio*dr.x;
            c1.r.y -= ratio*dr.y;
            c1.r.z -= ratio*dr.z;
            c2.r.x += ratio*dr.x;
            c2.r.y += ratio*dr.y;
            c2.r.z += ratio*dr.z;
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
