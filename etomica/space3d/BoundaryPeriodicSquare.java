package etomica.space3d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class BoundaryPeriodicSquare extends Boundary implements Boundary.Periodic  {
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx, double ly, double lz) {super(p);dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
        public BoundaryPeriodicSquare(double lx, double ly, double lz) {super();dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
        public etomica.space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        private final Vector temp = new Vector();
        private final Vector modShift = new Vector();//must be used only by centralImage and nearestImage methods
        protected final Vector dimensions = new Vector();
        protected final Vector dimensionsCopy = new Vector();
        protected final Vector dimensionsHalf = new Vector();
        public final etomica.space.Vector dimensions() {return dimensionsCopy;}
        public etomica.space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble();
            temp.y = dimensions.y*Simulation.random.nextDouble();
            temp.z = dimensions.z*Simulation.random.nextDouble();
            return temp;
        }
        private final void updateDimensions() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
        }
        public void nearestImage(etomica.space.Vector dr) {nearestImage((Vector) dr);}
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
        //Converts dr to its nearest-image value and returns in shift
        //the difference between the new and old values of dr
		public void nearestImage(Vector dr, Vector shift) {
			shift.EMod2Shift(dr, dimensionsHalf);
			if(!shift.isZero()) dr.PE(shift);
		}

        public boolean centralImage(Coordinate c) {
        	modShift.EModShift((Vector)c.position(), dimensions);
//			Vector r = (Space3D.Vector)c.position();
//        	modShift.E(r);
//			modShift.mod(dimensions);
//			modShift.ME(r);  //shift = (r mod dimensions) - r
        	if(modShift.isZero()) return false;
        	c.translateBy(modShift);
        	return true;
//        	return centralImage(c.r);
        }
 /*		public boolean centralImage(Coordinate c) {
 			Vector displace = centralImage(c.position());
 			c.translateBy(displace);
 		}
 */
        public boolean centralImage(etomica.space.Vector v) {return centralImage((Vector) v);}
        public boolean centralImage(Vector v) {
            temp.E(v);
            v.mod(dimensions);
            return temp.equals(v);
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
        public void inflate(etomica.space.Vector scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void setDimensions(etomica.space.Vector v) {dimensions.E(v); updateDimensions();}
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
            origins = new double[nImages][3];
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
        public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {
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
              shift = new float[1][3];
              shift[0][0] = (float)(shiftX*dimensions.x);
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[1][3];
              shift[0][1] = (float)(shiftY*dimensions.y);
            } else if((shiftX == 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[1][3];
              shift[0][2] = (float)(shiftZ*dimensions.z);
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ == 0)) {
              shift = new float[3][3];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][0] = shift[0][0];
              shift[2][1] = shift[1][1];
            } else if((shiftX != 0) && (shiftY == 0) && (shiftZ != 0)) {
              shift = new float[3][3];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][0] = shift[0][0];
              shift[2][2] = shift[1][2];
            } else if((shiftX == 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[3][3];
              shift[0][1] = (float)(shiftY*dimensions.y);
              shift[1][2] = (float)(shiftZ*dimensions.z);
              shift[2][1] = shift[0][1];
              shift[2][2] = shift[1][2];
            } else if((shiftX != 0) && (shiftY != 0) && (shiftZ != 0)) {
              shift = new float[7][3];
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
