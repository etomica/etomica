package etomica.lattice;
import etomica.Space;
import etomica.math.geometry.Polytope;

/**
 * Primitive group for a hexagonal system.  Primitive-vector angles
 * are (90,90,120) degrees and two vectors are of equal length.
 */
public class PrimitiveHexagonal extends Primitive implements Primitive3D {
    
    private double ab = 1.0, c = 1.0;
    private int ix = 0, iy = 1;
    private final double gamma = etomica.units.Degree.UNIT.toSim(120.);
//    private final double gamma = etomica.units.Degree.UNIT.toSim(60.);
    private final double cosGamma = Math.cos(gamma);
    private final double sinGamma = Math.sin(gamma);
    
    public PrimitiveHexagonal(Space space) {
        this(space, 1.0, 1.0);
    }
    public PrimitiveHexagonal(Space space, double ab, double c) {
        super(space);
        setAB(ab);
        setC(c);
    }
    
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveHexagonal(Space space, Primitive direct) {
        super(space, direct);
        ix = 1;
        iy = 0;
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveHexagonal(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveHexagonal)reciprocal()).setAB(2.0*Math.PI/(ab*sinGamma));
        ((PrimitiveHexagonal)reciprocal()).setC(2.0*Math.PI/c);
    }
    
    public void setA(double a) {setAB(a);}
    public double getA() {return ab;}
    
    public void setB(double b) {setAB(b);}
    public double getB() {return ab;}
    
    public void setAB(double ab) {
        if(immutable || ab <= 0.0) return;
        this.ab = ab;
        size[0] = size[1] = ab;
        
        //direct lattice (ix = 0, iy = 1)
        // v[0] = (1,0,0); v[1] = (s,c,0); v[2] = (0,0,1)  (times ab, c)
         
        //reciprocal lattice (ix = 1, iy = 0)
        // v[0] = (s,-c,0); v[1] = (0,1,0); v[2] = (0,0,1);  (times ab, c)
        latticeVectors[ix].setX(ix,ab);
        latticeVectors[iy].setX(ix,((ix==0)?+1:-1)*ab*cosGamma);
        latticeVectors[iy].setX(iy,ab*sinGamma);
        update();
    }
    
    public void setC(double c) {
        if(immutable || c <= 0.0) return;
        this.c = c;
        size[2] = c;
        latticeVectors[2].setX(2, c);
        update();
    }
    public double getC() {return c;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {}
    public double getBeta() {return rightAngle;}
    
    public void setGamma(double t) {}
    public double getGamma() {return gamma;}
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return false;}
    public boolean isEditableC() {return true;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return false;}
    public boolean isEditableGamma() {return false;}

    /**
     * Returns a new PrimitiveTetragonal with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveHexagonal(space, ab, c);
    }
    
    public void scaleSize(double scale) {
        setAB(ab*scale);
        setC(c*scale);
    }

    public int[] latticeIndex(Space.Vector q) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
   */ }
    
    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    */}
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method wignerSeitzCell not yet implemented");
    }
    
    public Polytope unitCell() {
        throw new RuntimeException("method unitCellFactory not yet implemented");
    //    return new UnitCellFactory(simulation);
    }
    
    public String toString() {return "Hexagonal";}
    
    /**
     * Main method to demonstrate use and to aid debugging
     * /
    public static void main(String[] args) {
        System.out.println("main method for PrimitiveCubic");
        Space space = new Space2D();
        Simulation sim = new Simulation(space);
        int D = space.D();
        PrimitiveCubic primitive = new PrimitiveCubic(sim);
        AtomFactory basis = primitive.unitCellFactory();
        final int nx = 4;
        final int ny = 5;
        BravaisLattice lattice = BravaisLattice.makeLattice(sim, 
                                basis, 
                                new int[] {nx,ny},
                                primitive);
        lattice.shiftFirstToOrigin();
        System.out.println("Total number of sites: "+lattice.siteList().size());
        System.out.println();
        System.out.println("Coordinate printout");
        AtomAction printSites = new AtomAction() {
            public void actionPerformed(Atom s) {
                System.out.print(s.coord.position().toString()+" ");
                if(((Site)s).latticeCoordinate()[1]==ny-1) System.out.println();
            }
        };
        AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
        iterator.allAtoms(printSites);
        System.out.println();
                
        AbstractCell testSite = (AbstractCell)lattice.site(new int[] {1,2});
        
        Space.Vector vector = space.makeVector(new double[] {1.5, 2.7});
        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
        int[] idx = primitive.latticeIndex(vector);
        System.out.print("; is in: ");
        for(int i=0; i<D; i++) System.out.print(idx[i]);
        System.out.println();
 
        vector = space.makeVector(new double[] {3.5, 5.1});
        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
        idx = primitive.latticeIndex(vector);
        System.out.print("; is in: ");
        for(int i=0; i<D; i++) System.out.print(idx[i]);
        System.out.println();
        
        System.out.println("cell volume: "+testSite.volume());
        System.out.println();
        System.out.println("setting size to 2.0");
        primitive.setSize(2.0);
        lattice.update();//recompute all lattice positions
        lattice.shiftFirstToOrigin();
        
        iterator.allAtoms(printSites);
        System.out.println();
        vector = space.makeVector(new double[] {1.5, 2.7});
        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
        idx = primitive.latticeIndex(vector);
        System.out.print("; is in: ");
        for(int i=0; i<D; i++) System.out.print(idx[i]);
        System.out.println();

        vector = space.makeVector(new double[] {3.5, 5.1});
        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
        idx = primitive.latticeIndex(vector);
        System.out.print("; is in: ");
        for(int i=0; i<D; i++) System.out.print(idx[i]);
        System.out.println();

        System.out.println("cell volume: "+testSite.volume());
  
        /*
        //write out vertices of some cells
        System.out.println();
        System.out.println("Unit-cell vertices of first cell");
        Space.Vector[] vertices = ((AbstractCell)lattice.siteList().getFirst()).vertex();
        for(int i=0; i<vertices.length; i++) {
            System.out.println(vertices[i].toString());
        }
        System.out.println();
        System.out.println("Unit-cell vertices of last cell");
        vertices = ((AbstractCell)lattice.siteList().getLast()).vertex();
        for(int i=0; i<vertices.length; i++) {
            System.out.println(vertices[i].toString());
        }
        * /
    }//end of main
  // */  
}//end of PrimitiveCubic
    
