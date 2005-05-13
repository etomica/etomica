package etomica.lattice.crystal;
import etomica.Space;
import etomica.lattice.Primitive;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Polytope;
import etomica.space.Vector;

/**
 * Primitive group for a tetragonal system.  All primitive
 * vectors orthogonal and two are of equal length.
 */
public class PrimitiveTetragonal extends Primitive implements Primitive3D {
    
    private double ab = 1.0, c = 1.0;
    
    public PrimitiveTetragonal(Space space) {
        this(space, 1.0, 1.0);
    }
    public PrimitiveTetragonal(Space space, double ab, double c) {
        super(space);//also makes reciprocal
        setAB(ab); //also sets reciprocal via update
        setC(c);
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveTetragonal(Space space, Primitive direct) {
        super(space, direct);
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveTetragonal(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveTetragonal)reciprocal()).setAB(2.0*Math.PI/ab);
        ((PrimitiveTetragonal)reciprocal()).setC(2.0*Math.PI/c);
    }
    
    public void setA(double a) {setAB(a);}
    public double getA() {return ab;}
    
    public void setB(double b) {setAB(b);}
    public double getB() {return ab;}
    
    public void setAB(double ab) {
        if(immutable || ab <= 0.0) return;
        this.ab = ab;
        size[0] = size[1] = ab;
        latticeVectors[0].setX(0,ab);
        latticeVectors[1].setX(1,ab);
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
    public double getGamma() {return rightAngle;}
    
    
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
        return new PrimitiveTetragonal(space, ab, c);
    }
    
    public void scaleSize(double scale) {
        setAB(ab*scale);
        setC(c*scale);
    }

    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
   */ }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
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
    
    /**
     * Returns a new Cuboid with edges given by the lengths of the
     * primitive vectors.
     */
    public Polytope unitCell() {
        return new Cuboid(space, ab, ab, c);
    }
    
    public String toString() {return "Tetragonal";}
    
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
    
