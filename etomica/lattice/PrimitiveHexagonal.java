package etomica.lattice;
import etomica.*;

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
    
    public PrimitiveHexagonal(Simulation sim) {
        this(sim, 1.0, 1.0);
    }
    public PrimitiveHexagonal(Simulation sim, double ab, double c) {
        super(sim);
        setAB(ab);
        setC(c);
    }
    
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveHexagonal(Simulation sim, Primitive direct) {
        super(sim, direct);
        ix = 1;
        iy = 0;
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveHexagonal(simulation, this);
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
        return new PrimitiveHexagonal(simulation, ab, c);
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
    
    public AtomFactory wignerSeitzCellFactory() {
        throw new RuntimeException("method wignerSeitzCell not yet implemented");
    }
    
    public AtomFactory unitCellFactory() {
        throw new RuntimeException("method unitCellFactory not yet implemented");
    //    return new UnitCellFactory(simulation);
    }
    
    public String toString() {return "Hexagonal";}
    
///////////////////////////////////////////////////////////////////////////////////////////

/*
public class UnitCellFactory extends AtomFactory {

    AtomType atomType;
    
    public UnitCellFactory(Simulation sim) {
        super(sim);
        setType(new AtomType(this));//default
    }
    
    public boolean isGroupFactory() {return false;}
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}

    /**
     * Builds a single unit cell.
     * /
    protected Atom build(AtomTreeNodeGroup parent) {
        return new UnitCell(space, atomType, parent);
    }
    public Atom build(Atom atom) {
        if(!(atom instanceof UnitCell)) throw new IllegalArgumentException("PrimitiveCubic.UnitCellFactory.build(Atom) attempted using an atom that is not an instance of UnitCell");
        return atom;
    }
}//end of UnitCellFactory

///////////////////////////////////////////////////////////////////////////////////////////

/**
 * A cubic unit cell.  Position of the cell is given by the vertex
 * in which each coordinate is minimized.
 */
/*public class UnitCell extends AbstractCell {
    
    private final Space.Vector delta;
    
    public UnitCell(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, parent);
        delta = space.makeVector();
    }
    /**
     * Dimension of the space occupied by the cell
     * /
     public int D() {return space.D();}
     
    /**
     * Returns the volume of the cubic cell.
     * /
    public double volume() {
        return space.powerD(size);
    }
    
    /**
     * Makes vertices for the unit cubic cell, positioned relative to the origin.
     * Must be scaled and translated to cell coordinate position to get true vertex locations.
     */
   /* private Space.Vector[] makeUnitVertices() {
        Space.Vector[] vertices = new Space.Vector[space.powerD(2)];//number of vertices is 2^D
        for(int i=0; i<vertices.length; i++) {
            vertices[i] = space.makeVector();
            int mask = 1;
            for(int j=0; j<D; j++) {
                //the bits of i indicate whether corresponding primitive vector is added
                //do bitwise "and" with mask (which has exactly 1 nonzero bit) to see if
                //each bit is 1 or 0
                if((i & mask) != 0) vertices[i].PE(latticeVectors[j]);//unit lattice vector
                mask *= 2;
            }
        }
        return vertices;
    }//end of vertex
    */
    
    /**
     * Returns the absolute positions of the vertices .
     * Note that vertices might be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values if using them often, 
     * but if doing so be careful to update them if any transformations are done to the lattice.
     * /
    public Space.Vector[] vertex() {
        Space.Vector[] vertices = new Space.Vector[space.powerD(2)];//number of vertices is 2^D
        for(int i=0; i<vertices.length; i++) {
            vertices[i] = space.makeVector();
            int mask = 1;
            for(int j=0; j<D; j++) {
                //the bits of i indicate whether corresponding primitive vector is added
                //do bitwise "and" with mask (which has exactly 1 nonzero bit) to see if
                //each bit is 1 or 0
                if((i & mask) != 0) vertices[i].PE(latticeVectors[j]);//unit vector
                mask *= 2;
            }
            vertices[i].TE(size);//scale to size (good for PrimitiveCubic only)
            vertices[i].PE(coord.position());//translate to position
  //          System.out.println(vertices[i].toString());
        }
        return vertices;
    }//end of vertex
    
    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     * /
    public boolean inCell(Space.Vector v) {
        delta.Ev1Mv2(v, coord.position());
        double x = size;
        switch(D()) {
            case 3: x = delta.x(2);
                    if(x < 0.0 || x > size) return false;
            case 2: x = delta.x(1);//fall through to check all dimensions
                    if(x < 0.0 || x > size) return false;
            case 1: x = delta.x(0);
                    if(x < 0.0 || x > size) return false;
                    break;
            default: throw new RuntimeException("PrimitiveCubic.UnitCell.inCell not implemented for given dimension");
        }
        return true;
    }
}//end of UnitCell

    /**
     * Main method to demonstrate use and to aid debugging
     * /
    public static void main(String[] args) {
        System.out.println("main method for PrimitiveCubic");
        Space space = new Space2D();
        Simulation sim = new Simulation(space);
        int D = space.D();
        PrimitiveCubic primitive = new PrimitiveCubic(sim);
        AtomFactory siteFactory = primitive.unitCellFactory();
        final int nx = 4;
        final int ny = 5;
        BravaisLattice lattice = BravaisLattice.makeLattice(sim, 
                                siteFactory, 
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
    
