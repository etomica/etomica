package etomica.lattice;
import etomica.*;

/**
 * Primitive group for a face-centered-cubic system.
 */
public class PrimitiveFcc extends Primitive implements Primitive3D {
    
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    //default size is 1.0
    private double size;
    private Space.Vector[] unitVectors;
    private static final double FCC_ANGLE = Math.acos(0.5);
    
    public PrimitiveFcc(Simulation sim) {
        this(sim, 1.0);
    }
    public PrimitiveFcc(Simulation sim, double size) {
        super(sim); //also makes reciprocal
        //set up orthogonal vectors of unit size
        unitVectors = new Space.Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = sim.space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
        setSize(size); //also sets reciprocal via update
    }
    /**
     * Constructor used by makeReciprocal method of PrimitiveBcc.
     */
    PrimitiveFcc(Simulation sim, Primitive direct) {
        super(sim, direct);
        unitVectors = new Space.Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = sim.space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveBcc(simulation, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveBcc)reciprocal()).setSize(4.0*Math.PI/size);
    }
    
    public void setA(double a) {setSize(a);}
    public double getA() {return size;}
    
    public void setB(double b) {setSize(b);}
    public double getB() {return size;}
        
    public void setC(double c) {setSize(c);}
    public double getC() {return size;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return FCC_ANGLE;}
    
    public void setBeta(double t) {}
    public double getBeta() {return FCC_ANGLE;}
    
    public void setGamma(double t) {}
    public double getGamma() {return FCC_ANGLE;}
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return false;}
    public boolean isEditableC() {return false;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return false;}
    public boolean isEditableGamma() {return false;}

    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveFcc(simulation, size);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        this.size = size;
        for(int i=0; i<D; i++) latticeVectors[i].Ea1Tv1(size,unitVectors[i]);
        update();
    }
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getSize() {return size;}
    
    /**
     * Multiplies the size of the current vectors by the given value.
     */
    public void scaleSize(double scale) {
        setSize(scale*size);
    }

    public int[] latticeIndex(Space.Vector q) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
/*        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
*/    }
    
    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
 /*       for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
 */   }
    
    public AtomFactory wignerSeitzCellFactory() {
        throw new RuntimeException("method PrimitiveFcc.wignerSeitzCell not yet implemented");
    }
    
    public AtomFactory unitCellFactory() {
        throw new RuntimeException("method unitCellFactory not yet implemented");
//        return new UnitCellFactory(simulation);
    }
    
    public String toString() {return "Fcc";}
    
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
 * /
public class UnitCell extends AbstractCell {
    
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
     * /
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
}//end of PrimitiveFcc
    
