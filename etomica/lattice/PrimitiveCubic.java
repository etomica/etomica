package etomica.lattice;
import etomica.*;

/**
 * Primitive group for a cubic system.  All primitive
 * vectors orthogonal and of equal length.
 */
public class PrimitiveCubic extends Primitive {
    
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    private double size;
    
    public PrimitiveCubic(Simulation sim) {
        super(sim);
        //set up orthogonal vectors of unit size
        for(int i=0; i<D; i++) latticeVectors[i].setX(i, 1.0);
    }
    
    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        PrimitiveCubic copy = new PrimitiveCubic(simulation);
        copy.setSize(size);
        return copy;
    }
    
    //override superclass method to scale copy-vectors to current size
    protected Space.Vector[] copyVectors() {
        for(int i=0; i<D; i++) {
            latticeVectorsCopy[i].E(latticeVectors[i]);
            latticeVectorsCopy[i].TE(size);
        }
        return latticeVectorsCopy;
    }
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
 //       for(int i=0; i<D; i++) latticeVectors[i].setx(i,size);
        this.size = size;
        if(lattice != null) lattice.update();
    }
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getSize() {return size;}
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.  All values in given array must be
     * equal to each other or else and IllegalArgumentException is thrown.
     */
    public void setSize(double[] size) {
        if(size.length != D) throw new IllegalArgumentException("Error in PrimitiveCubic.setSize: Number of sizes given is inconsistent with number of primitive vectors");
        //check that all values in size array are equal to each other
        double newSize = size[0];
        for(int i=1; i<D; i++) {
            if(size[i] != newSize) throw new IllegalArgumentException("Error in PrimitiveCubic.setSize(double[]): elements of size array are not identical in value");
        }
        setSize(newSize);
    }
    
    public void scaleSize(double scale) {
        setSize(scale*size);
    }

    public int[] latticeIndex(Space.Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }
    
    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }

    public Primitive reciprocal() {
        throw new RuntimeException("method PrimitiveCubic.reciprocal not yet implemented");
    }
    
    public AtomFactory wignerSeitzCellFactory() {
        throw new RuntimeException("method PrimitiveCubic.wignerSeitzCell not yet implemented");
    }
    
    public AtomFactory unitCellFactory() {
        return new UnitCellFactory(simulation);
    }
    
///////////////////////////////////////////////////////////////////////////////////////////

public class UnitCellFactory extends AtomFactory {

    AtomType atomType;
    
    public UnitCellFactory(Simulation sim) {
        super(sim);
        setType(new AtomType(this));//default
    }
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}

    /**
     * Builds a single unit cell.
     */
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
public class UnitCell extends AbstractCell {
    
    private final Space.Vector delta;
    
    public UnitCell(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, parent);
        delta = space.makeVector();
    }
    /**
     * Dimension of the space occupied by the cell
     */
     public int D() {return space.D();}
     
    /**
     * Returns the volume of the cubic cell.
     */
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
     */
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
     */
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
     */
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
        */
    }//end of main
    
}//end of PrimitiveCubic
    
