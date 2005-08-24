package etomica.lattice.crystal;

import etomica.lattice.Primitive;
import etomica.math.geometry.Cube;
import etomica.math.geometry.Polytope;
import etomica.math.geometry.Square;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Primitive group for a cubic system.  All primitive
 * vectors orthogonal and of equal length.
 */
public class PrimitiveCubic extends Primitive implements Primitive2D, Primitive3D {
    
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    //default size is 1.0
    private double cubicSize;
    
    public PrimitiveCubic(Space space) {
        this(space, 1.0);
    }
    public PrimitiveCubic(Space space, double latticeConstant) {
        super(space); //also makes reciprocal
        //set up orthogonal vectors of unit size
        for(int i=0; i<D; i++) latticeVectors[i].setX(i, 1.0);
        setCubicSize(latticeConstant); //also sets reciprocal via update
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveCubic(Space space, Primitive direct) {
        super(space, direct);
        for(int i=0; i<D; i++) latticeVectors[i].setX(i, 1.0);
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveCubic(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveCubic)reciprocal()).setCubicSize(2.0*Math.PI/cubicSize);
    }
    
    public void setA(double a) {setCubicSize(a);}
    public double getA() {return cubicSize;}
    
    public void setB(double b) {setCubicSize(b);}
    public double getB() {return cubicSize;}
        
    public void setC(double c) {setCubicSize(c);}
    public double getC() {return cubicSize;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {}
    public double getBeta() {return rightAngle;}
    
    public void setGamma(double t) {}
    public double getGamma() {return rightAngle;}
    
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
        return new PrimitiveCubic(space, cubicSize);
    }
    
    //override superclass method to scale copy-vectors to current size
    protected Vector[] copyVectors() {
        for(int i=0; i<D; i++) {
            latticeVectorsCopy[i].E(latticeVectors[i]);
            latticeVectorsCopy[i].TE(cubicSize);
        }
        return latticeVectorsCopy;
    }
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setCubicSize(double cubicSize) {
        if(immutable || cubicSize <= 0.0) return;
        this.cubicSize = cubicSize;
        for(int i=0; i<D; i++) {
            size[i] = cubicSize;
        }
        update();
    }
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getCubicSize() {return cubicSize;}
    
    public void scaleSize(double scale) {
        setCubicSize(scale*cubicSize);
    }

    public int[] latticeIndex(Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/cubicSize;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/cubicSize;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    /**
     * Returns a new Square (if primitive is 2D) or Cube (if 3D) with edges
     * given by the size of the primitive vectors.
     */
    public Polytope wignerSeitzCell() {
        return (D == 2) ? (Polytope)new Square(space,cubicSize) :  (Polytope)new Cube(space,cubicSize);
    }
    
    /**
     * Returns a new Square (if primitive is 2D) or Cube (if 3D) with edges
     * given by the size of the primitive vectors.
     */
    public Polytope unitCell() {
        return (D == 2) ? (Polytope)new Square(space,cubicSize) :  (Polytope)new Cube(space,cubicSize);
    }
    
    public String toString() {return "Cubic";}
    
//    /**
//     * Main method to demonstrate use and to aid debugging
//     */
//    public static void main(String[] args) {
//        System.out.println("main method for PrimitiveCubic");
//        Space space = new Space2D();
//        int D = space.D();
//        PrimitiveCubic primitive = new PrimitiveCubic(space);
//        AtomFactory siteFactory = primitive.unitCellFactory();
//        final int nx = 4;
//        final int ny = 5;
//        BravaisLattice lattice = BravaisLattice.makeLattice(space, 
//                                new int[] {nx,ny},
//                                primitive);
//        lattice.shiftFirstToOrigin();
//        System.out.println("Total number of sites: "+lattice.siteList().size());
//        System.out.println();
//        System.out.println("Coordinate printout");
//        AtomsetActionAdapter printSites = new AtomsetActionAdapter() {
//            public void actionPerformed(Atom[] s) {
//                System.out.print(s[0].coord.position().toString()+" ");
//                if(((Site)s[0]).latticeCoordinate()[1]==ny-1) System.out.println();
//            }
//        };
//        AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
//        iterator.allAtoms(printSites);
//        System.out.println();
//                
//        Site testSite = (Site)lattice.site(new int[] {1,2});
//        
//        Space.Vector vector = Space.makeVector(new double[] {1.5, 2.7});
//        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
//        int[] idx = primitive.latticeIndex(vector);
//        System.out.print("; is in: ");
//        for(int i=0; i<D; i++) System.out.print(idx[i]);
//        System.out.println();
// 
//        vector = Space.makeVector(new double[] {3.5, 5.1});
//        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
//        idx = primitive.latticeIndex(vector);
//        System.out.print("; is in: ");
//        for(int i=0; i<D; i++) System.out.print(idx[i]);
//        System.out.println();
//        
//        System.out.println("cell volume: "+testSite.volume());
//        System.out.println();
//        System.out.println("setting size to 2.0");
//        primitive.setSize(2.0);
//        lattice.update();//recompute all lattice positions
//        lattice.shiftFirstToOrigin();
//        
//        iterator.allAtoms(printSites);
//        System.out.println();
//        vector = Space.makeVector(new double[] {1.5, 2.7});
//        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
//        idx = primitive.latticeIndex(vector);
//        System.out.print("; is in: ");
//        for(int i=0; i<D; i++) System.out.print(idx[i]);
//        System.out.println();
//
//        vector = Space.makeVector(new double[] {3.5, 5.1});
//        System.out.print(vector.toString()+" in cell "+testSite.toString()+"? "+testSite.inCell(vector));
//        idx = primitive.latticeIndex(vector);
//        System.out.print("; is in: ");
//        for(int i=0; i<D; i++) System.out.print(idx[i]);
//        System.out.println();
//
//        System.out.println("cell volume: "+primitive.unitCell().volume());
//  
//        /*
//        //write out vertices of some cells
//        System.out.println();
//        System.out.println("Unit-cell vertices of first cell");
//        Space.Vector[] vertices = ((AbstractCell)lattice.siteList().getFirst()).vertex();
//        for(int i=0; i<vertices.length; i++) {
//            System.out.println(vertices[i].toString());
//        }
//        System.out.println();
//        System.out.println("Unit-cell vertices of last cell");
//        vertices = ((AbstractCell)lattice.siteList().getLast()).vertex();
//        for(int i=0; i<vertices.length; i++) {
//            System.out.println(vertices[i].toString());
//        }
//        */
//    }//end of main
    
}//end of PrimitiveCubic
    
