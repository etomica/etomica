package etomica.lattice;
import etomica.Space;
import etomica.Space2D;

/**
 * Square lattice in 2D.  Primitive lattice vectors are orthogonal and equal in size.
 */
public class SquareLattice extends BravaisLattice {
    
    private double primitiveVectorLength;
    public Site sites[][];
    private Space2D.Vector translation2D;
    
    /**
     * Constructs a "square" square lattice, i.e., one with the same number of sites in each dimension.
     */
    public SquareLattice(int dimension, SiteFactory factory, double pVectorLength) {
        this(new int[] {dimension, dimension}, factory, pVectorLength);
    }
    /**
     * Constructs a square (2D) lattice with lattice constant given by <code>pVectorLength</code>.
     * Number of cells in each dimension is given by the dimensions array, and the sites
     * of the lattice are constructed by the given SiteFactory.
     */
    public SquareLattice(int[] dimensions, SiteFactory factory, double pVectorLength) {
        super(dimensions, factory, 
                new Space.Vector[] {new Space2D.Vector(pVectorLength,0.0), new Space2D.Vector(0.0,pVectorLength)});
        if(dimensions.length != 2) { //should throw an exception
            System.out.println("Dimension conflict in constructor of SquareLattice");
            System.exit(1);
        }
        primitiveVectorLength = pVectorLength;
        translation2D = (Space2D.Vector)translation;
        //Create array of handles for easy reference to sites
        sites = new Site[dimensions[0]][dimensions[1]];
        for(int i=0; i<dimensions[0]; i++) {
            for(int j=0; j<dimensions[1]; j++) {
                sites[i][j] = site(new IntegerBravaisLattice.Coordinate(new int[] {i,j}));
            }
        }
    }
    public Site nearestSite(Space2D.Vector r, Space2D.Vector scale) {
        return sites[(int)((r.x/scale.x-translation2D.x)/primitiveVectorLength + 0.5)][(int)((r.y/scale.y-translation2D.y)/primitiveVectorLength + 0.5)];
    }        
    public Site nearestSite(Space2D.Vector r) {
        return sites[(int)((r.x-translation2D.x)/primitiveVectorLength + 0.5)][(int)((r.y-translation2D.y)/primitiveVectorLength + 0.5)];
//        int i = (int)((r.x-translation2D.x)/primitiveVectorLength + 0.5);
//        int j = (int)((r.y-translation2D.y)/primitiveVectorLength + 0.5);
//        System.out.println(i+" "+j);
//        return sites[i][j];
    }
    /**
     * Overrides superclass method so that it PERFORMS NO ACTION.
     * The only way to change the primitive vector is through the method taking the <code>double</code> argument.
     */
    public final void setPrimitiveVector(Space.Vector[] pVectors) {return;}
    /**
     * Overrides superclass method so that it PERFORMS NO ACTION.
     * The only way to change the primitive vector is through the method taking the <code>double</code> argument.
     */
    public final void setPrimitiveVector(int i, Space.Vector pVector) {return;}
    /**
     * Overrides superclass method so that it PERFORMS NO ACTION.
     * The only way to change the primitive vector is through the method taking the <code>double</code> argument.
     */
    public final void setPrimitiveVector(int vectorIndex, int xyzIndex, double value) {return;}
    /**
     * Changes the size of the primitive vector, or lattice constant, for the square lattice
     */
    public void setPrimitiveVectorLength(double length) {
        primitiveVector[0].setComponent(0,length); //sets the x component of the first vector
        primitiveVector[1].setComponent(1,length); //sets the y component of the second vector
        primitiveVectorLength = length;
        translation.E(0.);
        updateCoordinates();//this zeros any translations done to the cells
        notifyObservers();
    }
    public double getPrimitiveVectorLength() {return primitiveVectorLength;}
    
    /**
     * A cubic primitive cell.
     */
    public static class Cell extends AbstractCell {
        private SquareLattice squareLattice;
        //would prefer this constructor, but Cell is made by IntegerBravaisLattice, which inserts its Coordinate
        //it is subsequently replaced by a BravaisLattice.Coordinate
        public Cell(SquareLattice parent, BravaisLattice.Coordinate coord) {
            super(parent, coord);
            squareLattice = parent;
        }
//        public Cell(AbstractLattice.Coordinate coord) {
 //           super(SquareLattice.this, coord);
 //       }
        /**
        * Dimension of the space occupied by the cell
        */
        public int D() {return 2;}
         
        /**
        * Returns the volume of the cell.
        */
        public double volume() {return squareLattice.primitiveVectorLength*squareLattice.primitiveVectorLength;}
        /**
        * Returns the absolute positions of the vertices relative to the cell position.
        */
        public Space.Vector[] vertex() {
            double a = 0.5*squareLattice.primitiveVectorLength;
            Space.Vector r = ((AbstractLattice.PositionCoordinate)coordinate()).position();
            Space.Vector[] v = new Space.Vector[] {
            new Space2D.Vector(+a,+a),
            new Space2D.Vector(+a,-a),
            new Space2D.Vector(-a,-a),
            new Space2D.Vector(-a,+a)};
            for(int i=0; i<v.length; i++) {v[i].PE(r);}
            return v;
/*            return new Space.Vector[] {
            new Space2D.Vector(0.5,0.5),
            new Space2D.Vector(0.5,-0.5),
            new Space2D.Vector(-0.5,-0.5),
            new Space2D.Vector(-0.5,0.5)};*/
        }
        /**
        * Returns <code>true</code> if the given vector lies inside the cell, <code>false</code> otherwise.
        */
        public boolean inCell(Space.Vector v) {
            Space.Vector r = ((AbstractLattice.PositionCoordinate)coordinate()).position();
            double dx = v.component(0)-r.component(0);
            double dy = v.component(1)-r.component(1);
            double a = 0.5*squareLattice.primitiveVectorLength;
            return -0.99999*a <= dx && dx <= 1.000001*a && -0.99999*a <= dy && dy <= 1.00001*a;
        }
        
    }//end of SquareLattice.Cell
    
    /**
     * A factory that makes Sites of type SquareLattice.Cell
     */
    public static class CellFactory implements SiteFactory {
        public Site makeSite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
            if(!(parent instanceof SquareLattice)) {
                throw new IllegalArgumentException("SquareLattice.CellFactory: parent lattice must be of type SquareLattice");
            }
            if(!(coord instanceof BravaisLattice.Coordinate)) {
                throw new IllegalArgumentException("SquareLattice.Cell.Factory: coordinate must be of type BravaisLattice.Coordinate");
            }
            return new Cell((SquareLattice)parent, (BravaisLattice.Coordinate)coord);
//            return ((SquareLattice)parent).new Cell(coord);
        }
    }//end of SquareLattice.CellFactory

    /**
     * Main method to demonstrate use of class and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for SquareLattice");
        int D = 2;
        SquareLattice lattice = 
            new SquareLattice(new int[] {3,3}, new CellFactory(), 1.0);        
        SiteIterator iterator = lattice.iterator();
        System.out.println();
        System.out.println("Total number of sites: "+lattice.siteCount());
        System.out.println();
        System.out.println("Coordinate printout");
        SiteAction printSites = new SiteAction() {public void actionPerformed(Site s) {System.out.print(s.toString()+" ");}};
        iterator.allSites(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Translating lattice by (-1.0, 2.0)");
        lattice.translateBy(Space.makeVector(new double[] {-1.0, 2.0}));
        iterator.allSites(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Changing primitive vector to 2.0 and translating by (1,1)");
        lattice.setPrimitiveVectorLength(2.0);
        lattice.translateBy(new Space2D.Vector(1.0,1.0));
        iterator.allSites(printSites);
        System.out.println();
        System.out.println();
        
        System.out.print("Vertices for cell (1,1): ");
        Cell testSite = (Cell)lattice.site(new IntegerBravaisLattice.Coordinate(new int[] {1,1}));
        Space.Vector[] v = testSite.vertex();
        for(int i=0; i<v.length; i++) {
            System.out.print("("+((Space2D.Vector)v[i]).x+","+((Space2D.Vector)v[i]).y+") ");
        }
        System.out.println();
        System.out.println();
        System.out.println("(2.5,3.5) is in this region: "+testSite.inCell(new Space2D.Vector(2.5,3.5)));
        System.out.println("(0.5,2.5) is in this region: "+testSite.inCell(new Space2D.Vector(0.5,2.5)));
        System.out.println();
        System.out.println("Nearest site to (1.5,2.5): "+lattice.nearestSite(new Space2D.Vector(1.5,2.5)));
        System.out.println("Nearest site to (0.5,1.0): "+lattice.nearestSite(new Space2D.Vector(0.5,1.0)));
        System.out.println("Nearest site to (5.5,5.9): "+lattice.nearestSite(new Space2D.Vector(5.5,5.9)));
        System.out.println();
                
        System.out.print("A randomly selected site: ");
        System.out.println(lattice.randomSite().toString());
        
        System.out.println();
        System.out.println("All sites by looping through array:");
        for (int i=0; i<lattice.dimensions()[0]; i++) {
            for(int j=0; j<lattice.dimensions()[1]; j++) {
                System.out.print(lattice.sites[i][j].toString()+" ");
            }
            System.out.println();
        }
    }//end of main
        
   
}//end of SquareLattice