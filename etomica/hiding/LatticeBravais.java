package simulate;

//import java.util.Array;  //java 1.2  need for a (unimplemented) constructor
import java.util.Random;

/**
 * Arbitrary-dimension Bravais Lattice.  
 * Lattice is constructed recursively, so that at top level (dimension D) the lattice is
 * represented by an array of lattices each of dimension D-1.  This continues on down to
 * the zero-dimensional lattice, which contains a single Lattice.Site object.
 * Each dimension of the lattice may be of different length (i.e., the lattice can be rectangular and need not be strictly cubic)
 * 
 * The lattice occupies a D'-dimensional space.  Normally D' = D at the top level, but this need not be 
 * the case (e.g., we may have a 2-dimensional (planar) lattice given at some orientation in 3-d).
 * The placement of the sites in the D'-dimensional space is given by the primitiveVector array.  One
 * primitive vector is stored with each sublattice, describing how that sublattice was translated
 * when copied to form the D-dimensional lattice.
 */
public class LatticeBravais implements Lattice {

    /*  ****Compiler (sometimes) won't let me declare these final 
    private final Lattice.Site.Iterator iterator;
    private final int nRows, siteCount;
    public final int D;                              //dimension of lattice
    public final static Random random = new Random(); */   //random-number generator for selecting a random site
    private BravaisIterator iterator;
    private int nRows, siteCount;
    public int D;                                //dimension of lattice
    public static Random random;  //random-number generator for selecting a random site
    static {random = new Random();}
    
    private LatticeBravais[] rows;                       //array of D-1 dimensional lattices (generically called a "row" here)
    private LatticeBravais parentLattice;                //D+1 dimensional lattice in which this lattice is a row (null for highest-level lattice)
    private int index;                                 //index of row array for this lattice in parentLattice 
    
    private double[] primitiveVector;
    
    //modify to construct lattice from arbitrary site objects
    /**
     * Default constructor is a "zero-dimensional" lattice consisting of one site
     */
    public LatticeBravais() { 
        D = 0;
        iterator = new SingletIterator(new Site(this));
        siteCount = 1;
        nRows = 0;
    }
    public LatticeBravais(int[] dim, double[][] prim, final double r2Max) {    //constructor for top-level lattice
        this(dim.length,dim,prim);
        class Criterion implements Lattice.Site.NeighborIterator.Criterion {
            public boolean isNeighbor(Lattice.Site s1, Lattice.Site s2) {return r2((Site)s1,(Site)s2) < r2Max;}
        }
        Criterion nbrTest = new Criterion();
        Lattice.Site.Iterator iter = (Lattice.Site.Iterator)iterator.clone();
        iter.reset();   //need to pass iterator to setNeighbors, which changes its state, so we do outer loop with a clone of iterator
        while(iter.hasNext()) {
            iterator.reset();
            iter.next().adjacentIterator().setNeighbors(iterator,nbrTest);
        }
    }
    private LatticeBravais(int D, int[] dim, double[][] prim) {  //constructor for sub-lattices
        this.D = D;
        int d = D-1;
        primitiveVector = prim[d];
        nRows = dim[d];
        rows = new LatticeBravais[nRows];
        for(int i=0; i<nRows; i++) {
            rows[i] = (d>0) ? new LatticeBravais(d,dim, prim) : new LatticeBravais();
            rows[i].setParentLattice(this, i);
        }
        /*
        rows[0] = (d>0) ? new LatticeBravais(newDim, prim) : new LatticeBravais();
        rows[0].setParentLattice(this, 0);
        iterator = new LatticeIterator(rows);  //rows[0] needs instantiation before creating iterator
        for(int i=1; i<nRows; i++) {  //loop to set up neighbors in adjacent rows
            rows[i] = (d>0) ? new LatticeBravais(newDim, prim) : new LatticeBravais();
            rows[i].setParentLattice(this, i);
            rows[i-1].iterator().reset();
            rows[i].iterator().reset();
            while (rows[i].iterator().hasNext()) {
                Lattice.Site pSite = rows[i-1].iterator().next(); 
                Lattice.Site nSite = rows[i].iterator.next();
                ((AdjacentIterator)pSite.adjacentIterator()).add(nSite);
                ((AdjacentIterator)nSite.adjacentIterator()).add(pSite);
            }
        }
        //This does the periodic boundary
        rows[0].iterator().reset();
        rows[nRows-1].iterator().reset();
        while (rows[0].iterator().hasNext()) {
            Lattice.Site pSite = rows[0].iterator().next(); 
            Lattice.Site nSite = rows[nRows-1].iterator.next();
            ((AdjacentIterator)pSite.adjacentIterator()).add(nSite);
            ((AdjacentIterator)nSite.adjacentIterator()).add(pSite);
        }
        */
        iterator = new LatticeIterator(rows);
        siteCount = nRows*rows[0].siteCount();
    }
    /**
     * Create a square lattice of dimension "d" with "size" elements in each dimension
     * Need java 1.2 for Array.fill
     */
//    public LatticeBravais(int d, int size) {
//        this(construct array of int with d elements each having value size);
//          this(Array.fill(new int[d],size));
//    }

    //Lattice methods
    public int D() {return D;}
    public final int siteCount() {return siteCount;}  
    public int coordinationNumber() {return 2*D;}
    public Lattice.Site site(Lattice.Coordinate coord) {return site((Coordinate)coord);}
    public Lattice.Site site(Coordinate coord) {
        return (D==0) ? iterator.first() : rows[coord.index[D-1]].site(coord);
    }
        
    /**
     * Returns a random site in the lattice
     * Iteratively chooses a row at random and calls randomSite for that row
     */
    public Lattice.Site randomSite() {
        int i = (int)Math.floor(nRows*random.nextDouble());
        return (D==0) ? iterator.first() : rows[i].randomSite();}
        
    public Lattice.Site.Iterator iterator() {return iterator;}     //iterator for all sites in lattice
    
    
    //LatticeBravais methods
    public double[] primitiveVector(int i) {  //indexing starts at zero
        return (i==D-1) ? primitiveVector : rows[0].primitiveVector(i);
    }
    public double r2(Site s1, Site s2) {
        double[] deltaR = new double[primitiveVector.length];
        calculateDeltaR(deltaR, s1.coordinate.index, s2.coordinate.index);
        double r2 = 0.0;
        for(int i=0; i<deltaR.length; i++) {r2 += deltaR[i]*deltaR[i];}
        return r2;
    }
    private void calculateDeltaR(double[] deltaR, int[] index1, int[] index2) {
        int delta = index1[D-1] - index2[D-1];
        for(int i=0; i<primitiveVector.length; i++) {
            deltaR[i] += delta*primitiveVector[i];
        }
        if(D > 0) rows[0].calculateDeltaR(deltaR, index1, index2);
    }
        
    public LatticeBravais parentLattice() {return parentLattice;}
    public int index() {return index;}
    public void setParentLattice(LatticeBravais lat, int i) {parentLattice = lat; index = i;}
    
    public static class Coordinate implements Lattice.Coordinate{
        //first value (index[0]) indicates the site in the row
        //second value indicates which row in the plane, third value which plane in the cube, etc.
        public int[] index;
        public int D;
        public Coordinate(int D) {this.D = D; index = new int[D];}
        public Coordinate(int[] i) {setIndex(i);}
        public void setIndex(int[] i) {
            if(i.length == D) {index = i;}
            else {   //should throw exception
                index = new int[D];
                System.out.println("LatticeBravais.Coordinate instantiation error");
            }
        }
        public void setIndex(int d, int value) {  //doesn't check for correct range of d and value
            // require 0 <= d < D; 0 <= value < (size of row at dimension depth d)
            index[d] = value;
        }
        
    }
    
    /**
     * General iterator of sites on LatticeBravais.
     * SingletIterator is used for zero-D lattices
     */
     
     // Here's a recursive (and, unfortunately, hard to clone) version of the LatticeIterator
     
/*    private static final class LatticeIterator implements Lattice.Site.Iterator {  //might instead do this by creating a big array of all sites and loop through it
        private boolean hasNext;
        private int iRow;
        private Lattice.Site.Iterator current;
        private final Lattice[] rows;
        private final int nRows;
        public LatticeIterator(Lattice[] r) {rows = r; nRows = r.length; reset();}   //constructor
        public Lattice.Site first() {return rows[0].iterator().first();}
        public boolean hasNext() {return hasNext;}
        public void reset() {
            iRow = 0;
            current = rows[0].iterator();
            current.reset();
            hasNext = current.hasNext();
        }
        public Lattice.Site next() {
            Lattice.Site nextSite = current.next();
            if(!current.hasNext()) {  //no more in current row
                iRow++;
                if(iRow == nRows) {hasNext = false;} //that was the last row
                else {
                    current = rows[iRow].iterator();
                    current.reset();
                }
            }
            return nextSite;
        }
        public void allSites(Lattice.Site.Action act) {
            for(int i=0; i<nRows; i++) {rows[i].iterator().allSites(act);}
        }
        public Object clone() { 
            try {return super.clone();} catch(CloneNotSupportedException e) {return null;}
        }
    }
 */       
    /**
     * General iterator of sites on LatticeBravais.
     * SingletIterator is used for zero-D lattices
     */
     
     //Don't like this much, but it should work
    private interface BravaisIterator extends Lattice.Site.Iterator {
        public Lattice.Site.Linker firstLink();
        public Lattice.Site.Linker lastLink();
    }        
    private static final class LatticeIterator implements BravaisIterator {  //might instead do this by creating a big array of all sites and loop through it
        private boolean hasNext;
        private Lattice.Site.Linker firstLink, nextLink, lastLink;
        public LatticeIterator(LatticeBravais[] rows) {  //constructor
            int nRows = rows.length;
            firstLink = ((BravaisIterator)rows[0].iterator()).firstLink();
            lastLink = ((BravaisIterator)rows[nRows-1].iterator()).lastLink();
            for(int i=0; i<nRows-1; i++) {  //Loop through all rows, joining last link of each row with first link of next row
                ((BravaisIterator)rows[i].iterator()).lastLink().next = ((BravaisIterator)rows[i+1].iterator()).firstLink();
            }
            reset();
        }   
        public Lattice.Site first() {return firstLink.site;}
        public Lattice.Site.Linker firstLink() {return firstLink;}
        public Lattice.Site.Linker lastLink() {return lastLink;}
        public boolean hasNext() {return nextLink != null;}
        public void reset() {nextLink = firstLink;}
        public Lattice.Site next() {
            Lattice.Site nextSite = nextLink.site;
            nextLink = nextLink.next;
            return nextSite;
        }
        public void allSites(Lattice.Site.Action act) {
            for(Lattice.Site.Linker link=firstLink; link!=null; link=link.next) {act.action(link.site);}
        }
        public Object clone() { 
            try {return super.clone();} catch(CloneNotSupportedException e) {return null;}
        }
    }

    //Iterator for a lattice that contains only one site (zero-D lattice)
    private static final class SingletIterator implements BravaisIterator {
        public final Site site;
        private final Lattice.Site.Linker link;
        private boolean hasNext;
        public SingletIterator(Site s) {
            site = s; 
            hasNext = true; 
            link = new Lattice.Site.Linker(null,s);
        }
        public boolean hasNext() {return hasNext;}
        public Lattice.Site first() {return site;}
        public Lattice.Site.Linker firstLink() {return link;}
        public Lattice.Site.Linker lastLink() {return link;}
        public void reset() {hasNext = true;}
        public Lattice.Site next() {hasNext = false; return site;}
        public void allSites(Lattice.Site.Action act) {act.action(site);}
        public Object clone() { 
            try {return super.clone();} catch(CloneNotSupportedException e) {return null;}
        }
    }
    
    public static final class Site implements Lattice.Site {
        private final Lattice.Site.NeighborIterator adjacentIterator = new Lattice.Site.NeighborIterator(this);
        private final LatticeBravais lattice;
        final Coordinate coordinate;
        public Site(LatticeBravais p) {
            LatticeBravais lat = p;
            int i = 0;
            while(lat.parentLattice() != null) {lat = lat.parentLattice();}
            this.lattice = lat;
            coordinate = new Coordinate(lattice.D());
            lat = p;
            while(lat.parentLattice() != null) {coordinate.setIndex(i++, lat.index());}
        }
        public Lattice lattice() {return lattice;}  //returns the top-level lattice on which this site resides
        public Lattice.Site.NeighborIterator adjacentIterator() {return adjacentIterator;}
        public Lattice.Coordinate coordinate() {return coordinate;}
        public boolean isAdjacentTo(Lattice.Site s) {
            adjacentIterator.reset();
            while(adjacentIterator.hasNext()) {
                if(s == adjacentIterator.next()) {return true;}
            }
            return false;
        }
            

        /**
         * Returns the hierarchy of coordinates for this site
         */
/*        public String toString() {
            String string = "Site: ";
            for(LatticeBravais lat=parent; lat!=null; lat=lat.parentLattice()) {
                string += (String.valueOf(lat.index()) + ", ");
            }
            return string;
        }
        */    
    }
    
    /**
     * Main method to demonstrate use of LatticeBravais and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for LatticeBravais");
    }
}