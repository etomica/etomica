package simulate;

//import java.util.Array;  //java 1.2  need for a (unimplemented) constructor
import java.util.Random;

/**
 * Arbitrary-dimension Lattice with cubic coordination.  
 * Lattice is constructed recursively, so that at top level (dimension D) the lattice is
 * represented by an array of lattices each of dimension D-1.  This continues on down to
 * the zero-dimensional lattice, which contains a single Lattice.Site object.
 * Each dimension of the lattice may be of different length (i.e., the lattice can be rectangular and need not be strictly cubic)
 */
public class LatticeCubic implements Lattice {

    /*  ****Compiler (sometimes) won't let me declare these final 
    private final Lattice.Site.Iterator iterator;
    private final int nRows, siteCount;
    public final int D;                              //dimension of lattice
    public final static Random random = new Random(); */   //random-number generator for selecting a random site
    private Lattice.Site.Iterator iterator;
    private int nRows, siteCount;
    public int D;                                //dimension of lattice
    public static Random random = new Random();  //random-number generator for selecting a random site
    
    private LatticeCubic[] rows;                       //array of D-1 dimensional lattices
    private LatticeCubic parentLattice;                //D+1 dimensional lattice in which this lattice is a row (null for highest-level lattice)
    private int index;                                 //index of row array in parentLattice
    
    //modify to construct lattice from arbitrary site objects
    /**
     * Default constructor is a "zero-dimensional" lattice consisting of one site
     */
    public LatticeCubic() { 
        D = 0;
        iterator = new SingletIterator(new Site(this));
        siteCount = 1;
        nRows = 0;
    }
    public LatticeCubic(int[] dim) {
        D = dim.length;
        int d = D-1;
        int[] newdim = new int[d];
        if(d > 0) System.arraycopy(dim,0,newdim,0,d);  //create new dim array from all but last element of given dim array
        nRows = dim[d];
        rows = new LatticeCubic[nRows];
        rows[0] = (d>0) ? new LatticeCubic(newdim) : new LatticeCubic();
        rows[0].setParentLattice(this, 0);
        for(int i=1; i<nRows; i++) {  //loop to set up neighbors in adjacent rows
            rows[i] = (d>0) ? new LatticeCubic(newdim) : new LatticeCubic();
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
        iterator = new LatticeIterator(rows);
        siteCount = nRows*rows[0].siteCount();
    }
    /**
     * Create a square lattice of dimension "d" with "size" elements in each dimension
     * Need java 1.2 for Array.fill
     */
//    public LatticeCubic(int d, int size) {
//        this(construct array of int with d elements each having value size);
//          this(Array.fill(new int[d],size));
//    }

    public final int D() {return D;}
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
    
    public LatticeCubic parentLattice() {return parentLattice;}
    public int index() {return index;}
    public void setParentLattice(LatticeCubic lat, int i) {parentLattice = lat; index = i;}
    
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
                System.out.println("LatticeCubic.Coordinate instantiation error");
            }
        }
        public void setIndex(int d, int value) {  //doesn't check for correct range of d and value
            // require 0 <= d < D; 0 <= value < (size of row at dimension depth d)
            index[d] = value;
        }
        
    }
    
    /**
     * General iterator of sites on LatticeCubic.
     * SingletIterator is used for zero-D lattices
     */
    private static final class LatticeIterator implements Lattice.Site.Iterator {  //might instead do this by creating a big array of all sites and loop through it
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
    }
        
    //Iterator for a lattice that contains only one site (zero-D lattice)
    private static final class SingletIterator implements Lattice.Site.Iterator {
        public final Site site;
        private boolean hasNext;
        public SingletIterator(Site s) {site = s; hasNext = true;}
        public boolean hasNext() {return hasNext;}
        public Lattice.Site first() {return site;}
        public void reset() {hasNext = true;}
        public Lattice.Site next() {hasNext = false; return site;}
        public void allSites(Lattice.Site.Action act) {act.action(site);}
    }
    
    //Iterator to generate neighbors of site on LatticeCubic
    private static final class AdjacentIterator implements Lattice.Site.Iterator {
        private Lattice.Site.Linker first, nextLink;
        public void add(Lattice.Site s) {first = new Lattice.Site.Linker(first, s);}
        public boolean hasNext() {return nextLink != null;}
        public void reset() {nextLink = first;}
        public Lattice.Site first() {return first.site;}
        public Lattice.Site next() {
           Lattice.Site nextSite = nextLink.site;
           nextLink = nextLink.next;
           return nextSite;
        }
        public void allSites(Lattice.Site.Action act) {
            for(Lattice.Site.Linker l=first; l!=null; l=l.next) {act.action(l.site);}
        }
    }
    
    public static final class Site implements Lattice.Site {
        private final Lattice.Site.Iterator adjacentIterator = new AdjacentIterator();
        private final LatticeCubic parent;      //immediate parent lattice (not top-level) on which this site resides
        public Site(LatticeCubic p) {parent = p;}
        public final Lattice lattice() {  //returns the top-level lattice on which this site resides
            LatticeCubic lat = parent;
            while(lat.parentLattice() != null) {lat = lat.parentLattice();}
            return lat;
        }
        public LatticeCubic parent() {return parent;}
        public Lattice.Site.Iterator adjacentIterator() {return adjacentIterator;}
        public Lattice.Coordinate coordinate() {
            Coordinate coord = new LatticeCubic.Coordinate(lattice().D());
            int i=0;
            for(LatticeCubic p=parent(); p!=null; p=p.parentLattice()) {coord.setIndex(i++,p.index());}
            return coord;
        }
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
        public String toString() {
            String string = "Site: ";
            for(LatticeCubic lat=parent; lat!=null; lat=lat.parentLattice()) {
                string += (String.valueOf(lat.index()) + ", ");
            }
            return string;
        }    
    }
    
    /**
     * Main method to demonstrate use of LatticeCubic and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for LatticeCubic");
    }
}