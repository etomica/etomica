package simulate;

public class LatticeCubic {
    
    private final Lattice.Site.Iterator iterator;
    protected LatticeCubic[] rows;  //not private because of possible problem in sharing with Iterator inner class
    private final boolean singlet;
    private final int nRows, siteCount;
    
    //modify to construct lattice from arbitrary site objects
    public LatticeCubic(int[] dim) {
        if(dim.length==0) {  //recursion termination; make lattice of a single site
            singlet = true;
            iterator = new SingletIterator(new Site());
            siteCount = 1;
            nRows = 0;
        }
        else {  //make array of lower-dimensional lattices
            singlet = false;
            int n = dim.length-1;
            int[] newdim = new int[n];
            System.arraycopy(dim,0,newdim,0,n);  //create new dim array from all but last element of given dim array
            nRows = dim[n];
            rows = new LatticeCubic[nRows];
            rows[0] = new LatticeCubic(newdim);
            for(int i=1; i<nRows; i++) {
                rows[i] = new LatticeCubic(newdim);
                rows[i-1].iterator().reset();
                rows[i].iterator().reset();
                while (rows[i].iterator().hasNext()) {
                    Lattice.Site pSite = rows[i-1].iterator().next(); 
                    Lattice.Site nSite = rows[i].iterator.next();
                    ((NeighborIterator)pSite.neighborIterator()).setFirst(nSite);
                    ((NeighborIterator)nSite.neighborIterator()).setFirst(pSite);
                }
            }
            //This does the periodic boundary
            rows[0].iterator().reset();
            rows[nRows-1].iterator().reset();
            while (rows[0].iterator().hasNext()) {
                Lattice.Site pSite = rows[0].iterator().next(); 
                Lattice.Site nSite = rows[nRows-1].iterator.next();
                ((NeighborIterator)pSite.neighborIterator()).setFirst(nSite);
                ((NeighborIterator)nSite.neighborIterator()).setFirst(pSite);
            }
            iterator = new LatticeIterator();
            siteCount = nRows*rows[0].siteCount();
        }
    }
    public int siteCount() {return siteCount;}              
    public Site site(Coordinate coord) {return null;}  //not done
    public Site randomSite() {return null;}           //not done
    public Lattice.Site.Iterator iterator() {return iterator;}     //iterator for all sites in lattice
//    public abstract void draw(Graphics g, int[] origin, double scale);
    
    public class Coordinate {}  //not done
    
    public final class LatticeIterator implements Lattice.Site.Iterator {  //might instead do this by creating a big array of all sites and loop through it
        private boolean hasNext;
        private int iRow;
        private Lattice.Site.Iterator current;
        public LatticeIterator() {reset();}
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
    }
        
    //Iterator for a lattice that contains only one site (zero-D lattice)
    public final class SingletIterator implements Lattice.Site.Iterator {
        private final Site site;
        private boolean hasNext;
        public SingletIterator(Site s) {site = s; hasNext = true;}
        public boolean hasNext() {return hasNext;}
        public void reset() {hasNext = true;}
        public Lattice.Site next() {hasNext = false; return site;}
    }
    
    public final class NeighborIterator implements Lattice.Site.Iterator {
        private Site site;
        private boolean hasNext;
        public void setFirst(Lattice.Site s) {}  //not done
        public boolean hasNext() {return hasNext;}
        public void reset() {hasNext = true;}
        public Lattice.Site next() {hasNext = false; return site;}  //not done
    }
    
    public class Site implements Lattice.Site {
        private int neighborCount;
        private final Lattice.Site.Iterator neighborIterator = new NeighborIterator();
        public final Lattice.Site.Iterator neighborIterator() {return neighborIterator;}
        public final int neighborCount() {return neighborCount;}
    }
}