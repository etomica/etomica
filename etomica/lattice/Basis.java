package etomica.lattice;
import etomica.Space;
import java.util.Observer;
import java.util.Observable;

/**
 * Class that defines and creates the basis for a crystal lattice.  Each basis instance has features of both a 
 * site and a lattice---it is treated as a site on a Bravais lattice, but it also comprises the 
 * sites that make up the basis.  Thus the basis is instantiated as a SiteLattice object (an inner 
 * class of Basis), and is placed as a Site on a BravaisLattice.  In creating SiteLattice objects the
 * Basis class implements the SiteFactory interface.<br>
 * The constructor of the SiteLattice object creates the sites of the basis.  The number of sites created
 * and their position relative to a Bravais lattice site is defined by the basisPosition array.
 * The type of sites used to make up the basis is defined by a SiteFactory object, which defaults
 * to Site.Factory if one is not specified.<br>
 * A Lattice object combines a BravaisLattice and a Basis to form a crystal lattice.
 *
 * @author David Kofke
 */

public class Basis implements SiteFactory, java.io.Serializable {

    /** Vectors describing locations of the atoms relative to the bravais lattice site
     *  Each element of the array describes the position of one site of the basis.
     */
    private Space.Vector[] basisPosition;
    /**
     * Object used by this basis to create its sites
     */
    private SiteFactory siteFactory;
    /**
     * Dimension of the space in which this basis resides.
     */
    private int D;
        
    /**
    * Basis of one site (of type Site), placed at the origin in a space of dimension D
    */
    public Basis(int D) {
        this(new Space.Vector[] {Space.makeVector(D)});
    }
    /**
     * Basis of one site (created by given SiteFactory), placed at the origin in a space of dimension D.
     */
    public Basis(int D, SiteFactory factory) {
        this(new Space.Vector[] {Space.makeVector(D)},factory);
    }
    /**
     * Basis of several sites (each of type Site), placed as specified by the position array
     */
    public Basis(Space.Vector[] position) {
        this(position, new Site.Factory());
    }
    /**
     * Basis of several sites (created by the given SiteFactory), placed as specified by the position array
     */
    public Basis(Space.Vector[] position, SiteFactory factory) {
        D = position[0].length();
        basisPosition = position;
        siteFactory = factory;
    }
    
    /**
     * Number of sites in the basis.
     */
    public int siteCount() {return basisPosition.length;}
    
    /**
     * Specifies the SiteFactory object that will create the sites that form each instance of the basis.
     */
    public void setSiteFactory(SiteFactory factory) {siteFactory = factory;}
    
    /**
     * Observable anonymous inner class that notifies interested objects if the basis vectors change.
     * The observer would then update the sites created with this basis.  Used by the Lattice class.
     */
    private Observable monitor = new Observable() {
        public void notifyObservers() {this.notifyObservers(null);}
        public void notifyObservers(Object obj) {
            setChanged();
            super.notifyObservers(obj);
        }
        public void addObserver(Observer o) {if(o != null) super.addObserver(o);}
    };
    /**
     * Registers the given Observer so that it is notified any time the basis vectors are modified.
     */
    public void addObserver(Observer o) {monitor.addObserver(o);}
    /**
     * Notifies observers that the basis vectors are modified.
     */
    protected void notifyObservers() {monitor.notifyObservers();}
    
    /**
     * Sets the array of space vectors that describe the number and position of all sites in the basis.
     * Number of sites is given by the length of the array.
     */
    public void setBasisPosition(Space.Vector[] positions) {
        basisPosition = positions;
        notifyObservers();
    }
    /**
     * Sets the i<sup>th</sup> element of the basis position array to the given value.
     * 
     * @param i index of the position to be modified (with zero as first index). Must be less than position.length-1.
     * @param position the new position for vector i.
     */
    public void setBasisPosition(int i, Space.Vector position) {
        if(i < 0 || i >= basisPosition.length) throw new IllegalArgumentException("Out-of-bounds index in Basis.setBasisPosition");
        basisPosition[i] = position;
        notifyObservers();
    }
    /**
     * Sets the xyzIndex component of the positionIndex element of the basis position array to the given value.
     * 
     * @param positionIndex index of the position to be modified (with zero as first index). Must be less than position.length-1.
     * @param xyzIndex describing the spatial coordinate being modified (with zero as the first index).
     * @param value the new value of the coordinate.
     */
    public void setBasisPosition(int positionIndex, int xyzIndex, double value) {
        if(positionIndex < 0 || positionIndex >= basisPosition.length) throw new IllegalArgumentException("Out-of-bounds index in Basis.setBasisPosition");
        if(xyzIndex < 0 || xyzIndex >= D) throw new IllegalArgumentException("Out-of-bounds index in Basis.setBasisPosition");
        basisPosition[positionIndex].setComponent(xyzIndex,value);
        notifyObservers();
    }
    
    /**
    * Returns a new Basis.SiteLattice, 
    * which becomes a site on a BravaisLattice, and also forms a lattice of other sites that make up the basis.
    */
    public Site makeSite(AbstractLattice parent, AbstractLattice.Coordinate coord) {
        return new SiteLattice((BravaisLattice)parent, (BravaisLattice.Coordinate)coord);
    }

    /**
     * Class that forms an instance of a basis on a Bravais lattice.
     * Object of this type is a site of the Bravais lattice, but is also a lattice itself, 
     * holding the instance of the sites that form the basis on the Bravais lattice site.
     */
    public class SiteLattice extends Site implements AbstractLattice {
        
        /**
        * Coordinate object of bravais site hosting this basis instance.
        * Same as superclass coordinate field, but cast to BravaisLattice.Coordinate.
        */
        BravaisLattice.Coordinate bravaisCoordinate;
        /**
        * Iterator that gives sites in the basis.
        */
        SiteList basisList = new SiteList();
        SiteIteratorList basisIterator = new SiteIteratorList(basisList);
        
        /**
         * Constructor places instance as a site on a parent BravaisLattice and creates
         * sites that form the instance of the basis.
         */
        public SiteLattice(BravaisLattice parent, BravaisLattice.Coordinate coord) {
            super(parent, coord);
            bravaisCoordinate = coord;
            for(int i=0; i<basisPosition.length; i++) {
                Site newSite = siteFactory.makeSite(this,new Coordinate(i));
                basisList.add(newSite);
            }
       //     updateCoordinates();
        }
        //methods implementing AbstractLattice
        public int siteCount() {return basisIterator.size();}
        public Site site(AbstractLattice.Coordinate coord) {return basisIterator.list().get(((Coordinate)coord).basisIndex);}
        public int D() {return D;}
        public Site randomSite() {return basisIterator.list().getRandom();}         
        public SiteIterator iterator() {return basisIterator;} 
        
        /**
         * Updates all site coordinates for the instance of the basis.
         */
 /*       public void updateCoordinates() {
            basisIterator.reset();
            while(basisIterator.hasNext()) {
                ((Coordinate)basisIterator.next().coordinate()).update();
            }
        }
        /**
        * Coordinate for a site generated by the basis.
        * Position of the site is determined from the bravais-lattice coordinate and the basis.
        */
        public class Coordinate implements AbstractLattice.PositionCoordinate {
            final Space.Vector r;
            final int basisIndex;
            public Coordinate(int i) { //put coord in Basis
                basisIndex = i;
                r = Space.makeVector(D);
   //             update();
            }
            /**
             * Computes the position of the site using the bravais-lattice position and the appropriate basis vector
             */
            public final void update() {
                r.E(bravaisCoordinate.position());
                r.PE(basisPosition[basisIndex]);
            }
            /**
             * Spatial position of the site.
             */
            public Space.Vector position() {return r;}
            /**
             * Writes the site as a string that lists its spatial coordinates.
             */
            public String toString() {
                String value = "( ";
                for(int i=0; i<r.length(); i++) value += r.component(i) + " ";
                value += ") ";
                return value;
            }
        }//end of Basis.SiteLattice.Coordinate
    }//end of Basis.SiteLattice
}