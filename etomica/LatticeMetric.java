package simulate;

/**
 * Metric lattice, i.e., a lattice with distances defined between the sites
 */
public class LatticeMetric implements Lattice {
    
    protected double neighborIndexCutoff;
    
    public LatticeMetric(Lattice lat) {
        lattice = lat;
        D = lat.D();
        basis = new double[D][D];
    }
  //Lattice interface methods
  
    public int D() {return D;}           
    public int siteCount() {return lattice.siteCount();} 
    public int coordinationNumber() {return lattice.coordinationNumber();}    
    public Site site(Coordinate coord) {return }
    public Site randomSite() {return lattice.randomSite();}            
    public Site.Iterator iterator() {return lattice.iterator();}     
    public double r2(Site s1, Site s2) {
    }
    
    
    public final void setNeighborIndexCutoff(double c) {neighborIndexCutoff = c; setupNeighbors();}
    public final double getNeighborIndexCutoff() {return neighborIndexCutoff;}
    
    public abstract int[] dimensions();
    public abstract void setupNeighbors();
    
    public abstract double[][] getBasis();
    public abstract void setBasis(double[][] b);
        
    public interface Site {
        public Occupant first();
        public void setFirst(Occupant o);
        
//        public Linker firstUpNeighbor();
//        public Linker firstDownNeighbor();
        public int[] coordinate();
        public void setCoordinate(int[] i, double[][] basis);
        public double neighborIndex(Site s, int[] d);
    }
    
    public class Site implements Lattice.Site {
        public Lattice.Site nativeSite;
        public double[] position();
    }
    
    public interface Cell extends Point {
        public double[][] vertices();
        public double volume();
    }

}