package simulate;

/**
 * Metric lattice, i.e., a lattice with distances defined between the sites
 */
public abstract class LatticeMetric extends Lattice {
    
    protected double neighborIndexCutoff;
    
    public LatticeMetric() {
        setNeighborIndexCutoff(2.5);
    }
    
    public final void setNeighborIndexCutoff(double c) {neighborIndexCutoff = c; setupNeighbors();}
    public final double getNeighborIndexCutoff() {return neighborIndexCutoff;}
    
    public abstract int[] dimensions();
    public abstract int siteCount();
    public abstract void setupNeighbors();
    
    public abstract double[][] getBasis();
    public abstract void setBasis(double[][] b);
        
    public interface Occupant {
        public Site site();
    }
    
    public interface Site {
        public Occupant first();
        public void setFirst(Occupant o);
        
//        public Linker firstUpNeighbor();
//        public Linker firstDownNeighbor();
        public int[] coordinate();
        public void setCoordinate(int[] i, double[][] basis);
        public double neighborIndex(Site s, int[] d);
    }
    
    public interface Point extends Site {
        public double[] position();
    }
    
    public interface Cell extends Point {
        public double[][] vertices();
        public double volume();
    }

}