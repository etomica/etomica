package simulate;

public abstract class Lattice {
    
    protected double neighborIndexCutoff;
    
    public Lattice() {
        setNeighborIndexCutoff(1.0);
    }
    
    public final void setNeighborIndexCutoff(double c) {neighborIndexCutoff = c; setupNeighbors();}
    public final double getNeighborIndexCutoff() {return neighborIndexCutoff;}
    
    public abstract int[] dimensions();
    public abstract int siteCount();
    public abstract void setupNeighbors();
    
    public abstract double[][] getBasis();
    public abstract void setBasis(double[][] b);
        
    public interface Site {
        public Space.AtomCoordinate firstAtom();
        public void setFirstAtom(Space.AtomCoordinate c);
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