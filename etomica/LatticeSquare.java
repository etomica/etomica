package simulate;

public class LatticeSquare {
    
    Site origin;
    private int rows, columns;
    public static final double[] basis = {1.0, 1.0};
    
    public LatticeSquare() {
        this(20,20);
    }
    
    public LatticeSquare(int nx, int ny) {
        rows = ny;
        columns = nx;
        origin = new Site(0,0);
        makeRow(origin);
        Site firstXOld = origin;
        for(int iy=1; iy<rows; iy++) {
            Site firstX = new Site(0,iy);
            makeRow(firstX);
            joinRows(firstXOld, firstX);
            firstXOld = firstX;
        }
        joinRows(firstXOld, origin);  //periodic boundary
        
        for(Site s=origin; s!=null; s=s.nextSite) {s.makeNeighborList();}
    }
    
    public int nSites() {return rows*columns;}
    
    private void makeRow(Site first) {
        Site last = first;
        for(int i=1; i<columns; i++) {
            last.setE(new Site(i,first.iy));
            last.setNextSite(last.e);
            last = last.e;
        }
        last.setE(first);  //periodic boundary
    }
    
    private void joinRows(Site topFirst, Site bottomFirst) {
        topFirst.w.setNextSite(bottomFirst);//first site of bottom row is nextSite for last site of top row
        Site bottomLast = bottomFirst;
        Site topLast = topFirst;
        for (int i=0; i<columns; i++) {
            bottomLast.setN(topLast);
            bottomLast = bottomLast.e;
            topLast = topLast.e;
        }
    }
    
    private class Site implements LatticeSite {
        private Site n, e, w, s;  //neighbors
        private AtomP atom;
        private Site nextSite;
        private Site[] neighbors;
        private int ix, iy;       //indices describing "position", for painting image to screen
        Site(int i, int j) {ix = i; iy = j;}
        private void setN(Site k) {n = k; k.s = this;}
        private void setE(Site k) {e = k; k.w = this;}
        private void setW(Site k) {w = k; k.e = this;}
        private void setS(Site k) {s = k; k.n = this;}
        private void setNextSite(Site k) {nextSite = k;}
        private void makeNeighborList() {
            neighbors[0] = n;
            neighbors[1] = e;
            neighbors[2] = w;
            neighbors[3] = s;
        }
        public void putAtom(AtomP a) {atom = a;}
        public AtomP atom() {return atom;}
        public LatticeSite[] neighbors() {return neighbors;}
    }
}