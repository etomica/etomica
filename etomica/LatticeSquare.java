package simulate;

public class LatticeSquare extends Lattice {
    
    public Site origin;
    public double[][] basis = {{1.0,0.0}, {0.0,1.0}};
    private int columns, rows;
    private Site[][] sites;

    public LatticeSquare(Site site, int nx, int ny) {
        rows = ny;
        columns = nx;
        basis[0][0] = 1.0/(double)nx;    //square lattice
        basis[1][1] = 1.0/(double)ny;
        Site last = null;
        sites = new Site[nx][ny];
        for(int j=0; j<ny; j++) {
            for(int i=0; i<nx; i++) {
                sites[i][j] = (Site)site.makeSite(new int[] {i,j});
                if(last != null) last.setNextSite(sites[i][j]);
                last = sites[i][j];
            }
        }
        origin = sites[0][0];
    }
    
    public final double[][] getBasis() {return basis;}
    public final void setBasis(double[][] b) {basis = b;}
    
    public final int siteCount() {return rows*columns;}
    
    public final Site[][] sites() {return sites;}
    
    public void setupNeighbors() {
        for(Site s1=origin; s1!=null; s1=s1.nextSite()) {
            for(Site s2=s1.nextSite(); s2!=null; s2=s2.nextSite()) {
                if (s1.neighborIndex(s2) < neighborIndexCutoff) {
                     s1.firstUpNeighbor = new Linker(s1.firstUpNeighbor,s2);  //first upNeighbor of s1 now points to what was the first upneighbor, and represents s2
                     s2.firstDownNeighbor = new Linker(s2.firstDownNeighbor,s1); //likewise for first downNeighbor of s2;
                }
            }
        }
    }
    
    public final static class Linker {
        public final Linker next;
        public final Site site;
        public Linker(Linker n, Site s) {
            next = n;
            site = s;
        }
        public final Linker next() {return next;}
        public final Site site() {return site;}
    }
    
    public class Site implements Lattice.Site {
        private Space.Coordinate firstAtom;
        private Linker firstUpNeighbor;
        private Linker firstDownNeighbor;
        private Site nextSite;
        
        private int[] coordinate;
        public Site(int[] i) {setCoordinate(i);}
        public final Lattice.Site makeSite(int[] i) {return new Site(i);}
        public final Space.Coordinate firstAtom() {return firstAtom;}
        public final void setFirstAtom(Space.Coordinate a) {firstAtom = a;}
        public final Linker firstUpNeighbor() {return firstUpNeighbor;}
        public final Linker firstDownNeighbor() {return firstDownNeighbor;}
        public final int[] coordinate() {return coordinate;}
        public final void setCoordinate(int[] i) {coordinate = i;}
        public final Site nextSite() {return nextSite;}
        public final void setNextSite(Site s) {nextSite = s;}
        public double neighborIndex(simulate.Lattice.Site ss) {
            Site s = (Site)ss;
            int dx = Math.abs(s.coordinate()[0] - coordinate[0]);  
            int dy = Math.abs(s.coordinate()[1] - coordinate[1]);
            if(2*dx > columns) dx -= columns;
            if(2*dy > rows) dy -= rows;
            return (double)(dx*dx + dy*dy);
        }
    }
    
    public class Point extends Site implements Lattice.Point {
        private double[] position;
        public Point(int[] i) {
            super(i);
            position[0] = coordinate()[0]*basis[0][0] + coordinate()[1]*basis[1][0];
            position[1] = coordinate()[0]*basis[0][1] + coordinate()[1]*basis[1][1];
        }
        public double[] position() {return position;}
        public double neighborIndex(simulate.Lattice.Site s) {
            Point p = (Point)s;
            double dx = Math.abs(p.position()[0] - position[0]);
            double dy = Math.abs(p.position()[1] - position[1]);
            dx -= (dx > 0.0) ? Math.floor(dx+0.5) : Math.ceil(dx-0.5);  //PBC, unit-width lattice
            dy -= (dy > 0.0) ? Math.floor(dy+0.5) : Math.ceil(dy-0.5);
            return (double)(dx*dx + dy*dy);
        }
    }
    
    public final class Cell extends Point implements Lattice.Cell {
        private final double[][] vertices = {{0.5,0.5},{-0.5,0.5},{0.5,-0.5},{-0.5,-0.5}};
        private int nVertices = 4;
        public Cell(int[] i) {
            super(i);
            for(int k=0; k<nVertices; k++) {
                vertices[k][0] *= basis[0][0];
                vertices[k][1] *= basis[1][1];
                vertices[k][0] += position()[0];
                vertices[k][1] += position()[1];
            }
        }
        
        //Returns square distances between nearest vertices of the two cells
        public double neighborIndex(simulate.Lattice.Site s) {
            Cell c = (Cell)s;
            double r2Min = Double.MAX_VALUE;
            for(int k1=0; k1<nVertices; k1++) {
                for(int k2=0; k2<nVertices; k2++) {
                    double dx = Math.abs(c.vertices()[k1][0] - vertices[k2][0]);
                    double dy = Math.abs(c.vertices()[k1][1] - vertices[k2][1]);
                    dx -= (dx > 0.0) ? Math.floor(dx+0.5) : Math.ceil(dx-0.5);  //PBC, unit-width lattice
                    dy -= (dy > 0.0) ? Math.floor(dy+0.5) : Math.ceil(dy-0.5);
                    r2Min = Math.min(r2Min, dx*dx+dy*dy);
                }
            }
            return r2Min;
        }
        public double[][] vertices() {return vertices;}
        public double volume() {return basis[0][0]*basis[1][1];}   //square cell
    }
}
    
    