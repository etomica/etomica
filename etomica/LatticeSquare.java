package simulate;
import java.awt.Color;

public class LatticeSquare extends Lattice {
    
    public Site origin;
    public double[][] basis = {{1.0,0.0}, {0.0,1.0}};
    private int[] dimensions;
    private Site[][] sites;

    public LatticeSquare(Class siteType, int[] dim) {
        dimensions = dim;
        basis[0][0] = 1.0/(double)dimensions[0];    //square lattice
        basis[1][1] = 1.0/(double)dimensions[1];
        Site last = null;
        sites = new Site[dimensions[0]][dimensions[1]];
        for(int j=0; j<dimensions[1]; j++) {
            for(int i=0; i<dimensions[0]; i++) {
                try {sites[i][j] = (Site)siteType.newInstance();} 
                    catch (InstantiationException e) {}
                    catch (IllegalAccessException e) {}
                sites[i][j].setCoordinate(new int[] {i,j}, basis);
//                sites[i][j] = (Site)site.makeSite(new int[] {i,j});
                if(last != null) last.setNextSite(sites[i][j]);
                last = sites[i][j];
                if(i > 0) {last.setW(sites[i-1][j]);}
                if(i == dimensions[0]-1) {last.setE(sites[0][j]);}
                if(j > 0) {last.setN(sites[i][j-1]);}
                if(j == dimensions[1]-1) {last.setS(sites[i][0]);}
            }
        }
        origin = sites[0][0];
    }
    
    public final void clearCells() {
        for(Site s=origin; s!=null; s=s.nextSite()) {s.setFirst(null);}
    }
    
    public final Lattice.Occupant firstOccupant() {
        Site s = origin;
        Lattice.Occupant a = s.first();
        while(a == null) {
            s = s.nextSite();
            if(s == null) return null;
            a = s.first();
        }
        return a;
    }
    
    public final double[][] getBasis() {return basis;}
    public final void setBasis(double[][] b) {basis = b;}
    
    public final int[] dimensions() {return dimensions;}
    public final int siteCount() {return dimensions[0]*dimensions[1];}
    
    public final Site[][] sites() {return sites;}
    
    public void setupNeighbors() {
        for(Site s1=origin; s1!=null; s1=s1.nextSite()) {
            for(Site s2=s1.nextSite(); s2!=null; s2=s2.nextSite()) {
                if (s1.neighborIndex(s2,dimensions) < neighborIndexCutoff) {
                     s1.firstUpNeighbor = new Linker(s1.firstUpNeighbor,s2);  //first upNeighbor of s1 now points to what was the first upneighbor, and represents s2
                     s2.firstDownNeighbor = new Linker(s2.firstDownNeighbor,s1); //likewise for first downNeighbor of s2;
                }
            }
        }
    }
    
    public static final class Linker {
        public final Linker next;
        public final Site site;
        public Linker(Linker n, Site s) {
            next = n;
            site = s;
        }
        public final Linker next() {return next;}
        public final Site site() {return site;}
    }
    
    public static class Site implements Lattice.Site {
        public Lattice.Occupant first;
        private Linker firstUpNeighbor;
        private Linker firstDownNeighbor;
        public final Color color = Constants.RandomColor();
        private Site nextSite, previousSite;
        private int[] coordinate;
        Site n, e, s, w;
        public Site() {}
        public Site N() {return n;}
        public Site S() {return s;}
        public Site E() {return e;}
        public Site W() {return w;}
        protected void setN(Site k) {n = k; k.s = this;}
        protected void setE(Site k) {e = k; k.w = this;}
        protected void setW(Site k) {w = k; k.e = this;}
        protected void setS(Site k) {s = k; k.n = this;}
        public final Lattice.Occupant first() {return first;}
        public final void setFirst(Lattice.Occupant o) {first = o;}
        public final Linker firstUpNeighbor() {return firstUpNeighbor;}
        public final Linker firstDownNeighbor() {return firstDownNeighbor;}
        public final int[] coordinate() {return coordinate;}
        public void setCoordinate(int[] i, double[][] basis) {coordinate = i;}
        public final Site nextSite() {return nextSite;}
        public final Site previousSite() {return previousSite;}
        public final void setNextSite(Site s) {nextSite = s; if(s!=null) s.previousSite = this;}
        public double neighborIndex(simulate.Lattice.Site ss, int[] dimensions) {
            Site s = (Site)ss;
            int dx = Math.abs(s.coordinate()[0] - coordinate[0]);  
            int dy = Math.abs(s.coordinate()[1] - coordinate[1]);
            if(2*dx > dimensions[0]) dx -= dimensions[0];
            if(2*dy > dimensions[1]) dy -= dimensions[1];
            return (double)(dx*dx + dy*dy);
        }
    }
    
    public static class Point extends Site implements Lattice.Point {
        protected final double[] position = new double[2];
        public void setCoordinate(int[] i, double[][] basis) {
            super.setCoordinate(i,basis);
            position[0] = coordinate()[0]*basis[0][0] + coordinate()[1]*basis[1][0];
            position[1] = coordinate()[0]*basis[0][1] + coordinate()[1]*basis[1][1];
        }
        public final double[] position() {return position;}
        public double neighborIndex(simulate.Lattice.Site s, int[] dimensions) {
            Point p = (Point)s;
            double dx = Math.abs(p.position()[0] - position[0]);
            double dy = Math.abs(p.position()[1] - position[1]);
            dx -= (dx > 0.0) ? Math.floor(dx+0.5) : Math.ceil(dx-0.5);  //PBC, unit-width lattice
            dy -= (dy > 0.0) ? Math.floor(dy+0.5) : Math.ceil(dy-0.5);
            return (double)(dx*dx + dy*dy);
        }
    }
    
    public static final class Cell extends Point implements Lattice.Cell {
        private final double[][] vertices = {{0.5,0.5},{0.5,-0.5},{-0.5,-0.5},{-0.5,0.5}};
        private final int nVertices = 4;
        private double volume;
        public Cell() {}
        public void setCoordinate(int[] i, double[][] basis) {
            super.setCoordinate(i,basis);
            position[0] += 0.5*basis[0][0];     //shift position to center of cell
            position[1] += 0.5*basis[1][1];
            volume = basis[0][0]*basis[1][1];   //square cell
            for(int k=0; k<nVertices; k++) {
                vertices[k][0] *= basis[0][0];
                vertices[k][1] *= basis[1][1];
                vertices[k][0] += position[0];
                vertices[k][1] += position[1];
            }
        }
        
//        protected void setN(Site k) {super.setN(k); vertices[1]=((Cell)k).vertices[0]; vertices[2]=((Cell)k).vertices[3];}
//        protected void setE(Site k) {super.setE(k); vertices[0]=((Cell)k).vertices[3]; vertices[1]=((Cell)k).vertices[2];}
//        protected void setW(Site k) {super.setW(k); vertices[3]=((Cell)k).vertices[0]; vertices[2]=((Cell)k).vertices[1];}
//        protected void setS(Site k) {super.setS(k); vertices[3]=((Cell)k).vertices[2]; vertices[0]=((Cell)k).vertices[1];}
 
        public boolean inCell(Atom a) {
            double x = a.coordinate.position().component(0);
            double y = a.coordinate.position().component(1);
            return 0.99999*vertices[2][0] <= x && x <= 1.000001*vertices[0][0] && 0.99999*vertices[2][1] <= y && y <= 1.00001*vertices[0][1];
        }
        
        //Returns square distances between nearest vertices of the two cells
        public double neighborIndex(simulate.Lattice.Site s, int[] dimensions) {
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
        public double volume() {return volume;}   
    }
}
    
    