package etomica.lattice;

public class LatticeSC extends Lattice {
    
    private static double latticeConstant = 1.0;
    public LatticeSC(){
        super(3, 2, primaryVectors(), basis1());
        //set up neighbors
        SiteIterator.Neighbor iterator = new SiteIterator.Neighbor();
        iterator.reset();
        SiteIterator.Neighbor.Cursor cursor =  iterator.makeCursor();
        Criterion1 criterion1 = new Criterion1() ;
        while(iterator.hasNext()) {
            cursor.reset();
            Site site = iterator.next();
            site.adjacentIterator().setNeighbors(cursor, criterion1);
        }
    }
    
    public static Space3D.Vector [] primaryVectors(){
        Space3D.Vector[] r = new Space3D.Vector[3];
        for(int i=0; i<3; i++) {
            r[i] = new Space3D.Vector();
        }
            r[0].x = 1.0*latticeConstant;
            r[0].y = 0.0;
            r[0].z = 0.0;
            r[1].x = 0.0;
            r[1].y = 1.0*latticeConstant;
            r[1].z = 0.0;
            r[2].x = 0.0;
            r[2].y = 0.0;
            r[2].z = 1.0*latticeConstant;
            return r;
    }
    
    public static Space3D.Vector [] positions(){
        Space3D.Vector[] p = new Space3D.Vector[1];
        for(int i=0; i<1; i++) {
            p[i] = new Space3D.Vector();
        }
            p[0].x = 0.0;
            p[0].y = 0.0;
            p[0].z = 0.0;
            return p;
    }   
    
    
    public static Basis basis1(){
        Site.Factory factory = new Site.Factory();
        Basis b = new Basis(positions(),factory);
        return b;
    }
    
      public static class Criterion1 implements SiteIterator.Neighbor.Criterion{
                Space3D.BoundaryPeriodicSquare periodicBoundary = new Space3D.BoundaryPeriodicSquare(8.0,8.0,8.0);
                public boolean areNeighbors(Site site1,Site site2){
                    double r2 = Space3D.r2(
                        ((BravaisLattice.Coordinate)site1.coordinate()).position(),
                        ((BravaisLattice.Coordinate)site2.coordinate()).position(),
                        periodicBoundary);
                        if(r2 <= latticeConstant){
                            return true;
                        }
                        else
                            {return false;}
                            
                }
      }
    
    
}