package simulate;

public abstract class Lattice {
        
    public abstract int siteCount();              //total number of sites on lattice
    public abstract Site site(Coordinate coord);  //obtain a site given some specification in coordinate
    public abstract Site randomSite();            //get a random site
    public abstract Site.Iterator iterator();     //iterator for all sites in lattice
//    public abstract void draw(Graphics g, int[] origin, double scale);
    
    public static abstract class Coordinate {}
    
    public interface Site {
        public Iterator neighborIterator();
        
        public interface Iterator {
            public boolean hasNext();
            public void reset();
            public Site next();
        }
    }
}