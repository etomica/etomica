//see main method in IntegerBravaisLattice for examples of SiteIterator usage.
package simulate.lattice;
import java.util.Random;

public interface SiteIterator extends java.io.Serializable {
    public boolean hasNext();         //true if another iterate in list
    public void reset();              //put iterator back in beginning state
    public Site first();              //return first element of list without affecting state of iterator
    public Site next();               //advance iterator and return next element
    public void allSites(SiteAction act); //method to perform an action sequentially to every site given by iterator
  //consider adding these methods:
    //public int siteCount();
    //public Site randomSite();
    //public Site get(int i);
    public abstract class Cursor implements SiteIterator {
        public abstract boolean hasNext();         //true if another iterate in list
        public abstract void reset();              //put iterator back in beginning state
        public abstract Site first();              //return first element of list without affecting state of iterator
        public abstract Site next();               //advance iterator and return next element
        public void allSites(SiteAction act) {
            reset();
            while(hasNext()) {act.actionPerformed(next());}
        }
    }
        
    /**
     * Generic iterator that permits addition and removal of Sites.
     */
       public static final class List implements SiteIterator {
           public static Random random;  //random-number generator for selecting a random site
           static {random = new Random();}
           
           private Site.Linker first, next;
           private int siteCount;
           private boolean hasNext;
           public List() {hasNext = false;}
           /**
            * Makes a list iterator that duplicates the set of sites output by another iterator
            */
           public List(SiteIterator iterator) {
               iterator.reset();
               first = iterator.hasNext() ? new Site.Linker(iterator.next(),null) : null;
               Site.Linker last = first;
               while(iterator.hasNext()) {
                  last.next = new Site.Linker(iterator.next(), null);
                  last = last.next;
               }
               reset();
           }
               
           public boolean hasNext() {return hasNext;}
           public void reset() {
              next = first;
              hasNext = (next != null);
           }
           public Site first() {return (first == null) ? null : first.site;}
           public Site next() { //does not check that next is non-null
              Site site = next.site;
              next = next.next;
              hasNext = (next != null);
              return site;
           }
           public int siteCount() {return siteCount;}
            /**
            * Returns a random site in the lattice
            * Iteratively chooses a row at random and calls randomSite for that row
            */
            public Site randomSite() {
                int i = (int)Math.floor(siteCount*random.nextDouble());
                return get(i);
            }
            
            public Site get(int i) {
                if(i > siteCount-1) return null;
                Site.Linker link = first;
                for(int j=0; j<i; j++) {link = link.next;}
                return link.site;
            }

           public void allSites(SiteAction act) {
              for(Site.Linker link=first; link!=null; link=link.next) {
                  act.actionPerformed(link.site);
              }
           }
           /**
            * Adds an site to the set of sites given by this iterator
            */
           public void addSite(Site s) {
              for(Site.Linker link=first; link!=null; link=link.next) {if(link.site == s) return;} //check that site is not already in list
              first = new Site.Linker(s,first);
           }
           //will someday add a removeSite method
           
           public Cursor makeCursor() {return new Cursor();}
           
           public class Cursor extends SiteIterator.Cursor {
              private Site.Linker cursor;
              public boolean hasNext() {return cursor != null;}
              public Site next() {
                 if(cursor == null) return null;
                 Site site = cursor.site;
                 cursor = cursor.next;
                 return site;
              }
              public Site first() {return (first == null) ? null : first.site;}
              public void reset() {cursor = first;}
           }
           
       }//end of SiteIterator.List
     

    //Iterator to generate adjacent sites on LatticeBravais
    
    //rewrite this using SiteIterator.List
    public static class Neighbor implements SiteIterator {
        private Site.Linker firstUp, firstDown, nextLink;
        private boolean doBoth = false;
        private int count;
        private Site site;
        public Neighbor() {;}
        public Neighbor(Site s) {site = s;}
        public void setSite(Site s) {site = s;}
        public Site site() {return site;}
        public int neighborCount() {return count;}
        public boolean hasNext() {return nextLink != null;}
        public void clearAll() {firstUp = null; firstDown = null; nextLink = null; count = 0;}
        public void reset() {nextLink = firstUp; doBoth = true;}
        public void resetUp() {nextLink = firstUp; doBoth = false;} //put iterator in state ready to traverse up-neighbors
        public void resetDown() {nextLink = firstDown; doBoth = false;} //put iterator in state ready to traverse down-neighbors
        public Site first() {return firstUp.site;}
        public Site firstUp() {return firstUp.site;} //return first element of up-list without affecting state of iterator
        public Site firstDown() {return firstDown.site;} //return first element of down-list without affecting state of iterator
        public Site next() {
            Site nextSite = nextLink.site;
            nextLink = nextLink.next;
            if(nextLink==null && doBoth==true) {
                nextLink = firstDown; 
                doBoth = false;
            }
            return nextSite;
        }
        public void allSitesUp(SiteAction act) {  //method to perform an action sequentially to every up-site given by iterator
            for(Site.Linker l=firstUp; l!=null; l=l.next) {act.actionPerformed(l.site);}
        }
        public void allSitesDown(SiteAction act) { //method to perform an action sequentially to every down-site given by iterator
            for(Site.Linker l=firstDown; l!=null; l=l.next) {act.actionPerformed(l.site);}
        }
        public void allSites(SiteAction act) {
            for(Site.Linker l=firstUp; l!=null; l=l.next) {act.actionPerformed(l.site);}
            for(Site.Linker l=firstDown; l!=null; l=l.next) {act.actionPerformed(l.site);}
        }
        public void addUp(Site s) {
            for(Site.Linker l=firstUp; l!=null; l=l.next) {if(s==l.site) return;} //check that site isn't already in list
            firstUp = new Site.Linker(s, firstUp); count++;
        }
        public void addDown(Site s) {
            for(Site.Linker l=firstDown; l!=null; l=l.next) {if(s==l.site) return;} //check that site isn't already in list
            firstDown = new Site.Linker(s, firstDown); count++;
        }
        public void setNeighbors(SiteIterator iterator, Criterion criterion) {  //set up neighbors according to given criterion
            iterator.reset();
            boolean down = true;
            while(iterator.hasNext()) {              //begin outer loop
                Site s = iterator.next();
                if(s == site) {down = false;}     //subsequent neighbors go in up-list
                else if(criterion.areNeighbors(s,site)) {
                    if(down) {addDown(s);}
                    else     {addUp(s);}
                }
            }
        }//end of SiteIterator.Neighbor.setNeighbors
        
        public Cursor makeCursor() {return new Cursor();}
        
        public class Cursor extends SiteIterator.Cursor {
            private Site.Linker cursor = null;
            private boolean doBoth;
            public boolean hasNext() {return cursor != null;}
            public void reset() {cursor = firstUp; doBoth = true;}
            public void resetUp() {cursor = firstUp; doBoth = false;} //put iterator in state ready to traverse up-neighbors
            public void resetDown() {cursor = firstDown; doBoth = false;} //put iterator in state ready to traverse down-neighbors
            public Site first() {return firstUp.site;}
            public Site firstUp() {return firstUp.site;} //return first element of up-list without affecting state of iterator
            public Site firstDown() {return firstDown.site;} //return first element of down-list without affecting state of iterator
            public Site next() {
                Site nextSite = cursor.site;
                cursor = cursor.next;
                if(cursor==null && doBoth==true) {
                    cursor = firstDown; 
                    doBoth = false;
                }
                return nextSite;
            }
        }
        
        
        /**
         * Defines a criterion for specifying whether two sites on a lattice are to be designated as neighbors of each other.
         * Used by the SiteIterator.Neighbor to construct an iterator that will return all "neighbors" of a given site.
         * The "neighbors" returned by the iterator are those for which the areNeighbors method of this class returns 
         * <code>true</code>.
         */
        public interface Criterion {
            public boolean areNeighbors(Site s1, Site s2);
             
            /**
             * Criterion that defines all sites on the lattice to be neighbors of each other.
             */
            public class All implements Criterion {
                public boolean areNeighbors(Site s1, Site s2) {return s1 != s2;}
            }//end of Criterion.All
        }//end of Neighbor.Criterion
    }//end of SiteIterator.Neighbor
}//end of SiteIterator
                
