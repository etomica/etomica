//see main method in IntegerBravaisLattice for examples of SiteIterator usage.
package etomica.lattice;

public interface SiteIterator extends java.io.Serializable {
    public boolean hasNext();         //true if another iterate in list
    public void reset();              //put iterator back in beginning state
    public Site first();              //return first element of list without affecting state of iterator
    public Site next();               //advance iterator and return next element
    public void allSites(SiteAction act); //method to perform an action sequentially to every site given by iterator
    public int size();
  //consider adding these methods:
    //public Site randomSite();
    //public Site get(int i);
    
}//end of SiteIterator
                
