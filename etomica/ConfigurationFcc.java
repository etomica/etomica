package etomica;
import etomica.lattice.*;

public class ConfigurationFcc extends Configuration {
    
    public ConfigurationFcc(Space space) {
        super(space);
    }
    
    public void initializePositions(AtomIterator[] iterators){
        AtomIteratorCompound iterator = new AtomIteratorCompound(iterators);//lump 'em all together
        if(iterator == null || iterator.size() == 0) {return;}
        if(iterator.size() == 1) {
            iterator.reset();
            iterator.next().coord.translateTo(space.origin());
            return;
        }
    
    // Count number of molecules
        int sumOfMolecules = iterator.size();
        if(sumOfMolecules == 0) {return;}
        
        Space3D.Vector[] rLat = new Space3D.Vector[sumOfMolecules];
        rLat = lattice(sumOfMolecules);

   // Place molecules     
        int i = 0;
        iterator.reset();
        while(iterator.hasNext()) {
            iterator.next().coord.translateTo(rLat[i++]);
        }
    }//end of initializePositions
    
    public static Space3D.Vector[] lattice(int n) { 
        Space3D.Vector[] r = new Space3D.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space3D.Vector();}
        LatticeFCC fcc = new LatticeFCC(n, Default.BOX_SIZE);
        SiteIterator iteratorsites = fcc.iterator();
        iteratorsites.reset();
        int i = 0;
        while (iteratorsites.hasNext()&& i < n){
            Site site = iteratorsites.next();
            r[i].E(((AbstractLattice.PositionCoordinate)site.coordinate()).position());
            i++ ;
        }
        return r;
    }//end of lattice
}//end of ConfigurationFcc