package etomica;
import etomica.lattice.*;

public class ConfigurationFcc extends Configuration {
    
    private AtomIteratorSequential iterator = new AtomIteratorSequential();
    public ConfigurationFcc(Space space) {
        super(space);
    }
    
    //need to revise so that coordinates of given group are
    //initialized, instead of parentphase object
    public void initializeCoordinates(Atom group){
        if(group == null) {return;}
        Phase parentPhase = group.parentPhase();
    
    // Count number of molecules
        int sumOfMolecules = 0;
        for(SpeciesAgent s=parentPhase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            sumOfMolecules += s.moleculeCount();
        }
        if(sumOfMolecules == 0) {return;}
        
        Space3D.Vector[] rLat = new Space3D.Vector[sumOfMolecules];
           rLat = fccLattice(sumOfMolecules);

   // Place molecules     
        int i = 0;
        for(SpeciesAgent s=parentPhase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            iterator.setBasis(s);
            iterator.reset();
            while(iterator.hasNext()) {
                iterator.next().coord.translateTo(rLat[i]);
                i++;
            }
        }
        initializeMomenta(parentPhase.speciesMaster());
    }
    
 /*   public Space3D.Vector[] fccLattice(int n) { 
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
    }
    */
}