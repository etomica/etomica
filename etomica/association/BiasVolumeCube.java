package etomica.association;

import etomica.*;

public class BiasVolumeCube extends BiasVolume {
    
    private final Space.Vector dimensions;
    private final Space.Vector work;
    
    BiasVolumeCube(Space space){
        super(space);
        dimensions = space.makeVector();
        work = space.makeVector();
        dimensions.E(2.5);
    }
    
    public void setBiasCubeDimensions(Space.Vector dim) {
        dimensions.E(dim);
    }
    public void setBiasCubeDimensions(double d) {dimensions.E(d);}
    public Space.Vector getBiasCubeDimensions() {return dimensions;}
    
    public double biasVolume() {
        double prod = 1.0;
        for(int i=0; i<space.D(); i++) {
            prod *= dimensions.component(i);
        }
        return prod;
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(Atom atom1, Atom atom2) {
        work.setRandomCube();
        work.TE(dimensions);
        atom2.coord.translateTo(atom1.coord.position());
        atom2.coord.translateBy(work);
    }

    /**
     *Function to check for bonding
     */
    
    boolean isAssociated(Atom atom1, Atom atom2){
    
        work.E(atom2.coord.position());
        work.ME(atom1.coord.position());
        work.DE(dimensions);
        return work.abs().max() < 1.0;
    }
}
