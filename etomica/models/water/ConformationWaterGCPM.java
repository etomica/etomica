package etomica.models.water;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.config.Conformation;
import etomica.space.Space;

public class ConformationWaterGCPM extends Conformation {

    private double bondLengthOH = 0.9572;
    private double angleHOH = 104.52*Math.PI/180.;
    private double rOM=0.27;

    public ConformationWaterGCPM(Space space) {
        super(space);
    }
    
    public void initializePositions(IAtomSet list){
        
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        
        
        IAtomPositioned o = (IAtomPositioned)list.getAtom(SpeciesWater4P.indexO);
        o.getPosition().E(new double[] {x, y, z});
               
        IAtomPositioned h1 = (IAtomPositioned)list.getAtom(SpeciesWater4P.indexH1);
        h1.getPosition().E(new double[] {x-bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y, z+bondLengthOH*Math.sin((Math.PI-angleHOH)/2)});
                
        IAtomPositioned h2 = (IAtomPositioned)list.getAtom(SpeciesWater4P.indexH2);
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y, z+bondLengthOH*Math.sin((Math.PI-angleHOH)/2)});
        
        IAtomPositioned m = (IAtomPositioned)list.getAtom(SpeciesWater4P.indexM);
        m.getPosition().E(new double[] {x, y, z+rOM});

        
/*        Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x-bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y+bondLengthOH*Math.sin((Math.PI-angleHOH)/2), 0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y+bondLengthOH*Math.sin((Math.PI-angleHOH)/2), 0});
        
        Atom m = iterator.nextAtom();
        m.coord.position().E(new double[] {x, y+rOM, 0});
*/
        
 /*       Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});
        
        Atom m = iterator.nextAtom();
        m.coord.position().E(new double[] {x+rOM*Math.cos(angleHOH/2.0), y+rOM*Math.sin(angleHOH/2.0), 0.0});
*/
    }//end of initializePositions
    
    
}
