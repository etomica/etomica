
package etomica.chem.models.water;
import etomica.*;

public class ConfigurationTIP4PWater extends Configuration {

    private double bondLengthOH = 0.9572;
    private double angleHOH = 104.52*Math.PI/180.;
    private double bondLengthOcharge = 0.15;
    
    public ConfigurationTIP4PWater() {
        super();
    }
    
    public void initializePositions(AtomIterator[] iterators){
        if(iterators == null || iterators.length == 0) return;
        
        AtomIterator iterator = iterators[0];
        double x = 0.0;
        double y = 0.0;
        
        iterator.reset();
        
        Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

        Atom charge = iterator.nextAtom();
        charge.coord.position().E(new double[] {x+bondLengthOcharge*Math.cos(angleHOH/2), y+bondLengthOcharge*Math.sin(angleHOH/2), 0.0});
    }//end of initializePositions
    
    
}






