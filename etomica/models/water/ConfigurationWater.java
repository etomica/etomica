package etomica.models.water;
import etomica.*;

public class ConfigurationWater extends Configuration {

    private double bondLengthOH = 1.0;
    private double angleHOH = 109.5*Math.PI/180.;

    public ConfigurationWater(Simulation sim) {
        super(sim);
    }
    
    public void initializePositions(AtomIterator[] iterators){
        if(iterators == null || iterators.length == 0) return;
        
        AtomIterator iterator = iterators[0];
        double x = 0.0;
        double y = 0.0;
        
        iterator.reset();
        
        Atom o = iterator.next();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = iterator.next();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = iterator.next();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }//end of initializePositions
    
    
}
