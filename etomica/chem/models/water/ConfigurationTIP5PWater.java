package etomica.chem.models.water;
import etomica.Atom;
import etomica.Configuration;
import etomica.Space;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;

public class ConfigurationTIP5PWater extends Configuration {

	private double bondLengthOH = 0.9572;
	private double angleHOH = 104.52*Math.PI/180.;
	private double angleChargeOCharge = 109.47*Math.PI/180;
	private double bondLengthOcharge = 0.70;
    private final AtomIteratorListSimple iterator; 
    
	public ConfigurationTIP5PWater(Space space) {
		super(space);
        iterator = new AtomIteratorListSimple();
	}
    
	public void initializePositions(AtomList list){

		iterator.setList(list);
		double x = 0.0;
		double y = 0.0;
        double z = 0.0;
        
		iterator.reset();
        
		Atom o = iterator.nextAtom();
		o.coord.position().E(new double[] {x, y, 0.0});
               
		Atom h1 = iterator.nextAtom();
		h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
		Atom h2 = iterator.nextAtom();
		h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

		Atom charge1 = iterator.nextAtom();
		charge1.coord.position().E(new double[] {x+bondLengthOcharge*Math.cos(angleChargeOCharge/2)*Math.cos((2*Math.PI-angleHOH)/2), y+bondLengthOcharge*Math.cos(angleChargeOCharge/2)*Math.sin(2*Math.PI-(2*Math.PI-angleHOH)/2), z+bondLengthOcharge*Math.sin(angleChargeOCharge/2)});
	
		Atom charge2 = iterator.nextAtom();
		charge2.coord.position().E(new double[] {x+bondLengthOcharge*Math.cos(angleChargeOCharge/2)*Math.cos((2*Math.PI-angleHOH)/2), y+bondLengthOcharge*Math.cos(angleChargeOCharge/2)*Math.sin(2*Math.PI-(2*Math.PI-angleHOH)/2), z-bondLengthOcharge*Math.sin(angleChargeOCharge/2)});
	
	
	}//end of initializePositions
	
    
    
}

