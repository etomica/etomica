package simulate;
import java.util.*;

public class MoleculeAtomic extends Molecule {

    public MoleculeAtomic() {
        Space.uEa1(r,0.0);
        Space.uEa1(p,0.0);
        this.zeroForce();
        this.zeroDisplacement();
        neighbors = new Vector();
    }
    
    public MoleculeAtomic(int index) {
        this();
        this.speciesIndex = index;
    }
    
    public MoleculeAtomic(int index, double rm) {
        this();
        this.speciesIndex = index;
        this.rm = rm;
    }
    
    public void translate(double[] dr) {
        Space.uPEv1(r,dr);
        Space.uPEv1(displacement,dr);
    }
    
    public void accelerate(double[] dp) {Space.uPEv1(p,dp);}
    public void decelerate(double[] dp) {Space.uMEv1(p,dp);}
    
    public double getKineticEnergy() {return 0.5*rm*Space.v1Dv2(p,p);}
}