package simulate;
import java.awt.*;

//Extension of Atom that informs parent molecule of its movement by updating
//the molecule's center of mass

public class AtomInformer extends Atom {
    
    public AtomInformer(Molecule parent, int index) {
        super(parent, index);
    }
    
    public AtomInformer(Molecule parent, int index, double mass, double diameter) {
        super(parent,index,mass,diameter);
    }
    
    //Expect a problem with PBC for multi-atom molecule
    
    public final void setR(int i, double dr) {
        r[i] = dr;
        parentMolecule.updateCOM();
    }
    public final void translate(double[] dr) {
        Space.uPEv1(r,dr);
        parentMolecule.updateCOM(COMFraction, dr);
    }
    public final void translate(int i, double dr) {
        r[i] += dr;
        parentMolecule.updateCOM(COMFraction, i, dr);
    }
    
}