package simulate;
import java.awt.*;

public class SpeciesDisks extends Species {

    public AtomType.Disk protoType;
              
    /**
    * Default constructor.  Creates species containing 20 molecules, each with 1 disk atom.
    */
    public SpeciesDisks() {
        this(20,1);
    }
              
    public SpeciesDisks(int nM, int nA) {
        setSpeciesIndex(0);       //would like to have this set automatically, based on number of species added
        protoType = new AtomType.Disk(1.0, Color.black, 0.1);  // arguments are mass, color, diameter
        atomsPerMolecule = nA;
        setNMolecules(nM);
            
        colorScheme = new ColorSchemeNull();
        this.add(new ConfigurationMoleculeLinear());
    }
              
    public Molecule makeMolecule() {
        return new Molecule(this, protoType, atomsPerMolecule);
    } 
              
    // Exposed Properties
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}
}


