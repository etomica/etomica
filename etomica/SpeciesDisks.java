package simulate;
import java.awt.*;

/**
 * Species in which molecules are made of arbitrary number of disks (same number for all molecules, though) 
 * with each disk having the same mass and size (same type).
 */
public class SpeciesDisks extends Species {

//  The atomType is not declared final here becuase it makes setting up the constructors easier,
//  but effectively it cannot be changed once initialized; the instance is passed to the atoms, where
//  it is declared final
//  Note that the parameters of the type can be changed; only the instance of it is frozen once the atoms are made
//    (this is the same behavior as declaring it final)
    public AtomType.Disk protoType;
              
    public SpeciesDisks() {
        setSpeciesIndex(0);
        setAtomsPerMolecule(1);
        if(Simulation.D == 1) {
            protoType = new AtomType.Rod(1.0,Color.black,0.1);
        }
        else {
            protoType = new AtomType.Disk(1.0,Color.black,0.1);  //mass, color, diameter
        }
        colorScheme = new ColorSchemeByType();
        configurationMolecule = new ConfigurationMoleculeLinear();
        setNMolecules(20);
//        this.add(new ConfigurationMoleculeLinear());
    }
              
    public SpeciesDisks(Simulation ps) {
        this(ps,20,1);
    }

    public SpeciesDisks(Simulation ps, int nM, int nA) {
        this(ps, nM, nA, new AtomType.Disk(1.0, Color.black, 0.1));
    }
              
    public SpeciesDisks(Simulation ps, int nM, int nA, AtomType.Disk type) {
        parentSimulation = ps;
        setSpeciesIndex(0);       //would like to have this set automatically, based on number of species added
        protoType = type;
        atomsPerMolecule = nA;
        setNMolecules(nM);
            
        colorScheme = new ColorSchemeByType();
        this.add(new ConfigurationMoleculeLinear());
    }
    
            
    protected Molecule makeMolecule(Phase phase) {
        return new Molecule(this, phase, protoType, atomsPerMolecule);
    } 
              
    // Exposed Properties
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}
}


