package simulate;
import java.awt.*;
import simulate.units.Dimension;

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
        this(Simulation.instance);
    }
    public SpeciesDisks(Simulation sim) {
        super(sim);
        atomsPerMolecule = 1;
        if(sim.space().D() == 1) {
            protoType = new AtomType.Rod(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE);
        }
        else {
            protoType = new AtomType.Disk(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE);
        }
        nMolecules = Default.MOLECULE_COUNT;        
        moleculeConfiguration = new Molecule.Configuration.Linear(this);
    }
    
    public SpeciesDisks(int nM, int nA) {
        this(Simulation.instance, nM, nA);
    }
    public SpeciesDisks(Simulation sim, int nM, int nA) {
        this(sim, nM, nA, new AtomType.Disk(Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE));
    }
              
    public SpeciesDisks(Simulation sim, int nM, int nA, AtomType.Disk type) {
        super(sim);
        protoType = type;
        atomsPerMolecule = nA;
        nMolecules = nM;
            
        moleculeConfiguration = new Molecule.Configuration.Linear(this);
    }
    
            
    protected Molecule makeMolecule(Phase phase) {
        return new Molecule(this, phase, protoType, atomsPerMolecule);
    } 
              
    // Exposed Properties
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}
}


