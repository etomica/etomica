package etomica.gui;

import etomica.*;
import java.awt.*;

/**
 * Species in which molecules are made of arbitrary number of disks (same number for all molecules, though) 
 * with each disk having the same mass and size (same type).
 */
public class SpeciesUserDefinedDisks extends etomica.Species {
    public AtomType[] typeArray; 
//  The atomType is not declared final here becuase it makes setting up the constructors easier,
//  but effectively it cannot be changed once initialized; the instance is passed to the atoms, where
//  it is declared final
//  Note that the parameters of the type can be changed; only the instance of it is frozen once the atoms are made
//    (this is the same behavior as declaring it final)
    public AtomType.Disk currentType;
        
    public SpeciesUserDefinedDisks(AtomType[] types) {
        this(Simulation.instance, types);
    }
    
    public SpeciesUserDefinedDisks(Simulation sim, AtomType[] types) {
        super(sim);
        typeArray = ((AtomType[])types.clone());
        atomsPerMolecule = typeArray.length;
        currentType = (AtomType.Disk)typeArray[0];
        nMolecules = Default.MOLECULE_COUNT;        
        moleculeConfiguration = new Molecule.Configuration.Linear(this);
    }
              
    protected Molecule makeMolecule(Phase phase) {
        return new Molecule(this, phase, typeArray);
    } 
              
    // Exposed Properties
    public final double getMass() { return currentType.mass(); }
    public final void setMass(double mass) { currentType.setMass(mass); }
                
    public final double getDiameter() { return currentType.diameter(); }
    public void setDiameter(double d) { currentType.setDiameter(d); }
                    
    public final Color getColor() { return currentType.color(); }
    public final void setColor(Color c) { currentType.setColor(c); }
}


