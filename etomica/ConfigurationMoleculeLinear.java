package simulate;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeLinear extends ConfigurationMolecule {
    
    private double bondLength = 0.02;
    private PhaseSpace.Vector orientation;
    
    public ConfigurationMoleculeLinear(){
    }
  
    public void setBondLength(double b) {
        bondLength = b;
        computeDimensions();
    }
    public double getBondLength() {return bondLength;}
    
    public void setOrientation(PhaseSpace.Vector e) {orientation.E(e);}
  
  /**
   * Sets all atoms coordinates to lie on a straight line along the x-axis, with the
   * center of mass unchanged from the value before method was called
   */
    public void initializeCoordinates(Molecule m) {
        PhaseSpace.Vector OldCOM = m.parentSpecies.parentPhaseSpace.makeVector();
        OldCOM.E(m.COM());
        double xNext = 0.0;
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
            a.coordinate.translateTo(OldCOM);  //put all atoms at same point
            a.coordinate.translateToward(orientation,xNext);  //move xNext distance in direction orientation
            xNext += bondLength;
        }
        m.coordinate.translateTo(OldCOM);  //shift molecule to original COM
    }
    
    protected void computeDimensions() {
        dim[1] = 0.0;
        Molecule m = parentSpecies().firstMolecule();  //a typical molecule
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
            dim[1] = Math.max(dim[1], ((AtomType.Disk)a.type).diameter());  //width is that of largest atom
        }
        dim[0] = 0.5*(((AtomType.Disk)m.firstAtom().type).diameter() + ((AtomType.Disk)m.lastAtom().type).diameter()) + (m.nAtoms-1) * bondLength;
    }

  public void setParentSpecies(Species s) {
    parentSpecies = s;
    orientation = (PhaseSpace2D.Vector)s.parentPhaseSpace.makeVector();   //temporary
    ((PhaseSpace2D.Vector)orientation).x = 1.0;    //temporary;  fix orientation to x-axis of 2-D space
  }


}
