package simulate;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeLinear extends ConfigurationMolecule {
    
    private double bondLength = 0.02;
    private Space.Vector orientation;
    
    public ConfigurationMoleculeLinear(){
    }
  
    public void setBondLength(double b) {
        bondLength = b;
        computeDimensions();
    }
    public double getBondLength() {return bondLength;}
    
    public void setOrientation(Space.Vector e) {orientation.E(e);}
  
  /**
   * Sets all atoms coordinates to lie on a straight line along the x-axis, with the
   * center of mass unchanged from the value before method was called
   */
    public void initializeCoordinates(Molecule m) {
        Space.Vector OldCOM = m.parentSpecies.parentSimulation.space.makeVector();
        OldCOM.E(m.COM());
        double xNext = 0.0;
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
            a.translateTo(OldCOM);  //put all atoms at same point
            a.translateBy(xNext,orientation);  //move xNext distance in direction orientation
            xNext += bondLength;
        }
        m.translateTo(OldCOM);  //shift molecule to original COM
    }
    
    protected void computeDimensions() {
        dim[1] = 0.0;
        Molecule m = parentSpecies().makeMolecule();  //a typical molecule
        initializeCoordinates(m);
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
            dim[1] = Math.max(dim[1], ((AtomType.Disk)a.type).diameter());  //width is that of largest atom
        }
        dim[0] = 0.5*(((AtomType.Disk)m.firstAtom().type).diameter() + ((AtomType.Disk)m.lastAtom().type).diameter()) + (m.nAtoms-1) * bondLength;
    }

  public void setParentSpecies(Species s) {
    parentSpecies = s;
    orientation = (Space2DCell.Vector)s.parentSimulation.space.makeVector();   //temporary
    orientation.setComponent(1,1.0);   //fix orientation to x-axis of 2-D space
  }


}
