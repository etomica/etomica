package simulate;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeLinear extends ConfigurationMolecule {
    
    private double bondLength = 0.02;
    
    public ConfigurationMoleculeLinear(){
    }
  
    public void setBondLength(double b) {
        bondLength = b;
        computeDimensions();
    }
    public double getBondLength() {return bondLength;}
  
  /**
   * Sets all atoms coordinates to lie on a straight line along the x-axis, with the
   * center of mass unchanged from the value before method was called
   */
    public void initializeCoordinates(Molecule m) {
        double[] OldCOM = new double[Space.D];
        double[] NewCOM = new double[Space.D];
        Space.uEv1(OldCOM,m.COM());
        Space.uEa1(NewCOM,0.0);
        double xNext = 0.0;
        for(AtomC a=(AtomC)m.firstAtom(); a!=m.terminationAtom(); a=(AtomC)a.getNextAtom()) {
            Space.uEa1(a.r,0.0);   //zero all coordinates
            a.r[0] = xNext;
            xNext += bondLength;
            NewCOM[0] += a.r[0] * a.mass;
        }
        
        Space.uDEa1(NewCOM,(double)m.nAtoms);
        
        for(AtomC a=(AtomC)m.firstAtom(); a!=m.terminationAtom(); a=(AtomC)a.getNextAtom()) {
            Space.uMEv1(a.r,NewCOM);  //zero the molecule center-of-mass
            Space.uPEv1(a.r,OldCOM);  //move com to original position
        }
    }
    
    protected void computeDimensions() {
        dim[1] = 0.0;
        Molecule m = parentSpecies.firstMolecule();  //a typical molecule
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
            dim[1] = Math.max(dim[1], ((AtomDisk)a).getDiameter());  //width is that of largest atom
        }
        dim[0] = 0.5*(((AtomDisk)m.firstAtom()).getDiameter() + ((AtomDisk)m.lastAtom()).getDiameter()) + (m.nAtoms-1) * bondLength;
    }
}
