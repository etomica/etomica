package simulate;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeLinear extends ConfigurationMolecule {
    
    private double bondLength = 0.02;
    private Space.Vector orientation;
    private double[] angle = new double[Simulation.D];
//    private double theta = 45.;
    
    public ConfigurationMoleculeLinear(){
    }
  
    public void setBondLength(double b) {
        bondLength = b;
        computeDimensions();
    }
    public double getBondLength() {return bondLength;}
    
    public void setAngle(int i, double t) {
        angle[i] = Math.PI*t/180.;
        switch(Simulation.D) {
            case 1:
                return;
            case 2:
                setOrientation(new Space2D.Vector(Math.cos(angle[0]),Math.sin(angle[0])));
                return;
//            case 3:
//                setOrientation(new Space3D.Vector(Math.sin(angle[1])*Math.cos(angle[0]),
//                                                  Math.sin(angle[1])*Math.sin(angle[0]),
//                                                  Math.cos(angle[1])));
//                return;
        }
    }
    public double getAngle(int i) {return angle[i];}
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
/*        if(parentSpecies()==null) return;
        dim[1] = 0.0;
        Molecule m = parentSpecies().getMolecule();  //a typical molecule
        initializeCoordinates(m);
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
            dim[1] = Math.max(dim[1], ((AtomType.Disk)a.type).diameter());  //width is that of largest atom
        }
        dim[0] = 0.5*(((AtomType.Disk)m.firstAtom().type).diameter() + ((AtomType.Disk)m.lastAtom().type).diameter()) + (m.atomCount-1) * bondLength;
*/
    }

  public void setParentSpecies(Species s) {
    parentSpecies = s;
    orientation = s.parentSimulation.space.makeVector();   //temporary
//    setAngle(theta);   //fix orientation to x-axis of 2-D space
//    computeDimensions();
  }


}
