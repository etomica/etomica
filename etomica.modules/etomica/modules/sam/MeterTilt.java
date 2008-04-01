package etomica.modules.sam;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.Angle;

/**
 * Meter that measures the average tilt angle (not the angle of average tilt!)
 *
 * @author Andrew Schultz
 */
public class MeterTilt extends DataSourceScalar {

    public MeterTilt(Space space, ISpecies species) {
        super("Tilt", Angle.DIMENSION);
        this.species = species;
        dr = space.makeVector();
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }

    public double getDataAsScalar() {
        IAtomSet molecules = box.getMoleculeList(species);
        int nMolecules = molecules.getAtomCount();
        int leafCount = species.getNumLeafAtoms();
        double thetasum = 0;
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(i);
            dr.E(((IAtomPositioned)molecule.getChildList().getAtom(leafCount-1)).getPosition());
            dr.ME(((IAtomPositioned)molecule.getChildList().getAtom(1)).getPosition());
            double costheta = dr.x(1)/Math.sqrt(dr.squared());
            thetasum += Math.acos(costheta);
        }
        return thetasum / nMolecules;
    }

    private static final long serialVersionUID = 1L;
    protected final ISpecies species;
    protected IBox box;
    protected final IVector dr;
}
