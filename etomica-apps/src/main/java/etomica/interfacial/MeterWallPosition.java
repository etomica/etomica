package etomica.interfacial;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;

public class MeterWallPosition extends DataSourceScalar {

    protected final Box box;
    protected final ISpecies wallSpecies;
    protected final double zShift;
    
    public MeterWallPosition(Box box, ISpecies wallSpecies, double zShift) {
        super("Wall position", Length.DIMENSION);
        this.box = box;
        this.wallSpecies = wallSpecies;
        this.zShift = zShift;
    }

    public double getDataAsScalar() {
        return box.getMoleculeList(wallSpecies).getMolecule(0).getChildList().getAtom(0).getPosition().getX(2) - zShift;
    }

}
