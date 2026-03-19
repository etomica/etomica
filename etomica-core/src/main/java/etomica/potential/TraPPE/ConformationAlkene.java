package etomica.potential.TraPPE;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.models.co2.SpeciesTraPPECO2;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Electron;

public class ConformationAlkene implements IConformation, java.io.Serializable{

    public ConformationAlkene(Space space){
        this.space = space;
        vector = space.makeVector();
    }

    public void initializePositions(IAtomList atomList) {

        IAtom n1 = atomList.get(SpeciesTraPPECO2.indexC);
        n1.getPosition().E(new double[] {0, 0, 0});

        IAtom n2 = atomList.get(SpeciesTraPPECO2.indexOleft);
        n2.getPosition().E(new double[] {-bondlength, 0, 0});

        IAtom n3 = atomList.get(SpeciesTraPPECO2.indexOright);
        n3.getPosition().E(new double[] {bondlength, 0, 0});


    }

    public final static double [] Echarge = new double [3];
    static {

        //add + charge to C, - charge to O, what is the unit of the pt charge?

        etomica.models.co2.ConformationCO2.Echarge[SpeciesTraPPECO2.indexC] = Electron.UNIT.toSim( 0.70);

        etomica.models.co2.ConformationCO2.Echarge[SpeciesTraPPECO2.indexOleft] = Electron.UNIT.toSim( -0.75);
        etomica.models.co2.ConformationCO2.Echarge[SpeciesTraPPECO2.indexOright] = Electron.UNIT.toSim( -0.75);

    }

    protected final Space space;
    protected static final double bondlength = 1.16;


    protected Vector vector;

    private static final long serialVersionUID = 1L;
}
