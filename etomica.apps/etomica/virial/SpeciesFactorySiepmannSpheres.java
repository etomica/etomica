package etomica.virial;

import etomica.api.IVector;
import etomica.config.ConformationChainZigZag;
import etomica.api.ISimulation;
import etomica.space.ISpace;
import etomica.api.ISpecies;

/**
 * SpeciesFactory that makes Siepmann's alkane model.
 */
public class SpeciesFactorySiepmannSpheres implements SpeciesFactory, java.io.Serializable {

    public SpeciesFactorySiepmannSpheres(ISpace space, int nA) {
        this(space, nA, nominalBondL, nominalBondTheta);
    }
    
    public SpeciesFactorySiepmannSpheres(ISpace space, int nA, double bondL, double bondTheta) {
        this.nA = nA;
        this.bondL = bondL;
        this.bondTheta = bondTheta;
        this.space = space;
        init();
    }
    
    public void setBondL(double newBondL) {
        bondL = newBondL;
        init();
    }
    
    public double getBondL() {
        return bondL;
    }
    
    public void setBondTheta(double newBondTheta) {
        bondTheta = newBondTheta;
        init();
    }
    
    public double getBondTheta() {
        return bondTheta;
    }
    
    public void init() {
        IVector vector1 = space.makeVector();
        vector1.setX(0, bondL);
        IVector vector2 = space.makeVector();
        vector2.setX(0, bondL*Math.cos(bondTheta));
        vector2.setX(1, bondL*Math.sin(bondTheta));
        conformation = new ConformationChainZigZag(space, vector1, vector2);
    }
    
    public ISpecies makeSpecies(ISimulation sim, ISpace _space) {
        SpeciesAlkane species = new SpeciesAlkane(sim, _space, nA);
        species.setConformation(conformation);
        return species;
    }
    
    private static final long serialVersionUID = 1L;
    protected static final double nominalBondL = 1.54;
    protected static final double nominalBondTheta = Math.PI*114/180;
    protected final ISpace space;
    protected double bondL;
    protected double bondTheta;
    private final int nA;
    private ConformationChainZigZag conformation;
}
