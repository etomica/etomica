package etomica.normalmode;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomPositioned;
import etomica.data.DataSourceScalar;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space1d.Space1D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Energy;

/**
 * Meter that calculates the harmonic energy of a configuration given
 * eigenvectors and omegas corresponding to wave vectors.
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergy extends DataSourceScalar {

    public MeterHarmonicEnergy(CoordinateDefinition coordinateDefinition, NormalModes normalModes) {
        super("Harmonic Energy", Energy.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public double getDataAsScalar() {
        double energySum = 0;
        
        boolean borked = false;
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            int coordinateDim = coordinateDefinition.getCoordinateDim();
            for (int i=0; i<coordinateDim; i++) {
                if (Double.isInfinite(omegaSquared[iVector][i])) {
                    continue;
                }
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenvectors[iVector][j][i] * realT[j];
                    imaginaryCoord += eigenvectors[iVector][j][i] * imaginaryT[j];
                }
                double normalCoord = realCoord*realCoord + imaginaryCoord*imaginaryCoord;
                energySum += waveVectorCoefficients[iVector] * normalCoord * omegaSquared[iVector][i];
                if (!borked && energySum > Math.log(1e200)) {
                    borked = true;
                    System.out.println(iVector+" "+i+" "+realCoord+" "+imaginaryCoord+" "+(waveVectorCoefficients[iVector] * normalCoord * omegaSquared[iVector][i]));
                }
            }
        }
        return energySum;//don't multiply by 1/2 because we're summing over only half of the wave vectors
    }

    public Phase getPhase() {
        return coordinateDefinition.getPhase();
    }

    public void setPhase(Phase newPhase) {

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];

        normalModes.getWaveVectorFactory().makeWaveVectors(newPhase);
        setWaveVectors(normalModes.getWaveVectorFactory().getWaveVectors(),normalModes.getWaveVectorFactory().getCoefficients());
        setEigenvectors(normalModes.getEigenvectors(newPhase));
        setOmegaSquared(normalModes.getOmegaSquared(newPhase));
    }
    
    protected void setWaveVectors(IVector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    protected void setEigenvectors(double[][][] eigenvectors) {
        this.eigenvectors = (double[][][])eigenvectors.clone();
    }
    
    protected void setOmegaSquared(double[][] omega2) {
        omegaSquared = new double[omega2.length][omega2[0].length];
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = omega2[i][j];
            }
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected double[] realT, imaginaryT;
    protected IVector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenvectors;
    protected double[][] omegaSquared;
    protected NormalModes normalModes;
    
    public static void main(String[] args) {
        
        int numAtoms = 8;
        double L = 10;
        
        Simulation sim = new Simulation(Space1D.getInstance(), true, new PotentialMaster(Space1D.getInstance()));

        sim.getDefaults().makeLJDefaults();
        sim.getDefaults().atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);

        Phase phase = new Phase(new BoundaryRectangularPeriodic(sim.getSpace(), sim.getRandom(), L));
        sim.addPhase(phase);
        phase.getAgent(species).setNMolecules(numAtoms);

        AtomArrayList atoms = phase.getSpeciesMaster().getLeafList();
        
        Primitive primitive = new PrimitiveCubic(sim.getSpace());

        CoordinateDefinition coordinateDefinition = new CoordinateDefinitionLeaf(phase, primitive);
        coordinateDefinition.initializeCoordinates(new int[]{numAtoms});

        for(int i=0; i<numAtoms; i++) {
            ((IAtomPositioned)atoms.get(i)).getPosition().E((i+0.5)*L/numAtoms - 0.5*L);
            System.out.println(((IAtomPositioned)atoms.get(i)).getPosition().x(0));
        }
        
        NormalModes normalModes = new NormalModes1DHR();

        MeterHarmonicEnergy meter = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        meter.setPhase(phase);
        ((IAtomPositioned)atoms.get(1)).getPosition().PE(0.5);
        ((IAtomPositioned)atoms.get(6)).getPosition().PE(-1.5);
        System.out.println("Harmonic energy: "+meter.getDataAsScalar());
    }
}
