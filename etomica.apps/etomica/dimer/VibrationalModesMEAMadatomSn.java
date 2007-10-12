package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class VibrationalModesMEAMadatomSn extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "VibrationalModesMEAMadatomSn";
    public final PotentialMaster potentialMaster;
    public Box box;
    public SpeciesSpheresMono sn, snFix, snAdatom, movable;
    public PotentialMEAM potential;
    CalcGradientDifferentiable calcGradientDifferentiable;
    CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    
    public static void main(String[] args){
        final String APP_NAME = "SimDimerMinMEAMadatomSn";
        final VibrationalModesMEAMadatomSn sim = new VibrationalModesMEAMadatomSn();

        sim.getController().actionPerformed();
        
        
    }
    
    public VibrationalModesMEAMadatomSn() {
        
        super(Space3D.getInstance(), true);
        
        potentialMaster = new PotentialMaster(space);
                
        // Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
        
        snFix = new SpeciesSpheresMono(this, tinFixed);
        sn = new SpeciesSpheresMono(this, tinFixed);
        snAdatom = new SpeciesSpheresMono(this, Tin.INSTANCE);
        movable = new SpeciesSpheresMono(this, Tin.INSTANCE);
        
        getSpeciesManager().addSpecies(snFix);
        getSpeciesManager().addSpecies(sn);
        getSpeciesManager().addSpecies(snAdatom);
        getSpeciesManager().addSpecies(movable);
        
        ((AtomTypeSphere)snFix.getMoleculeType()).setDiameter(3.022); 
        ((AtomTypeSphere)sn.getMoleculeType()).setDiameter(3.022); 
        ((AtomTypeSphere)snAdatom.getMoleculeType()).setDiameter(3.022);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(3.022);
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        // Sn
        box.setNMolecules(snFix, 72);
        box.setNMolecules(sn, 144); 
        box.setNMolecules(snAdatom, 0);
        
        potential = new PotentialMEAM(space);
        
        potential.setParameters(snFix, ParameterSetMEAM.Sn);
        potential.setParameters(sn, ParameterSetMEAM.Sn);
        potential.setParameters(snAdatom, ParameterSetMEAM.Sn);
        potential.setParameters(movable, ParameterSetMEAM.Sn);
        
        this.potentialMaster.addPotential(potential, new Species[]{sn, snFix, snAdatom, movable});

        box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
     
        Configuration config = new ConfigurationLattice(crystal);
        config.initializeCoordinates(box); 
        
        IAtom iAtom = snAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(snAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 10.0);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.1);
        ((IAtomPositioned)iAtom).getPosition().setX(2, -0.1);
        
        String fTail = "-6091374.806734505";
        
        ConfigurationFile configFile = new ConfigurationFile(fTail+"_minimum");
        configFile.initializeCoordinates(box);
        System.out.println("Reading in system coordinates...");
        
        calcGradientDifferentiable = new CalcGradientDifferentiable(box, potentialMaster, 216, 216);
        d = new int[3];
        positions = new double[d.length];
        dForces = new double[positions.length][positions.length];
        
        // setup position array
        for(int i=0; i<positions.length/3; i++){
            for(int j=0; j<3; j++){
                positions[(3*i)+j] = ((IAtomPositioned)box.getLeafList().getAtom(216)).getPosition().x(j);
            }
        }
        
        // fill dForces array
        for(int l=0; l<d.length; l++){
            d[l] = 1;                
            System.arraycopy(calcGradientDifferentiable.df2(d, positions), 0, dForces[l], 0, d.length);
            System.out.println("Calculating force constant row "+l+"...");
            d[l] = 0;
        }
        
        calcVibrationalModes = new CalcVibrationalModes(dForces);
        modeSigns = new int[3];
    
        // calculate vibrational modes and frequencies
        System.out.println("Calculating lambdas...");
        lambdas = calcVibrationalModes.getLambdas();
        System.out.println("Calculating frequencies...");
        frequencies = calcVibrationalModes.getFrequencies();
        modeSigns = calcVibrationalModes.getModeSigns();
        
        System.out.println("Writing data...");
        // output data
        FileWriter writer;
        
        //LAMBDAS
        try { 
            writer = new FileWriter(fTail+"_lambdas");
            for(int i=0; i<lambdas.length; i++){
                writer.write(lambdas[i]+"\n");
            }
            writer.close();
        }catch(IOException e) {
            System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            return;
        }
        
        //FREQUENCIES
        try { 
            writer = new FileWriter(fTail+"_frequencies");
            for(int i=0; i<frequencies.length; i++){
                writer.write(frequencies[i]+"\n");
            }
            writer.close();
        }catch(IOException e) {
            System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            return;
        }
        
        //MODE INFO
        try { 
            writer = new FileWriter(fTail+"_modeSigns");
            writer.write(modeSigns[0]+" positive modes"+"\n");
            writer.write(modeSigns[1]+" negative modes"+"\n");
            writer.write(modeSigns[2]+" total modes"+"\n");
            writer.close();
        }catch(IOException e) {
            System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            return;
        }
        
        System.out.println("Done.");
    }
    
}
