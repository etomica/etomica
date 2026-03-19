package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.util.ParameterBase;

import java.util.HashMap;
import java.util.Map;

public class GCMCMOPCombines extends Simulation {
    ISpecies speciesOne, speciesTwo, speciesMOP;
    Box box;
    public GCMCMOPCombines(String mopGeom, String confNameOne, String confNameTwo) {
        super(Space3D.getInstance());
        PDBReaderMOPCombine pdbReaderMOPOne = new PDBReaderMOPCombine();
        PDBReaderMOPCombine pdbReaderMOPTwo = new PDBReaderMOPCombine();
        Map<String, AtomType> typeMapNew = new HashMap<>();
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();

     //   speciesMOP = autoMOPCombines.makeSpeciesMOP( box, mopGeom, confNameOne, confNameTwo,pdbReaderMOPOne, pdbReaderMOPTwo, speciesOne, speciesTwo, typeMapNew);
        autoMOPCombines.makeMoleculeSideLength( "trunTetra", confNameOne, confNameTwo, speciesOne, speciesTwo);
        speciesOne = autoMOPCombines.getSpeciesOne();
        speciesTwo = autoMOPCombines.getSpeciesTwo();
        addSpecies(speciesOne);
        addSpecies(speciesTwo);
        box = this.makeBox();
        box.addNewMolecule(speciesOne);
        autoMOPCombines.getSideLength(box, speciesOne, autoMOPCombines);
    }
    public static void main(String[] args) {
        GCMCMOPCombineParams params = new GCMCMOPCombineParams();
        String mopGeom = params.mopGeom;
        String confNameOne = params.confNameOne;
        String confNameTwo = params.confNameTwo;
        GCMCMOPCombines gcmcMOP = new GCMCMOPCombines(mopGeom, confNameOne, confNameTwo);
    }
    
    public static class GCMCMOPCombineParams extends ParameterBase{
        public String mopGeom = "3py42L6";
        public String confNameOne = "D:\\Sem-IX\\pdbs\\3py\\Zr3O";
        public String confNameTwo = "D:\\Sem-IX\\pdbs\\3py\\L1";
    }
}
