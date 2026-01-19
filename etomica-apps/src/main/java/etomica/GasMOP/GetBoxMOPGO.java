package etomica.GasMOP;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.potential.UFF.GeneralGrapheneReader;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;

import java.util.ArrayList;
import java.util.List;

public class GetBoxMOPGO {
    public Box box;
    public SpeciesManager sm;
    public ISpecies speciesMOP;
    public List<String> goFileName;
    public List<Vector> vectorList;
    public ISpecies speciesGO1, speciesGO2, speciesGO3, speciesGO4, speciesGO5, speciesGO6, speciesGO7, speciesGO8, speciesGO9, speciesGO10, speciesGO11, speciesGO12, speciesGO13, speciesGO14, speciesGO15, speciesGO16, speciesGO17, speciesGO18, speciesGO19, speciesGO20;
    public Vector com1, com2, com3, com4, com5, com6, com7, com8, com9, com10, com11, com12, com13, com14, com15, com16, com17, com18, com19, com20;
    public List<ISpecies> speciesGO;
    public List<Vector> speciesGOCOM;
    public SpeciesManager getSpeciesManagerMOPGO(Space space, SpeciesManager sm, GetBoxMOPGO getBoxMOPGO){
        List<String> goFileNames = getBoxMOPGO.getGoFileName();
        List<ISpecies> speciesList = new ArrayList<>();
        List<Vector> comFiles = getBoxMOPGO.getListCOM();
        GeneralGrapheneReader grapheneReader = new GeneralGrapheneReader();
        for (int i = 0; i < goFileNames.size(); i ++){
            Vector vecI = comFiles.get(i);
            String fileName = goFileNames.get(i);
            ISpecies species = grapheneReader.getSpecies(fileName, new Vector3D(0.0, 0.0, 0.0), true);
            speciesList.add(species);
            box.addNewMolecule(species);
        }
        int speciesLargest = 0;
        int moleculeNum  =0;
        for (int i = 0; i < speciesList.size(); i ++){
            int size  = speciesList.get(i).makeMolecule().getChildList().size();
            if (size > speciesLargest){
                speciesLargest = size;
                moleculeNum = i;
            }
        }
        List<Vector> oldPositions = new ArrayList<>();
        IMolecule moleculeLargest = box.getMoleculeList().get(moleculeNum);
        while (oldPositions.size() < moleculeLargest.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        for (int i = 0; i < speciesList.size(); i++){



        }
        return sm;
    }

    public void getListString(List<String> goFileNames){
        this.goFileName = goFileNames;
    }
    public List<String> getGoFileName(){return  goFileName;}
    public void getListCOM(List<Vector> comList){
        this.vectorList = comList;
    }
    public List<Vector> getListCOM(){return vectorList;}
    public void getspeciesMOP(ISpecies speciesMOP){
        this.speciesMOP = speciesMOP;
    }
    public ISpecies getSpeciesMOP(){return speciesMOP;}
}
