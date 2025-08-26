package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.potential.TraPPE.SpeciesGasTraPPE;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.potential.UFF.PDBReaderReplica;
import etomica.potential.UFF.cifReader;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static etomica.potential.TraPPE.SpeciesGasTraPPE.ChemForm.*;

public class GetSpeciesMOP {
    public Box box;
    public SpeciesManager sm;
    ISpecies speciesLigand;
    public ISpecies speciesGas;
    public void setSpeciesGas (ISpecies speciesGas){
        this.speciesGas = speciesGas;
    }
    public ISpecies getSpeciesGas(){
        return speciesGas;
    }

    public ISpecies getSpeciesMOP (PDBReaderMOP pdbReaderMOP, PDBReaderReplica pdbReaderReplica, SpeciesManager sm, String conf1, String conf2, String struc, String structName, ArrayList<ArrayList<Integer>> connectivityLigand, Map<Integer, String> mapLigand, ArrayList<ArrayList<Integer>> connectivityMOP, Map<Integer, String> mapMOP, boolean ifVirial, boolean setDynamic, boolean ifMOPpresent, boolean ifAutomMOP, boolean ifGasTrappe){
        AutoMOP autoMOP = new AutoMOP();
        GetSpeciesMOP getSpeciesMOP = new GetSpeciesMOP();
       // GetSpeciesMOP getSpeciesMOPNew = new GetSpeciesMOP(Space3D.getInstance());
       /* if (false) {
            if (ifMultipleGraphenePresent) {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                addSpecies(speciesGrapheneOne);
            } else {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                addSpecies(speciesGrapheneOne);
            }
        }*/
        Simulation sim = new Simulation(Space3D.getInstance());
        if (ifMOPpresent) {
            speciesLigand = pdbReaderMOP.getSpeciesMOP(conf1, false, new Vector3D(0.0, 0.0, 0.0), false);
            cifReader cifReader = new cifReader();
            sim.addSpecies(speciesLigand);
        }else if(ifAutomMOP){
            speciesLigand = pdbReaderMOP.getSpeciesMOP(conf1, false, new Vector3D(0.0, 0.0, 0.0), false);
            sim.addSpecies(speciesLigand);
        }


        SpeciesGasTraPPE speciesGasTraPPE = new SpeciesGasTraPPE();
        if (ifGasTrappe) {
            if (conf2.equals("F://Avagadro//molecule//ch4")) {
                SpeciesGasTraPPE.ChemForm = CH4;
            } else if (conf2.equals("F://Avagadro//molecule//ethane")) {
                SpeciesGasTraPPE.ChemForm = C2H6;
            } else if (conf2.equals("F://Avagadro//molecule//propane")) {
                SpeciesGasTraPPE.ChemForm = C3H8;
            } else if (conf2.equals("F://Avagadro//molecule//ethene")) {
                SpeciesGasTraPPE.ChemForm = C2H4;
            } else if (conf2.equals("F://Avagadro//molecule//propene")) {
                SpeciesGasTraPPE.ChemForm = C3H6;
            } else if (conf2.equals("F://Avagadro//molecule//co2")) {
                SpeciesGasTraPPE.ChemForm = CO2;
            } else if (conf2.equals("F://Avagadro//molecule//o2")) {
                SpeciesGasTraPPE.ChemForm = O2;
            } else if (conf2.equals("F://Avagadro//molecule//n2")) {
                SpeciesGasTraPPE.ChemForm = N2;
            }else if (conf2.equals("F://Avagadro//molecule//nh3")) {
                SpeciesGasTraPPE.ChemForm = NH3;
            }
            setSpeciesGas(speciesGas);
            sim.addSpecies(speciesGas);
        }

        //speciesMOP = cifReader.speciesCIF(confNameOne,false);

        //species Gas
       /* if (isGasCOMPASS) {
            //  speciesGas = pdbReaderCOMPASS.getSpecies(confNameGasOne, false, true, new Vector3D(0,0,0));
        } else if (isGasTraPPE) {
            speciesGas = speciesGasTraPPE.speciesGasTraPPE(Space3D.getInstance(), SpeciesGasTraPPE.ChemForm, false);
        } else {
            speciesGas = pdbReaderReplica.getSpecies(confNameGasOne, true, new Vector3D(0, 0, 0), false);
        }*/

        getSpeciesMOP.setSpeciesGas(speciesGas);
        Map<Integer, Integer> oldnewAtomNumMap = new HashMap<>();
        List<Integer[]> metalAtomsNew = new ArrayList<>();
        IMolecule molecule = new Molecule(speciesLigand, speciesLigand.getLeafAtomCount());
        System.out.println(molecule);
        Map<Integer, List<Integer[]>> metalAtomsConnect = new HashMap<>();
        Map<String, AtomType> typeMapNew = new HashMap<>();
        ISpecies species = null;
        connectivityLigand = pdbReaderMOP.connectivityModified;
       // mapLigand = pdbReaderMOP.getAtomMapModified();
        mapLigand  = pdbReaderMOP.getAtomIdentifierMapMod();
      //  sm =  new SpeciesManager.Builder().addSpecies(speciesLigand).build();
        box = new Box(Space3D.getInstance());
        sim.addBox(box);
        box.getBoundary().setBoxSize(new Vector3D(200,200,200));
        mapMOP = new HashMap<>();
        connectivityMOP = new ArrayList<>();
       // box= autoMOP.getAutoMOPBox(space, struc, structName, speciesLigand, box, boxNew, connectivityLigand, mapLigand, ifVirial, setDynamic, 1);
       // System.out.println(box.getMoleculeList());
       // if (ifAutomMOP){
            SpeciesBuilder speciesBuilder =  autoMOP.getAutoMOPBox(Space3D.getInstance(), struc, structName, speciesLigand, box, connectivityLigand, mapLigand, ifVirial, setDynamic);
            ISpecies species1 = speciesBuilder.setDynamic(true).build();
       /* } else {
            return speciesLigand;
        }*/
        return species1;
    }
}
