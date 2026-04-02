package etomica.potential.OPLS_AA;

import etomica.atom.AtomType;
import etomica.chem.elements.*;
import etomica.potential.TraPPE.SpeciesGasTraPPE;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SpeciesGasOPLSAA {
    public static ChemForm ChemForm;
    protected boolean polar;
    protected static Element elementM = new ElementSimple("M", 0.0);
    Map<String, AtomType> typeMapNew = new HashMap<>();
    protected ChemForm chemForm;
    protected Space space;
    public AtomType[] atomTypes;
    public double[] sigma;
    public double[] epsilon;
    public double[] charge;
    public Map<Integer, Double> chargeMap;
    public ISpecies species, species1;
    public Map<Integer, String> atomMap;

    public ISpecies speciesGasOPLSAA(Space space, String confName, ChemForm chemForm, boolean isDynamic){
        this.space = space;
        AtomType typeNew;
        this.chemForm = chemForm;
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        double nMass = Nitrogen.INSTANCE.getMass();
        PDBDataExtracterOPLS pdbDataExtracterOPLS = new PDBDataExtracterOPLS();
        species = pdbDataExtracterOPLS.getSpeciesNew(confName);
        AtomType C3 = new AtomType(Carbon.INSTANCE, "C3");
        AtomType C1 = new AtomType(Carbon.INSTANCE, "C1");
        AtomType O = new AtomType(Oxygen.INSTANCE, "O");
        AtomType N1 = new AtomType(Nitrogen.INSTANCE, "N1");
        AtomType HA = new AtomType(Hydrogen.INSTANCE, "HA");
        AtomType HC = new AtomType(Hydrogen.INSTANCE, "HC");
        Map<Integer, AtomType> atomTypeMap = new HashMap<>();
        if (chemForm == ChemForm.CH4){
            atomTypeMap.put(0, C3);
            atomTypeMap.put(1, HC);
            atomTypeMap.put(2, HC);
            atomTypeMap.put(3, HC);
            atomTypeMap.put(4, HC);
        } else if (chemForm == ChemForm.C2H4) {
            atomTypeMap.put(0, C3);
            atomTypeMap.put(1, C3);
            atomTypeMap.put(2, HC);
            atomTypeMap.put(3, HC);
            atomTypeMap.put(4, HC);
            atomTypeMap.put(5, HC);
            atomTypeMap.put(6, HC);
            atomTypeMap.put(7, HC);
        } else if (chemForm == ChemForm.CO2) {
            atomTypeMap.put(0, C1);
            atomTypeMap.put(1, O);
            atomTypeMap.put(2, O);
        } else if (chemForm == ChemForm.N2) {
            atomTypeMap.put(0, N1);
            atomTypeMap.put(1, N1);
        }else if (chemForm == ChemForm.H2) {
            atomTypeMap.put(0, HA);
            atomTypeMap.put(1, HA);
        }
        Map<Integer, String> atomIdentifierMap = pdbDataExtracterOPLS.getAtomIdentifierMapModified();
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        Map<Integer, Vector> positions = pdbDataExtracterOPLS.getPositions();
        for(int i = 0; i < atomIdentifierMap.size(); i++) {
            AtomType newName = atomTypeMap.get(i);
            Vector position = positions.get(i);
            speciesBuilderNew.addAtom(newName, position,  "");
        }
        species= speciesBuilderNew.build();
        return species;
    }

    public enum ChemForm{
        CH4, C2H6, C3H8, C2H4, C3H6, CO2, N2, O2, NH3, H2
    }


    public double[] atomicPot (String atomType){
        HashMap<String, double[]> atomicConstant = new HashMap<>();



        if (!atomicConstant.containsKey(atomType)) {
            // If atomtype doesn't exist, throw a RuntimeException with a custom message
            throw new RuntimeException("Atom type '" + atomType + "' does not exist in atomic constants.");
        }

        // If atomtype exists, return the corresponding values
        return atomicConstant.get(atomType);
    }
}
