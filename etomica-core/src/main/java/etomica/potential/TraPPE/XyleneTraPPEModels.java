package etomica.potential.TraPPE;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Degree;

import java.util.Map;

import static etomica.potential.TraPPE.SpeciesGasTraPPE.ChemForm.*;

public class XyleneTraPPEModels {
    protected boolean polar;
    protected static Element elementM = new ElementSimple("M", 0.0);
    protected ChemForm ChemForm;
    protected Space space;
    public enum ChemForm{
        oC8H10, pC8H10, mC8H10, eC8H10;
    }
    public AtomType[] atomTypes;
    public double[] sigma;
    public double[] epsilon;
    public double[] charge;
    public Map<Integer, Double> chargeMap;
    public ISpecies species;
    public Map<Integer, String> atomMap;
    double sigmaC = 3.88, sigmaCH = 3.695, sigmaCH3 = 3.75;
    double epsilonC = 21.0, epsilonCH = 50.5, epsilonCH3 = 98;
    public ISpecies selectSpeciesModule (int moduleNum, String atomName, boolean isDynamic){
        XyleneTraPPEModels xyleneTraPPEModels = new XyleneTraPPEModels();
        if (atomName.equals("F://Avagadro//molecule//oxylene")) {
            xyleneTraPPEModels.ChemForm = xyleneTraPPEModels.ChemForm.oC8H10;
        } else if (atomName.equals("F://Avagadro//molecule//pxylene")) {
            xyleneTraPPEModels.ChemForm = xyleneTraPPEModels.ChemForm.pC8H10;
        } else if (atomName.equals("F://Avagadro//molecule//mxylene")) {
            xyleneTraPPEModels.ChemForm = xyleneTraPPEModels.ChemForm.mC8H10;
        } else if (atomName.equals("F://Avagadro//molecule//ethylbenzene")) {
            xyleneTraPPEModels.ChemForm = xyleneTraPPEModels.ChemForm.eC8H10;
        }else {
            throw new RuntimeException("Error in molecule selected");
        }
        if (moduleNum == 1){
            throw new RuntimeException("To be coded");
        } else if (moduleNum == 2) {
            species = xyleneTraPPEModels.speciesXylene( space, ChemForm, isDynamic);
        }
        return species;
    }

    //pure TraPPE
    public ISpecies speciesXylene (Space space, ChemForm chemForm, boolean isDynamic){
        this.space = space;
        this.ChemForm = chemForm;
        XyleneTraPPEModels xyleneTraPPE = new XyleneTraPPEModels();
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        AtomType typeCH3 = AtomType.simple("ch3", cMass + 3*hMass);
        AtomType typeCH = AtomType.simple("ch", cMass + hMass);
        AtomType typeC = AtomType.simple("c", cMass);
        atomTypes = new AtomType[]{typeCH3, typeCH, typeC};
        double p = 1.54, r = 1.4;
        double sinTheta = Math.sin(Degree.UNIT.toSim(60));
        double cosTheta = Math.cos(Degree.UNIT.toSim(60));
        double sr = sinTheta * r;
        double cr = cosTheta * r;
        double sp = sinTheta * p;
        double cp = cosTheta * p;
        Vector3D vec1 = new Vector3D(new double[]{0,0,0});
        Vector3D vec2 = new Vector3D(new double[]{-sr, -cr, 0});
        Vector3D vec3 = new Vector3D(new double[]{-sr, -cr -r, 0});
        Vector3D vec4 = new Vector3D(new double[]{ 0, -r-(2*cr), 0});
        Vector3D vec5 = new Vector3D(new double[]{ sr, -cr, 0});
        Vector3D vec6 = new Vector3D(new double[]{ sr, -cr-r, 0});
        if (chemForm == ChemForm.oC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ -sr -sp, -cr + cp, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .build();
        }else if (chemForm == ChemForm.mC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{-sr-sp, -cr-r- cp, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .build();
        } else if (chemForm == ChemForm.pC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ 0, -r-2*cr-p, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .build();
        }else if (chemForm == ChemForm.eC8H10) {
            double alpha = Degree.UNIT.toSim(66);
            double sine = Math.sin(alpha);
            double cosine = Math.cos(alpha);
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ p * sine, p + p*cosine, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .build();
        }
        return species;
    }


    public ISpecies speciesXyleneModified(Space space, ChemForm chemForm, boolean isDynamic){
        this.space = space;
        this.ChemForm = chemForm;
        XyleneTraPPEModels xyleneTraPPE = new XyleneTraPPEModels();
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        AtomType typeCH3 = AtomType.simple("ch3", cMass + 3*hMass);
        AtomType typeCH = AtomType.simple("ch", cMass + hMass);
        AtomType typeC = AtomType.simple("c", cMass);
        AtomType typeM = AtomType.simple("m", 0.0);
        atomTypes = new AtomType[]{typeCH3, typeCH, typeC, typeM};
        double p = 1.54, r = 1.4;
        double sinTheta = Math.sin(Degree.UNIT.toSim(60));
        double cosTheta = Math.cos(Degree.UNIT.toSim(60));
        double sr = sinTheta * r;
        double cr = cosTheta * r;
        double sp = sinTheta * p;
        double cp = cosTheta * p;
        Vector3D vec1 = new Vector3D(new double[]{0,0,0});
        Vector3D vec2 = new Vector3D(new double[]{-sr, -cr, 0});
        Vector3D vec3 = new Vector3D(new double[]{-sr, -cr -r, 0});
        Vector3D vec4 = new Vector3D(new double[]{ 0, -r-(2*cr), 0});
        Vector3D vec5 = new Vector3D(new double[]{ sr, -cr, 0});
        Vector3D vec6 = new Vector3D(new double[]{ sr, -cr-r, 0});
        Vector3D vecCentre = new Vector3D();
        vecCentre.PE(vec1);
        vecCentre.PE(vec2);
        vecCentre.PE(vec3);
        vecCentre.PE(vec4);
        vecCentre.PE(vec5);
        vecCentre.PE(vec6);
        vecCentre.TE(1.0/6);
        Vector3D vec9 = new Vector3D();
        Vector3D vec10 = new Vector3D();
        vec9 = new Vector3D(new double[]{vecCentre.getX(0),vecCentre.getX(1), vecCentre.getX(2) + 0.785});
        vec10 = new Vector3D(new double[]{vecCentre.getX(0),vecCentre.getX(1), vecCentre.getX(2) - 0.785});
        if (chemForm == ChemForm.oC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ -sr -sp, -cr + cp, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .addAtom(typeM, vec9, "M")
                    .addAtom(typeM, vec10, "M")
                    .build();
        }else if (chemForm == ChemForm.mC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{-sr-sp, -cr-r- cp, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .addAtom(typeM, vec9, "M")
                    .addAtom(typeM, vec10, "M")
                    .build();
        } else if (chemForm == ChemForm.pC8H10) {
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ 0, -r-2*cr-p, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .addAtom(typeM, vec9, "M")
                    .addAtom(typeM, vec10, "M")
                    .build();
        }else if (chemForm == ChemForm.eC8H10) {
            double alpha = Degree.UNIT.toSim(66);
            double sine = Math.sin(alpha);
            double cosine = Math.cos(alpha);
            Vector3D vec7 = new Vector3D(new double[]{0, p, 0});
            Vector3D vec8 = new Vector3D(new double[]{ p * sine, p + p*cosine, 0});
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, vec1, "C")
                    .addAtom(typeC, vec2, "C")
                    .addAtom(typeCH, vec3, "CH")
                    .addAtom(typeCH, vec4, "CH")
                    .addAtom(typeCH, vec5, "CH")
                    .addAtom(typeCH, vec6, "CH")
                    .addAtom(typeCH3, vec7, "CH3")
                    .addAtom(typeCH3, vec8, "CH3")
                    .addAtom(typeM, vec9, "M")
                    .addAtom(typeM, vec10, "M")
                    .build();
        }
        return species;
    }

    public ISpecies speciesXyleneOPLS(ISpecies speciesXy){


        return speciesXy;
    }


}
