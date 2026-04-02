package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.chem.elements.*;
import etomica.potential.PotentialMasterBonding;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class COMPASSParams {
    public ISpecies species;
    public ISpecies speciesNoble(Space space, String atomType, boolean isgasMassInfinite){
        Vector3D posM = new Vector3D(new double[]{0, 0, 0});
        AtomType atomType1 = new AtomType(Helium.INSTANCE);
        String symbol = "";
        if(isgasMassInfinite){
            if(atomType.equals("He")){
                atomType1 = new AtomType(new ElementSimple("He", Double.POSITIVE_INFINITY), "He");
                symbol = "He";
            }else if (atomType.equals("Ar")){
                atomType1 = new AtomType(new ElementSimple("Ar", Double.POSITIVE_INFINITY), "Ar");
                symbol = "Ar";
            } else if (atomType.equals("Kr")){
                atomType1 = new AtomType(new ElementSimple("Kr", Double.POSITIVE_INFINITY), "Kr");
                symbol = "Kr";
            }else if (atomType.equals("Ne")){
                atomType1 = new AtomType(new ElementSimple("Ne", Double.POSITIVE_INFINITY), "Ne");
                symbol = "Ne";
            }
        }else{
            if(atomType.equals("He")){
                atomType1 = new AtomType(Helium.INSTANCE);
                symbol = "He";
            }else if (atomType.equals("Ar")){
                atomType1 = new AtomType(Argon.INSTANCE);
                symbol = "Ar";
            } else if (atomType.equals("Kr")){
                atomType1 = new AtomType(Krypton.INSTANCE);
                symbol = "Kr";
            }else if (atomType.equals("Ne")){
                atomType1 = new AtomType(Neon.INSTANCE);
                symbol = "Ne";
            }
        }
        species = new SpeciesBuilder(space).addAtom(atomType1, posM, symbol ).build();
        return species;
    }

    public ISpecies speciesAlkaneAlkene(Space space, String gastype, List<int[]> bonds, List<Double[]> nonBondLJ){
        SetPotential setPotential = new SetPotential();
        if(gastype.equals("methane")){
            Vector posnC =new Vector3D(new double[]{  -3.133 ,  1.888 , -0.720  });
            Vector posnH1 =new Vector3D(new double[]{ -3.840 ,  1.696 ,  0.059   });
            Vector posnH2 =new Vector3D(new double[]{-2.763 ,  2.888  ,-0.627 });
            Vector posnH3 =new Vector3D(new double[]{  -2.319 ,  1.198,  -0.642 });
            Vector posnH4 = new Vector3D(new double[ ]{ -3.609,   1.772 , -1.671 });
            species = new SpeciesBuilder(space)
                    .addAtom(new AtomType(Carbon.INSTANCE), posnC, "C")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH1, "H1")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH2, "H2")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH3, "H3")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH4, "H4").build();

        }else if(gastype.equals("ethane")){
            Vector posnC1 = new Vector3D(new double[]{ -7.839,   0.664,  -0.045});
            Vector posnC2 = new Vector3D(new double[]{ -7.059,   1.851,  -0.583 });
            Vector posnH1 = new Vector3D(new double[]{ -6.704,   2.481,   0.259});
            Vector posnH2 = new Vector3D(new double[]{-7.187,   0.054,   0.615});
            Vector posnH3 = new Vector3D(new double[]{-8.714,   1.020,   0.537});
            Vector posnH4 = new Vector3D(new double[]{ -6.184,   1.495,  -1.166});
            Vector posnH5 = new Vector3D(new double[]{ -7.711,   2.461,  -1.244});
            Vector posnH6 = new Vector3D(new double[]{-8.194,   0.034 , -0.888});
            species = new SpeciesBuilder(space)
                    .addAtom(new AtomType(Carbon.INSTANCE), posnC1, "C1")
                    .addAtom(new AtomType(Carbon.INSTANCE), posnC2, "C2")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH2, "H1")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH3, "H2")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH6, "H3")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH1, "H4")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH4, "H5")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH5, "H6").build();
        }else if (gastype.equals("ethene")){
            Vector posnC1 = new Vector3D(new double[]{-9.343,   0.704,  -0.016  });
            Vector posnC2 = new Vector3D(new double[]{-8.634 ,  1.730 , -0.482  });
            Vector posnH1 = new Vector3D(new double[]{-8.872,  -0.047 ,  0.609 });
            Vector posnH2 = new Vector3D(new double[]{-10.397,   0.613,  -0.259 });
            Vector posnH3 = new Vector3D(new double[]{-7.580,   1.822 , -0.239 });
            Vector posnH4 = new Vector3D(new double[]{-9.105,   2.482,  -1.107 });
            species = new SpeciesBuilder(space)
                    .addAtom(new AtomType(Carbon.INSTANCE), posnC1, "C1")
                    .addAtom(new AtomType(Carbon.INSTANCE), posnC2, "C2")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH1, "H1")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH2, "H2")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH3, "H3")
                    .addAtom(new AtomType(Hydrogen.INSTANCE), posnH4, "H4").build();
        }
    //    setPotential.setBondStretchCOMPASS(species, bonds, pmBonding, gastype);
        getLJ(species, nonBondLJ);
        return species;
    }


    public void getLJ(ISpecies species, List<Double[]> allAtomNonBond){
        List<AtomType> atomTypes = new ArrayList<>();
        atomTypes = species.getUniqueAtomTypes();
        Double[] nonBond = new Double[2];
        Map<String, Double[]> ljBondMap = new HashMap<>();
        ljBondMap.put("C", new Double[]{3.8540,0.0620});
        ljBondMap.put("H", new Double[]{2.8780,0.0230});
        ljBondMap.put("Ne", new Double[]{3.2000,0.0550});
        ljBondMap.put("He", new Double[]{2.9000,0.0050});
        ljBondMap.put("Ar", new Double[]{3.8800,0.2000});
        ljBondMap.put("Kr", new Double[]{4.3000,0.2800});

        for (int i = 0; i < atomTypes.size(); i++){
            String atomName = String.valueOf(atomTypes.get(i));
            int startIndex = atomName.indexOf("[") + 1;
            int endIndex = atomName.indexOf("]");
            String nameNew = atomName.substring(startIndex, endIndex);
            nonBond = ljBondMap.get(nameNew);
            allAtomNonBond.add(nonBond);
        }
    }

    public static void main(String[] args) {
        COMPASSParams params = new COMPASSParams();
        ISpecies species1 = params.speciesNoble(Space3D.getInstance(),"Ar", true );
        System.out.println(species1);
        System.out.println(species1.getUniqueAtomTypes());
        List<int[]> bonds = new ArrayList<>();
        List<Double[]> nonBondLJ = new ArrayList<>();
        species1 = params.speciesAlkaneAlkene(Space3D.getInstance(), "ethene", bonds, nonBondLJ);
        System.out.println(species1.getUniqueAtomTypes());
    }

    //E = K2 * (R - R0)^2  +  K3 * (R - R0)^3  +  K4 * (R - R0)^4
    public double[] quadraticBonds (String atomOne, String atomTwo){
        double[] vals = new double[4];
        if (atomOne.equals("c4") && atomTwo.equals("c4")){
            vals = new double[]{1.5300,    299.6700,   -501.7700,    679.8100};
        } else if (atomOne.equals("c4") && atomTwo.equals("h1")) {
            vals = new double[]{1.1010,    345.0000,   -691.8900,    844.6000};
        } else if (atomOne.equals("c3A") && atomTwo.equals("c3A")) {
            vals = new double[]{1.4170, 470.8361, -627.6179, 1327.6345 };
        }else if (atomOne.equals("c3A") && atomTwo.equals("c4")) {
            vals = new double[]{1.5010, 321.9021, -521.8208, 572.1628};
        }else if (atomOne.equals("c3A") && atomTwo.equals("h1")) {
            vals = new double[]{1.0982, 372.8251, -803.4526, 894.3173};
        }else if (atomOne.equals("c3A") && atomTwo.equals("o2e")) {
            vals = new double[]{1.3768, 428.8798, -738.2351, 1114.9655};
        } else if (atomOne.equals("c4") && atomTwo.equals("o2e")) {
            vals = new double[]{1.4200, 400.3954, -835.1951, 1313.0142};
        } else if (atomOne.equals("c4") && atomTwo.equals("o2h")) {
            vals = new double[]{1.4200, 400.3954, -835.1951, 1313.0142};
        } else if (atomOne.equals("h1") && atomTwo.equals("o2h")) {
            vals = new double[]{0.9494, 540.3633, -1311.8663, 2132.4446};
        } else if (atomOne.equals("c3'") && atomTwo.equals("o2e")) {
            vals = new double[]{1.3750, 368.7309, -832.4784, 1274.0231};
        }else if (atomOne.equals("c3'") && atomTwo.equals("c4")) {
            vals = new double[]{1.5140, 312.3719, -465.8290, 473.8300};
        }else if (atomOne.equals("c3'") && atomTwo.equals("o1=")) {
            vals = new double[]{1.2160, 823.7948, -1878.7939, 2303.5310};
        }else if (atomOne.equals("c3'") && atomTwo.equals("c3a")) {
            vals = new double[]{1.4890, 339.3574, -655.7236, 670.2362};
        }
        return vals;
    }

    //Delta = Theta - Theta0
    // E = K2 * Delta^2  +  K3 * Delta^3  +  K4 * Delta^4
    public double[] quadraticAngle(String atomOne, String atomTwo, String atomThree){
        double[] vals = new double[2];
        if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("c4")){
            vals = new double[]{112.6700,     39.5160,     -7.4430,     -9.5583};
        }else if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){
            vals = new double[]{110.7700,     41.4530,    -10.6040,      5.1290};
        }else if (atomOne.equals("h1") && atomTwo.equals("c4") && atomThree.equals("h1")){
            vals = new double[]{107.6600,     39.6410,    -12.9210,     -2.4318};
        }
        return vals;
    }


    public double bondBondDemo(String atomOne, String atomTwo, String atomThree){
        double val = 0.0;
        if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){

        }else if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){

        }else if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){

        }
        return val;
    }

    //E = K(b,b') * (R - R0) * (R' - R0')
    public double bondBond(String atomOne, String atomTwo, String atomThree){
        double val = 0.0;
        if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){
            val =  3.3872;
        }else if (atomOne.equals("h1") && atomTwo.equals("c4") && atomThree.equals("h1")){
            val = 5.3316;
        }
        return val;
    }

    //Not required
    //E = K(b,b') * (R - R0) * (R' - R0')
    public double bondBondAlt(String atomOne, String atomTwo, String atomThree){
        double val = 0.0;
        return val;
    }

    //E = K * (R - R0) * (Theta - Theta0)
    public double[] bondAngle(String atomOne, String atomTwo, String atomThree){
        double val[] = new double[2];
        if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("c4")){
            val = new double[]{8.0160};
        }else if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("h1")){
            val = new double[]{ 20.7540,     11.4210};
        }else if (atomOne.equals("h1") && atomTwo.equals("c4") && atomThree.equals("h1")){
            val = new double[]{18.1030};
        }
        return val;
    }

   // E = SUM(n=1,3) { V(n) * [ 1 - cos(n*Phi - Phi0(n)) ] }
    public double[] torsion(String atomOne, String atomTwo, String atomThree, String atomFour){
        double val[] =new double[6];
        if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("c4") && atomFour.equals("c4")){
            val = new double[]{0.0000 ,   0.0 ,     0.0514  ,  0.0  ,   -0.1430 ,   0.0};
        }else if (atomOne.equals("c4") && atomTwo.equals("c4") && atomThree.equals("c4")&& atomFour.equals("h1")){
            val = new double[]{ 0.0000 ,   0.0 ,     0.0316 ,   0.0,     -0.1681 ,   0.0};
        }else if (atomOne.equals("h1") && atomTwo.equals("c4") && atomThree.equals("h1") && atomFour.equals("h1")){
            val = new double[]{-0.1432,    0.0 ,     0.0617 ,   0.0,     -0.1530  ,  0.0};
        }
        return val;
    }

    public double[] LJ (String atomOne){
        double val[] = new double[2];
        if (atomOne.equals("c4")){
            val = new double[]{  3.8540,  0.0620};
        }else if (atomOne.equals("h1")){
            val = new double[]{ 2.8780,   0.0230};
        }
        return val;
    }

}
