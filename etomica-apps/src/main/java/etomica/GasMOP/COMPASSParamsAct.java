package etomica.GasMOP;

import java.util.HashMap;
import java.util.Map;

public class COMPASSParamsAct {
    enum AtomType {
        c4, h1, c3A, c3D, o2e, o2h, h1o, o2;
    }

    private String p2Key(AtomType a, AtomType b) {
        return (a.ordinal() <= b.ordinal()) ? a + "-" + b : b + "-" + a;
    }

    private String p3Key(AtomType a, AtomType b, AtomType c) {
        return (a.ordinal() <= b.ordinal()) ? a + "-" + b + "-" +c : c+"-"+b + "-" + a;
    }

    private String p4Key(AtomType a, AtomType b, AtomType c, AtomType d){
        return (a.ordinal() <= b.ordinal()) ? a + "-" + b + "-" +c +"-" +d : d+"-"+c+"-"+b + "-" + a;
    }

    private static final Map<String, double[]> BOND_PARAMS = new HashMap<>();

    static {
        BOND_PARAMS.put("c4-c4", new double[]{1.5300, 299.6700, -501.7700, 679.8100});
        BOND_PARAMS.put("c4-h1", new double[]{1.1010, 345.0000, -691.8900, 844.6000});
        BOND_PARAMS.put("c3A-c3A", new double[]{1.4040, 469.0000, -627.6000, 860.1000});
    }


    public double[] quadraticBonds(AtomType atomOne, AtomType atomTwo) {
        String key = p2Key(atomOne, atomTwo);
        double[] params = BOND_PARAMS.get(key);

        if (params == null) {
            throw new IllegalArgumentException(
                    "No COMPASS bond parameters for " + atomOne + " - " + atomTwo
            );
        }
        return params;
    }

    public static void main(String[] args) {
        COMPASSParamsAct compassParamsAct = new COMPASSParamsAct();
        double[] bondParams = compassParamsAct.quadraticBonds(AtomType.c3A, AtomType.c3A);
        System.out.println(bondParams);
    }




}
