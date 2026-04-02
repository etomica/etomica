package etomica.potential.TraPPE;

public class PDBReaderTraPPE {

    public double setBondLength (String atomType){
        double r0 = 0;
        switch (atomType){
            case "C_3":
                r0 = 1.54;
                break;
            case "C_Ar":
                r0 = 1.40;
                break;
            case "C_2":
                r0 = 1.33;
                break;
        }
        return r0;
    }

    public double[] setBondAngle (String atomOne, String atomTwo, String atomThree, int numHydrogen){
        double[] bondAngle = new double[2];
        if(atomOne.equals("C_3") && atomThree.equals("C_3") && atomTwo.equals("C_3") && numHydrogen ==2){
            bondAngle = new double[]{114, 62500};
        } else if (atomOne.equals("C_3") && atomThree.equals("C_3") && atomTwo.equals("C_3") && numHydrogen ==1) {
            bondAngle = new double[]{112, 62500};
        } else if ((atomOne.equals("C_2") && atomThree.equals("C_3") && atomTwo.equals("C_2") && numHydrogen == 2) ||(atomThree.equals("C_2") && atomOne.equals("C_3") && atomTwo.equals("C_2") && numHydrogen == 2)) {
            bondAngle = new double[]{119.7, 70420};
        } else if ((atomOne.equals("C_2") && atomThree.equals("C_3") && atomTwo.equals("C_2") && numHydrogen == 1) ||(atomThree.equals("C_2") && atomOne.equals("C_3") && atomTwo.equals("C_2") && numHydrogen == 1)) {
            bondAngle = new double[]{119.7, 70420};
        } else if (atomOne.equals("C_3") && atomThree.equals("C_3") && atomTwo.equals("C_2")) {
            bondAngle = new double[]{119.7, 70420};
        }
        return bondAngle;
    }

    public double[] setBondTorsionQuad(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] bondTorsion = new double[4];
        if(atomOne.equals("C_3") && atomThree.equals("C_3") && atomTwo.equals("C_3") && atomFour.equals("C_3")){
            bondTorsion = new double[]{0.0, 355.03, -68.19, 791.32};
        } else if ((atomThree.equals("C_2")&& atomFour.equals("C_2"))||(atomOne.equals("C_2")&&atomTwo.equals("C_2"))) {
            bondTorsion = new double[]{688.5, 86.36, -109.77, -282.24};
        }
        return bondTorsion;
    }

    public double[] setBondTorsionDouble (boolean ifCis){
        double[] bondTorsionDouble = new double[2];
        if(ifCis) {
            bondTorsionDouble = new double[]{12400, Math.PI};
        }else {
            bondTorsionDouble = new double[]{13400, 0.0};
        }
        return bondTorsionDouble;
    }

    public double[] setLJ (String atomOne, int numHydrogen, boolean ifbothAtomsSP2){
        double[] paramsLJ = new double[3];
        if (atomOne.equals("C_4")){
            paramsLJ = new double[]{148, 3.73, 0.0};
        } else if (atomOne.equals("C_3") && numHydrogen == 3) {
            paramsLJ = new double[]{98, 3.75, 0.0};
        } else if (atomOne.equals("C_3") && numHydrogen ==2) {
            paramsLJ = new double[]{46, 3.95, 0.0};
        } else if (atomOne.equals("C_3") && numHydrogen ==0) {
            paramsLJ = new double[]{0.5, 6.4,0.0};
        }else if (atomOne.equals("C_3") && numHydrogen ==1) {
            paramsLJ = new double[]{10, 4.68,0.0};
        }else if (atomOne.equals("C_2") && numHydrogen ==2) {
            paramsLJ = new double[]{85, 3.675,0.0};
        }else if (atomOne.equals("C_2") && numHydrogen ==1 && !ifbothAtomsSP2) {
            paramsLJ = new double[]{47, 3.73,0.0};
        } else if (atomOne.equals("C_2") && numHydrogen ==1 && ifbothAtomsSP2) {
            paramsLJ = new double[]{52, 3.710};
        }

        return paramsLJ;
    }
}

