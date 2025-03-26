package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.IConformation;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.units.Degree;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;


public class DistCalc {
    double dist=0;
    public double phi = (1+ Math.sqrt(5))/2;
    public double phi2 = phi*phi;
    List<Vector> listVecPoly = new ArrayList<>();
    RotationTensor3D rotationTensor;
    RotationTensor3D rotationalTensor1, rotationalTensor2;
    public List<Integer[]> getStrucList(String struc){
        Map<Double, List<Integer[]>> distMap = new HashMap<>();
        List<Integer[]>distEqual = new ArrayList<>();
        DistCalc distCalc = new DistCalc();
        if(struc.equals("tetra")){
            List<Vector> listVectorTetra = new ArrayList<>();
            Vector vec1 = Vector.of(1, 1, 1);
            Vector vec2 = Vector.of(-1, -1, 1);
            Vector vec3 = Vector.of(-1, 1, -1);
            Vector vec4 = Vector.of(1, -1, -1);
            listVectorTetra.add(vec1);
            listVectorTetra.add(vec2);
            listVectorTetra.add(vec3);
            listVectorTetra.add(vec4);
            distEqual = distCalc.edgeBetweenPoints(distMap, listVectorTetra);
            setVecPositionList(listVectorTetra);
        } else if (struc.equals("dodeca")) {
            List<Vector> listVecDodeca = new ArrayList<>();
            Vector vec1 = Vector.of(1, 0, phi2);
            Vector vec2 = Vector.of(1, 0, -phi2);
            Vector vec3 = Vector.of(-1, 0, phi2);
            Vector vec4 = Vector.of(-1, 0, -phi2);
            Vector vec5 = Vector.of(phi, phi, phi);
            Vector vec6 = Vector.of(-phi, -phi, -phi);
            Vector vec7 = Vector.of(phi, -phi, -phi);
            Vector vec8 = Vector.of(phi, phi, -phi);
            Vector vec9 = Vector.of(-phi, phi, -phi);
            Vector vec10 = Vector.of(-phi, phi, phi);
            Vector vec11 = Vector.of(-phi, -phi, phi);
            Vector vec12 = Vector.of(phi, -phi, phi);
            Vector vec13 = Vector.of(0, phi2, 1);
            Vector vec14 = Vector.of(0, phi2, -1);
            Vector vec15 = Vector.of(0, -phi2, 1);
            Vector vec16 = Vector.of(0, -phi2, -1);
            Vector vec17 = Vector.of(phi2, 1, 0);
            Vector vec18 = Vector.of(phi2, -1, 0);
            Vector vec19 = Vector.of(-phi2, 1, 0);
            Vector vec20 = Vector.of(-phi2, -1, 0);
            listVecDodeca.add(vec1);
            listVecDodeca.add(vec2);
            listVecDodeca.add(vec3);
            listVecDodeca.add(vec4);
            listVecDodeca.add(vec5);
            listVecDodeca.add(vec6);
            listVecDodeca.add(vec7);
            listVecDodeca.add(vec8);
            listVecDodeca.add(vec9);
            listVecDodeca.add(vec10);
            listVecDodeca.add(vec11);
            listVecDodeca.add(vec12);
            listVecDodeca.add(vec13);
            listVecDodeca.add(vec14);
            listVecDodeca.add(vec15);
            listVecDodeca.add(vec16);
            listVecDodeca.add(vec17);
            listVecDodeca.add(vec18);
            listVecDodeca.add(vec19);
            listVecDodeca.add(vec20);
            distEqual = distCalc.edgeBetweenPoints(distMap, listVecDodeca);
            setVecPositionList(listVecDodeca);
        } else if (struc.equals("cube")) {
            List<Vector> listVecCubic = new ArrayList<>();
            Vector vec1 = Vector.of(1, 1, 1);
            Vector vec2 = Vector.of(-1, -1, -1);
            Vector vec3 = Vector.of(-1, 1, 1);
            Vector vec4 = Vector.of(1, -1, 1);
            Vector vec5 = Vector.of(1, 1, -1);
            Vector vec6 = Vector.of(-1, -1, 1);
            Vector vec7 = Vector.of(1, -1, -1);
            Vector vec8 = Vector.of(-1, 1, -1);
            listVecCubic.add(vec1);
            listVecCubic.add(vec2);
            listVecCubic.add(vec3);
            listVecCubic.add(vec4);
            listVecCubic.add(vec5);
            listVecCubic.add(vec6);
            listVecCubic.add(vec7);
            listVecCubic.add(vec8);
            distEqual = distCalc.edgeBetweenPoints(distMap, listVecCubic);
            setVecPositionList(listVecCubic);
        } else if (struc.equals("icosa")) {
            List<Vector> listVecIcosa = new ArrayList<>();
            Vector vec1 = Vector.of(0, 1, phi);
            Vector vec2 = Vector.of(0, 1, -phi);
            Vector vec3 = Vector.of(0, -1, phi);
            Vector vec4 = Vector.of(0, -1, -phi);
            Vector vec5 = Vector.of(1, phi, 0);
            Vector vec6 = Vector.of(1, -phi, 0);
            Vector vec7 = Vector.of(-1, phi, 0);
            Vector vec8 = Vector.of(-1, -phi, 0);
            Vector vec9 = Vector.of(phi, 0, 1);
            Vector vec10 = Vector.of(-phi, 0, 1);
            Vector vec11 = Vector.of(phi, 0, -1);
            Vector vec12 = Vector.of(-phi, 0, -1);
            listVecIcosa.add(vec1);
            listVecIcosa.add(vec2);
            listVecIcosa.add(vec3);
            listVecIcosa.add(vec4);
            listVecIcosa.add(vec5);
            listVecIcosa.add(vec6);
            listVecIcosa.add(vec7);
            listVecIcosa.add(vec8);
            listVecIcosa.add(vec9);
            listVecIcosa.add(vec10);
            listVecIcosa.add(vec11);
            listVecIcosa.add(vec12);
            distEqual = distCalc.edgeBetweenPoints(distMap, listVecIcosa);
            setVecPositionList(listVecIcosa);
        } else if (struc.equals("octa")) {
            List<Vector> listVecRhombic = new ArrayList<>();
            Vector vec1 = Vector.of(1, 0, 0);
            Vector vec2 = Vector.of(-1, 0, 0);
            Vector vec3 = Vector.of(0, 1, 0);
            Vector vec4 = Vector.of(0, -1, 0);
            Vector vec5 = Vector.of(0, 0, 1);
            Vector vec6 = Vector.of(0, 0, -1);
            listVecRhombic.add(vec1);
            listVecRhombic.add(vec2);
            listVecRhombic.add(vec3);
            listVecRhombic.add(vec4);
            listVecRhombic.add(vec5);
            listVecRhombic.add(vec6);
            distEqual = distCalc.edgeBetweenPoints(distMap, listVecRhombic);
            setVecPositionList(listVecRhombic);
        }
        return distEqual;
    }

    public void reOrientMolecule(Space space, Box box,ISpecies species, ArrayList<Integer> arrayListLinkerCarb, ArrayList<Integer> arrayListLinkerOxy, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map ){
        if(arrayListLinkerCarb.size() == 3){
            int carbA = arrayListLinkerCarb.get(0);
            int carbB = arrayListLinkerCarb.get(1);
            int carbC = arrayListLinkerCarb.get(2);
            Vector carbAPosn = species.makeMolecule().getChildList().get(carbA).getPosition();
            Vector carbBPosn = species.makeMolecule().getChildList().get(carbB).getPosition();
            Vector carbCPosn = species.makeMolecule().getChildList().get(carbC).getPosition();
            IMolecule molecule = species.makeMolecule();
            //rotateAB into same plane
             Vector vA = new Vector3D(), vB = new Vector3D(), vC = new Vector3D(), carbAB = new Vector3D(), vAB = new Vector3D(), vABCopy = new Vector3D(), carbABCopy = new Vector3D(), finalvecAB = new Vector3D(), finalcarbAB = new Vector3D(), vAC = new Vector3D();
             vA.E(carbAPosn);
             //B if in same plane as A
             vB = Vector.of(carbBPosn.getX(0), carbBPosn.getX(1), carbAPosn.getX(2) );
             vAB.E(vA);
             vAB.ME(vB);
             finalvecAB.E(vA);
             vABCopy.E(vAB);
             vAB.normalize();
             carbAB.E(carbA);
             carbAB.ME(carbBPosn);
             finalcarbAB.E(carbAB);
             carbABCopy.E(carbAB);
             carbAB.normalize();

             double prodDotABRotate = vAB.dot(carbAB);
             double thetaRotate = Math.acos(prodDotABRotate);

             vABCopy.XE(carbABCopy);
             vABCopy.normalize();

            rotationTensor.setRotationAxis(vABCopy, thetaRotate);
            molecule.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListLinkerCarb.get(0)) {
                    return;
                }

                atom.getPosition().ME(carbAPosn);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(carbAPosn);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });


             //rotate AC into same plane as AB
            Vector vACcopy = new Vector3D(), carbACCopy = new Vector3D(), finalvecAC = new Vector3D(), finalcarbAC = new Vector3D(), carbAC = new Vector3D();
            vC = Vector.of(carbCPosn.getX(0), carbCPosn.getX(1), carbAPosn.getX(2) );
            vAC.E(carbA);
            vAC.ME(vC);
            vACcopy.E(vAC);
            vAC.normalize();
            finalvecAC.E(vAC);
            carbAC.E(carbAPosn);
            carbAC.ME(carbCPosn);
            carbACCopy.E(carbAC);
            finalcarbAC.E(carbAC);
            carbAC.normalize();

            double prodDotACRotate = vAC.dot(carbAC);
            double deltaRotate = Math.acos(prodDotACRotate);

            vACcopy.XE(carbACCopy);
            vACcopy.normalize();
            rotationTensor.setRotationAxis(vACcopy, deltaRotate);
            molecule.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListLinkerCarb.get(0)) {
                    return;
                }

                atom.getPosition().ME(carbAPosn);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(carbAPosn);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });





            //to be changed
            Vector vecCact = new Vector3D();
            vecCact = idealStruc(molecule, arrayListLinkerCarb);
            Vector vecAB = new Vector3D();
            vecAB.equals(carbAPosn);
            vecAB.ME(carbBPosn);
            double distAB = Math.sqrt(vecAB.squared());
            Vector unitAB = new Vector3D();
            unitAB.E(vecAB);
            unitAB.normalize();




            Vector vecC = new Vector3D();
            Vector vecCCact = new Vector3D();
            vecCCact.E(vecC);
            vecCCact.ME(vecCact);
            Vector unitCCact = new Vector3D();
             unitCCact.E(vecCCact);
            unitCCact.normalize();

            double distCCact = Math.sqrt(vecCCact.squared());

            IMolecule moleculeMOPOne = species.makeMolecule();

            List<Vector> oldPosition = new ArrayList<>();
            IMolecule moleculeMOP = box.getMoleculeList().get(0);
            while (oldPosition.size() < moleculeMOP.getChildList().size()) {
                oldPosition.add(space.makeVector());
            }

            Vector vecAC = new Vector3D();
            vecAC.equals(carbAPosn);
            Vector vecABCCross = new Vector3D();
            vecAC.ME(carbCPosn);
            vecABCCross.E(vecAC);
            vecABCCross.XE(vecAB);
            double distABC = Math.sqrt(vecABCCross.squared()) / Math.sqrt(vecAB.squared());
            double denominator = Math.sqrt(vecAB.squared());

            moleculeMOPOne.getChildList().forEach(atom -> {
                Vector vecAM = new Vector3D();
                Vector vecABMCross = new Vector3D();
                oldPosition.get(atom.getIndex()).E(atom.getPosition());
                vecAM.equals(carbAPosn);
                vecAM.ME(atom.getPosition());
                vecABMCross.E(vecAM);
                vecABMCross.XE(vecAB);
                double distABM = Math.sqrt(vecABMCross.squared()) / denominator;

                Vector vecCCMod = new Vector3D();
                vecCCMod.E(vecCCact);
                vecCCMod.TE(distABM/distABC);
                atom.getPosition().PE(vecCCMod);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }
    }


    public void reOrientMoleculeMethodTwo(Space space, Box box,ISpecies species, ArrayList<Integer> arrayListLinkerCarb, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Vector> moleculeProjection){
        if(arrayListLinkerCarb.size() == 3){
            //allign along Z axis
            int carbA = arrayListLinkerCarb.get(0);
            int carbB = arrayListLinkerCarb.get(1);
            int carbC = arrayListLinkerCarb.get(2);
            Vector carbAPosn = species.makeMolecule().getChildList().get(carbA).getPosition();
            Vector carbBPosn = species.makeMolecule().getChildList().get(carbB).getPosition();
            Vector carbCPosn = species.makeMolecule().getChildList().get(carbC).getPosition();
            //System.out.println("\n "+carbAPosn + " "+ carbBPosn);
            //System.out.println(carbCPosn +"\n");
            IMolecule molecule = species.makeMolecule();

            //Normal to the plane
            Vector vecAB = new Vector3D(), vecBC = new Vector3D(), vecCross = new Vector3D(), vecNormal = new Vector3D(), vecNormalCopy = new Vector3D(), vecCrossPrd = new Vector3D();
            vecAB.E(carbAPosn);
            vecAB.ME(carbBPosn);
            vecBC.E(carbCPosn);
            vecBC.ME(carbBPosn);
            vecCross.E(vecAB);
            vecCross.XE(vecBC);
            vecNormal.E(vecCross);
            vecNormal.normalize();
            Vector vecZNormal = Vector.of(0,0,1);
            vecNormalCopy.E(vecNormal);
            double productDot = vecNormal.dot(vecZNormal);
            double theta = Math.acos(productDot);
            vecCrossPrd.E(vecNormal);
            vecCrossPrd.XE(vecZNormal);
            vecCrossPrd.normalize();
            rotationTensor.setRotationAxis(vecCrossPrd, theta);
            molecule.getChildList().forEach(atom -> {
                rotationTensor.transform(atom.getPosition());
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            // rearrange atoms in the molecule so as to align iun form of equilateral

            Vector vecCact = new Vector3D();
            carbAPosn = molecule.getChildList().get(carbA).getPosition();
            carbBPosn = molecule.getChildList().get(carbB).getPosition();
            carbCPosn = molecule.getChildList().get(carbC).getPosition();
            //System.out.println("\n "+carbAPosn + " "+ carbBPosn);
            //System.out.println(carbCPosn +"\n");

            Vector carbAPosnCopy = new Vector3D(), carbBPosnCopy = new Vector3D(), carbCPosnCopy = new Vector3D(), vecAC = new Vector3D(), vecACNorm = new Vector3D(), vecABNorm = new Vector3D(), vecCrossRotate = new Vector3D(), vecCC = new Vector3D(), vecCCNorm = new Vector3D();
            carbAPosnCopy.E(carbAPosn);
            carbBPosnCopy.E(carbBPosn);
            carbCPosnCopy.E(carbCPosn);
            vecAB.E(carbAPosnCopy);
            vecAB.ME(carbBPosnCopy);
            vecABNorm.E(vecAB);
            vecABNorm.normalize();
            vecAC.E(carbAPosnCopy);
            vecAC.ME(carbCPosnCopy);
            vecACNorm.E(vecAC);
            vecACNorm.normalize();

            //movement of every atom
            vecCrossRotate.E(vecAC);
            vecCrossRotate.XE(vecABNorm);
            double distABC = Math.sqrt(vecCrossRotate.squared()) / Math.sqrt(vecACNorm.squared());

            vecCact = idealStruc(molecule, arrayListLinkerCarb);
            vecCC.E(carbCPosn);
            vecCC.ME(vecCact);
            vecCCNorm.E(vecCC);
            vecCCNorm.normalize();
            double distCCact = Math.sqrt(vecCC.squared());

            Vector vecAM = new Vector3D();
            Vector vecABMCross = new Vector3D();
            Vector vecABMCrossNorm = new Vector3D();
            Vector vecCCMove = new Vector3D();
            double[] distABM = new double[1];
            double[] multiplier = new double[1];
            Vector finalCarbAPosn = carbAPosn;
            AtomicInteger i = new AtomicInteger();
            molecule.getChildList().forEach(atom -> {
                Vector vecProj = moleculeProjection.get(i.get());
                vecAM.E(finalCarbAPosn);
                vecAM.ME(vecProj);
                vecABMCross.E(vecAM);
                vecABMCross.XE(vecAB);
                vecABMCrossNorm.E(vecABMCross);
                vecABMCrossNorm.normalize();
                distABM[0] = Math.sqrt(vecABMCrossNorm.squared())/Math.sqrt(vecABNorm.squared());
                multiplier[0] = distABM[0]/distABC;
                vecCCMove.E(vecCC);
                vecCCMove.TE(multiplier[0]);
                i.getAndIncrement();
            });
        }
    }

    public Vector idealStruc (IMolecule molecule, ArrayList<Integer> listExtremeties ){
        Vector vecReturn = new Vector3D();
        if(listExtremeties.size() ==3) {
            Vector v1 = molecule.getChildList().get(listExtremeties.get(0)).getPosition();
            Vector v2 = molecule.getChildList().get(listExtremeties.get(1)).getPosition();
            Vector v3 = molecule.getChildList().get(listExtremeties.get(2)).getPosition();
            double angleOne = 60;
            double angleTwo = -60;

            double stheta = Math.sin(Degree.UNIT.toSim(angleOne));
            double ctheta = Math.cos(Degree.UNIT.toSim(angleTwo));

            double cxOne =  ctheta * (v1.getX(0) - v2.getX(0)) - stheta * (v1.getX(1) - v2.getX(1));
            double cyOne =  stheta * (v1.getX(0) - v2.getX(0)) + ctheta * (v1.getX(1) - v2.getX(1));
            double czOne = (v1.getX(3) +  v2.getX(3)) / 2;

            double cxTwo =  ctheta * (v1.getX(0) - v2.getX(0)) + stheta * (v1.getX(1) - v2.getX(1));
            double cyTw0 = -stheta * (v1.getX(0) - v2.getX(0)) + ctheta * (v1.getX(1) - v2.getX(1));
            double czTwo = (v1.getX(3) +  v2.getX(3)) / 2;

            Vector vactOne = Vector.of(cxOne, cyOne, czOne);
            Vector vactTwo = Vector.of(cxTwo, cyTw0, czTwo);


            double distOne = Math.sqrt(vactOne.Mv1Squared(v3));
            double distTwo = Math.sqrt(vactTwo.Mv1Squared(v3));

            if ( distOne < distTwo){
                vecReturn.E(vactOne);
            }else {
                vecReturn.E(vactTwo);
            }
        } else if (listExtremeties.size() == 4) {

        } else {
            throw new RuntimeException("Incorrect struct");
        }
        return vecReturn;
    }
    public void getSpecies(ISpecies speciesMOP, Space space,  String struc, String structName, ISpecies speciesLigand, Box box, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map){
        Box box1 = new Box(space);
        box1 = getAutoMOPBox(space, struc, structName, speciesLigand, box, connectivity, map);
        IMoleculeList molecules = new MoleculeArrayList();
        molecules = box.getMoleculeList(speciesLigand);
        for (int i =0; i< molecules.size(); i++){
            System.out.println(molecules.get(i));
        }
        System.exit(1);

    }
    public void makespeciesMOP(Space space, String struc, String structName, ISpecies speciesLigand, ISpecies speciesMOP, Box box, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map, ArrayList<ArrayList<Integer>> connectivityMOP, Map<Integer, String> mapMOP){
        DistCalc distCalc = new DistCalc();
        box= distCalc.getAutoMOPBox(space, struc, structName, speciesLigand, box, connectivity, map);
        System.out.println(box.getMoleculeList());
        System.exit(1);

    }

    public Box getAutoMOPBox(Space space, String struc, String structName, ISpecies species, Box box, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map){
       // System.out.println("inside automop");
        rotationTensor = new RotationTensor3D();
        List<Integer[]> linkedVertices = getStrucList(structName);
        ArrayList<Integer> arrayListMetals = new ArrayList<>();
        ArrayList<Integer> arrayListLinkerOxy = new ArrayList<>();
        ArrayList<Integer> arrayListLinkerCarb = new ArrayList<>();
        calcExtremities(connectivity, map, arrayListMetals);
        locateLinkerOxy(species, connectivity, map, arrayListLinkerOxy, arrayListLinkerCarb);
        ArrayList<Vector> moleculeProjection = new ArrayList<>();
        if(struc.equals("face") || struc.equals("bent")){
            moleculeProjection = getMoleculeProjection(space, box, species, arrayListLinkerCarb, arrayListLinkerOxy, connectivity, map );
            reOrientMoleculeMethodTwo(space, box, species, arrayListLinkerCarb, arrayListLinkerOxy,moleculeProjection);
        }
        if(struc.equals("face")) {
            dist = calcDistFace(species, arrayListLinkerOxy);
        }
        int numEdges = 0;
        double coordinateMultiplier = 0; // converts edge length to circumsphere radius for modifying vertices
        double sqrt2 = Math.sqrt(2);
        if(structName.equals("tetra")){
            numEdges = 6;
            coordinateMultiplier = dist/(2*sqrt2);
        }  else if (structName.equals("icosa")) {
            numEdges = 30;
            // coordinateMultiplier = phi2 * dist / (2* sqrt3);
            coordinateMultiplier = 0.5 * dist;
        } else if (structName.equals("dodeca")) {
            numEdges = 30;
            //  coordinateMultiplier = phi* sqrt3/2 * dist;
            coordinateMultiplier = 0.5 * dist;
        } else if (structName.equals("octa")) {
            numEdges = 12;
            coordinateMultiplier =  dist/sqrt2;
        } else if (structName.equals("cube")) {
            numEdges = 12;
              //coordinateMultiplier = Math.sqrt(3)/2 * dist;
            coordinateMultiplier = 0.5 * dist;
        }
        coordinateMultiplier = 0.0;
        int atomOne = arrayListMetals.get(0);
        int atomTwo = arrayListMetals.get(1);
        IMolecule molecule = species.makeMolecule();
        Vector posnOne = molecule.getChildList().get(atomOne).getPosition();
        Vector posntTwo = molecule.getChildList().get(atomTwo).getPosition();
        coordinateMultiplier= Math.sqrt(posntTwo.Mv1Squared(posnOne));
       /// System.out.println(structName);
        if(struc.equals("edge")){
            getMOPEdge(box, structName,  species,linkedVertices, arrayListMetals, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, coordinateMultiplier);
        } else if (struc.equals("face") || struc.equals("bent")) {
            getMOPFaceMod(box,struc, structName, species, connectivity, map,linkedVertices, arrayListMetals, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, dist, moleculeProjection );
        }/*else if (struc.equals("bent")) {
            getMOPBent(box, structName, species, connectivity, map, linkedVertices, arrayListMetals, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, dist, moleculeProjection);
        }*/else {
            box.setNMolecules(species, 1);
        }

       // box.setNMolecules(species, 1);
        //System.exit(1);
        return box;
    }


    public Box getAutoMOPComplex(Space space, String struc1, String struct2, String StrucGeom1, ISpecies species1, ISpecies species2, Box box, ArrayList<ArrayList<Integer>> connectivity1, ArrayList<ArrayList<Integer>> connectivity2, Map<Integer, String> map1, Map<Integer, String> map2){
        rotationalTensor1 = new RotationTensor3D();
        rotationalTensor2 = new RotationTensor3D();
        ArrayList<Integer> arrayListLinkerOxy = new ArrayList<>();
        ArrayList<Integer> arrayListLinkerCarb = new ArrayList<>();
        List<Integer[]> linkedVertices = getStrucList(StrucGeom1);
        ArrayList<Integer> metalList1 = new ArrayList<>(), oxyList1 = new ArrayList<>(), metalList2 = new ArrayList<>(), oxyList2 = new ArrayList<>();
        Map<Integer, Integer[]> oxyConnectedCommonCarbon1 = new HashMap<>(), oxyConnectedCommonCarbon2 = new HashMap<>();
        getOxyList(species1.makeMolecule(), connectivity1, map1, metalList1, oxyList1, oxyConnectedCommonCarbon1);
        getOxyList(species2.makeMolecule(), connectivity2, map2, metalList2, oxyList2, oxyConnectedCommonCarbon1);
        if(metalList1.size() != 0 || oxyConnectedCommonCarbon1.size() != 0){
            parallelizeOxygenOrganic( box, species1.makeMolecule(), oxyConnectedCommonCarbon1);
        }
        locateLinkerOxy(species1, connectivity1, map1, arrayListLinkerOxy, arrayListLinkerCarb);
        parallelizeOxygenOrganic( box, species1.makeMolecule(), oxyConnectedCommonCarbon1);
        if(metalList2.size() != 0 || oxyConnectedCommonCarbon2.size() != 0){
            parallelizeOxygenOrganic( box, species2.makeMolecule(), oxyConnectedCommonCarbon2);
        }
        ArrayList<Vector> moleculeProjection = new ArrayList<>();
        if(struc1.equals("face") ){
            moleculeProjection = getMoleculeProjection(space, box, species1, arrayListLinkerCarb, arrayListLinkerOxy, connectivity1, map1 );
            reOrientMoleculeMethodTwo(space, box, species1, arrayListLinkerCarb, arrayListLinkerOxy,moleculeProjection);
        }
        Vector centre1 = new Vector3D();
        if(StrucGeom1.equals("edge")){
           // getMOPEdge(box, s,linkedVertices, arrayListMetals, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, coordinateMultiplier);
        } else if (StrucGeom1.equals("face")) {
            getFaceMOPComplex ( box, StrucGeom1,species1, species2, linkedVertices, map1, map2, metalList1, centre1,  metalList2,  oxyList1,  oxyList2,  moleculeProjection);
        }else if (StrucGeom1.equals("bent")) {
            //getMOPBent(box, structName, species, connectivity, map, linkedVertices, arrayListMetals, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, dist, moleculeProjection);
        }else {
         box.setNMolecules(species1, 1);
         box.setNMolecules(species2, 1);
        }
        return box;
    }

    public void getFaceMOPComplex (Box box, String strucName, ISpecies species1, ISpecies species2, List<Integer[]> linkedVertices, Map<Integer, String> map1, Map<Integer, String> map2, ArrayList<Integer> arrayListMetal1, Vector centre1, ArrayList<Integer> arrayListMetal2, ArrayList<Integer> arrayListLinkerOxy1, ArrayList<Integer> arrayListLinkerOxy2, ArrayList<Vector> moleculeProjection){
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        int numFace, numVertices;
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();

        //gets coordinates depending on struct
        verticesComboFaces(strucName, vertices, faces, faceEdges );

        //get dist bw two extremities
        Vector vecVert = vertices.get(0);
        Vector vecVertTwo = vertices.get(1);
        double faceLen = Math.sqrt(vecVert.Mv1Squared(vecVertTwo));

        Vector vecOne = species1.makeMolecule().getChildList().get(arrayListLinkerOxy1.get(0)).getPosition();
        Vector vecTwo = species1.makeMolecule().getChildList().get(arrayListLinkerOxy1.get(1)).getPosition();
        double lenLig1 = Math.sqrt(vecOne.Mv1Squared(vecTwo));

        Vector vecThree = species2.makeMolecule().getChildList().get(arrayListMetal1.get(0)).getPosition();
        Vector vecFour = species2.makeMolecule().getChildList().get(arrayListMetal1.get(1)).getPosition();
        double lenLig2 = Math.sqrt(vecOne.Mv1Squared(vecTwo));

        double atomOneLJ = pdbReaderMOP.atomicPot(map1.get(arrayListLinkerOxy1.get(0)))[0];
        double atomTwoLJ = pdbReaderMOP.atomicPot(map2.get(arrayListMetal2.get(1)))[0];
        double coordinateMultiplier = (lenLig1 +  2 * lenLig2 + atomTwoLJ+ atomOneLJ) / faceLen;

        if(strucName.equals("tetra")){
            numFace = 4;
            numVertices = 4;
        }else if (strucName.equals("cube")) {
            numFace = 6;
            numVertices = 8;
        } else if (strucName.equals("octa")) {
            numFace = 8;
            numVertices = 6;
        } else if (strucName.equals("icosa")) {
            numFace = 20;
            numVertices = 12;
        } else if (strucName.equals("dodeca")) {
            numFace = 12;
            numVertices = 20;
        }else {
            numFace = 1;
            numVertices = 1;
            throw  new RuntimeException("Shape not defined in face MOP");
        }

        //Move all ligands to the zeroth coordinate
        Space space = Space3D.getInstance();
        box.setNMolecules(species1, numFace);
        box.setNMolecules(species2, numVertices);
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }

        List<Vector> verticesModified = new ArrayList<>();
        //Multiply unit coordinates to coordinates of Polyhedra acc to Ligand size

        for(int i = 0; i< vertices.size(); i++){
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vertices.get(i);
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
            verticesModified.add(i, copyOriginOne);
        }
        Map<Integer, Vector> centreBentLigand = new HashMap<>();

        // For every face  of MOP
        Map<Integer, List<List<Vector>>> vectorListMap = new HashMap<>();  // vector positions
        Map<Integer, List<List<Integer>>> atomListMap = new HashMap<>();   // atom numbers at the positions
        Map<Integer, Vector> centreMap = new HashMap<>();
        Vector centre = new Vector3D();
        for (int i =0; i< numVertices; i++){
            centreMap.put(i, centre);
        }
        //Move every ligand to first vertice of the face
        for (int i = 0; i < numVertices; i++){
            Integer[] vecArray = faces.get(i);
            moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne =verticesModified.get(vecArray[0]);
            copyOriginOne.E(originOne);
            Vector moleculeStart = moleculeMOPOne.getChildList().get(arrayListMetal2.get(0)).getPosition();
            copyOriginOne.ME(moleculeStart);
            Vector centreMod = centreMap.get(i);
            moleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().ME(copyOriginOne);
                atom.getPosition().PE(shift);
            });
            centreMod.ME(copyOriginOne);
            Vector shift = box.getBoundary().centralImage(centreMod);
            centreMod.PE(copyOriginOne);
            centreMod.PE(shift);
            //rotate vertice ligand
            // translate or rotate it as is so as to make sense
        }
        //add face ligand



    }


    public List<Integer> getOxyList(IMolecule molecule, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map, ArrayList<Integer> metalList, ArrayList<Integer> oxyList, Map<Integer, Integer[]> oxyConnectedCommonCarbon){
        //check if metal is present
        for(int i =0; i<connectivity.size(); i++){
            int numZero = connectivity.get(i).get(0);
            String atomZero = map.get(numZero);
            if(atomZero.equals("FE+2")|| atomZero.equals("FE+3")||atomZero.equals("CR")||atomZero.equals("MO")||atomZero.equals("AL")||atomZero.equals("PD")||atomZero.equals("PT")||atomZero.equals("CU")||atomZero.equals("ZN")||atomZero.equals("ZR")||atomZero.equals("RH")||atomZero.equals("RU")){
                metalList.add(numZero);
            } else if (atomZero.equals("NI") || atomZero.equals("CO") ||atomZero.equals("PO") ||atomZero.equals("MO") ||atomZero.equals("MG") ||atomZero.equals("V") ||atomZero.equals("W") ||atomZero.equals("IR") ||atomZero.equals("TI") ||atomZero.equals("MN")) {
                metalList.add(numZero);
            }
        }

        //check if metal atom is connected to two oxygen atoms
        Integer[] oxyConnected = new Integer[2];
        for(int i =0; i< metalList.size(); i++){
            ArrayList<Integer> connectivityMetal = connectivity.get(metalList.get(i));
            if(connectivityMetal.size() == 3){
                for (int j =0; j<2; j++){
                    oxyConnected = new Integer[2];
                    int numOne = connectivityMetal.get(j);
                    if(map.get(numOne).equals("O")){
                        if(connectivity.get(numOne).size() == 2){
                            oxyList.add(numOne);
                            oxyConnected[j] = numOne;
                        }
                    }
                }
            }
            oxyConnectedCommonCarbon.put(i, oxyConnected);
        }
        return oxyList;
    }

    public void parallelizeOxygenOrganic (Box box, IMolecule molecule, Map<Integer, Integer[]> oxyConnectedCommonCarbon){
        int oxyOne = oxyConnectedCommonCarbon.get(0)[0];
        int oxyTwo = oxyConnectedCommonCarbon.get(0)[1];
        int oxyThree = oxyConnectedCommonCarbon.get(1)[0];
        int oxyFour = oxyConnectedCommonCarbon.get(1)[1];
        Vector oxyOnePosn = molecule.getChildList().get(oxyOne).getPosition();
        Vector oxyTwoPosn = molecule.getChildList().get(oxyTwo).getPosition();
        Vector oxyThreePosn = molecule.getChildList().get(oxyThree).getPosition();
        Vector oxyFourPosn = molecule.getChildList().get(oxyFour).getPosition();
        Vector oxyAB = new Vector3D(), oxyABCopy = new Vector3D();
        Vector oxyCD = new Vector3D(), oxyCDCopy = new Vector3D();
        oxyAB.E(oxyOnePosn);
        oxyAB.ME(oxyTwoPosn);
        oxyABCopy.E(oxyAB);
        oxyAB.normalize();
        oxyCD.E(oxyThreePosn);
        oxyCD.ME(oxyFourPosn);
        oxyCDCopy.E(oxyCD);
        oxyCD.normalize();

        double prdDot = oxyAB.dot(oxyCD);
        double theta = Math.acos(prdDot);
        Vector vecABCDCross = new Vector3D();
        vecABCDCross.E(oxyAB);
        vecABCDCross.XE(oxyCD);
        vecABCDCross.normalize();

        rotationalTensor1.setRotationAxis(vecABCDCross, theta);
        molecule.getChildList().forEach(atom -> {
            if(atom.getIndex() == oxyThree || atom.getIndex() == oxyFour) {
               // atom.getPosition().ME(vecANew);
                rotationTensor.transform(atom.getPosition());
               // atom.getPosition().PE(vecANew);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            }
        });
    }

    public ArrayList<Vector> getMoleculeProjection(Space space, Box box,ISpecies species, ArrayList<Integer> arrayListLinkerCarb, ArrayList<Integer> arrayListLinkerOxy, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map ){
        ArrayList<Vector> atomProj = new ArrayList<>();
        if(arrayListLinkerCarb.size() == 3){
            //allign along Z axis
            int carbA = arrayListLinkerCarb.get(0);
            int carbB = arrayListLinkerCarb.get(1);
            int carbC = arrayListLinkerCarb.get(2);
            Vector carbAPosn = species.makeMolecule().getChildList().get(carbA).getPosition();
            Vector carbBPosn = species.makeMolecule().getChildList().get(carbB).getPosition();
            Vector carbCPosn = species.makeMolecule().getChildList().get(carbC).getPosition();
           /// System.out.println("\n "+carbAPosn + " "+ carbBPosn);
            //System.out.println(carbCPosn +"\n");
            IMolecule molecule = species.makeMolecule();

            //Normal to the plane
            Vector vecAB = new Vector3D(), vecBC = new Vector3D(), vecCross = new Vector3D(), vecNormal = new Vector3D(), vecNormalCopy = new Vector3D(), vecCrossPrd = new Vector3D();
            vecAB.E(carbAPosn);
            vecAB.ME(carbBPosn);
            vecBC.E(carbCPosn);
            vecBC.ME(carbBPosn);
            vecCross.E(vecAB);
            vecCross.XE(vecBC);
            vecNormal.E(vecCross);
            vecNormal.normalize();
            double A = vecCross.getX(0);
            double B = vecCross.getX(1);
            double C = vecCross.getX(2);
            //System.out.println(molecule.getChildList().size());

            molecule.getChildList().forEach(atom -> {
                Vector atomPosn = atom.getPosition();
                double xA = atomPosn.getX(0);
                double yB = atomPosn.getX(1);
                double zC = atomPosn.getX(2);
                double t =  - ((A * xA)  + (B *yB) + (C * zC) )/((A*A)+(B*B)+(C*C));
                Vector vecProj = Vector.of((xA+(t*A)), (yB+(t*C)), (zC+(t*C)));
                atomProj.add(vecProj);
               // System.out.println(t);
              //  System.out.println(atomPosn +" "+  vecProj);
            });

        }else {
            throw  new RuntimeException("Code not written " + arrayListLinkerCarb.size());
        }
        return atomProj;
    }
    public Box getMOPEdge(Box box, String strucName,  ISpecies species,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double coordinateMultiplier){
        box.setNMolecules(species, numEdges);
        List<Vector> vecListPoly = getListVecPoly();
        Space space = Space3D.getInstance();
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOP = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOP.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }
      //  System.out.println("inside mopedgfe");
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();
        //gets coordinates depending on struct
        verticesComboFaces(strucName, vertices, faces, faceEdges);
        List<Integer[]> faceOne = faceEdges.get(0);
        double actMult = 100000;
        for (int i =0; i< faceOne.size(); i++){
            //System.out.println(Arrays.toString(faceOne.get(i)));
            Vector vertOne = vertices.get(faceOne.get(i)[0]);
            Vector vertTwo = vertices.get(faceOne.get(i)[1]);
            double distVal = Math.sqrt(vertOne.Mv1Squared(vertTwo));
          //  System.out.println(distVal  +" "+ actMult);
            if (distVal < actMult){
                actMult = distVal;
               // System.out.println("changed " + distVal +" "+ actMult);
            }
        }
     //   System.exit(1);
       // Vector vecOne = vecListPoly.get(0);
        //Vector vecTwo = vecListPoly.get(1);
        //actMult = Math.sqrt(vecOne.Mv1Squared(vecTwo));
        coordinateMultiplier =  coordinateMultiplier/actMult;
        for(int i=0; i<numEdges; i++) {
            //ideal vertice
            Integer[] vecArray = linkedVertices.get(i); //vertices that form an edge
            IMolecule moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vecListPoly.get(vecArray[0] ); // first end of vertice
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne); //convert actual to real vertice
            Vector bAct = new Vector3D();
            bAct.E(vecListPoly.get(vecArray[1] ));
            bAct.TE(coordinateMultiplier);
           // System.out.println("Expected posnt " + copyOriginOne +" "+ bAct);
           // System.out.println(originOne +" "+ Arrays.toString(vecArray));
            //move molecule first extremity to the real vertice

            IMolecule finalMoleculeMOPOne = moleculeMOPOne;
            Vector originStart = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            copyOriginOne.ME(originStart);
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

           // System.out.println("final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition() +" final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition());

            //find angle between vector connecting extremities and vertices that need to be connected
            Vector vecAId = vecListPoly.get(vecArray[0]);
        //    System.out.println("vecAID " +copyOriginOne);
            Vector vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
            Vector vecBId = vecListPoly.get(vecArray[1]);
            Vector vecBIdCopy = new Vector3D();
            vecBIdCopy.Ev1Mv2(vecBId, vecAId);
            //System.out.println(vecAId +" "+vecBId);
            vecBIdCopy.normalize();
            //System.out.println("vecAID norm " + vecBIdCopy);
            Vector vecANew = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vecBNew = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition();
           // System.out.println("vecA "+ vecANew);
            //System.out.println("vecB "+vecBNew);
            Vector vecBCopy = new Vector3D();
            vecBCopy.Ev1Mv2(vecBNew, vecANew);
            vecBCopy.normalize();

            //angle between molecule extremity vector and edge vector
            double prodDot = vecBCopy.dot(vecBIdCopy);
            double theta = Math.acos(prodDot);
            Vector vecBCross = new Vector3D();
            vecBCross.E(vecBCopy);
            vecBCross.XE(vecBIdCopy);
            vecBCross.normalize();

            rotationTensor.setRotationAxis(vecBCross, theta);
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(vecANew);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecANew);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
           // System.out.println("after moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition() +" final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition() +"\n");
            vecArray = linkedVertices.get(i);
            vecAId = vecListPoly.get(vecArray[0]);
            vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

            //edges
            Integer[] vecBArray = linkedVertices.get(i);
            vecBId = vecListPoly.get(vecBArray[1]);
            Vector vecBIdVerify = new Vector3D();
            vecBIdVerify.Ea1Tv1(coordinateMultiplier,vecBId );
            Vector vecABId = new Vector3D();
            Vector vecABVerify = new Vector3D();
            vecABId.Ev1Mv2(vecAIdVerify,vecBIdVerify);
            vecABVerify.E(vecABId);
            vecABVerify.normalize();

            Vector vecCId = new Vector3D();
            Vector vecCVerify = new Vector3D();
            vecCId.Ev1Pv2(vecAIdVerify, vecBIdVerify);
            vecCVerify.E(vecCId);
            vecCId.TE(0.5);
            vecCId.normalize();

            //double cPerp = vecCId.dot(vecABVerify);
            Vector vectorOxyoneCopy = new Vector3D();
            Vector vectorCarCopy = new Vector3D();
            Vector vectorOxyone = finalMoleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
            Vector vectorCarbone = finalMoleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
            vectorOxyoneCopy.E(vectorOxyone);
            vectorCarCopy.E(vectorCarbone);
            vectorOxyoneCopy.ME(vectorCarCopy);

            vectorOxyoneCopy.normalize();
            vectorCarCopy.normalize();
            double project = vectorOxyoneCopy.dot(vecABVerify);
            Vector projection = new Vector3D();
            projection.Ea1Tv1(project, vecABVerify);

            vectorOxyoneCopy.ME(projection);
            vectorOxyoneCopy.normalize();

            double prodDotRotate = vectorOxyoneCopy.dot(vecCId);
            double thetaRotate = Math.acos(prodDotRotate);

            vectorOxyoneCopy.XE(vecCId);
            vectorOxyoneCopy.normalize();

            rotationTensor.setRotationAxis(vectorOxyoneCopy, thetaRotate);
            Vector finalVecAIdVerify = vecAIdVerify;
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(finalVecAIdVerify);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(finalVecAIdVerify);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }
        return box;
    }


    /*



        public Box getMOPEdge(Box box, ISpecies species,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double coordinateMultiplier){
        box.setNMolecules(species, numEdges);
        List<Vector> vecListPoly = getListVecPoly();
        Space space = Space3D.getInstance();
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOP = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOP.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }

        for(int i=0; i<linkedVertices.size(); i++) {
            //ideal vertice
            Integer[] vecArray = linkedVertices.get(i);
            IMolecule moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vecListPoly.get(vecArray[0]);
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne); //convert actual to real vertice

            //move molecule first extremity to the real vertice
            IMolecule finalMoleculeMOPOne = moleculeMOPOne;
            Vector originStart = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            copyOriginOne.ME(originStart);
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            //find angle between vector connecting extremities and vertices that need to be connected
            Vector vecAId = vecListPoly.get(vecArray[0]);
            System.out.println("vecAID " +vecAId);
            Vector vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
            Vector vecBId = vecListPoly.get(vecArray[1]);
            System.out.println("vecBID" +vecBId);
            Vector vecBIdCopy = new Vector3D();
            vecBIdCopy.Ev1Mv2(vecBId, vecAId);
            vecBIdCopy.normalize();
            Vector vecANew = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vecBNew = finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition();
            System.out.println("vecA "+ vecANew);
            System.out.println("vecB "+vecBNew);
            Vector vecBCopy = new Vector3D();
            vecBCopy.Ev1Mv2(vecBNew, vecANew);
            vecBCopy.normalize();

            //angle between molecule extremity vector and edge vector
            double prodDot = vecBCopy.dot(vecBIdCopy);
            double theta = Math.acos(prodDot);
            Vector vecBCross = new Vector3D();
            vecBCross.E(vecBCopy);
            vecBCross.XE(vecBIdCopy);
            vecBCross.normalize();

            rotationTensor.setRotationAxis(vecBCross, theta);
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(vecANew);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecANew);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

          /*  vecArray = linkedVertices.get(i);
            vecAId = vecListPoly.get(vecArray[0]);
            vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

            //edges
            Integer[] vecBArray = linkedVertices.get(i);
            vecBId = vecListPoly.get(vecBArray[1]);
            Vector vecBIdVerify = new Vector3D();
            vecBIdVerify.Ea1Tv1(coordinateMultiplier,vecBId );
            Vector vecABId = new Vector3D();
            Vector vecABVerify = new Vector3D();
            vecABId.Ev1Mv2(vecAIdVerify,vecBIdVerify);
            vecABVerify.E(vecABId);
            vecABVerify.normalize();

            Vector vecCId = new Vector3D();
            Vector vecCVerify = new Vector3D();
            vecCId.Ev1Pv2(vecAIdVerify, vecBIdVerify);
            vecCVerify.E(vecCId);
            vecCId.TE(0.5);
            vecCId.normalize();

            //double cPerp = vecCId.dot(vecABVerify);
            Vector vectorOxyoneCopy = new Vector3D();
            Vector vectorCarCopy = new Vector3D();
            Vector vectorOxyone = finalMoleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
            Vector vectorCarbone = finalMoleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
            vectorOxyoneCopy.E(vectorOxyone);
            vectorCarCopy.E(vectorCarbone);
            vectorOxyoneCopy.ME(vectorCarCopy);

            vectorOxyoneCopy.normalize();
            vectorCarCopy.normalize();
            double project = vectorOxyoneCopy.dot(vecABVerify);
            Vector projection = new Vector3D();
            projection.Ea1Tv1(project, vecABVerify);

            vectorOxyoneCopy.ME(projection);
            vectorOxyoneCopy.normalize();

            double prodDotRotate = vectorOxyoneCopy.dot(vecCId);
            double thetaRotate = Math.acos(prodDotRotate);

            vectorOxyoneCopy.XE(vecCId);
            vectorOxyoneCopy.normalize();

            rotationTensor.setRotationAxis(vectorOxyoneCopy, thetaRotate);
            Vector finalVecAIdVerify = vecAIdVerify;
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(finalVecAIdVerify);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(finalVecAIdVerify);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

}
        System.exit(1);
                return box;
                }
     */

    public Box getMOPFaceMod(Box box,String struc, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, ArrayList<Vector> moleculeProjection){
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();
        Map<Integer, List<Vector>> midPoints = new HashMap<>();
        //gets coordinates depending on struct
        verticesComboFaces(strucName, vertices, faces, faceEdges);
        for (int i =0; i< faceEdges.size(); i++){
            List<Integer[]> facesInOne = faceEdges.get(i);
            List<Vector> midPointPerFace = new ArrayList<>();
            for (int j =0; j< faceEdges.get(i).size(); j++){
                Integer[] verts = facesInOne.get(j);
                int numZero = verts[0];
                int numOne = verts[1];
                Vector vecZero = vertices.get(numZero);
                Vector vecOne = vertices.get(numOne);
                Vector vecMid = new Vector3D();
                vecMid.E(vecZero);
                vecMid.PE(vecOne);
                vecMid.TE(0.5);
                midPointPerFace.add(vecMid);
            }
            midPoints.put(i, midPointPerFace);
        }
        Vector vecOne = midPoints.get(0).get(0);
        Vector vecTwo = midPoints.get(0).get(1);
        double distNominal = Math.sqrt(vecOne.Mv1Squared(vecTwo));
        Vector bentEnd =  new Vector3D();
        double length = bentLength(species.makeMolecule(), arrayListExtremities, bentEnd);
        double coordinateMultiplier = length/distNominal;
        Map<Integer, List<Vector>> midPointsMod = new HashMap<>();
        for (int i =0; i<midPoints.size(); i++){
            List<Vector> facepointMod = new ArrayList<>();
            for (int j =0; j<midPoints.get(i).size(); j++){
                Vector vec = midPoints.get(i).get(j);
                vec.TE(coordinateMultiplier);
                facepointMod.add(vec);
            }
            midPointsMod.put(i, facepointMod);
        }

        for(int i =0; i< midPointsMod.size(); i++){
            System.out.println(Arrays.deepToString(midPointsMod.get(i).toArray()));
        }
        System.exit(1);
        return box;
    }

    public Box getMOPFace(Box box,String struc, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, ArrayList<Vector> moleculeProjection){
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();

        //gets coordinates depending on struct
        verticesComboFaces(strucName, vertices, faces, faceEdges );

        //get dist bw two extremities
        Vector vecVert = vertices.get(0);
        Vector vecVertTwo = vertices.get(1);
        double edgeLen = Math.sqrt(vecVert.Mv1Squared(vecVertTwo));

        Vector vecOne = species.makeMolecule().getChildList().get(arrayListExtremities.get(0)).getPosition();
        Vector vecTwo = species.makeMolecule().getChildList().get(arrayListExtremities.get(1)).getPosition();
        double lenLig = Math.sqrt(vecOne.Mv1Squared(vecTwo));

       //add ligands depending on MOP geometry and modify ideal coordinates
        double coordinateMultiplier = lenLig/ edgeLen;
        if(strucName.equals("tetra")){
            numEdges = 4;
        }else if (strucName.equals("cube")) {
            numEdges = 6;
        } else if (strucName.equals("octa")) {
            numEdges = 8;
        } else if (strucName.equals("icosa")) {
            numEdges = 20;
        } else if (strucName.equals("dodeca")) {
            numEdges = 12;
        }
        //Move all ligands to the zeroth coordinate
        Space space = Space3D.getInstance();
        box.setNMolecules(species, numEdges);
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }
        List<Vector> verticesModified = new ArrayList<>();
        Map<Integer, List<Vector>> bentVertices = new HashMap<>();

        // get vertices for Bent MOP (as centre of edges)
        if(struc.equals("bent")){
            for(int i = 0; i< vertices.size(); i++){
                Vector originOne, copyOriginOne = new Vector3D();
                originOne = vertices.get(i);
                copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
                verticesModified.add(i, copyOriginOne);
            }

            Vector bentEnd =  new Vector3D();
            double length = bentLength(species.makeMolecule(), arrayListExtremities, bentEnd);
            double lengthStandard = calcDist( vertices, faces, faceEdges );
            if(strucName.equals("tetra")){
                coordinateMultiplier = 2 * length/lengthStandard;
            } else {
                throw new RuntimeException("Conversion factor not written");
            }

            //converts vertices of tetra MOP to final position where they need to be translated
            for (int i =0; i<faceEdges.size(); i++){
                List<Vector> newList = new ArrayList<>();
                for (int j= 0; j<faceEdges.get(i).size(); j++){
                    int vertIndOne  = faceEdges.get(i).get(j)[0];
                    int vertIndTwo  = faceEdges.get(i).get(j)[1];
                    Vector vertVecOne = verticesModified.get(vertIndOne);
                    Vector vertVecTwo = verticesModified.get(vertIndTwo);
                    Vector vecMod = new Vector3D();
                    vecMod.E(vertVecOne);
                    vecMod.PE(vertVecTwo);
                    vecMod.TE(0.5);
                    newList.add(vecMod);
                }
                bentVertices.put(i, newList);
            }
            System.out.println("\n");
        }
        // For every face  of MOP
        Map<Integer, List<List<Vector>>> vectorListMap = new HashMap<>();  // vector positions
        Map<Integer, List<List<Integer>>> atomListMap = new HashMap<>();   // atom numbers at the positions

        //Move zeroth extremities atom to the zeroth position of the MOP vertices
        if(struc.equals("face")){
            //Multiply unit coordinates to coordinates of Polyhedra acc to Ligand size
            for(int i = 0; i< vertices.size(); i++){
                Vector originOne, copyOriginOne = new Vector3D();
                originOne = vertices.get(i);
                copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
                verticesModified.add(i, copyOriginOne);
            }

            for(int i = 0; i<numEdges; i++){
                Integer[] vecArray = faces.get(i);
                moleculeMOPOne = box.getMoleculeList().get(i);
                Vector originOne, copyOriginOne = new Vector3D();
                originOne =verticesModified.get(vecArray[0]);
                copyOriginOne.E(originOne);
                Vector moleculeStart = moleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
                copyOriginOne.ME(moleculeStart);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(copyOriginOne);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                everyMolecule(i, faces.get(i), moleculeMOPOne, verticesModified, arrayListExtremities, vectorListMap, atomListMap, false );
            }
        }else{

            for(int i = 0; i<numEdges; i++){
                IMolecule finalmoleculeMOPOne;
                List<Vector> faceVect = bentVertices.get(i);
                finalmoleculeMOPOne = box.getMoleculeList().get(i);
                Vector originOne, copyOriginOne = new Vector3D();
                originOne =faceVect.get(0);
                copyOriginOne.E(originOne);
                // System.out.println("Before "+finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition());
                Vector moleculeStart = finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
                copyOriginOne.ME(moleculeStart);
                finalmoleculeMOPOne.getChildList().forEach(atom -> {
                    oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(copyOriginOne);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                // System.out.println("After "+finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition()+ "\n");
                // finds which extremity of Ligand is closest to which vertice of MOP face
                everyMolecule(i, faces.get(i),  finalmoleculeMOPOne, faceVect, arrayListExtremities, vectorListMap, atomListMap, true);
            }
        }

        for(int i=0; i<4; i++){
            IMolecule moleculeMOPOneNew = box.getMoleculeList().get(i);
            //points on faces List not facelist
            List<List<Vector>> faceList = vectorListMap.get(i);
            List<List<Integer>> ligandList = atomListMap.get(0);

            List<Integer> moleculeOne = ligandList.get(0);
            List<Vector> faceOne = faceList.get(0);
            Vector vecA = faceOne.get(0); // first position of face (already translated)
            Vector vecP = moleculeMOPOneNew.getChildList().get(moleculeOne.get(1)).getPosition(); // first position of molecule (already translated)
            Vector vecACopy = new Vector3D(), vecPCopy = new Vector3D();
            vecACopy.E(vecA);
            vecPCopy.E(vecP);
          //  for(int j = 1; j<faceList.size(); j++){
            int j = 1;
           // if(int j = 1){
                Vector vecB = faceList.get(j).get(0); // second vertice of the face (to be rotated) / nth vertice
                Vector vecQ = moleculeMOPOneNew.getChildList().get(ligandList.get(j).get(1)).getPosition(); // second endpoint of molecule (to be rotated)
                Vector vecAB = new Vector3D(), vecPQ = new Vector3D();
                vecAB.Ev1Mv2(vecB, vecACopy);
                //System.out.println("vecAB : "+vecAB );
                vecAB.normalize();
                vecPQ.Ev1Mv2(vecQ, vecPCopy);
               // System.out.println("vecPQ : "+vecPQ );
                vecPQ.normalize();
                double prodDotRotate = vecAB.dot(vecPQ);
                double thetaRotate = Math.acos(prodDotRotate);
                Vector vecABPQCross = new Vector3D();
                vecABPQCross.E(vecPQ);
                vecABPQCross.XE(vecAB);
                vecABPQCross.normalize();
                rotationTensor.setRotationAxis(vecABPQCross,thetaRotate);
                moleculeMOPOneNew.getChildList().forEach(atom -> {

                  /*  if(atom.getIndex() == arrayListExtremities.get(0)) {
                        return;
                    }*/
                       atom.getPosition().ME(vecACopy);
                       rotationTensor.transform(atom.getPosition());
                       atom.getPosition().PE(vecACopy);
                       Vector shift = box.getBoundary().centralImage(atom.getPosition());
                       atom.getPosition().PE(shift);
                });
            //System.out.println("molB "+ vectorListMap.get(i).get(1).get(0));
            //System.out.println("PosnB "+ moleculeMOPOneNew.getChildList().get(ligandList.get(1).get(1)).getPosition()+"\n");
            // }
           /* System.out.println("molA "+ vectorListMap.get(i).get(0).get(0));
            System.out.println("PosnA "+ moleculeMOPOneNew.getChildList().get(ligandList.get(0).get(1)).getPosition());
            System.out.println("molB "+ vectorListMap.get(i).get(1).get(0));
            System.out.println("PosnB "+ moleculeMOPOneNew.getChildList().get(ligandList.get(1).get(1)).getPosition());
            System.out.println("molC "+ vectorListMap.get(i).get(2).get(0));
            System.out.println("PosnC "+ moleculeMOPOneNew.getChildList().get(ligandList.get(2).get(1)).getPosition());
            System.out.println("\n");*/

           // if(true){
                vecP = moleculeMOPOneNew.getChildList().get(ligandList.get(0).get(1)).getPosition(); // first translated atom
                vecQ = moleculeMOPOneNew.getChildList().get(ligandList.get(1).get(1)).getPosition(); // first rotated atom
                Vector vecR = moleculeMOPOneNew.getChildList().get(ligandList.get(2).get(1)).getPosition();  // third atom not reached final position
                Vector vecC = faceList.get(2).get(0); //final position of third atom
                Vector midPQ = new Vector3D(), midPQCopy = new Vector3D(), vecPQC = new Vector3D(), vecPQR = new Vector3D(), vecPQCNorm = new Vector3D(), vecPQRNorm = new Vector3D();

                //find angle  between plane before (PQR)  and after rotation (PQC)
                //line connecting midpoint of P & Q and initial and final destion of R
                midPQ.E(vecP);
                midPQ.PE(vecQ);
                midPQ.TE(0.5);
                midPQCopy.E(midPQ);
                vecPQC.E(midPQ);
                vecPQC.ME(vecC);
                vecPQCNorm.E(vecPQC);
                vecPQCNorm.normalize();

                vecPQR.E(midPQ);
                vecPQR.ME(vecR);
                vecPQRNorm.E(vecPQR);
                vecPQRNorm.normalize();
                //System.out.println(vecPQCNorm + " "+vecPQRNorm);
                double prodplane = vecPQCNorm.dot(vecPQRNorm);
                double anglePlane = Math.acos(prodplane);
                Vector vecPQRCCross = new Vector3D();
                vecPQRCCross.E(vecPQRNorm);
                vecPQRCCross.XE(vecPQCNorm);
                vecPQRCCross.normalize();

                rotationTensor.setRotationAxis(vecPQRCCross,anglePlane);
                Vector finalVecQ = vecQ;
                moleculeMOPOneNew.getChildList().forEach(atom -> {

                    /*if(atom.getIndex() == ligandList.get(0).get(1) || atom.getIndex() == ligandList.get(1).get(1)) {
                        return;
                    }*/
                    atom.getPosition().ME(vecA);
                    rotationTensor.transform(atom.getPosition());
                    atom.getPosition().PE(vecA);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                //System.out.println(vecPQCNorm + " "+vecPQRNorm);

                /*System.out.println( "after "+ vecR);
                System.out.println("molA "+ vectorListMap.get(i).get(0).get(0));
                System.out.println("PosnA "+ moleculeMOPOneNew.getChildList().get(ligandList.get(0).get(1)).getPosition());
                System.out.println("molB "+ vectorListMap.get(i).get(1).get(0));
                System.out.println("PosnB "+ moleculeMOPOneNew.getChildList().get(ligandList.get(1).get(1)).getPosition());
                System.out.println("molC "+ vectorListMap.get(i).get(2).get(0));
                System.out.println("PosnC "+ moleculeMOPOneNew.getChildList().get(ligandList.get(2).get(1)).getPosition());*/
           // }

        }
       //
        return box;
    }

    public Box getMOPBent(Box box, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, ArrayList<Vector> moleculeProj){
        double coordinateMultiplier = 0;
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();
        verticesComboFaces(strucName, vertices, faces, faceEdges );

        Vector bentEnd =  new Vector3D();
        double length = bentLength(species.makeMolecule(), arrayListExtremities, bentEnd);
        double lengthStandard = calcDist( vertices, faces, faceEdges );
        if(strucName.equals("tetra")){
            numEdges = 4;
            coordinateMultiplier = 2 * length/lengthStandard;
        } else {
            throw new RuntimeException("Conversion factor not written");
        }
       // System.out.println(length +" "+ lengthStandard);
       // System.out.println(coordinateMultiplier);
        // rotateAroundAxis(vecListPoly,coordinateMultiplier);
        //   distFace = calcDistFace(species, arrayListLinkerOxy);

       /* if(strucName.equals("tetra")){
            numEdges = 4;
            dist = 3.0/2 * bentEndLength;
        }else if (strucName.equals("cube")) {
            numEdges = 6;
            dist = 2 * bentEndLength;
        } else if (strucName.equals("octa")) {
            numEdges = 8;
            dist = 3.0/2 * bentEndLength;
        } else if (strucName.equals("icosa")) {
            numEdges = 20;
            dist = 3.0/2 * bentEndLength;
        } else if (strucName.equals("dodeca")) {
            numEdges = 12;
            dist = 2 * Math.sin(Degree.UNIT.toSim(36)) * bentEndLength;
        }*/
        Space space = Space3D.getInstance();
        box.setNMolecules(species, numEdges);
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }

        List<Vector> verticesModified = new ArrayList<>();
        for(int i = 0; i< vertices.size(); i++){
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vertices.get(i);
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
            verticesModified.add(i, copyOriginOne);
        }

        //converts vertices of tetra MOP to final position where they need to be translated
        Map<Integer, List<Vector>> bentVertices = new HashMap<>();
        for (int i =0; i<faceEdges.size(); i++){
            List<Vector> newList = new ArrayList<>();
            for (int j= 0; j<faceEdges.get(i).size(); j++){
                int vertIndOne  = faceEdges.get(i).get(j)[0];
                int vertIndTwo  = faceEdges.get(i).get(j)[1];
                Vector vertVecOne = verticesModified.get(vertIndOne);
                Vector vertVecTwo = verticesModified.get(vertIndTwo);
                Vector vecMod = new Vector3D();
                vecMod.E(vertVecOne);
                vecMod.PE(vertVecTwo);
                vecMod.TE(0.5);
                newList.add(vecMod);
            }
            bentVertices.put(i, newList);
        }
        System.out.println("\n");

        // print all midpoints which are final positions of extremeties
      /*  for(int i = 0; i<faceEdges.size(); i++){
           // System.out.println(Arrays.deepToString(faceEdges.get(i).toArray()));
            System.out.println(Arrays.deepToString(bentVertices.get(i).toArray()));
        }*/

        Map<Integer, List<List<Vector>>> vectorListMap = new HashMap<>();
        Map<Integer, List<List<Integer>>> atomListMap = new HashMap<>();
        for(int i = 0; i<faceEdges.size(); i++){
            IMolecule finalmoleculeMOPOne;
            List<Vector> faceVect = bentVertices.get(i);
            finalmoleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne =faceVect.get(0);
            copyOriginOne.E(originOne);
           // System.out.println("Before "+finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition());
            Vector moleculeStart = finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            copyOriginOne.ME(moleculeStart);
            finalmoleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
           // System.out.println("After "+finalmoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition()+ "\n");
            // finds which extremity of Ligand is closest to which vertice of MOP face
            everyMolecule(i, faces.get(i),  finalmoleculeMOPOne, faceVect, arrayListExtremities, vectorListMap, atomListMap, true);
        }

        for(int i=0; i<faceEdges.size(); i++){
            IMolecule finalmoleculeMOPOne = box.getMoleculeList().get(i);
            // Integer[] face = faces.get(i);
            List<List<Vector>> faceList = vectorListMap.get(i);
            List<Vector> faceOne = faceList.get(0);
            List<List<Integer>> atomList = atomListMap.get(i);
            Vector vecA = faceOne.get(0);
            Vector vecP = faceOne.get(1);
            Vector vecACopy = new Vector3D(), vecPCopy = new Vector3D();
            vecACopy.E(vecA);
            vecPCopy.E(vecP);

            for(int j = 1; j< faceList.size(); j++){
                int numAtom = atomList.get(j).get(1);
                Vector vecB = faceList.get(j).get(0);
                System.out.println(vecB);
                Vector vecQ = faceList.get(j).get(1);
                Vector vecAB = new Vector3D(), vecPQ = new Vector3D();
                vecAB.Ev1Mv2(vecACopy, vecB);// ideal coordinates vector
                vecAB.normalize();
                vecPQ.Ev1Mv2(vecPCopy, vecQ); // ligand vector
                vecPQ.normalize();

                System.out.println("Before " + numAtom + " " + finalmoleculeMOPOne.getChildList().get(numAtom).getPosition() );
                double prodDotRotate = vecPQ.dot(vecAB);
                double thetaRotate = Math.acos(prodDotRotate);// angle of rotation
                Vector vecABPQCross = new Vector3D();
                vecABPQCross.E(vecPQ);
                vecABPQCross.XE(vecAB);
                vecABPQCross.normalize();// axis of rotation

                rotationTensor.setRotationAxis(vecABPQCross, thetaRotate);
               // int atomZero = atomList.get(0).get(0);
               // int atomOne = atomList.get(1).get(0);
               // int finalJ = j;
                finalmoleculeMOPOne.getChildList().forEach(atom -> {
                    /*if(atom.getIndex() == atomZero) {
                        return;
                    }
                    if(finalJ == 2 &&  atom.getIndex() == atomOne){
                        return;
                    }*/
                    //atom.getPosition().ME(vecP);
                    rotationTensor.transform(atom.getPosition());
                    //atom.getPosition().PE(vecP);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                System.out.println("After " + numAtom + " " + finalmoleculeMOPOne.getChildList().get(numAtom).getPosition());

            }
        }
        return box;
    }

    private void actDist(List<Vector> vertices, List<Integer[]> faces, Map<Integer, List<Integer[]>> faceEdges, double divider) {
        Integer[] faceZero = faces.get(0);
        List<Integer[]> faceZeroEdges = faceEdges.get(0);

        System.out.println(Arrays.toString(faceZero));
     //   for (int i =0; i< faceZeroEdges.size(); i++){
            Integer[] faceZeroEdgeZero = faceZeroEdges.get(0);
            System.out.println(Arrays.toString(faceZeroEdgeZero));
            Vector vecVertZero = vertices.get(faceZeroEdgeZero[0]);
            Vector vecVertOne = vertices.get(faceZeroEdgeZero[1]);
            Vector dist = new Vector3D();
            dist.Ev1Mv2(vecVertZero, vecVertOne);
            double distVert = Math.sqrt(dist.squared());
        System.out.println(distVert);
    //    }
      //  System.out.println(faceZeroEdges);
    }

    /*private void stretchMolecule(ISpecies species, ArrayList<Integer> arrayListLinkerOxy) {
        IMolecule molecule = species.makeMolecule();
        Map<Double,Integer[]> distanceMap = new HashMap<>();
        for(int i =0; i< arrayListLinkerOxy.size(); i++){
            for(int j=i+1; j<arrayListLinkerOxy.size(); j++){
                Vector vectorOne = molecule.getChildList().get(arrayListLinkerOxy.get(i)).getPosition();
                Vector vectorTwo = molecule.getChildList().get(arrayListLinkerOxy.get(j)).getPosition();
                double distance = Math.sqrt(vectorOne.Mv1Squared(vectorTwo));
             //   System.out.println( i + " " +j +" " + vectorOne+ " "+ vectorTwo + distance);
                distanceMap.put(distance, new Integer[]{i, j});
            }
        }
        List<Integer[]> adjactntPoints = new ArrayList<>();
        List<Double> sortedKeys = new ArrayList<>(distanceMap.keySet());
        Collections.sort(sortedKeys);
        int n =10;
        Map<Double,Integer[]> smallestDistanceMap = new HashMap<>();
        // Print the n smallest keys and their respective values
        for (int i = 0; i < Math.min(arrayListLinkerOxy.size(), sortedKeys.size()); i++) {
            Double key = sortedKeys.get(i);
            Integer[] value = distanceMap.get(key);
            smallestDistanceMap.put(key, value);
         //   System.out.println("Key: " + key + " Value: " + Arrays.toString(value));
        }
        double largestKey = Integer.MIN_VALUE;
      //  System.out.println(smallestDistanceMap);
        for (Map.Entry<Double, Integer[]> entry : smallestDistanceMap.entrySet()) {
            System.out.println("Key: " + entry.getKey() + " Value: " + Arrays.toString(entry.getValue()));
        }
        Double smallestKey = 0.0;

        double smallestSize = Integer.MIN_VALUE;

        for (Map.Entry<Double, Integer[]> entry : smallestDistanceMap.entrySet()) {
            double currentSize = entry.getKey();
            // Check for smaller size or the same size with a smaller key
            if (smallestSize < currentSize || (currentSize == smallestSize && (smallestKey == null || entry.getKey() > smallestKey))) {
                smallestSize = currentSize;
                smallestKey = entry.getKey();
            }
        }
        double largestKeyFurther = smallestKey;
        Integer[] valLargestKey = smallestDistanceMap.get(largestKeyFurther);
        int zeroLargest = valLargestKey[0];
        int oneLargest = valLargestKey[1];
        System.out.println(Arrays.toString(valLargestKey));
        List<Integer> toCompare = new ArrayList<>();
        List<int[]> toCompareNew = new ArrayList<>();
        System.out.println(smallestKey + " " + Arrays.toString(smallestDistanceMap.get(smallestKey)));

        for(Map.Entry<Double, Integer[]> entry : smallestDistanceMap.entrySet()){
            double currentSize = entry.getKey();
            double diff = Math.abs(currentSize - largestKeyFurther);
            if(diff > 1e-9){

                Integer[] intNums = smallestDistanceMap.get(currentSize);
                int zeroNum = intNums[0];
                int oneNum = intNums[1];
                System.out.println("num " + Arrays.toString(intNums) + " largest "+ Arrays.toString(valLargestKey));

                if(zeroLargest == zeroNum){
                    toCompare.add(oneNum);
                }else if(zeroLargest == oneNum){
                    toCompare.add(zeroNum);
                }else if(oneLargest == zeroNum){
                    toCompare.add(oneNum);
                }else if(oneLargest == oneNum){
                    toCompare.add(zeroNum);
                }
            }
        }
        System.out.println(toCompare);
        Vector vectorA = molecule.getChildList().get(zeroLargest).getPosition();
        Vector vectorB = molecule.getChildList().get(oneLargest).getPosition();
        Vector vectorC = molecule.getChildList().get(toCompare.get(0)).getPosition();
        Vector vectorD = molecule.getChildList().get(toCompare.get(1)).getPosition();
        Vector vectorAB = new Vector3D();
        vectorAB.Ev1Mv2(vectorA, vectorB);
        vectorAB.normalize();
        vectorC.normalize();
        vectorD.normalize();
        double cosThetaCAB = vectorAB.dot(vectorC);
        double cosThetaDAB = vectorAB.dot(vectorD);
      //  System.out.println(toCompareNew);
        System.exit(1);
    }*/



    private double calcDist(List<Vector> vertices, List<Integer[]> faces, Map<Integer, List<Integer[]>> faceEdges) {
        double dist ;
        Integer[] faceOne = faceEdges.get(0).get(0);
        int vertOne = faceOne[0];
        int vertTwo = faceOne[1];
        Vector vecOne = vertices.get(vertOne);
        Vector vecTwo = vertices.get(vertTwo);
        dist = vecTwo.Mv1Squared(vecOne);
        return Math.sqrt(dist);
    }

    public double bentLength(IMolecule molecule, List<Integer> extremeties, Vector bentLength) {
        Vector vecOne = molecule.getChildList().get(extremeties.get(0)).getPosition();
        Vector vecTwo = molecule.getChildList().get(extremeties.get(1)).getPosition();
        return Math.sqrt(vecOne.Mv1Squared(vecTwo));
    }


    private double calcDistFace(ISpecies species, ArrayList<Integer> arrayListLinkerOxy) {
        List<Vector> listVect = new ArrayList<>();
        IMolecule molecule = species.makeMolecule();
        Vector vectMean = new Vector3D();
        for(int i = 0; i< arrayListLinkerOxy.size(); i++){
            Vector vect = molecule.getChildList().get(arrayListLinkerOxy.get(i)).getPosition();
            Vector vectCopy = new Vector3D();
            vectCopy.E(vect);
            vectMean.PE(vectCopy);
            listVect.add(vect);
        }
        //System.out.println(vectMean);
        vectMean.TE(1.0/arrayListLinkerOxy.size());
     //   System.out.println(listVect);
      //  System.out.println(vectMean);
        Vector vectCopy = molecule.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
        double distFace = vectCopy.Mv1Squared(vectMean);
       // System.out.println(distFace);
        return Math.sqrt(distFace);
    }

    private void calcFaceCentroidFace (List<Integer[]>faces,List<Vector> vertices, List<Vector> centres){
        for(int i =0; i< faces.size(); i++){
            Integer[] vecArray = faces.get(i);
            Vector sum= new Vector3D();
            for(int j=0; j<vecArray.length; j++){
                Vector vecOne = new Vector3D();
              //  System.out.println(vecArray[j]);
                vecOne = vertices.get(vecArray[j]);
                sum.PE(vecOne);
            }
            sum.TE(1.0/vecArray.length);
            centres.add(sum);
        }
    }

    public void getAdjacentDistane(int i, List<Vector> vertices, Integer[] face,  Map<Integer, Map<Integer, Double>> largerDistMap, Map<Integer, Map<Integer, Double>> largerAtomMap, ArrayList<Integer> arrayListExtremities, IMolecule molecule){
            Map<Integer, Double> distanceMap = new HashMap<>();
            int vertOne = face[0];
            Vector vertOnePosn = vertices.get(vertOne);
            Vector vertOnePosnCopy = new Vector3D();
            vertOnePosnCopy.E(vertOnePosn);
            for(int j =1; j<face.length; j++){
                int vertTwo = face[j];
                Vector vertTwoPosn = vertices.get(vertTwo);
                Vector vertTwoPosnCopy = new Vector3D();
                vertTwoPosnCopy.E(vertTwoPosn);
                double distEdges = vertTwoPosnCopy.Mv1Squared(vertOnePosnCopy);
                distanceMap.put(vertTwo, distEdges);
              //  System.out.println(distanceMap);
            }
            largerDistMap.put(i, distanceMap);
            Map<Integer, Double> distanceAtomMap = new HashMap<>();
            int numOne = arrayListExtremities.get(0);
            Vector atomOnePosn = molecule.getChildList().get(numOne).getPosition();
            Vector atomOnePosnCopy = new Vector3D();
            atomOnePosnCopy.E(atomOnePosn);
            for(int k = 1; k < arrayListExtremities.size(); k++){
                int numTwo = arrayListExtremities.get(k);
                Vector atomTwoPosn = molecule.getChildList().get(numTwo).getPosition();
                Vector vertTwoPosnCopy = new Vector3D();
                vertTwoPosnCopy.E(atomTwoPosn);
                double distVertAtoms = vertTwoPosnCopy.Mv1Squared(atomOnePosnCopy);
                distanceAtomMap.put(arrayListExtremities.get(k), distVertAtoms);
            }
            largerAtomMap.put(i, distanceAtomMap);
       // System.exit(1);
    }

    public void everyMolecule (int i, Integer[] face, IMolecule molecule, List<Vector> modifiedVertices,ArrayList<Integer> arrayListExtremities, Map<Integer, List<List<Vector>>> map ,   Map<Integer, List<List<Integer>>>  atomListmap, boolean ifBent ){
        List<List<Integer>> faceLigand = new ArrayList<>();
        List<List<Vector>> faceList = new ArrayList<>();
        //ABCDE represent  face alphabets and PQRST represent molecule alphabets
        if (face.length == 3){
            Vector vectorA = new Vector3D(), vectorB= new Vector3D(), vectorC = new Vector3D();
            if(ifBent){
                vectorA = modifiedVertices.get(0);
                vectorB = modifiedVertices.get(1);
                vectorC = modifiedVertices.get(2);
            }else {
                vectorA = modifiedVertices.get(face[0]);
                vectorB = modifiedVertices.get(face[1]);
                vectorC = modifiedVertices.get(face[2]);
            }

            Vector vectorP = molecule.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vectorQ = molecule.getChildList().get(arrayListExtremities.get(1)).getPosition();
            Vector vectorR = molecule.getChildList().get(arrayListExtremities.get(2)).getPosition();
            Vector vectorQCopy = new Vector3D(), vectorBCopy = new Vector3D(), vectorCCopy = new Vector3D(), vectorBQ = new Vector3D(), vectorCQ = new Vector3D();
            vectorBCopy.E(vectorB);
            vectorCCopy.E(vectorC);
            vectorQCopy.E(vectorQ);
            double doubleBQ = vectorBCopy.Mv1Squared(vectorQCopy);
            double doubleCQ = vectorCCopy.Mv1Squared(vectorQCopy);
            List<Vector> listA = new ArrayList<>();
            List<Vector> listB = new ArrayList<>();
            List<Vector> listC = new ArrayList<>();
            listA.add(vectorA);
            listA.add(vectorP);
            listB.add(vectorB);
            listC.add(vectorC);

            List<Integer> listvA = new ArrayList<>();
            List<Integer> listvB = new ArrayList<>();
            List<Integer> listvC = new ArrayList<>();
            listvA.add(face[0]);
            listvA.add(arrayListExtremities.get(0));
            listvB.add(face[1]);
            listvC.add(face[2]);

            if(doubleBQ<doubleCQ){
                listB.add(vectorQ);
                listC.add(vectorR);

                listvB.add(arrayListExtremities.get(1));
                listvC.add(arrayListExtremities.get(2));
            }else {
                listB.add(vectorR);
                listC.add(vectorQ);

                listvB.add(arrayListExtremities.get(2));
                listvC.add(arrayListExtremities.get(1));
            }
            faceList.add(listA);
            faceList.add(listB);
            faceList.add(listC);

            faceLigand.add(listvA);
            faceLigand.add(listvB);
            faceLigand.add(listvC);

            map.put(i, faceList);
            atomListmap.put(i, faceLigand);
        } else if (face.length == 4) {
            //molecule Math
            Vector vectorP = molecule.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vectorQ = molecule.getChildList().get(arrayListExtremities.get(1)).getPosition();
            Vector vectorR= molecule.getChildList().get(arrayListExtremities.get(2)).getPosition();
            Vector vectorS = molecule.getChildList().get(arrayListExtremities.get(3)).getPosition();
            Vector vectorPCopy = new Vector3D();
            vectorPCopy.E(vectorP);
          /*  double doublePQ = vectorPCopy.Mv1Squared(vectorQ);
            double doublePR = vectorPCopy.Mv1Squared(vectorR);
            double doublePS = vectorPCopy.Mv1Squared(vectorS);*/


            //faceMath
            Vector vectorC = new Vector3D(), vectorCCopy = new Vector3D(), vectorD = new Vector3D(), vectorA = new Vector3D(), vectorB = new Vector3D();
            vectorA = modifiedVertices.get(0);
            vectorB = modifiedVertices.get(1);
            vectorC = modifiedVertices.get(2);
            vectorD = modifiedVertices.get(3);

            double doubleBS = vectorB.Mv1Squared(vectorS);
            double doubleBQ = vectorB.Mv1Squared(vectorQ);

            List<Vector> listA = new ArrayList<>();
            List<Vector> listB = new ArrayList<>();
            List<Vector> listC = new ArrayList<>();
            List<Vector> listD = new ArrayList<>();
            listA.add(vectorA);
            listA.add(vectorP);
            listB.add(vectorB);
            listC.add(vectorC);
            listC.add(vectorR);
            listD.add(vectorD);

            if(doubleBQ > doubleBS ){
                listB.add(vectorS);
                listD.add(vectorQ);
                //System.out.println("AP " + " BS " + " DQ " + " CR");
            } else {
                listD.add(vectorS);
                listB.add(vectorQ);
               // System.out.println("AP " + " BQ " + " DS " +" CR");
            }
            faceList.add(listA);
            faceList.add(listB);
            faceList.add(listC);
            faceList.add(listD);
            map.put(i, faceList);

        } else if (face.length == 5) {
            Vector vectorP = molecule.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vectorQ = molecule.getChildList().get(arrayListExtremities.get(1)).getPosition();
            Vector vectorR= molecule.getChildList().get(arrayListExtremities.get(4)).getPosition();
            Vector vectorS = molecule.getChildList().get(arrayListExtremities.get(3)).getPosition();
            Vector vectorT = molecule.getChildList().get(arrayListExtremities.get(2)).getPosition();
            Vector vectorPCopy = new Vector3D();
            vectorPCopy.E(vectorP);
            double doublePQ = vectorPCopy.Mv1Squared(vectorQ);
            double doublePR = vectorPCopy.Mv1Squared(vectorR);
            double doublePS = vectorPCopy.Mv1Squared(vectorS);
            double doublePT = vectorPCopy.Mv1Squared(vectorT);


            //faceMath
            Vector vectorC = new Vector3D(), vectorCCopy = new Vector3D(), vectorD = new Vector3D(), vectorA = new Vector3D(), vectorB = new Vector3D(), vectorE = new Vector3D();
            vectorA = modifiedVertices.get(0);
            vectorB = modifiedVertices.get(2);
            vectorE = modifiedVertices.get(3);
            vectorC = modifiedVertices.get(4);
            vectorD = modifiedVertices.get(5);

            double doubleBT = vectorB.Mv1Squared(vectorT);
            double doubleBQ = vectorB.Mv1Squared(vectorQ);
            double doubleCR = vectorC.Mv1Squared(vectorR);
            double doubleCS = vectorC.Mv1Squared(vectorS);

            List<Vector> listA = new ArrayList<>();
            List<Vector> listB = new ArrayList<>();
            List<Vector> listC = new ArrayList<>();
            List<Vector> listD = new ArrayList<>();
            List<Vector> listE = new ArrayList<>();
            listA.add(vectorA);
            listA.add(vectorP);
            listB.add(vectorB);
            listC.add(vectorC);
            listD.add(vectorD);
            listE.add(vectorE);

            if(doubleBT > doubleBQ ){
                listB.add(vectorQ);
                listE.add(vectorT);
                if(doubleCR > doubleCS){
                    listC.add(vectorS);
                    listD.add(vectorR);
                   // System.out.println("AP " + " BQ " + " CS " + " DR " +" ET ");
                } else {
                    listC.add(vectorR);
                    listD.add(vectorS);
                   // System.out.println("AP " + " BQ " + " CR " + " DS " +" ET ");
                }
            }else {
                listE.add(vectorQ);
                listB.add(vectorT);
                if(doubleCR > doubleCS){
                    listC.add(vectorS);
                    listD.add(vectorR);
                  //  System.out.println("AP " + " BT " + " CS " + " DR " +" EQ");
                } else {
                    listC.add(vectorR);
                    listD.add(vectorS);
                   // System.out.println("AP " + " BT " + " CR " + " DS " +" EQ ");
                }
            }
            faceList.add(listA);
            faceList.add(listB);
            faceList.add(listC);
            faceList.add(listD);
            faceList.add(listE);
            map.put(i, faceList);
        }
    }

    public Map<Integer, Double> sortIntegers(Map<Integer, Double> map){
        List<Map.Entry<Integer, Double>> list = new LinkedList<>(map.entrySet());

        // Sort the list based on values, handling potential floating-point precision issues
        Collections.sort(list, Comparator.comparingDouble(Map.Entry::getValue));

        // Create a new LinkedHashMap to preserve insertion order
        Map<Integer, Double> sortedMap = new LinkedHashMap<>();
        for (Map.Entry<Integer, Double> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }

    public Box getMOPVertice(Box box, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, double coordinateMultiplier){
        double distFace = 0;
        box.setNMolecules(species, numEdges);
        List<Vector> vecListPoly = getListVecPoly();
        rotateAroundAxis(vecListPoly,coordinateMultiplier);
        //   System.exit(1);
        if(strucName.equals("tetra")){
            numEdges = 4;
            dist = 2 * Math.sqrt(3) * distFace;
        } else if (strucName.equals("cube")) {
            numEdges = 6;
            dist = 2 * distFace;
        }else if (strucName.equals("octa")) {
            numEdges = 8;
            dist = 2 * Math.sqrt(3) * distFace;
        }  else if (strucName.equals("icosa")) {
            numEdges = 20;
            dist = 2 * Math.sqrt(3) * distFace;
        } else if (strucName.equals("dodeca")) {
            numEdges = 12;
            dist =  2 * distFace / (Math.tan(Degree.UNIT.toSim(36)));
        }

        Space space = Space3D.getInstance();
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }

        return box;
    }

    public void verticesComboFaces(String strucName, List<Vector> vertices,  List<Integer[]> faces, Map<Integer, List<Integer[]>> faceEdges){
        double phi = ( 1 + Math.sqrt(5)) / 2;
        if(strucName.equals("tetra")){
           // System.out.println("tetra");
            Vector vec0 = Vector.of(1, 1, 1);
            Vector vec1 = Vector.of(1, -1, -1);
            Vector vec2 = Vector.of(-1, 1, -1);
            Vector vec3 = Vector.of(-1, -1, 1);


            vertices.add(vec0);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);

            faces.add(new Integer[]{ 0, 1, 3});
            faces.add(new Integer[]{ 0, 1, 2});
            faces.add(new Integer[]{ 0, 2, 3});
            faces.add(new Integer[]{ 1, 3, 2});

            List<Integer[]> faceZero = new ArrayList<>();
            List<Integer[]> faceOne = new ArrayList<>();
            List<Integer[]> faceTwo = new ArrayList<>();
            List<Integer[]> faceThree = new ArrayList<>();

            faceZero.add(new Integer[]{0,3});
            faceZero.add(new Integer[]{2,0});
            faceZero.add(new Integer[]{2,3});
            faceEdges.put(0, faceZero);

            faceOne.add(new Integer[]{0,1});
            faceOne.add(new Integer[]{0,2});
            faceOne.add(new Integer[]{1,2});
            faceEdges.put(1, faceOne);

            faceTwo.add(new Integer[]{1,2});
            faceTwo.add(new Integer[]{1,3});
            faceTwo.add(new Integer[]{2,3});
            faceEdges.put(2, faceTwo);

            faceThree.add(new Integer[]{0,1});
            faceThree.add(new Integer[]{0,3});
            faceThree.add(new Integer[]{1,3});
            faceEdges.put(3, faceThree);


         /*   faces.add(new Integer[]{ 1, 2, 4});
            faces.add(new Integer[]{ 1, 2, 3});
            faces.add(new Integer[]{ 1, 3, 4});
            faces.add(new Integer[]{ 2, 4, 3});*/

        } else if (strucName.equals("cube") ){
            Vector vec0 = Vector.of(1, 1, 1);
            Vector vec7 = Vector.of(-1, -1, -1);
            Vector vec4 = Vector.of(-1, 1, 1);
            Vector vec2 = Vector.of(1, -1, 1);
            Vector vec1 = Vector.of(1, 1, -1);
            Vector vec6 = Vector.of(-1, -1, 1);
            Vector vec3 = Vector.of(1, -1, -1);
            Vector vec5 = Vector.of(-1, 1, -1);
            vertices.add(vec0);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);
            vertices.add(vec6);
            vertices.add(vec7);


           /* faces.add(new Integer[]{ 0, 2, 3, 5});
            faces.add(new Integer[]{ 0, 3, 4, 6});
            faces.add(new Integer[]{ 1, 6, 7, 4});
            faces.add(new Integer[]{ 1, 2, 5, 7});
            faces.add(new Integer[]{ 1, 5, 3, 6});*/

            faces.add(new Integer[]{ 0, 3, 5, 2});
            faces.add(new Integer[]{ 0, 3, 4, 6});
            faces.add(new Integer[]{ 4, 6, 1, 7});
            faces.add(new Integer[]{ 5, 2, 7, 1});
            faces.add(new Integer[]{ 3, 5, 1, 6});
            faces.add(new Integer[]{ 0, 2, 7, 4});

            List<Integer[]> faceZero = new ArrayList<>();
            List<Integer[]> faceOne = new ArrayList<>();
            List<Integer[]> faceTwo = new ArrayList<>();
            List<Integer[]> faceThree = new ArrayList<>();
            List<Integer[]> faceFour = new ArrayList<>();
            List<Integer[]> faceFive = new ArrayList<>();

            faceZero.add(new Integer[]{ 0,3 });
            faceZero.add(new Integer[]{ 3,5 });
            faceZero.add(new Integer[]{ 5,2 });
            faceZero.add(new Integer[]{0 ,2 });
            faceEdges.put(0, faceZero);

            faceTwo.add(new Integer[]{ 1,6 });
            faceTwo.add(new Integer[]{ 1,7 });
            faceTwo.add(new Integer[]{ 4,6 });
            faceTwo.add(new Integer[]{ 4,7 });
            faceEdges.put(1, faceOne);

            faceTwo.add(new Integer[]{ 0, 2});
            faceTwo.add(new Integer[]{2 , 7});
            faceTwo.add(new Integer[]{ 4, 7});
            faceTwo.add(new Integer[]{ 0, 4});
            faceEdges.put(2, faceTwo);

            faceThree.add(new Integer[]{ 5,2 });
            faceThree.add(new Integer[]{ 2,7 });
            faceThree.add(new Integer[]{ 7,1 });
            faceThree.add(new Integer[]{ 1, 5});
            faceEdges.put(3, faceThree);

            faceFour.add(new Integer[]{ 3,5 });
            faceFour.add(new Integer[]{ 5,1 });
            faceFour.add(new Integer[]{ 1, 6});
            faceFour.add(new Integer[]{ 6, 3});
            faceEdges.put(4, faceFour);

            faceFive.add(new Integer[]{ 0, 2});
            faceFive.add(new Integer[]{ 2, 7});
            faceFive.add(new Integer[]{ 4, 7});
            faceFive.add(new Integer[]{ 0, 4});
            faceEdges.put(5, faceFive);


         /*   faces.add(new Integer[]{ 1, 3, 5, 8});
            faces.add(new Integer[]{ 1, 3, 4, 6});
            faces.add(new Integer[]{ 1, 4, 5, 7});
            faces.add(new Integer[]{ 2, 7, 8, 5});
            faces.add(new Integer[]{ 2, 3, 6, 8});
            faces.add(new Integer[]{ 2, 6, 4, 7});*/


        } else if (strucName.equals("octa")){
            Vector vec0 = Vector.of(1,0 , 0);
            Vector vec1 = Vector.of(-1,0 ,0 );
            Vector vec2 = Vector.of(0,1 , 0);
            Vector vec3 = Vector.of(0, -1, 0);
            Vector vec4 = Vector.of(0,0 , 1);
            Vector vec5 = Vector.of(0,0 , -1);
            vertices.add(vec0);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);

            faces.add(new Integer[]{ 4, 2, 1});
            faces.add(new Integer[]{ 0, 2, 4 });
            faces.add(new Integer[]{ 1, 4, 3});
            faces.add(new Integer[]{ 0, 4, 3});
            faces.add(new Integer[]{ 0, 3, 5});
            faces.add(new Integer[]{ 0, 2, 5});
            faces.add(new Integer[]{ 2, 1, 5});
            faces.add(new Integer[]{ 1, 3, 5});

            List<Integer[]> faceZero = new ArrayList<>();
            List<Integer[]> faceOne = new ArrayList<>();
            List<Integer[]> faceTwo = new ArrayList<>();
            List<Integer[]> faceThree = new ArrayList<>();
            List<Integer[]> faceFour = new ArrayList<>();
            List<Integer[]> faceFive = new ArrayList<>();
            List<Integer[]> faceSix = new ArrayList<>();
            List<Integer[]> faceSeven = new ArrayList<>();

            faceZero.add(new Integer[]{ 4,2 });
            faceZero.add(new Integer[]{ 2, 1});
            faceZero.add(new Integer[]{ 4,1 });
            faceEdges.put(0, faceZero);

            faceOne.add(new Integer[]{ 0,2 });
            faceOne.add(new Integer[]{ 2, 4});
            faceOne.add(new Integer[]{ 0,4 });
            faceOne.add(new Integer[]{ , });
            faceEdges.put(1, faceOne);

            faceTwo.add(new Integer[]{1 , 4});
            faceTwo.add(new Integer[]{ 4, 3});
            faceTwo.add(new Integer[]{ 1, 3});
            faceTwo.add(new Integer[]{ , });
            faceEdges.put(2, faceTwo);

            faceThree.add(new Integer[]{ 0, 4});
            faceThree.add(new Integer[]{ 4, 3});
            faceThree.add(new Integer[]{ 0, 3});
            faceThree.add(new Integer[]{ , });
            faceEdges.put(3, faceThree);

            faceFour.add(new Integer[]{ 0, 3});
            faceFour.add(new Integer[]{ 3, 5});
            faceFour.add(new Integer[]{ 0,5 });
            faceFour.add(new Integer[]{ , });
            faceEdges.put(4, faceFour);

            faceFive.add(new Integer[]{ 0, 2});
            faceFive.add(new Integer[]{ 2,5 });
            faceFive.add(new Integer[]{ 0, 5});
            faceFive.add(new Integer[]{ , });
            faceEdges.put(5, faceFive);

            faceSix.add(new Integer[]{ 2, 1});
            faceSix.add(new Integer[]{ 1,5 });
            faceSix.add(new Integer[]{ 2, 5});
            faceSix.add(new Integer[]{ , });
            faceEdges.put(6, faceSix);

            faceSeven.add(new Integer[]{ 1, 3});
            faceSeven.add(new Integer[]{ 3, 5});
            faceSeven.add(new Integer[]{ 1, 5});
            faceEdges.put(5, faceSeven);


          /*  faces.add(new Integer[]{ 5, 3, 2});
            faces.add(new Integer[]{ 1, 3, 5 });
            faces.add(new Integer[]{ 2, 5, 4});
            faces.add(new Integer[]{ 1, 5, 4});
            faces.add(new Integer[]{ 1, 4, 6});
            faces.add(new Integer[]{ 1, 3, 6});
            faces.add(new Integer[]{ 3, 2, 6});
            faces.add(new Integer[]{ 2, 4, 6});*/

        } else if (strucName.equals("icosa")) {
            Vector vec0 = Vector.of(0,1 , phi);
            Vector vec1 = Vector.of(0,1,-phi );
            Vector vec2 = Vector.of(0,-1 , phi);
            Vector vec3 = Vector.of(0,-1,-phi );
            Vector vec4 = Vector.of(phi, 0, 1);
            Vector vec5 = Vector.of(phi, 0, -1);
            Vector vec6 = Vector.of(-phi, 0, 1);
            Vector vec7 = Vector.of(-phi, 0, -1);
            Vector vec8 = Vector.of(1, phi, 0);
            Vector vec9 = Vector.of(1, -phi, 0);
            Vector vec10 = Vector.of(-1, phi, 0);
            Vector vec11 = Vector.of(-1, -phi, 0);
            vertices.add(vec0);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);
            vertices.add(vec6);
            vertices.add(vec7);
            vertices.add(vec8);
            vertices.add(vec9);
            vertices.add(vec10);
            vertices.add(vec11);

          /*  faces.add(new Integer[]{ 0, 4, 6});
            faces.add(new Integer[]{ 0, 4, 8});
            faces.add(new Integer[]{ 0, 8, 2});
            faces.add(new Integer[]{ 0, 2, 9});
            faces.add(new Integer[]{ 0, 6, 9});
            faces.add(new Integer[]{ 1, 4, 6});
            faces.add(new Integer[]{ 1, 4, 10});
            faces.add(new Integer[]{ 1, 3, 10});
            faces.add(new Integer[]{ 1, 3, 11});
            faces.add(new Integer[]{ 1, 6, 11});
            faces.add(new Integer[]{ 2, 9, 7});
            faces.add(new Integer[]{ 2, 7, 5});
            faces.add(new Integer[]{ 2, 5, 8});
            faces.add(new Integer[]{ 3, 7, 11});
            faces.add(new Integer[]{ 3, 7, 5});
            faces.add(new Integer[]{ 3, 5, 10});
            faces.add(new Integer[]{ 8, 10, 4});
            faces.add(new Integer[]{ 5, 8, 10});
            faces.add(new Integer[]{ 6, 9, 11});
            faces.add(new Integer[]{ 7, 9, 11});


            faces.add(new Integer[]{ 1,  5, 7});
            faces.add(new Integer[]{ 1,  5, 9});
            faces.add(new Integer[]{ 1,  9, 3});
            faces.add(new Integer[]{ 1,  3, 10});
            faces.add(new Integer[]{ 1,  7, 10});
            faces.add(new Integer[]{ 2,  5, 7});
            faces.add(new Integer[]{ 2,  5, 11});
            faces.add(new Integer[]{ 2,  4, 11});
            faces.add(new Integer[]{ 2,  4, 12});
            faces.add(new Integer[]{ 2,  7, 12});
            faces.add(new Integer[]{ 3,  10, 8});
            faces.add(new Integer[]{ 3, 8, 6});
            faces.add(new Integer[]{ 3, 6, 9});
            faces.add(new Integer[]{ 4, 8, 12});
            faces.add(new Integer[]{ 4, 8, 6});
            faces.add(new Integer[]{ 4, 6, 11});
            faces.add(new Integer[]{ 9, 11, 5});
            faces.add(new Integer[]{ 6, 9, 11});
            faces.add(new Integer[]{ 7,  10, 12});
            faces.add(new Integer[]{ 8,  10, 12});*/
            faces.add(new Integer[]{2,0,8 });
            faces.add(new Integer[]{ 8,0,4});
            faces.add(new Integer[]{4,0,6 });
            faces.add(new Integer[]{6,0,9});
            faces.add(new Integer[]{9,0,2});
            faces.add(new Integer[]{ 1,4,10});
            faces.add(new Integer[]{ 6,11,9});
            faces.add(new Integer[]{ 1,11,3});
            faces.add(new Integer[]{7,3,5 });
            faces.add(new Integer[]{1,6,11 });
            faces.add(new Integer[]{2,9,7 });
            faces.add(new Integer[]{ 5,2,7});
            faces.add(new Integer[]{ 8,2,5});
            faces.add(new Integer[]{ 5,8,10});
            faces.add(new Integer[]{8,4,10});
            faces.add(new Integer[]{1,4,6 });
            faces.add(new Integer[]{1,3,10});
            faces.add(new Integer[]{3,7,11 });
            faces.add(new Integer[]{5,3,10 });
            faces.add(new Integer[]{9,7,11});

            List<Integer[]> faceZero = new ArrayList<>();
            List<Integer[]> faceOne = new ArrayList<>();
            List<Integer[]> faceTwo = new ArrayList<>();
            List<Integer[]> faceThree = new ArrayList<>();
            List<Integer[]> faceFour = new ArrayList<>();
            List<Integer[]> faceFive = new ArrayList<>();
            List<Integer[]> faceSix = new ArrayList<>();
            List<Integer[]> faceSeven = new ArrayList<>();
            List<Integer[]> faceEight = new ArrayList<>();
            List<Integer[]> faceNine = new ArrayList<>();
            List<Integer[]> faceTen = new ArrayList<>();
            List<Integer[]> faceEleven = new ArrayList<>();
            List<Integer[]> faceTwelve = new ArrayList<>();
            List<Integer[]> faceThirteen = new ArrayList<>();
            List<Integer[]> faceFourteen = new ArrayList<>();
            List<Integer[]> faceFifteen = new ArrayList<>();
            List<Integer[]> faceSixteen = new ArrayList<>();
            List<Integer[]> faceSeventeen = new ArrayList<>();
            List<Integer[]> faceEighteen = new ArrayList<>();
            List<Integer[]> faceNineteen = new ArrayList<>();

            faceZero.add(new Integer[]{ 2,0 });
            faceZero.add(new Integer[]{ 0,8 });
            faceZero.add(new Integer[]{2 ,8 });
            faceEdges.put(0, faceZero);

            faceOne.add(new Integer[]{ 8, 0});
            faceOne.add(new Integer[]{ 0, 4});
            faceOne.add(new Integer[]{ 8,4 });
            faceEdges.put(1, faceOne);

            faceTwo.add(new Integer[]{ 4,0 });
            faceTwo.add(new Integer[]{ 4,6 });
            faceTwo.add(new Integer[]{6 , 0});
            faceEdges.put(2, faceTwo);

            faceThree.add(new Integer[]{ 6, 0});
            faceThree.add(new Integer[]{ 0,9 });
            faceThree.add(new Integer[]{ 6, 9});
            faceEdges.put(3, faceThree);

            faceFour.add(new Integer[]{ 9, 0});
            faceFour.add(new Integer[]{ 0, 2});
            faceFour.add(new Integer[]{ 9, 2});
            faceEdges.put(4, faceFour);

            faceFive.add(new Integer[]{ 1, 4});
            faceFive.add(new Integer[]{ 4, 10});
            faceFive.add(new Integer[]{ 1, 10});
            faceEdges.put(5, faceFive);

            faceSix.add(new Integer[]{ 6, 11});
            faceSix.add(new Integer[]{ 11,9 });
            faceSix.add(new Integer[]{ 6, 9});
            faceEdges.put(6, faceSix);

            faceSeven.add(new Integer[]{ 1, 11});
            faceSeven.add(new Integer[]{ 11, 3});
            faceSeven.add(new Integer[]{ 1, 3});
            faceEdges.put(7, faceSeven);

            faceEight.add(new Integer[]{ 7, 3});
            faceEight.add(new Integer[]{ 3, 5});
            faceEight.add(new Integer[]{ 7, 5});
            faceEdges.put(8, faceEight);

            faceNine.add(new Integer[]{ 1, 6});
            faceNine.add(new Integer[]{ 6, 11});
            faceNine.add(new Integer[]{1 ,11 });
            faceEdges.put(9, faceNine);

            faceTen.add(new Integer[]{ 2, 9});
            faceTen.add(new Integer[]{ 9, 7});
            faceTen.add(new Integer[]{ 2,7 });
            faceEdges.put(10, faceTen);

            faceEleven.add(new Integer[]{ 5, 2});
            faceEleven.add(new Integer[]{ 2, 7});
            faceEleven.add(new Integer[]{ 5, 7});
            faceEdges.put(11, faceEleven);

            faceTwelve.add(new Integer[]{ 8,2 });
            faceTwelve.add(new Integer[]{2 , 5});
            faceTwelve.add(new Integer[]{ 8,5 });
            faceEdges.put(12, faceTwelve);

            faceThirteen.add(new Integer[]{ 5, 8});
            faceThirteen.add(new Integer[]{ 8,10 });
            faceThirteen.add(new Integer[]{ 5,10 });
            faceEdges.put(13, faceThirteen);

            faceFourteen.add(new Integer[]{ 8, 4});
            faceFourteen.add(new Integer[]{ 4, 10});
            faceFourteen.add(new Integer[]{ 8,10 });
            faceEdges.put(14, faceFourteen);

            faceFifteen.add(new Integer[]{ 1, 4});
            faceFifteen.add(new Integer[]{ 4, 6});
            faceFifteen.add(new Integer[]{ 1, 6});
            faceEdges.put(15, faceFifteen);

            faceSixteen.add(new Integer[]{ 1,3 });
            faceSixteen.add(new Integer[]{ 3, 10});
            faceSixteen.add(new Integer[]{ 1, 10});
            faceEdges.put(16, faceSixteen);

            faceSeventeen.add(new Integer[]{ 3, 7});
            faceSeventeen.add(new Integer[]{ 7, 11});
            faceSeventeen.add(new Integer[]{ 3, 11});
            faceEdges.put(17, faceSeventeen);

            faceEighteen.add(new Integer[]{ 5, 3});
            faceEighteen.add(new Integer[]{ 3, 10});
            faceEighteen.add(new Integer[]{ 5, 10});
            faceEdges.put(18, faceEighteen);

            faceNineteen.add(new Integer[]{ 9,7 });
            faceNineteen.add(new Integer[]{ 7, 11});
            faceNineteen.add(new Integer[]{ 9,11 });
            faceEdges.put(19, faceNineteen);



        } else if (strucName.equals("dodeca")) {
            Vector vec0 = Vector.of(1, 0, phi2);
            Vector vec1 = Vector.of(1, 0, -phi2);
            Vector vec2 = Vector.of(-1, 0, phi2);
            Vector vec3 = Vector.of(-1, 0, -phi2);
            Vector vec4 = Vector.of(phi, phi,phi );
            Vector vec5 = Vector.of(-phi, -phi,-phi );
            Vector vec6 = Vector.of(phi, -phi, -phi);
            Vector vec7 = Vector.of(phi, phi, -phi);
            Vector vec8 = Vector.of(-phi, phi,-phi );
            Vector vec9 = Vector.of(-phi,phi ,phi );
            Vector vec10 = Vector.of(-phi,-phi , phi);
            Vector vec11 = Vector.of(phi,-phi ,phi);
            Vector vec12 = Vector.of(0, phi2,1 );
            Vector vec13 = Vector.of(0, phi2, -1);
            Vector vec14 = Vector.of(0,-phi2 ,1 );
            Vector vec15 = Vector.of(0, -phi2, -1);
            Vector vec16 = Vector.of(phi2, 1, 0);
            Vector vec17 = Vector.of(phi2, -1, 0);
            Vector vec18 = Vector.of(-phi2, 1, 0);
            Vector vec19 = Vector.of(-phi2, -1, 0);
            vertices.add(vec0);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);
            vertices.add(vec6);
            vertices.add(vec7);
            vertices.add(vec8);
            vertices.add(vec9);
            vertices.add(vec10);
            vertices.add(vec11);
            vertices.add(vec12);
            vertices.add(vec13);
            vertices.add(vec14);
            vertices.add(vec15);
            vertices.add(vec16);
            vertices.add(vec17);
            vertices.add(vec18);
            vertices.add(vec19);
          /*  faces.add(new Integer[]{ 10, 14, 19 , 15 ,5 });
            faces.add(new Integer[]{ 10, 19,  2,  9, 18  });
            faces.add(new Integer[]{ 10, 14, 11 , 0, 2 });
            faces.add(new Integer[]{ 11, 17,  6, 15, 14 });
            faces.add(new Integer[]{ 16, 4 , 12, 7, 13  });
            faces.add(new Integer[]{ 17, 16,  6, 7, 1 });
            faces.add(new Integer[]{ 13,  7,  1 , 3,8  });
            faces.add(new Integer[]{ 9,  12,  8, 18,13 });
            faces.add(new Integer[]{ 8,   3,  5, 18, 19  });
            faces.add(new Integer[]{ 0,   2,  4, 9, 12});
            faces.add(new Integer[]{ 0,  11, 16 ,17 , 4  });
            faces.add(new Integer[]{ 1,   3, 5, 15, 6  });*/

            faces.add(new Integer[]{0,11,17,16,4 });
            faces.add(new Integer[]{0,4,12,9,2 });
            faces.add(new Integer[]{0,11,14,10,2 });
            faces.add(new Integer[]{4,16,7,13,12 });
            faces.add(new Integer[]{ 2,9,18,19,10 });
            faces.add(new Integer[]{5,19, 10,14,15});
            faces.add(new Integer[]{ 11,14,15,6,17});
            faces.add(new Integer[]{1,3,8,13,7});
            faces.add(new Integer[]{ 1,7,16,17,6});
            faces.add(new Integer[]{ 1,6,15,5,3});
            faces.add(new Integer[]{8,13,12,9,18});
            faces.add(new Integer[]{ 8,3,5,19,18});

            List<Integer[]> faceZero = new ArrayList<>();
            List<Integer[]> faceOne = new ArrayList<>();
            List<Integer[]> faceTwo = new ArrayList<>();
            List<Integer[]> faceThree = new ArrayList<>();
            List<Integer[]> faceFour = new ArrayList<>();
            List<Integer[]> faceFive = new ArrayList<>();
            List<Integer[]> faceSix = new ArrayList<>();
            List<Integer[]> faceSeven = new ArrayList<>();
            List<Integer[]> faceEight = new ArrayList<>();
            List<Integer[]> faceNine = new ArrayList<>();
            List<Integer[]> faceTen = new ArrayList<>();
            List<Integer[]> faceEleven = new ArrayList<>();

            faceZero.add(new Integer[]{ 0,11 });
            faceZero.add(new Integer[]{ 11,17 });
            faceZero.add(new Integer[]{ 17,16 });
            faceZero.add(new Integer[]{16 ,4 });
            faceZero.add(new Integer[]{0 , 4});
            faceEdges.put(0, faceZero);

            faceOne.add(new Integer[]{0 ,4 });
            faceOne.add(new Integer[]{ 4,12 });
            faceOne.add(new Integer[]{ 12,9 });
            faceOne.add(new Integer[]{ 9,2 });
            faceOne.add(new Integer[]{ 2,0 });
            faceEdges.put(1, faceOne);

            faceTwo.add(new Integer[]{0 ,11 });
            faceTwo.add(new Integer[]{ 11,14 });
            faceTwo.add(new Integer[]{ 14,10 });
            faceTwo.add(new Integer[]{ 10,2 });
            faceTwo.add(new Integer[]{ 0,2 });
            faceEdges.put(2, faceTwo);

            faceThree.add(new Integer[]{ 4,16 });
            faceThree.add(new Integer[]{ 16,7 });
            faceThree.add(new Integer[]{ 7,13 });
            faceThree.add(new Integer[]{ 13,12 });
            faceThree.add(new Integer[]{ 4,12 });
            faceEdges.put(3, faceThree);

            faceFour.add(new Integer[]{ 2,9 });
            faceFour.add(new Integer[]{ 9,18 });
            faceFour.add(new Integer[]{ 18,19 });
            faceFour.add(new Integer[]{ 19,10 });
            faceFour.add(new Integer[]{ 2,10 });
            faceEdges.put(4, faceFour);

            faceFive.add(new Integer[]{ 5,19 });
            faceFive.add(new Integer[]{ 19,10 });
            faceFive.add(new Integer[]{ 10,14 });
            faceFive.add(new Integer[]{14 ,15 });
            faceFive.add(new Integer[]{15 ,5 });
            faceEdges.put(5, faceFive);

            faceSix.add(new Integer[]{ 11,14 });
            faceSix.add(new Integer[]{ 14,15 });
            faceSix.add(new Integer[]{ 15,6 });
            faceSix.add(new Integer[]{ 6,17 });
            faceSix.add(new Integer[]{ 17, 11});
            faceEdges.put(6, faceSix);

            faceSeven.add(new Integer[]{ 1,3 });
            faceSeven.add(new Integer[]{ 3,8 });
            faceSeven.add(new Integer[]{ 8,13 });
            faceSeven.add(new Integer[]{ 13,7 });
            faceSeven.add(new Integer[]{ 7,1 });
            faceEdges.put(7, faceSeven);

            faceEight.add(new Integer[]{1 ,7 });
            faceEight.add(new Integer[]{ 7,16 });
            faceEight.add(new Integer[]{ 16,17 });
            faceEight.add(new Integer[]{ 17,6 });
            faceEight.add(new Integer[]{ 6, 1});
            faceEdges.put(8, faceEight);

            faceNine.add(new Integer[]{ 1,6 });
            faceNine.add(new Integer[]{ 6,15 });
            faceNine.add(new Integer[]{ 15,5 });
            faceNine.add(new Integer[]{ 5,3 });
            faceNine.add(new Integer[]{3 , 1});
            faceEdges.put(9, faceNine);

            faceTen.add(new Integer[]{ 8,13 });
            faceTen.add(new Integer[]{ 13,12 });
            faceTen.add(new Integer[]{ 12,9 });
            faceTen.add(new Integer[]{ 9,18 });
            faceTen.add(new Integer[]{18 ,8 });
            faceEdges.put(10, faceTen);

            faceEleven.add(new Integer[]{ 8, 3});
            faceEleven.add(new Integer[]{ 3, 5});
            faceEleven.add(new Integer[]{ 5,19 });
            faceEleven.add(new Integer[]{19 ,18 });
            faceEleven.add(new Integer[]{ 8, 18});
            faceEdges.put(11, faceEleven);

        /*    faces.add(new Integer[]{ 11, 15 , 20 , 16 ,6 });
            faces.add(new Integer[]{ 11,  20,  3,10, 19  });
            faces.add(new Integer[]{ 11, 15 , 12 , 1, 3 });
            faces.add(new Integer[]{ 12, 18 ,  7, 16, 15 });
            faces.add(new Integer[]{ 17, 5 ,13 ,8, 14  });
            faces.add(new Integer[]{ 18, 17 , 7 , 8, 2 });
            faces.add(new Integer[]{ 14,  8, 2 , 4,9  });
            faces.add(new Integer[]{ 10,  13,  9, 19,14 });
            faces.add(new Integer[]{ 9,  4,  6, 19, 20  });
            faces.add(new Integer[]{ 1,  3,  5, 10, 13 });
            faces.add(new Integer[]{ 1,  12, 17 ,18 , 5  });
            faces.add(new Integer[]{ 2, 4 ,6  ,16, 7  });*/


        }
    }
    //finds which edges are in a MOP
   /* public Map<Double, List<Integer[]>> edgeBetweenPoints ( Map<Double, List<Integer[]>> distMap, List<Vector> listVect){
        List<Integer[]> doubleList = new ArrayList<>();
        Integer[] points = new Integer[2];
        for (int i = 0; i < listVect.size(); i++) {
            Vector pointA = listVect.get(i);
            for (int j = i + 1; j < listVect.size(); j++) { // Start from i + 1 to avoid duplicates
                Vector pointB = listVect.get(j);

                double distance = calculateDistance(pointA, pointB);
                doubleList = distMap.get(distance);
                if(doubleList == null){
                    doubleList = new ArrayList<>();
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }else {
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }
            }
        }
        return distMap;
    }*/

    public List<Integer[]> edgeBetweenPoints(Map<Double, List<Integer[]>> distMap,  List<Vector> listVecRhombic){
        Double smallestKey = null;
        List<Integer[]> doubleList = new ArrayList<>();
        Integer[] points = new Integer[2];
        for (int i = 0; i < listVecRhombic.size(); i++) {
            Vector pointA = listVecRhombic.get(i);
            for (int j = i + 1; j < listVecRhombic.size(); j++) { // Start from i + 1 to avoid duplicates
                Vector pointB = listVecRhombic.get(j);

                double distance = calculateDistance(pointA, pointB);
                doubleList = distMap.get(distance);
                if(doubleList == null){
                    doubleList = new ArrayList<>();
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }else {
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }
            }
        }

        double smallestSize = Integer.MAX_VALUE;

        for (Map.Entry<Double, List<Integer[]>> entry : distMap.entrySet()) {
            double currentSize = entry.getKey();
            // Check for smaller size or the same size with a smaller key
            if (currentSize < smallestSize || (currentSize == smallestSize && (smallestKey == null || entry.getKey() < smallestKey))) {
                smallestSize = currentSize;
                smallestKey = entry.getKey();
            }
        }

        return  distMap.get(smallestKey);
    }

    public double calculateDistance(Vector pointA, Vector pointB){
        return Math.sqrt(pointA.Mv1Squared(pointB));
    }


    public void setVecPositionList(List<Vector> vectPositions){
        this.listVecPoly  = vectPositions;
    }

    public void calcCentroid(List<Vector> polyVec, Vector centroid){
        Vector sum = new Vector3D();
        Vector polyVecEqual = new Vector3D();
        for(int i = 0; i < polyVec.size(); i++){
            polyVecEqual.E(polyVec.get(i));
            centroid.PE(polyVecEqual);
        }
        int size = polyVec.size();
        centroid.TE(1.0 / size);
    }

    public void rotateAroundAxis(List<Vector> listVecPoly, double coordinateMultiplier){
        List<Vector> listMultiplied = new ArrayList<>();
        for(int i = 0; i < listVecPoly.size() ; i++){
            Vector multiplied = new Vector3D();
            multiplied.Ea1Tv1(coordinateMultiplier, listVecPoly.get(i));
            listMultiplied.add(multiplied);
        }
        Vector centroid = new Vector3D();
        calcCentroid(listMultiplied, centroid);
        System.out.println(centroid);
    }

    public List<Vector> getListVecPoly(){
        return listVecPoly;
    }

    // presently considering simple single atom cation
    public void calcExtremities (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMap, ArrayList<Integer> atomNamesExtreme){
        //check for metal atom
        String[] metalNames = {"FE+2", "FE+3", "CR", "MO", "AL",
                "PD", "PT", "CU", "ZN", "ZR", "RH","RU", "NI","CO","PO","MO","MN","MG","V","W","IR","TI"};
        String[] nonMetalNames = {"N_1", "N_2", "N_3", "O_1", "O_2",
                "O_3", "S_2", "S_3", "P_3", "P_5", "O", "N", "S", "P"};
        // also verify it is connected to just one Nitrogen or Oxygen
        for(int i=0; i< connectivityModified.size(); i++){
            ArrayList<Integer> intenalArray = connectivityModified.get(i);
            String atomName = atomMap.get(intenalArray.get(0));
            boolean found = false;
            boolean connectedSingle = false;
            boolean nonMetalConnected = false;
            for (String name : metalNames) {
                if (name.equals(atomName)) {
                    found = true;
                    if(intenalArray.size() ==2){
                        connectedSingle = true;
                        String atomNameConnected = atomMap.get(intenalArray.get(1));
                        for (String nonName : nonMetalNames) {
                            if (nonName.equals(atomNameConnected)) {
                                nonMetalConnected = true;
                                atomNamesExtreme.add(i);
                                break; // Stop searching if a match is found
                            }
                        }
                    }
                    break; // Stop searching if a match is found
                }
            }
        }
    }

    public void locateLinkerOxy (ISpecies species, ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMap, ArrayList<Integer> atomNamesOxygenExtreme, ArrayList<Integer> atomNamesCarbonExtreme){
        String[] metalNames = {"FE+2", "FE+3", "CR", "MO", "AL",
                "PD", "PT", "CU", "ZN", "ZR", "RH","RU", "NI","CO","PO","MO","MN","MG","V","W","IR","TI"};
        int p =0;
        IMolecule molecule = species.makeMolecule();
        ArrayList<Integer> oxyNums = new ArrayList<>();
        for(int i=0; i< connectivityModified.size(); i++){
            ArrayList<Integer> arrayMetal = connectivityModified.get(i);
            String atomName = atomMap.get(arrayMetal.get(0));
            for (String name : metalNames) {
                if (name.equals(atomName)) {
                  //  metalFound = true;
                    //System.out.println(name);
                    if(arrayMetal.size() ==2){
                        int connectivityAdjOxy = arrayMetal.get(1);
                        ArrayList<Integer> arrayAdjOxy = connectivityModified.get(connectivityAdjOxy);
                       // adjOxyFound = true ;
                       // System.out.println("Array Oxy " + arrayAdjOxy + "  array Metal "+  arrayMetal );
                        int connectivityCarb = (arrayMetal.get(0) != arrayAdjOxy.get(1)) ? arrayAdjOxy.get(1) : arrayAdjOxy.get(2);
                        atomNamesCarbonExtreme.add(connectivityCarb);
                        atomNamesOxygenExtreme.add(arrayMetal.get(1));
                       // linkedCFound = true;
                        //ArrayList<Integer> arrayCarb = connectivityModified.get(connectivityCarb);
                        //int doubleBondOxy = ( arrayCarb.get(1) != connectivityAdjOxy) ? arrayCarb.get(1) : arrayCarb.get(2) ;
                       // atomNamesOxygenExtreme.add(doubleBondOxy);
                    }
                    break;
                }/* else if (atomName.equals("O")) {
                    if(arrayMetal.size() == 2){
                        System.out.println("O present "+ p + " "+ i + " "+ arrayMetal.size());
                       // System.out.println("size " );
                        oxyNums.add(i);
                    }
                    p++;
                    break;
                }*/
            }
        }
        List<AtomType> atomTypes = species.getAtomTypes();
        IConformation conformation = species.getConformation();
       /* for ( Integer atomNum : oxyNums){
            //AtomType name = atomTypes.get(atomNum);
           // System.out.println(name);
        }*/
      //  System.exit(1);
    }
}
