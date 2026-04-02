package etomica.GasMOP;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.units.Degree;

import java.util.*;

public class DistCalcCopy {
    double dist=0;
    public double phi = (1+ Math.sqrt(5))/2;
    public double phi2 = phi*phi;
    List<Vector> listVecPoly = new ArrayList<>();
    RotationTensor3D rotationTensor;

    public static void main(String[] args) {
        Map<Double, List<Integer[]>> distMap = new HashMap<>();
        List<Double> slopeXY = new ArrayList<>();
        List<Double> angleZ = new ArrayList<>();

        DistCalcCopy distCalc = new DistCalcCopy();
        // double slp = distCalc.slp(1,1,-1,-1);
        // System.out.println(slp);
        // System.exit(1);
        List<Vector> listVectorTetra = new ArrayList<>();
        Vector vec1 = Vector.of(1, 1, 1);
        Vector vec2 = Vector.of(-1, -1, 1);
        Vector vec3 = Vector.of(-1, 1, -1);
        Vector vec4 = Vector.of(1, -1, -1);
        listVectorTetra.add(vec1);
        listVectorTetra.add(vec2);
        listVectorTetra.add(vec3);
        listVectorTetra.add(vec4);
        distCalc.slopeFinder(distMap, listVectorTetra);
        List<Integer[]>distEqual = distCalc.printResults(distMap, listVectorTetra);
        for(int i=0; i<distEqual.size(); i++){
            System.out.println(Arrays.toString(distEqual.get(i)) + " " + distEqual.get(i)[0] + " "+ distEqual.get(i)[1]   +" " + listVectorTetra.get(distEqual.get(i)[0]) + " "+ listVectorTetra.get(distEqual.get(i)[1]));
        }
        // System.exit(1);
        //  distMap = distCalc.slopeFinder(distMap, listVectorTetra);
        //   distCalc.printResults(distMap);
        // System.exit(1);

        List<Vector> listVecDodeca = new ArrayList<>();
        vec1 = Vector.of(1, 0, distCalc.phi2);
        vec2 = Vector.of(1, 0, -distCalc.phi2);
        vec3 = Vector.of(-1, 0, distCalc.phi2);
        vec4 = Vector.of(-1, 0, -distCalc.phi2);
        Vector vec5 = Vector.of(distCalc.phi, distCalc.phi,distCalc.phi );
        Vector vec6 = Vector.of(-distCalc.phi, -distCalc.phi,-distCalc.phi );
        Vector vec7 = Vector.of(distCalc.phi, -distCalc.phi, -distCalc.phi);
        Vector vec8 = Vector.of(distCalc.phi, distCalc.phi, -distCalc.phi);
        Vector vec9 = Vector.of(-distCalc.phi, distCalc.phi,-distCalc.phi );
        Vector vec10 = Vector.of(-distCalc.phi,distCalc.phi ,distCalc.phi );
        Vector vec11 = Vector.of(-distCalc.phi,-distCalc.phi , distCalc.phi);
        Vector vec12 = Vector.of(distCalc.phi,-distCalc.phi ,distCalc.phi);
        Vector vec13 = Vector.of(0, distCalc.phi2,1 );
        Vector vec14 = Vector.of(0, distCalc.phi2, -1);
        Vector vec15 = Vector.of(0,-distCalc.phi2 ,1 );
        Vector vec16 = Vector.of(0, -distCalc.phi2, -1);
        Vector vec17 = Vector.of(distCalc.phi2, 1, 0);
        Vector vec18 = Vector.of(distCalc.phi2, -1, 0);
        Vector vec19 = Vector.of(-distCalc.phi2, 1, 0);
        Vector vec20 = Vector.of(-distCalc.phi2, -1, 0);
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
        //  System.out.println(listVecDodeca);
        //  System.exit(2);
        // distCalc.slopeFinder(distMap, listVecDodeca);
        // System.out.println(distMap.size());
        //System.exit(1);
        //  distCalc.printResults(distMap, listVecDodeca,  slopeXY, angleZ);
        // System.exit(1);


        List<Vector> listVecCubic = new ArrayList<>();
        vec1 = Vector.of(1, 1, 1);
        vec2 = Vector.of(-1, -1, -1);
        vec3 = Vector.of(-1, 1, 1);
        vec4 = Vector.of(1, -1, 1);
        vec5 = Vector.of(1, 1, -1);
        vec6 = Vector.of(-1, -1, 1);
        vec7 = Vector.of(1, -1, -1);
        vec8 = Vector.of(-1, 1, -1);
        listVecCubic.add(vec1);
        listVecCubic.add(vec2);
        listVecCubic.add(vec3);
        listVecCubic.add(vec4);
        listVecCubic.add(vec5);
        listVecCubic.add(vec6);
        listVecCubic.add(vec7);
        listVecCubic.add(vec8);
        //  distCalc.slopeFinder(distMap, listVecCubic);
        // distCalc.printResults(distMap, listVecCubic);
        //    System.exit(1);


        List<Vector> listVecRhombic = new ArrayList<>();
        vec1 = Vector.of(1,0 , 0);
        vec2 = Vector.of(-1,0 ,0 );
        vec3 = Vector.of(0,1 , 0);
        vec4 = Vector.of(0, -1, 0);
        vec5 = Vector.of(0,0 , 1);
        vec6 = Vector.of(0,0 , -1);
        listVecRhombic.add(vec1);
        listVecRhombic.add(vec2);
        listVecRhombic.add(vec3);
        listVecRhombic.add(vec4);
        listVecRhombic.add(vec5);
        listVecRhombic.add(vec6);
        //  System.out.println(listVecRhombic);
        //distCalc.slopeFinder(distMap, listVecRhombic);
        // distCalc.printResults(distMap, listVecRhombic,  slopeXY, angleZ);
        // System.exit(1);


        List<Vector> listVecIcosa = new ArrayList<>();
        vec1 = Vector.of(0,1 , distCalc.phi);
        vec2 = Vector.of(0, 1, -distCalc.phi);
        vec3 = Vector.of(0, -1, distCalc.phi);
        vec4 = Vector.of(0, -1, -distCalc.phi);
        vec5 = Vector.of(1,distCalc.phi ,0 );
        vec6 = Vector.of(1,-distCalc.phi,0 );
        vec7 = Vector.of(-1, distCalc.phi,0 );
        vec8 = Vector.of(-1,-distCalc.phi,0 );
        vec9 = Vector.of(distCalc.phi, 0,1 );
        vec10 = Vector.of(-distCalc.phi,0 ,1 );
        vec11 = Vector.of(distCalc.phi,0 ,-1 );
        vec12 = Vector.of(-distCalc.phi, 0, -1);
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
        System.out.println(listVecIcosa);
        distCalc.slopeFinder(distMap, listVecIcosa);
        // distCalc.printResults(distMap, listVecIcosa,  slopeXY, angleZ);
        // System.exit(1);
    }

    public List<Integer[]> getStrucList(String struc){
        Map<Double, List<Integer[]>> distMap = new HashMap<>();
        List<Double> slopeXY = new ArrayList<>();
        List<Double> angleZ = new ArrayList<>();
        List<Integer[]>distEqual = new ArrayList<>();
        DistCalcCopy distCalc = new DistCalcCopy();

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
            distCalc.slopeFinder(distMap, listVectorTetra);
            distEqual = distCalc.printResults(distMap, listVectorTetra);
            setVecPositionList(listVectorTetra);
        } else if(struc.equals("dodeca")){
            List<Vector> listVecDodeca = new ArrayList<>();
            Vector vec1 = Vector.of(1, 0, phi2);
            Vector vec2 = Vector.of(1, 0, -phi2);
            Vector vec3 = Vector.of(-1, 0, phi2);
            Vector vec4 = Vector.of(-1, 0, -phi2);
            Vector vec5 = Vector.of(phi, phi,phi );
            Vector vec6 = Vector.of(-phi, -phi,-phi );
            Vector vec7 = Vector.of(phi, -phi, -phi);
            Vector vec8 = Vector.of(phi, phi, -phi);
            Vector vec9 = Vector.of(-phi, phi,-phi );
            Vector vec10 = Vector.of(-phi,phi ,phi );
            Vector vec11 = Vector.of(-phi,-phi , phi);
            Vector vec12 = Vector.of(phi,-phi ,phi);
            Vector vec13 = Vector.of(0, phi2,1 );
            Vector vec14 = Vector.of(0, phi2, -1);
            Vector vec15 = Vector.of(0,-phi2 ,1 );
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
            distCalc.slopeFinder(distMap, listVecDodeca);
            distEqual = distCalc.printResults(distMap, listVecDodeca);
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
            distCalc.slopeFinder(distMap, listVecCubic);
            distEqual = distCalc.printResults(distMap, listVecCubic);
            setVecPositionList(listVecCubic);
        } else if (struc.equals("octa")) {
            List<Vector> listVecRhombic = new ArrayList<>();
            Vector vec1 = Vector.of(1,0 , 0);
            Vector vec2 = Vector.of(-1,0 ,0 );
            Vector vec3 = Vector.of(0,1 , 0);
            Vector vec4 = Vector.of(0, -1, 0);
            Vector vec5 = Vector.of(0,0 , 1);
            Vector vec6 = Vector.of(0,0 , -1);
            listVecRhombic.add(vec1);
            listVecRhombic.add(vec2);
            listVecRhombic.add(vec3);
            listVecRhombic.add(vec4);
            listVecRhombic.add(vec5);
            listVecRhombic.add(vec6);
            distCalc.slopeFinder(distMap, listVecRhombic);
            distEqual = distCalc.printResults(distMap, listVecRhombic);
            setVecPositionList(listVecRhombic);
        } else if (struc.equals("icosa")) {
            List<Vector> listVecIcosa = new ArrayList<>();
            Vector vec1 = Vector.of(0,1 , phi);
            Vector vec2 = Vector.of(0, 1, -phi);
            Vector vec3 = Vector.of(0, -1, phi);
            Vector vec4 = Vector.of(0, -1, -phi);
            Vector vec5 = Vector.of(1,phi ,0 );
            Vector vec6 = Vector.of(1,-phi,0 );
            Vector vec7 = Vector.of(-1, phi,0 );
            Vector vec8 = Vector.of(-1,-phi,0 );
            Vector vec9 = Vector.of(phi, 0,1 );
            Vector vec10 = Vector.of(-phi,0 ,1 );
            Vector vec11 = Vector.of(phi,0 ,-1 );
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
            distCalc.slopeFinder(distMap, listVecIcosa);
            distEqual = distCalc.printResults(distMap, listVecIcosa);
            setVecPositionList(listVecIcosa);
        }
        return distEqual;
    }

    public Box getAutoMOPBox(String struc, String structName, ISpecies species, Box box, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map){
        rotationTensor = new RotationTensor3D();
        List<Integer[]> linkedVertices = getStrucList(structName);
        List<Vector> vecListPolyInit = getListVecPoly();
        //double calcDistance = calculateDistance(vecListPolyInit.get(linkedVertices.get(0)[0]),vecListPolyInit.get(linkedVertices.get(0)[1]));
        // System.out.println(calcDistance);
        //  System.exit(1);
        ArrayList<Integer> arrayListExtremities = new ArrayList<>();
        ArrayList<Integer> arrayListLinkerOxy = new ArrayList<>();
        ArrayList<Integer> arrayListLinkerCarb = new ArrayList<>();
        IMolecule molecule = species.makeMolecule();
        calcMinMaxEdge(connectivity, map, arrayListExtremities);
        locateLinkerOxy(connectivity, map, arrayListLinkerOxy, arrayListLinkerCarb);
        System.out.println(arrayListLinkerOxy);
        System.out.println(arrayListExtremities);
        System.out.println(arrayListLinkerCarb);

        // System.exit(1);
        horzMinMax(molecule, arrayListExtremities);
        dist = getDist();
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
            //  coordinateMultiplier = sqrt3/2 * dist;
            coordinateMultiplier = 0.5 * dist;
        }

        if(struc.equals("edge")){
            getMOPEdge(box,structName, species, connectivity, map,linkedVertices, arrayListExtremities, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, dist, coordinateMultiplier);
        } else if (struc.equals("face")) {
            getMOPFace(box,structName, species, connectivity, map,linkedVertices, arrayListExtremities, arrayListLinkerOxy, arrayListLinkerCarb, numEdges, dist, coordinateMultiplier );
        }else if (struc.equals("bent")) {

        }






     /*   if(false){
            for(int i=0; i<linkedVertices.size(); i++) {
                //edge
                Integer[] vecArray = linkedVertices.get(i);
                Vector vecAId = vecListPoly.get(vecArray[0]);
                Vector vecAIdVerify = new Vector3D();
                vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
                Vector vecBId = vecListPoly.get(vecArray[1]);
                Vector vecBIdCopy = new Vector3D();
                vecBIdCopy.Ev1Mv2(vecBId, vecAId);
                vecBIdCopy.normalize();

                //oxygen vectors
                Vector vectorOxyone = new Vector3D();
                Vector vectorCarbone = new Vector3D();
                Vector vectorOxyoneCopy = new Vector3D();
                Vector vectorCarCopy = new Vector3D();
                moleculeMOPOne = box.getMoleculeList().get(0);
                vectorOxyone = moleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
                vectorCarbone = moleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
                vectorOxyoneCopy.E(vectorOxyone);
                vectorCarCopy.E(vectorCarbone);
                vectorOxyoneCopy.ME(vectorCarCopy);
                vectorOxyoneCopy.normalize();


                double prodDotRotate = vectorOxyoneCopy.dot(vecBIdCopy);
                double thetaRotate = Math.acos(prodDotRotate);
                Vector vecOxyCross = new Vector3D();
                vecOxyCross.E(vectorOxyoneCopy);
                vecOxyCross.XE(vecBIdCopy);
                vecOxyCross.normalize();
                rotationTensor.setRotationAxis(vecOxyCross, thetaRotate);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    if(atom.getIndex() == arrayListExtremities.get(0)) {
                        return;
                    }
                    // atom.getPosition().ME(vecANew);
                    rotationTensor.transform(atom.getPosition());
                    // atom.getPosition().PE(vecANew);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
            }
        }

        if(false){
            for(int i=0; i<linkedVertices.size(); i++) {

                //equation of plane
                Integer[] vecArray = linkedVertices.get(i);
                Vector vecAId = vecListPoly.get(vecArray[0]);
                Vector vecAIdVerify = new Vector3D();
                vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

                //edges
                Integer[] vecBArray = linkedVertices.get(i);
                Vector vecBId = vecListPoly.get(vecBArray[1]);
                Vector vecBIdVerify = new Vector3D();
                vecBIdVerify.Ea1Tv1(coordinateMultiplier,vecBId );
                Vector vecABNormal = new Vector3D();
                vecABNormal.E(vecAIdVerify);

                Vector vecCId = new Vector3D(0,0,0);

                vecABNormal.XE(vecAIdVerify);
                Vector vecBIdC = new Vector3D(), vecBIdC2 = new Vector3D(), vecBIdC3 = new Vector3D();

                vecBIdC.E(vecABNormal);
              //  vecBIdC.ME(vecCId);
                vecBIdC.normalize();

                //oxygen lines
                Vector vectorOxyoneCopy = new Vector3D();
                Vector vectorCarCopy = new Vector3D();
                moleculeMOPOne = box.getMoleculeList().get(0);
                Vector vectorOxyone = moleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
                Vector vectorCarbone = moleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
                vectorOxyoneCopy.E(vectorOxyone);
                vectorCarCopy.E(vectorCarbone);
                vectorOxyoneCopy.ME(vectorCarCopy);
                vectorOxyoneCopy.normalize();

                double prodDotRotate = vectorOxyoneCopy.dot(vecBIdC);
                double thetaRotate = Math.acos(prodDotRotate);
                Vector vecOxyCross = new Vector3D();
                vecOxyCross.E(vectorOxyoneCopy);
                vecOxyCross.XE(vecBIdC);
                vecOxyCross.normalize();
                rotationTensor.setRotationAxis(vecOxyCross, thetaRotate);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    if(atom.getIndex() == arrayListExtremities.get(0)) {
                        return;
                    }
                    // atom.getPosition().ME(vecANew);
                    rotationTensor.transform(atom.getPosition());
                    // atom.getPosition().PE(vecANew);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });

            }

        }

        if(false){
            for(int i=0; i<linkedVertices.size(); i++) {
                Integer[] vecArray = linkedVertices.get(i);
                Vector vecAId = vecListPoly.get(vecArray[0]);
                Vector vecAIdVerify = new Vector3D();
                vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

                //edges
                Integer[] vecBArray = linkedVertices.get(i);
                Vector vecBId = vecListPoly.get(vecBArray[1]);
                Vector vecBIdVerify = new Vector3D();
                vecBIdVerify.Ea1Tv1(coordinateMultiplier,vecBId );
                Vector vecABId = new Vector3D();
                Vector vecABVerify = new Vector3D();
                vecABId.Ev1Mv2(vecAIdVerify,vecBIdVerify);
                vecABVerify.E(vecABId);

                Vector vecCId = new Vector3D();
                Vector vecCVerify = new Vector3D();
                vecCId.Ev1Pv2(vecAIdVerify, vecBIdVerify);
                vecCVerify.E(vecCId);
                vecCId.TE(0.5);
                System.out.println(vecABVerify + " " + vecCVerify +" " + vecCId);


                Vector vectorOxyoneCopy = new Vector3D();
                Vector vectorCarCopy = new Vector3D();
                moleculeMOPOne = box.getMoleculeList().get(0);
                Vector vectorOxyone = moleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
                Vector vectorCarbone = moleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
                vectorOxyoneCopy.E(vectorOxyone);
                vectorCarCopy.E(vectorCarbone);
                vectorOxyoneCopy.ME(vectorCarCopy);
                vectorOxyoneCopy.normalize();

                double prodDotRotate = vectorOxyoneCopy.dot(vecCId);
                double thetaRotate = Math.acos(prodDotRotate);
                Vector vecOxyCross = new Vector3D();
                vecOxyCross.E(vectorOxyoneCopy);
                vecOxyCross.XE(vecCId);
                vecOxyCross.normalize();
                rotationTensor.setRotationAxis(vecOxyCross, thetaRotate);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    if(atom.getIndex() == arrayListExtremities.get(0)) {
                        return;
                    }
                    // atom.getPosition().ME(vecANew);
                    rotationTensor.transform(atom.getPosition());
                    // atom.getPosition().PE(vecANew);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });

            }
           // System.exit(1);

        }*/
        //    System.exit(1);

        //System.exit(1);
        return box;
    }



    public Box getMOPEdge(Box box, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map, List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, double coordinateMultiplier){
        box.setNMolecules(species, numEdges);
        List<Vector> vecListPoly = getListVecPoly();
        // rotateAroundAxis(vecListPoly,coordinateMultiplier);
        Space space = Space3D.getInstance();
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }

        for(int i=0; i<linkedVertices.size(); i++) {
            Integer[] vecArray = linkedVertices.get(i);
            moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vecListPoly.get(vecArray[0]);
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
            Vector originStart = moleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            copyOriginOne.ME(originStart);
            moleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            Vector vecAId = vecListPoly.get(vecArray[0]);
            Vector vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
            Vector vecBId = vecListPoly.get(vecArray[1]);
            Vector vecBIdCopy = new Vector3D();
            vecBIdCopy.Ev1Mv2(vecBId, vecAId);
            vecBIdCopy.normalize();

            moleculeMOPOne = box.getMoleculeList().get(i);
            Vector vecANew = moleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            Vector vecBNew = moleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition();
            Vector vecBCopy = new Vector3D();
            vecBCopy.Ev1Mv2(vecBNew, vecANew);
            vecBCopy.normalize();

            double prodDot = vecBCopy.dot(vecBIdCopy);
            double theta = Math.acos(prodDot);
            Vector vecBCross = new Vector3D();
            vecBCross.E(vecBCopy);
            vecBCross.XE(vecBIdCopy);
            vecBCross.normalize();
            rotationTensor.setRotationAxis(vecBCross, theta);
            moleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(vecANew);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecANew);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

        }

        for(int i=0; i<linkedVertices.size(); i++) {
            Integer[] vecArray = linkedVertices.get(i);
            Vector vecAId = vecListPoly.get(vecArray[0]);
            Vector vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

            //edges
            Integer[] vecBArray = linkedVertices.get(i);
            Vector vecBId = vecListPoly.get(vecBArray[1]);
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

            // System.out.println(vecABVerify + " " + vecCVerify +" " + vecCId);

            Vector vectorOxyoneCopy = new Vector3D();
            Vector vectorCarCopy = new Vector3D();
            moleculeMOPOne = box.getMoleculeList().get(0);
            Vector vectorOxyone = moleculeMOPOne.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
            Vector vectorCarbone = moleculeMOPOne.getChildList().get(arrayListLinkerCarb.get(0)).getPosition();
            vectorOxyoneCopy.E(vectorOxyone);
            vectorCarCopy.E(vectorCarbone);
            vectorOxyoneCopy.ME(vectorCarCopy);
            vectorOxyoneCopy.normalize();

            double prodDotRotate = vectorOxyoneCopy.dot(vecCId);
            double thetaRotate = Math.acos(prodDotRotate);
            Vector vecOxyCross = new Vector3D();
            vecOxyCross.E(vectorOxyoneCopy);
            vecOxyCross.XE(vecCId);
            vecOxyCross.normalize();
            rotationTensor.setRotationAxis(vecABVerify, thetaRotate);

            moleculeMOPOne.getChildList().forEach(atom -> {
                if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }
                atom.getPosition().ME(vecAIdVerify);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecAIdVerify);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }
        return box;
    }

    public void getAdjacentDistane(List<Vector> vertices, List<Integer[]> faces,  Map<Integer, Double> distanceMap, Map<Integer, Double> distanceAtomMap, ArrayList<Integer> arrayListExtremities, IMolecule molecule){
        for(int i = 0; i < faces.size(); i++){
            Integer[] face = faces.get(i);
            int vertOne = face[0];
            Vector vertOnePosn = vertices.get(vertOne);
            Vector vertOnePosnCopy = new Vector3D();
            vertOnePosnCopy.E(vertOnePosn);
            for(int j =1; j<faces.size(); j++){
                int vertTwo = face[j];
                Vector vertTwoPosn = vertices.get(vertTwo);
                Vector vertTwoPosnCopy = new Vector3D();
                vertTwoPosnCopy.E(vertTwoPosn);
                double distEdges = vertTwoPosnCopy.Mv1Squared(vertOnePosnCopy);
                distanceMap.put(vertTwo, distEdges);
            }
        }
        int numOne = arrayListExtremities.get(0);
        Vector atomOnePosn = molecule.getChildList().get(numOne).getPosition();
        Vector atomOnePosnCopy = new Vector3D();
        atomOnePosnCopy.E(atomOnePosn);
        for(int i = 1; i < arrayListExtremities.size(); i++){
            Vector atomTwoPosn = molecule.getChildList().get(numOne).getPosition();
            Vector vertTwoPosnCopy = new Vector3D();
            vertTwoPosnCopy.E(atomTwoPosn);
            double distVertAtoms = vertTwoPosnCopy.Mv1Squared(atomOnePosnCopy);
            distanceAtomMap.put(arrayListExtremities.get(i), distVertAtoms);
        }
    }


    public Box getMOPFace(Box box, String strucName, ISpecies species, ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> map,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double dist, double coordinateMultiplier){
        double distFace = 0;
        List<Integer[]> faces = new ArrayList<>();
        List<Vector> vertices = new ArrayList<>();
        verticesComboFaces(strucName, vertices, faces );
        ArrayList<Vector> meanVectorList = new ArrayList<>();
        // rotateAroundAxis(vecListPoly,coordinateMultiplier);
        distFace = calcDistFace(species, arrayListLinkerOxy);

        if(strucName.equals("tetra")){
            numEdges = 4;
            dist = 2 * Math.sqrt(2.0/3.0) * distFace;
        }else if (strucName.equals("cube")) {
            numEdges = 6;
            dist = Math.sqrt(2) * distFace;
        } else if (strucName.equals("octa")) {
            numEdges = 8;
            dist = 2 * Math.sqrt(2.0/3.0) * distFace;
        } else if (strucName.equals("icosa")) {
            numEdges = 20;
            dist = 2 * Math.sqrt(2.0/3.0) * distFace;
        } else if (strucName.equals("dodeca")) {
            numEdges = 12;
            dist = 2 * Math.sin(Degree.UNIT.toSim(36)) * distFace;
        }
        Space space = Space3D.getInstance();
        box.setNMolecules(species, numEdges);
        List<Vector> oldPositionsTwo = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsTwo.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsTwo.add(space.makeVector());
        }
        List<Vector> centre = new ArrayList<>();
        calcFaceCentroidFace(faces, vertices, centre);
        for(int i = 0; i<faces.size(); i++){
            Integer[] vecArray = linkedVertices.get(i);
            moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vertices.get(vecArray[0]);
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne);
            Vector originStart = moleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition();
            copyOriginOne.ME(originStart);
            moleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            Map<Integer, Double> distanceMap = new HashMap<>();
            Map<Integer, Double> distanceAtomMap = new HashMap<>();
            getAdjacentDistane(vertices, faces, distanceMap, distanceAtomMap, arrayListExtremities, moleculeMOPOne );

            for(int j = 0; j < map.size(); j++ ){
                Set<Integer> keys = distanceMap.keySet();
                int keyOne = (int) keys.toArray()[j];
            }

        }
        return box;
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
        System.out.println(vectMean);
        vectMean.TE(1.0/arrayListLinkerOxy.size());
        System.out.println(listVect);
        System.out.println(vectMean);
        Vector vectCopy = molecule.getChildList().get(arrayListLinkerOxy.get(0)).getPosition();
        double distFace = vectCopy.Mv1Squared(vectMean);
        System.out.println(distFace);
        return Math.sqrt(distFace);
    }

    private void calcFaceCentroidFace (List<Integer[]>faces,List<Vector> vertices, List<Vector> cetres){
        for(int i =0; i< faces.size(); i++){
            Integer[] vecArray = faces.get(i);
            Vector sum= new Vector3D();
            for(int j=0; j<vecArray.length; j++){
                Vector vecOne = new Vector3D();
                //  System.out.println(vecArray[j]);
                vecOne = vertices.get(vecArray[j]-1);
                sum.PE(vecOne);
            }
            sum.TE(1.0/vecArray.length);
            cetres.add(sum);
        }
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



    public void verticesComboFaces(String strucName, List<Vector> vertices,  List<Integer[]> faces){
        double phi = ( 1 + Math.sqrt(5)) / 2;
        if(strucName.equals("tetra")){
            Vector vec1 = Vector.of(1, 1, 1);
            Vector vec4 = Vector.of(-1, -1, 1);
            Vector vec3 = Vector.of(-1, 1, -1);
            Vector vec2 = Vector.of(1, -1, -1);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            faces.add(new Integer[]{1, 2, 3});
            faces.add(new Integer[]{1, 3, 4});
            faces.add(new Integer[]{1, 2, 4});
            faces.add(new Integer[]{2, 3, 4});
        } else if (strucName.equals("cube") ){
            Vector vec1 = Vector.of(1, 1, 1);
            Vector vec8 = Vector.of(-1, -1, -1);
            Vector vec5 = Vector.of(-1, 1, 1);
            Vector vec3 = Vector.of(1, -1, 1);
            Vector vec2 = Vector.of(1, 1, -1);
            Vector vec7 = Vector.of(-1, -1, 1);
            Vector vec4 = Vector.of(1, -1, -1);
            Vector vec6 = Vector.of(-1, 1, -1);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);
            vertices.add(vec6);
            vertices.add(vec7);
            vertices.add(vec8);
            faces.add(new Integer[]{1, 2, 4, 3});
            faces.add(new Integer[]{5, 6, 8, 7});
            faces.add(new Integer[]{1, 5, 7, 3});
            faces.add(new Integer[]{2, 6, 8, 4});
            faces.add(new Integer[]{1, 5, 6, 2});
            faces.add(new Integer[]{3, 7, 8, 4});

        } else if (strucName.equals("octa")){
            Vector vec1 = Vector.of(1,0 , 0);
            Vector vec2 = Vector.of(-1,0 ,0 );
            Vector vec3 = Vector.of(0,1 , 0);
            Vector vec4 = Vector.of(0, -1, 0);
            Vector vec5 = Vector.of(0,0 , 1);
            Vector vec6 = Vector.of(0,0 , -1);
            vertices.add(vec1);
            vertices.add(vec2);
            vertices.add(vec3);
            vertices.add(vec4);
            vertices.add(vec5);
            vertices.add(vec6);
            faces.add(new Integer[]{1, 3, 5});
            faces.add(new Integer[]{1, 3, 6});
            faces.add(new Integer[]{1, 4, 5});
            faces.add(new Integer[]{1, 4, 6});
            faces.add(new Integer[]{2, 3, 5});
            faces.add(new Integer[]{2, 3, 6});
            faces.add(new Integer[]{2, 4, 5});
            faces.add(new Integer[]{2, 4, 6});

        } /*else if (strucName.equals("icosa")) {

            vertices.add(new Double[]{1, phi, 0});
            vertices.add(new Double[]{-1, (int) Math.round(1.618), 0});
            vertices.add(new Double[]{1, (int) Math.round(-1.618), 0});
            vertices.add(new Double[]{-1, (int) Math.round(-1.618), 0});
            vertices.add(new Double[]{(int) Math.round(1.618), 0, 1});
            vertices.add(new Double[]{(int) Math.round(1.618), 0, -1});
            vertices.add(new Double[]{(int) Math.round(-1.618), 0, 1});
            vertices.add(new Double[]{(int) Math.round(-1.618), 0, -1});
            vertices.add(new Double[]{0, 1, (int) Math.round(1.618)});
            vertices.add(new Double[]{0, -1, (int) Math.round(1.618)});
            vertices.add(new Double[]{0, 1, (int) Math.round(-1.618)});
            vertices.add(new Double[]{0, -1, (int) Math.round(-1.618)});

            faces.add(new Integer[]{1, 9, 5});
            faces.add(new Integer[]{1, 5, 3});
            faces.add(new Integer[]{1, 3, 7});
            faces.add(new Integer[]{1, 7, 9});
            faces.add(new Integer[]{1, 9, 5});
            faces.add(new Integer[]{5, 7, 9});
            faces.add(new Integer[]{8, 4, 11});
            faces.add(new Integer[]{6, 12, 7});
            faces.add(new Integer[]{3, 7, 12});
            faces.add(new Integer[]{3, 12, 11});
            faces.add(new Integer[]{2, 6, 8});
            faces.add(new Integer[]{8, 4, 11});
            faces.add(new Integer[]{6, 8, 7});
            faces.add(new Integer[]{6, 7, 12});
            faces.add(new Integer[]{4, 3, 7});
            faces.add(new Integer[]{12, 7, 8});
            faces.add(new Integer[]{12, 3, 7});
            faces.add(new Integer[]{9, 11, 5});
            faces.add(new Integer[]{12, 1, 3});
            faces.add(new Integer[]{1, 9, 5});
        } else if (strucName.equals("dodeca")) {

        }*/
    }
    public Map<Double, List<Integer[]>> slopeFinder ( Map<Double, List<Integer[]>> distMap, List<Vector> listVect){
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
                //   System.out.println("Distance between point " + i + " and point " + j + ": " + distance);
            }
        }
        return distMap;
    }

    public List<Integer[]> printResults(Map<Double, List<Integer[]>> distMap,  List<Vector> listVecRhombic){
        Double smallestKey = null;
        System.out.println(  distMap.keySet() );
        //  System.exit(1);
        double smallestSize = Integer.MAX_VALUE;

        for (Map.Entry<Double, List<Integer[]>> entry : distMap.entrySet()) {
            double currentSize = entry.getKey();
            // Check for smaller size or the same size with a smaller key
            if (currentSize < smallestSize || (currentSize == smallestSize && (smallestKey == null || entry.getKey() < smallestKey))) {
                smallestSize = currentSize;
                smallestKey = entry.getKey();
            }
        }

        if (smallestKey != null) {
            System.out.println("Key with the smallest list size: " + smallestKey);
        } else {
            System.out.println("Map is empty.");
        }
        List<Integer[]> valueList = distMap.get(smallestKey);
        return  valueList;
    }
    public void printVal( List<Integer[]> valueList,  List<Vector> listVecRhombic, List<Double> slopeXY, List<Double> angleZ){
        for(int i=0; i<valueList.size(); i++){
            System.out.println(Arrays.toString(valueList.get(i)) + " " + valueList.get(i)[0] + " "+ valueList.get(i)[1]   +" " + listVecRhombic.get(valueList.get(i)[0]) + " "+ listVecRhombic.get(valueList.get(i)[1]));
        }
        //    System.exit(1);
    }
    public double slp(double x1, double y1, double x2, double y2){
        return (y2-y1)/(x2-x1);
    }



    public Vector boxSizeVector (Map<Integer, Vector> positions){
        double minX = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY;
        double maxZ = Double.NEGATIVE_INFINITY;
        Vector boxSize = new Vector3D();
        for (Vector vec : positions.values()) {
            if (vec.getX(0) < minX) minX = vec.getX(0) ;
            if (vec.getX(0)  > maxX) maxX = vec.getX(0) ;
            if (vec.getX(1) < minY) minY = vec.getX(1);
            if (vec.getX(1) > maxY) maxY = vec.getX(1);
            if (vec.getX(2) < minZ) minZ = vec.getX(2);
            if (vec.getX(2) > maxZ) maxZ = vec.getX(2);
        }
        double valX = (maxX - minX) * 1.1;
        double valY = (maxY - minY) * 1.1;
        double valZ = (maxZ - minZ) * 1.1;
        boxSize = Vector.of(valX, valY, valZ);
        return boxSize;
    }

    public List<Vector> horzMinMax(IMolecule molecule, List<Integer> arrayListExtremeties){
        List<Vector> vecList = new ArrayList<>();
        List<Integer> vecNum = new ArrayList<>();
        Vector vecNew = molecule.getChildList().get(arrayListExtremeties.get(0)).getPosition();
        Vector maxVec = vecNew;
        int numMax = arrayListExtremeties.get(0);
        int numMin = arrayListExtremeties.get(1);
        Vector minVec =  molecule.getChildList().get(arrayListExtremeties.get(1)).getPosition();
        vecList.add(minVec);
        vecList.add(maxVec);
        vecNum.add(numMin);
        vecNum.add(numMax);
        setDist(minVec, maxVec);
        return vecList;
    }

    public void setDist(Vector minVec, Vector maxVec){
        minVec.ME(maxVec);
        this.dist = Math.sqrt(minVec.squared());
    }

    public double getDist(){return dist;}

    public double truncRad(Vector dist) {
        double maxValue = Math.max(dist.getX(0), Math.max(dist.getX(1), dist.getX(2)));
        return maxValue/2;
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


    public void setMapLinkers (Map<Integer, ArrayList<Integer>>map, int keyOne , int keyTwo){
        ArrayList<Integer> two = map.get(keyOne);
        two.add(keyTwo);
        ArrayList<Integer> one = map.get(keyTwo);
        one.add(keyOne);
    }

    // presently cconsidering simple single atom cation
    public void calcMinMaxEdge (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMap, ArrayList<Integer> atomNamesExtreme){
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

    public void locateLinkerOxy (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMap, ArrayList<Integer> atomNamesOxygenExtreme, ArrayList<Integer> atomNamesCarbonExtreme){
        String[] metalNames = {"FE+2", "FE+3", "CR", "MO", "AL",
                "PD", "PT", "CU", "ZN", "ZR", "RH","RU", "NI","CO","PO","MO","MN","MG","V","W","IR","TI"};
        for(int i=0; i< connectivityModified.size(); i++){
            ArrayList<Integer> arrayMetal = connectivityModified.get(i);
            String atomName = atomMap.get(arrayMetal.get(0));
            //  boolean metalFound = false;
            //  boolean adjOxyFound = false;
            // boolean linkedCFound = false;
            //   boolean doubBondOxyFound = false;
            for (String name : metalNames) {
                if (name.equals(atomName)) {
                    //  metalFound = true;
                    if(arrayMetal.size() ==2){
                        int connectivityAdjOxy = arrayMetal.get(1);
                        ArrayList<Integer> arrayAdjOxy = connectivityModified.get(connectivityAdjOxy);
                        // adjOxyFound = true ;
                        int connectivityCarb = (arrayMetal.get(0) != arrayAdjOxy.get(1)) ? arrayAdjOxy.get(1) : arrayAdjOxy.get(2);
                        atomNamesCarbonExtreme.add(connectivityCarb);
                        // linkedCFound = true;
                        ArrayList<Integer> arrayCarb = connectivityModified.get(connectivityCarb);
                        int doubleBondOxy = ( arrayCarb.get(1) != connectivityAdjOxy) ? arrayCarb.get(1) : arrayCarb.get(2) ;
                        atomNamesOxygenExtreme.add(doubleBondOxy);
                    }
                    break;
                }
            }
        }
    }
}
