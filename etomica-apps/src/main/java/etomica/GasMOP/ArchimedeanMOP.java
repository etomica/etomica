package etomica.GasMOP;

import etomica.space.Vector;

import java.util.*;

public class ArchimedeanMOP {
    public static void main(String[] args) {
        AutoMOP autoMOP = new AutoMOP();
        List<Integer[]> distEqual = new ArrayList<>();
        List<List<Integer>> faceList = new ArrayList<>();
        Map<Integer, Vector > verticeMap = new HashMap<>();
        Map<Double, List<Integer[]>> distMap = new HashMap<>();
        Map<Integer, Vector> vectorListMap = new HashMap<>();  Map<Integer, Integer[]> midPointEdgeMap = new HashMap<>(); Map<Integer, List<Integer>> getMidpointVector = new HashMap<>();Map<String, List<Integer>> faceEdgeMap = new HashMap<>();
        autoMOP.verticesComboFacesArchimedian("tetra", distEqual, faceList, verticeMap, distMap,vectorListMap, midPointEdgeMap);
        faceEdgeMap = autoMOP.getFaceEdgeMap();
        for (Map.Entry<String, List<Integer>> entry : faceEdgeMap.entrySet()) {
            System.out.println("Edge " + entry.getKey() + " found in faces: " + entry.getValue());
        }
        Map<Integer, List<Integer>> midpointVector = autoMOP.getMidpointVector();
        for (Map.Entry<Integer, List<Integer>> entry : midpointVector.entrySet()) {
            Integer key = entry.getKey();
            List<Integer> value = entry.getValue();
            System.out.println("Midpoint ID " + key + " -> Face Indices " + value);
        }
        System.out.println(vectorListMap);
        for (Map.Entry<Integer, Integer[]> entry: midPointEdgeMap.entrySet()){
            System.out.println("Edge " + entry.getKey() + " out " + Arrays.toString(entry.getValue()));
        }
    }

       /* public Box getMOPEdge(Box box, String strucName,  ISpecies species,List<Integer[]> linkedVertices, ArrayList<Integer> arrayListExtremities, ArrayList<Integer> arrayListLinkerOxy, ArrayList<Integer> arrayListLinkerCarb, int numEdges, double coordinateMultiplier, boolean ifPlatonicGeom){
        box.setNMolecules(species, numEdges);
        List<Vector> vecVertList = getListVecPoly();
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
        AutoMOP autoMOP = new AutoMOP();
        List<Integer[]> listEdges = new ArrayList<>(); List<List<Integer>> faceList = new ArrayList<>(); Map<Integer, Vector> verticeMap = new HashMap<>(); Map<Double, List<Integer[]>> distMap = new HashMap<>(); Map<Integer, Vector> vectorListMap = new HashMap<>();  Map<Integer, Integer[]> midPointEdgeMap = new HashMap<>();
        //gets coordinates depending on struct
        if (ifPlatonicGeom){
            verticesComboFaces(strucName, vertices, faces, faceEdges);
        }else {
            autoMOP.verticesComboFacesArchimedian(strucName, listEdges, faceList, verticeMap, distMap,vectorListMap, midPointEdgeMap );
        }
        List<Integer[]> faceOne = new ArrayList<>();
        if (ifPlatonicGeom){
            faceOne = faceEdges.get(0);
        }else {
            faceOne = listEdges;
            linkedVertices = listEdges;

        }

        //System.out.println(vecVertList);
        for (int i =0; i < verticeMap.size(); i++){
            vertices.add(verticeMap.get(i));
        }
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
        //System.exit(1);
        // Vector vecOne = vecListPoly.get(0);
        //Vector vecTwo = vecListPoly.get(1);
        //actMult = Math.sqrt(vecOne.Mv1Squared(vecTwo));
        coordinateMultiplier =  coordinateMultiplier/actMult;
        List<Vector> vecListPoly = new ArrayList<>();
        if (ifPlatonicGeom){

            for (int i =0; i < vecVertList.size(); i++){
                Vector vecOne = vecVertList.get(i);
                //Vector vecTwo = new Vector3D();
                //vecTwo.Ea1Tv1(coordinateMultiplier, vecOne);
                vecListPoly.add(vecOne);

            }
        }/*else {
            for (int i =0; i < verticeMap.size(); i++){
                vertices.add(verticeMap.get(i));
                Vector vecOne = verticeMap.get(i);
                Vector vecTwo = new Vector3D();
                vecTwo.Ea1Tv1(coordinateMultiplier, vecOne);
                vecListPoly.add(vecTwo);
            }
        }
        System.out.println(vecListPoly);
        for(int i=0; i<linkedVertices.size(); i++) {
            //ideal vertice
            Integer[] vecArray = linkedVertices.get(i); //vertices that form an edge
            IMolecule moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vecListPoly.get(vecArray[0] ); // first end of vertice

            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne); //convert actual to real vertice
            Vector bAct = new Vector3D();
            bAct.E(vecListPoly.get(vecArray[1] ));
            bAct.TE(coordinateMultiplier);
            //System.out.println(i + " 0 " +vecListPoly.get(4));
            //System.out.println("Expected basic "+ vecListPoly.get(vecArray[0]) +  " "+ vecListPoly.get(vecArray[1]));
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

            //System.out.println("final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition() +" final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition());

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
            //System.out.println(i + " 1 " +vecListPoly.get(4));
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
            //System.out.println("after moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition() +" final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition() +"\n");
            vecArray = linkedVertices.get(i);
            vecAId = vecListPoly.get(vecArray[0]);
            vecAIdVerify = new Vector3D();
            vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
            //System.out.println(i + " 2 " +vecListPoly.get(4));
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
            //System.out.println(i + " 3 " +vecListPoly.get(4));
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
            //System.out.println(i + " 4 " +vecListPoly.get(4));
            //System.out.println("done with Ligand \n");
        }
        //System.exit(1);
        return box;
    }*/
}
