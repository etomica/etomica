package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.potential.UFF.UFF;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AutoMOPCombines {
    int numVert = 0, numEdge =0;
    double distAD = 0.0, distAC = 0.0;
    String geomName="";
    private RotationTensor3D rotationTensor;
    ISpecies speciesOne, speciesTwo;
    Vector vecCentrePl;
    List<Integer> fluorideList;
    double tolerance = 1E-2;
    double distLen;
    public ISpecies makeSpeciesMOP(Box box,String mopGeom, String confNameOne, String confNameTwo, PDBReaderMOPCombine pdbReaderMOPOne, PDBReaderMOPCombine pdbReaderMOPTwo, ISpecies speciesOne, ISpecies speciesTwo, Map<String, AtomType> typeMapNew){
        ISpecies speciesMOP;
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        PDBReaderMOP pdbReaderMOP1 = new PDBReaderMOP();
        PDBReaderMOP pdbReaderMOP2 = new PDBReaderMOP();
        speciesOne = pdbReaderMOPOne.getSpecies(pdbReaderMOP1, confNameOne, true, new Vector3D(0,0,0), false, false);
        speciesTwo = pdbReaderMOPTwo.getSpecies(pdbReaderMOP2, confNameTwo, true, new Vector3D(0,0,0), false, false);
        int numVert = 0, numEdge = 0;
        String strucAct ="";
        List<Integer[]> vertCoordinates = autoMOPCombines.getVertices(mopGeom, numVert, numEdge, strucAct);
        List<Integer> flurideOne = pdbReaderMOPOne.getSpeciesConnectList();
        List<Integer> fluorideTwo = pdbReaderMOPTwo.getSpeciesConnectList();
        //get centre for Vertice Ligand
        Vector vecCentrePyPl = autoMOPCombines.getCentreVertice(flurideOne, speciesOne);
        double geomSide = autoMOPCombines.getGeomSide(speciesOne, speciesTwo, vecCentrePyPl, flurideOne, fluorideTwo);
        autoMOPCombines.makeSpeciesCombined(autoMOPCombines,box, speciesOne, speciesTwo, vertCoordinates, geomSide, flurideOne, fluorideTwo, vecCentrePyPl);
        IMolecule moleculeI = null;
        AtomType atomI = null;
        AtomType typeNew;
        int alpha = 0;
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        Map<Integer, String> atomMapMOPPDB = new HashMap<>();
        Map<Integer, Vector3D> atomMapVector = new HashMap<>();
        for (int i  = 0; i < box.getMoleculeList().size(); i++){
            moleculeI = box.getMoleculeList().get(i);
            for (int j =0; j < moleculeI.getChildList().size(); j++){
                atomI = moleculeI.getChildList().get(j).getType();
                String atomIString = String.valueOf(atomI);
                Vector atomPosn = moleculeI.getChildList().get(j).getPosition();
                // System.out.println(atomI);
                int startIndex = atomIString.indexOf("[") + 1;
                int endIndex = atomIString.indexOf("]");
                String nameNew = atomIString.substring(startIndex, endIndex);
                AtomType newName = pdbReaderMOP.returnElement(nameNew, false);
                // System.out.println(newName);
                // AtomType newNameNew = returnElement(nameNew);
                if (typeMapNew.containsKey(nameNew)) {
                    typeNew = typeMapNew.get(nameNew);
                } else {
                    typeNew = newName;
                    typeMapNew.put(nameNew, typeNew);
                }
                String atomPDB = nameNew.replaceAll("_\\d+$", "");
                atomMapMOPPDB.put(alpha, atomPDB);
                atomMapVector.put(alpha, (Vector3D) atomPosn);
                speciesBuilder.addAtom(typeNew, atomPosn,  "");
                alpha++;
            }
        }
        speciesMOP = speciesBuilder.setDynamic(true).build();
        return speciesMOP;
    }

    //designed for basic edgeTypeMOPs
    public List<Integer[]> getVertices(String geomName, int numVert, int numEdge, String strucAct){
        AutoMOP autoMOP = new AutoMOP();
        List<Integer[]> vertCoordinates = new ArrayList<>();
        if (geomName.equals("3py4_2L6")){
            vertCoordinates = autoMOP.getStrucList("tetra");
            numVert = 4;
            numEdge = 6;
            strucAct = "tetra";
        } else if (geomName.equals("4py6_2L12")) {
            vertCoordinates = autoMOP.getStrucList("octa");
            numVert = 6;
            numEdge = 12;
            strucAct = "octa";
        } else if (geomName.equals("5py12_2L30")){
            vertCoordinates = autoMOP.getStrucList("icosa");
            numVert = 12;
            numEdge = 30;
            strucAct = "icosa";
        }
        setEdge(numEdge);
        setVert(numVert);
        setGeom(strucAct);
        return vertCoordinates;
    }
    public void setVert(int vert){this.numVert = vert;}
    public void setEdge(int edge){this.numEdge = edge;}
    public void setGeom(String geom){this.geomName = geom;}
    public int getVert(){return numVert;}
    public int getEdge(){return numEdge;}
    public String getGeom(){return geomName;}

    public Map<Integer, Integer[]> pointClosestVertMap(String geomName){
        Map<Integer, Integer[]> pointsMap = new HashMap<>();
        pointsMap.put(0, new Integer[]{1, 2, 3});
        pointsMap.put(1, new Integer[]{0, 2, 3});
        pointsMap.put(2, new Integer[]{0, 1, 3});
        pointsMap.put(3, new Integer[]{0, 1, 2});
        return pointsMap;
    }

    public Box makeSpeciesCombined(AutoMOPCombines autoMOPCombines,Box box, ISpecies speciesOne, ISpecies speciesTwo, List<Integer[]> vectorVertices, double coordinateMultiplier, List<Integer> fluorideOne, List<Integer> fluoridetwo, Vector vecCentrePl){
        ISpecies speciesMOP = null;
        List<Integer> arrayListLinkerOxy = new ArrayList<>();
        List<Integer> arrayListLinkerCarb = new ArrayList<>();
        int numEdge = autoMOPCombines.getEdge();
        int numVert = autoMOPCombines.getVert();
        String strucAct = autoMOPCombines.getGeom();
        box.setNMolecules(speciesOne, numEdge);
        box.setNMolecules(speciesTwo, numVert);
        Space space = Space3D.getInstance();
        IMolecule moleculeEdge = box.getMoleculeList().get(0);
        IMolecule moleculeVert = box.getMoleculeList().get(numEdge);
        List<Vector> oldPositionEdge = new ArrayList<>();
        List<Vector> oldPositionVert = new ArrayList<>();
        while (oldPositionVert.size() < moleculeVert.getChildList().size()) {
            oldPositionVert.add(space.makeVector());
        }
        while (oldPositionEdge.size() < moleculeEdge.getChildList().size()) {
            oldPositionEdge.add(space.makeVector());
        }
        List<Vector> vertices = new ArrayList<>();
        List<Integer[]> faces = new ArrayList<>();
        Map<Integer, List<Integer[]>> faceEdges = new HashMap<>();
        AutoMOP autoMOP = new AutoMOP();
        autoMOP.verticesComboFaces(strucAct, vertices, faces, faceEdges );
        List<Integer[]> faceOne = faceEdges.get(0);
        double actMult = 1000000;
        for (int i =0; i< faceOne.size(); i++){
            Vector vertOne = vertices.get(faceOne.get(i)[0]);
            Vector vertTwo = vertices.get(faceOne.get(i)[1]);
            double distVal = Math.sqrt(vertOne.Mv1Squared(vertTwo));
            if (distVal < actMult){
                actMult = distVal;
            }
        }
        coordinateMultiplier =  coordinateMultiplier/actMult;
        double distAC = autoMOPCombines.getDistAC();
        double disAD = autoMOPCombines.getDistAD();
        List<List<Vector>> edgeVectors = autoMOPCombines.getInternalCoordinatesEdges(vectorVertices, vertices, distAC, disAD);
        Map<Integer, Integer[]> verticeConnectionMap = autoMOPCombines.pointClosestVertMap(geomName);
        // rotate vertices
        for (int i = 0; i<vertices.size(); i++){
            IMolecule moleculeVertN = box.getMoleculeList(speciesOne).get(i);
            Vector vertice =vertices.get(i);
            Vector vecCopy = new Vector3D();
            vecCopy.E(vertice);
            vecCopy.ME(vecCentrePl);

            //move molecule centre to vertice
            final Vector vecCopyTranslate = vecCopy;
            moleculeVertN.getChildList().forEach(atom -> {
                oldPositionVert.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecCopyTranslate);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            Vector vecA = moleculeVert.getChildList().get(fluorideOne.get(0)).getPosition();
            Vector vecACopy = new Vector3D();
            vecACopy.E(vecA);
            vecCopy = new Vector3D();
            vecCopy.E(vecCentrePl);

            vecACopy.ME(vecCopy);
            vecACopy.normalize();
            Integer[] vertConnected = verticeConnectionMap.get(i);
            Vector vertAdj = vertices.get(vertConnected[0]);
            vertAdj.ME(vecCopy);
            vertAdj.normalize();

            double prodDot = vecACopy.dot(vertAdj);
            double theta = Math.acos(prodDot);
            Vector vecBCross = new Vector3D();
            vecBCross.E(vecACopy);
            vecBCross.XE(vertAdj);
            vecBCross.normalize();

            final Vector vecCopyFinal = vecCopy;
            rotationTensor.setRotationAxis(vecBCross, theta);
            moleculeVertN.getChildList().forEach(atom -> {
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
                atom.getPosition().ME(vecCopyFinal);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecCopyFinal);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });



            vecA = moleculeVert.getChildList().get(fluorideOne.get(1)).getPosition();
            vecACopy = new Vector3D();
            vecACopy.E(vecA);
            vecCopy = new Vector3D();
            vecCopy.E(vecCentrePl);

            vecACopy.ME(vecCopy);
            vecACopy.normalize();
            vertAdj = vertices.get(vertConnected[1]);
            vertAdj.ME(vecCopy);
            vertAdj.normalize();

            prodDot = vecACopy.dot(vertAdj);
            theta = Math.acos(prodDot);
            vecBCross = new Vector3D();
            vecBCross.E(vecACopy);
            vecBCross.XE(vertAdj);
            vecBCross.normalize();

            final Vector vecCopyFinal2 = vecCopy;
            rotationTensor.setRotationAxis(vecBCross, theta);
            moleculeVertN.getChildList().forEach(atom -> {
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
                atom.getPosition().ME(vecCopyFinal2);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(vecCopyFinal2);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }

        //rotate edges
        List<Integer[]> linkedVertices = vectorVertices;
        List<Vector> vecListPoly = vertices;


        for(int i=0; i<numEdge; i++) {
            //ideal vertice
            Integer[] vecArray = linkedVertices.get(i); //vertices that form an edge
            IMolecule moleculeMOPOne = box.getMoleculeList().get(i);
            Vector originOne, copyOriginOne = new Vector3D();
            originOne = vecListPoly.get(vecArray[0] ); // first end of vertice
            copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne); //convert actual to real vertice
            Vector bAct = new Vector3D();
            bAct.E(vecListPoly.get(vecArray[1] ));
            bAct.TE(coordinateMultiplier);
            //move molecule first extremity to the real vertice
            IMolecule finalMoleculeMOPOne = moleculeMOPOne;
            Vector originStart = finalMoleculeMOPOne.getChildList().get(fluorideOne.get(0)).getPosition();
            copyOriginOne.ME(originStart);
            finalMoleculeMOPOne.getChildList().forEach(atom -> {
                oldPositionEdge.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(copyOriginOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            //find angle between vector connecting extremities and vertices that need to be connected
            Vector vecAId = vecListPoly.get(vecArray[0]);
            //    System.out.println("vecAID " +copyOriginOne);
            Vector vecAIdVerify = new Vector3D();
            vecAIdVerify.E(vecAId);
            // vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
            Vector vecBId = vecListPoly.get(vecArray[1]);
            Vector vecBIdCopy = new Vector3D();
            vecBIdCopy.Ev1Mv2(vecBId, vecAId);
            //System.out.println(vecAId +" "+vecBId);
            vecBIdCopy.normalize();
            //System.out.println("vecAID norm " + vecBIdCopy);
            Vector vecANew = moleculeEdge.getChildList().get(fluoridetwo.get(0)).getPosition();
            Vector vecBNew = moleculeEdge.getChildList().get(fluoridetwo.get(1)).getPosition();
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
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
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
              /*  if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
                atom.getPosition().ME(finalVecAIdVerify);
                rotationTensor.transform(atom.getPosition());
                atom.getPosition().PE(finalVecAIdVerify);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }

        // make species MOP
        box.addMolecule(speciesMOP.makeMolecule());
        return box;
    }

    // get Centre of Pyramidal or planar ligand which is used to keep at vertice of the polyhedra
    public Vector getCentreVertice(List<Integer> fluorideConnect, ISpecies speciesOne){
        IMolecule molecule = speciesOne.makeMolecule();
        Vector vecSum = new Vector3D();
        Vector vecOne = new Vector3D();
        for (int i =0; i < fluorideConnect.size(); i++){
            vecOne = molecule.getChildList().get(fluorideConnect.get(i)).getPosition();
            vecSum.PE(vecOne);
        }
        double teVal =  1.0 / fluorideConnect.size();
        vecSum.TE(teVal);
        if (fluorideConnect.size() == 3){
            Vector vecA = molecule.getChildList().get(fluorideConnect.get(0)).getPosition();
            Vector vecB = molecule.getChildList().get(fluorideConnect.get(0)).getPosition();
            double length = Math.sqrt(vecB.Mv1Squared(vecA));
            setLenDist(length);
        }
        return vecSum;
    }
    public void setLenDist(double lenDist){
        this.distLen = lenDist;
    }
    public double getDistLen(){return distLen;}
    public void makeMoleculeSideLength(String mopGeom, String confNameOne, String confNameTwo, ISpecies speciesOne, ISpecies speciesTwo){
        Space space = Space3D.getInstance();
      /*  IMolecule moleculeVert = box.getMoleculeList().get(numEdge);
        List<Vector> oldPositionVert = new ArrayList<>();
        while (oldPositionVert.size() < moleculeVert.getChildList().size()) {
            oldPositionVert.add(space.makeVector());
        }*/
        PDBReaderMOPCombine pdbReaderMOP1 = new PDBReaderMOPCombine();
        PDBReaderMOPCombine pdbReaderMOP2 = new PDBReaderMOPCombine();
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        PDBReaderMOP pdbReaderMOPO1 = new PDBReaderMOP();
        PDBReaderMOP pdbReaderMOPO2 = new PDBReaderMOP();
        speciesOne = pdbReaderMOP1.getSpecies(pdbReaderMOPO1, confNameOne, false, new Vector3D(0.0, 0.0, 0.0), false, true);
        speciesTwo = pdbReaderMOP2.getSpecies(pdbReaderMOPO2, confNameTwo, false, new Vector3D(0.0, 0.0, 0.0), false, true);
        System.out.println(speciesOne);
        System.out.println(speciesTwo);
        ArrayList<ArrayList<Integer>> connectivityListOne = pdbReaderMOP1.getConnectivity();
        Map<Integer, String> atomMapOne = pdbReaderMOP1.getAtomMap();
        Map<Integer, String> atomMapModifiedOne = pdbReaderMOPO1.getAtomMapModified();
        //AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        double sideLength = autoMOPCombines.getDistLen();
        double actMult = 1000000;
        AutoMOP autoMOP = new AutoMOP();


       /* for (int i =0; i< faceOne.size(); i++){
            Vector vertOne = vertices.get(faceOne.get(i)[0]);
            Vector vertTwo = vertices.get(faceOne.get(i)[1]);
            double distVal = Math.sqrt(vertOne.Mv1Squared(vertTwo));
            if (distVal < actMult){
                actMult = distVal;
            }
        }*/
        List<Integer> fluorideList = pdbReaderMOP1.settingUpFluorideList(connectivityListOne, atomMapModifiedOne, speciesOne);
        Vector vecCentrePl = getCentreVertice(fluorideList, speciesOne);
        setSpeciesOne(speciesOne);
        setSpeciesTwo(speciesTwo);
        setVectorPl(vecCentrePl);
        setFluorideList(fluorideList);
    }

    private void setFluorideList(List<Integer> fluorideList) {this.fluorideList = fluorideList;}
    private List<Integer> getFluorideList() {return fluorideList;}

    private void setVectorPl(Vector vecCentrePl) {this.vecCentrePl = vecCentrePl;}

    private Vector getVectorPl() {return  vecCentrePl;}

    private void setSpeciesTwo(ISpecies speciesTwo) {this.speciesTwo = speciesTwo;}
    public ISpecies getSpeciesTwo() {return  speciesTwo;}
    private void setSpeciesOne(ISpecies speciesOne) {this.speciesOne = speciesOne;}
    public ISpecies getSpeciesOne() {return  speciesOne;}

    public void getSideLength(Box box, ISpecies speciesOne, AutoMOPCombines autoMOPCombines ){
        rotationTensor = new RotationTensor3D();
        AutoMOP autoMOP = new AutoMOP();
        Space space = Space3D.getInstance();
        List<Integer[]> listEdges = new ArrayList<>(); List<List<Integer>> faceList = new ArrayList<>(); Map<Integer, Vector> verticeMap = new HashMap<>(); Map<Double, List<Integer[]>> distMap = new HashMap<>(); Map<Integer, Vector> vectorListMap = new HashMap<>();  Map<Integer, Integer[]> midPointEdgeMap = new HashMap<>();
        //gets coordinates depending on struct
        autoMOP.verticesComboFacesArchimedian("tetra", listEdges, faceList, verticeMap, distMap,vectorListMap, midPointEdgeMap );
        IMolecule moleculeVertN = box.getMoleculeList(speciesOne).get(0);
        IMolecule moleculeVert = moleculeVertN;
        Vector verticeOne = verticeMap.get(0);
        Vector vertOneCopy = new Vector3D();
        vertOneCopy.E(verticeOne);
        vertOneCopy.ME(vecCentrePl);
        List<Vector> oldPositionVert = new ArrayList<>();
        while (oldPositionVert.size() < moleculeVert.getChildList().size()) {
            oldPositionVert.add(space.makeVector());
        }
        //move molecule centre to vertice
        final Vector vecCopyTranslate = vertOneCopy;
        moleculeVertN.getChildList().forEach(atom -> {
            oldPositionVert.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(vecCopyTranslate);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });


        //rotate fluoride to align it with edges of tetrahedron
        //Angle between (FluorideOne and vertTwo)  and (vertOne and vertTwo)
        Vector vecFluorOne = moleculeVert.getChildList().get(fluorideList.get(0)).getPosition();
        Vector vecFluorOneCopy = new Vector3D();
        vecFluorOneCopy.E(vecFluorOne);
        Vector verticeTwo = verticeMap.get(1);
        Vector vertTwoCopy = new Vector3D();
        vertTwoCopy.E(verticeTwo);

        Vector vertTwoFluor = new Vector3D();
        Vector vertOneTwo = new Vector3D();
        Vector vertTwoFluorCopy = new Vector3D();
        Vector vertOneTwoCopy = new Vector3D();

        vertTwoFluor.E(verticeTwo);
        vertTwoFluor.ME(vecFluorOne);
        vertTwoFluorCopy.E(vertTwoFluor);
        vertTwoFluorCopy.normalize();

        vertOneTwo.E(verticeOne);
        vertOneTwo.ME(verticeTwo);
        vertOneTwoCopy.E(vertOneTwo);
        vertOneTwoCopy.normalize();

        double prodDot = vertTwoFluorCopy.dot(vertOneTwoCopy);
        double theta = Math.acos(prodDot);
        vertTwoFluorCopy.XE(vertOneTwoCopy);


        final Vector vecCopyFinal = vertTwoCopy;
        rotationTensor.setRotationAxis(vertTwoFluorCopy, theta);
        moleculeVertN.getChildList().forEach(atom -> {
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
            atom.getPosition().ME(vecCopyFinal);
            rotationTensor.transform(atom.getPosition());
            atom.getPosition().PE(vecCopyFinal);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        Vector fluorOne = moleculeVertN.getChildList().get(fluorideList.get(0)).getPosition();
        Vector vecFluorVertOne = new Vector3D();
        vecFluorVertOne.E(fluorOne);
        vecFluorVertOne.ME(verticeOne);
        Vector fluorTwo = moleculeVertN.getChildList().get(fluorideList.get(1)).getPosition();
        Vector vecFluorVertTwo = new Vector3D();
        vecFluorVertTwo.E(fluorTwo);
        vecFluorVertTwo.ME(verticeTwo);
        autoMOPCombines.ifIntersects(vecFluorVertOne, verticeOne, vecFluorVertTwo, verticeTwo, 0.01);
    }
    public void ifIntersects (Vector vecOne , Vector pointA,  Vector vecTwo, Vector pointB, double tolerance){
        Line lineA = new Line(pointA, vecOne);
        Line lineB = new Line(pointB, vecTwo);
        LineIntersection lineIntersection = new LineIntersection();
        Vector intersection = lineIntersection.findIntersection(lineA, lineB, tolerance);
        if (intersection != null) {
            System.out.println("The intersection point is: " + intersection);
        }
    }

    // presently just for Pyramidal and Linear Combination
    public double getGeomSide(ISpecies speciesOne, ISpecies speciesTwo, Vector vecCentre, List<Integer> lenFluorideOne, List<Integer> lenFluorideTwo){
        double length = 0.0;
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        IMolecule moleculeOne = speciesOne.makeMolecule();
        IMolecule moleculeTwo = speciesTwo.makeMolecule();
        String symbol = moleculeOne.getChildList().get(lenFluorideOne.get(0)).getType().toString();
        String atomOne = autoMOPCombines.getAtomNameFromAtomType(symbol);
        symbol = moleculeTwo.getChildList().get(lenFluorideTwo.get(0)).getType().toString();
        String atomTwo = autoMOPCombines.getAtomNameFromAtomType(symbol);
        double[] potentialOne = pdbReaderMOP.atomicPot(atomOne);
        double[] potentialTwo = pdbReaderMOP.atomicPot(atomTwo);
        if (lenFluorideTwo.size() != 2){
            throw new RuntimeException("error geometry Combinations");
        }
        double dist = 0.0;
        Vector vecMolTwoAtomOne = moleculeTwo.getChildList().get(lenFluorideTwo.get(0)).getPosition();
        Vector vecMolTwoAtomTwo = moleculeTwo.getChildList().get(lenFluorideTwo.get(0)).getPosition();
        dist = Math.sqrt(vecMolTwoAtomOne.Mv1Squared(vecMolTwoAtomTwo));
        double[] bondParamsArray= UFF.bondUFF (potentialOne[0],  potentialTwo[0], potentialOne[5],  potentialTwo[5], potentialOne[6], potentialTwo[6], 1.0);
        double bondLength = bondParamsArray[1];
        // distance between centre of vertice molecule to the attachment part.
        double projectionLength = autoMOPCombines.getProjectionLength(moleculeTwo);
        dist = dist + 2 * bondLength + 2 * projectionLength;
        setDistAC(bondLength + projectionLength);
        setDistAD(dist + bondLength + projectionLength);
        return dist;
    }
    public void setDistAC(double distAC){this.distAC = distAC;}
    public void setDistAD(double distAD){this.distAD = distAD;}
    public double getDistAC(){return distAC;}
    public double getDistAD(){return distAD;}

    //projection of moleculeTwo from its centre to metal connector
    public double getProjectionLength(IMolecule molecule){
        double projectionLength = 0.0;
        return projectionLength;
    }


    public List<List<Vector>> getInternalCoordinatesEdges(List<Integer[]>vertConnected, List<Vector> vertices, double distAC, double distAD){
        List<List<Vector>> vectEdges = new ArrayList<>();
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        List<Vector> vecList = new ArrayList<>();
        for (int i = 0; i < vertConnected.size(); i++ ){
            Integer[] edgeI = vertConnected.get(i);
            Vector vertZero = vertices.get(edgeI[0]);
            Vector vertOne = vertices.get(edgeI[1]);
            Vector vecC = autoMOPCombines.edgeMath(vertZero, vertOne, distAC);
            Vector vecD = autoMOPCombines.edgeMath(vertZero, vertOne, distAD);
            vecList = new ArrayList<>();
            vecList.add(vecC);
            vecList.add(vecD);
            vectEdges.add(vecList);
        }
        return vectEdges;
    }
    public static SpeciesManager boxSpecies (SpeciesManager speciesManager, Box box, ISpecies speciesVertice, ISpecies speciesEdge){
        Vector vectorOne, vectorTwo = new Vector3D();
        vectorOne = new Vector3D(0,0,0);
        vectorTwo = new Vector3D(2,1,1);

        List<Vector> oldPositionVert = new ArrayList<>();
        Space space = Space3D.getInstance();
        IMolecule moleculeVert = speciesVertice.makeMolecule();
        while (oldPositionVert.size() < moleculeVert.getChildList().size()) {
            oldPositionVert.add(space.makeVector());
        }
        //bond length of atoms
        double bond = 1.23;
        //place vertice at a location

        //rotate edges of vertice

        //place and translate edge length

        //rotate edges along structure edge

        //translate second vertex length

        //second edge translate

        //calculate length of edge projection and length of edge

        //modify the  tetrahedron size and its coordinates and get where fluoride atoms lie

        // determine the fluoride atoms of the edge and vertice where it starts

        // translate vertice edges

        // vertice ligands rotate it according to new polyhedra

        // edge translate to bondlengths and rotate

      //  rotateEdge(box, vectorOne, vectorTwo, bond);
     //   rotateVertice(box, vectorOne);
        return speciesManager;
    }


    public void rotateVertice(Box box, List<Vector> oldPositionVert,  List<Vector> vertices,   Map<Integer, Integer[]> verticeConnectionMap, Vector vecCentrePl, List<Integer> fluorideOne, Vector vectorAtomOnePosn, Vector vectorAtomTwoPosn, Vector vectorAtomThreePosn, int moleculeNum, ISpecies speciesOne) {
        RotationTensor3D rotationTensor = new RotationTensor3D();
        IMolecule moleculeVertN = box.getMoleculeList(speciesOne).get(moleculeNum);
        Vector vertice =vectorAtomOnePosn;
        Vector vecCopy = new Vector3D();
        vecCopy.E(vertice);
        vecCopy.ME(vecCentrePl);

        //move molecule centre to vertice
        final Vector vecCopyTranslate = vecCopy;
        moleculeVertN.getChildList().forEach(atom -> {
            oldPositionVert.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(vecCopyTranslate);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });

        Vector vecA = moleculeVertN.getChildList().get(fluorideOne.get(0)).getPosition();
        Vector vecACopy = new Vector3D();
        vecACopy.E(vecA);
        vecCopy = new Vector3D();
        vecCopy.E(vecCentrePl);

        vecACopy.ME(vecCopy);
        vecACopy.normalize();
        Integer[] vertConnected = verticeConnectionMap.get(moleculeNum);
        Vector vertAdj = vertices.get(vertConnected[0]);
        vertAdj.ME(vecCopy);
        vertAdj.normalize();

        double prodDot = vecACopy.dot(vertAdj);
        double theta = Math.acos(prodDot);
        Vector vecBCross = new Vector3D();
        vecBCross.E(vecACopy);
        vecBCross.XE(vertAdj);
        vecBCross.normalize();

        final Vector vecCopyFinal = vecCopy;
        rotationTensor.setRotationAxis(vecBCross, theta);
        moleculeVertN.getChildList().forEach(atom -> {
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
            atom.getPosition().ME(vecCopyFinal);
            rotationTensor.transform(atom.getPosition());
            atom.getPosition().PE(vecCopyFinal);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });

        vecA = moleculeVertN.getChildList().get(fluorideOne.get(1)).getPosition();
        vecACopy = new Vector3D();
        vecACopy.E(vecA);
        vecCopy = new Vector3D();
        vecCopy.E(vecCentrePl);

        vecACopy.ME(vecCopy);
        vecACopy.normalize();
        vertAdj = vertices.get(vertConnected[1]);
        vertAdj.ME(vecCopy);
        vertAdj.normalize();

        prodDot = vecACopy.dot(vertAdj);
        theta = Math.acos(prodDot);
        vecBCross = new Vector3D();
        vecBCross.E(vecACopy);
        vecBCross.XE(vertAdj);
        vecBCross.normalize();

        final Vector vecCopyFinal2 = vecCopy;
        rotationTensor.setRotationAxis(vecBCross, theta);
        moleculeVertN.getChildList().forEach(atom -> {
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
            atom.getPosition().ME(vecCopyFinal2);
            rotationTensor.transform(atom.getPosition());
            atom.getPosition().PE(vecCopyFinal2);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
    }

    private static void rotateEdge(Box box, List<Integer> arrayListLinkerOxy, List<Integer> arrayListLinkerCarb, List<Vector> oldPositionEdge, List<Integer> fluorideOne, List<Integer> fluoridetwo, double coordinateMultiplier,  List<Vector> vecListPoly, List<Integer[]> linkedVertices, Integer moleculeNum, Vector vectorOne, Vector vectorTwo, double bond) {
        RotationTensor3D rotationTensor = new RotationTensor3D();
        //ideal vertice
        Integer[] vecArray = linkedVertices.get(moleculeNum); //vertices that form an edge
        IMolecule moleculeMOPOne = box.getMoleculeList().get(moleculeNum);
        Vector originOne, copyOriginOne = new Vector3D();
        originOne = vecListPoly.get(vecArray[0] ); // first end of vertice
        copyOriginOne.Ea1Tv1(coordinateMultiplier, originOne); //convert actual to real vertice
        Vector bAct = new Vector3D();
        bAct.E(vecListPoly.get(vecArray[1] ));
        bAct.TE(coordinateMultiplier);
        //move molecule first extremity to the real vertice
        IMolecule finalMoleculeMOPOne = moleculeMOPOne;
        IMolecule moleculeEdge = moleculeMOPOne;
        Vector originStart = finalMoleculeMOPOne.getChildList().get(fluorideOne.get(0)).getPosition();
        copyOriginOne.ME(originStart);
        finalMoleculeMOPOne.getChildList().forEach(atom -> {
            oldPositionEdge.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(copyOriginOne);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });

        //find angle between vector connecting extremities and vertices that need to be connected
        Vector vecAId = vecListPoly.get(vecArray[0]);
        //    System.out.println("vecAID " +copyOriginOne);
        Vector vecAIdVerify = new Vector3D();
        vecAIdVerify.E(vecAId);
        // vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );
        Vector vecBId = vecListPoly.get(vecArray[1]);
        Vector vecBIdCopy = new Vector3D();
        vecBIdCopy.Ev1Mv2(vecBId, vecAId);
        //System.out.println(vecAId +" "+vecBId);
        vecBIdCopy.normalize();
        //System.out.println("vecAID norm " + vecBIdCopy);
        Vector vecANew = moleculeEdge.getChildList().get(fluoridetwo.get(0)).getPosition();
        Vector vecBNew = moleculeEdge.getChildList().get(fluoridetwo.get(1)).getPosition();
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
                /*if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
            atom.getPosition().ME(vecANew);
            rotationTensor.transform(atom.getPosition());
            atom.getPosition().PE(vecANew);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        // System.out.println("after moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(0)).getPosition() +" final moved "+finalMoleculeMOPOne.getChildList().get(arrayListExtremities.get(1)).getPosition() +"\n");
        vecArray = linkedVertices.get(moleculeNum);
        vecAId = vecListPoly.get(vecArray[0]);
        vecAIdVerify = new Vector3D();
        vecAIdVerify.Ea1Tv1(coordinateMultiplier,vecAId );

        //edges
        Integer[] vecBArray = linkedVertices.get(moleculeNum);
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
              /*  if(atom.getIndex() == arrayListExtremities.get(0)) {
                    return;
                }*/
            atom.getPosition().ME(finalVecAIdVerify);
            rotationTensor.transform(atom.getPosition());
            atom.getPosition().PE(finalVecAIdVerify);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
    }

    public Vector edgeMath(Vector vecA, Vector vecB, double distPoint){
        double vectorAB_x = vecB.getX(0) - vecA.getX(0);
        double vectorAB_y = vecB.getX(1) - vecA.getX(1);
        double vectorAB_z = vecB.getX(2) - vecA.getX(2);

        double magnitudeAB = Math.sqrt(Math.pow(vectorAB_x, 2) + Math.pow(vectorAB_y, 2) + Math.pow(vectorAB_z, 2));

        double unitVector_x = vectorAB_x / magnitudeAB;
        double unitVector_y = vectorAB_y / magnitudeAB;
        double unitVector_z = vectorAB_z / magnitudeAB;

        double D_x = vecA.getX(0) + distPoint * unitVector_x;
        double D_y = vecA.getX(1) + distPoint * unitVector_y;
        double D_z = vecA.getX(2) + distPoint * unitVector_z;
        return new Vector3D(D_x, D_y, D_z);
    }


    public String getAtomNameFromAtomType (String symbol){
        int startIndex = symbol.indexOf("[") + 1;
        int endIndex = symbol.indexOf("]");
        return symbol.substring(startIndex, endIndex);
    }
}