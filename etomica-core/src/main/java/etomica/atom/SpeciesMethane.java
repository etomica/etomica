package etomica.atom;


import etomica.molecule.IMolecule;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpBuilder;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class SpeciesMethane {
    private static SpeciesBuilder speciesBuilder; //Etomica defined class
    //static String confName = "F:/ethane";
    public static ISpecies buildMethane(boolean isDynamic, String confName){
       // String confName = getConfName();
        ISpecies species = SpBuilder.buildSpecies(confName);
        IMolecule molecule = species.makeMolecule();
        speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        for(IAtom atom: molecule.getChildList()){
            System.out.println(atom.getType()+" "+ atom.getPosition());
        }

        for (IAtom atom: molecule.getChildList()){
            speciesBuilder.addAtom(atom.getType(), atom.getPosition());
            //speciesBuilder.withConformation(new ConfirmationMethane(Space3D.getInstance()));

        }
        speciesBuilder.setDynamic(isDynamic);
        return speciesBuilder.build();

    }

  /*  public static ArrayList<ArrayList<Integer>> getConnect(){
        ArrayList<ArrayList<Integer>> connect = SpBuilder.getConnectivity(confName);
        System.out.println(connect + "Connect");
        return connect;
    }

    public static ArrayList<ArrayList<Integer>> getConnectNew(){
        ArrayList<ArrayList<Integer>> connect = getConnect();
        ArrayList<ArrayList<Integer>> connectNew = SpBuilder.getconnectivityModified(connect);
        System.out.println(connectNew + "ConnectNew");
        return connectNew;
    }

    public static Map<Integer, String> getAtomLayout(){
        ArrayList<ArrayList<Integer>> connect = getConnect();
        Map<Integer, String> atomStyle = SpBuilder.getAtomMap(confName, connect);
        System.out.println(atomStyle + "atomStyle");
        return atomStyle;
    }

    public static HashMap<Integer, String> getAtomLayoutModified(){
        ArrayList<ArrayList<Integer>> connect = getConnect();
        Map<Integer, String> atomStyle = SpBuilder.getAtomMap(confName, connect);
        HashMap<Integer, String>  AtomLayoutModified = SpBuilder.getatomMapModified(atomStyle);
        System.out.println(AtomLayoutModified + "atomMapModified");
        return AtomLayoutModified;
    }

    public static void main(String[] args) {
        getConnect();
        getConnectNew();
        getAtomLayout();
        getAtomLayoutModified();

    }*/

}
