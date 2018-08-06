package etomica.modules;

public final class EtomicaModuleInfo {
    public final String moduleName;
    public final String description;
    public final Class<?> moduleClass;

    public EtomicaModuleInfo(String moduleName, String description, Class<?> moduleClass) {
        this.moduleName = moduleName;
        this.description = description;
        this.moduleClass = moduleClass;
    }

    public static EtomicaModuleInfo mod(String moduleName, String description, Class<?> moduleClass) {
        return new EtomicaModuleInfo(moduleName, description, moduleClass);
    }

    public String toString() {
        return this.moduleName;
    }

    public static final EtomicaModuleInfo[] ETOMICA_MODULES = {
            mod("stats","No description yet",etomica.modules.statistics.StatisticsMCGraphic.class),
            mod("adsorption","No description yet",etomica.modules.adsorption.AdsorptionGraphic.class),
            mod("catalysis","No description yet",etomica.modules.catalysis.CatalysisGraphic.class),
            mod("chainequilibrium","No description yet",etomica.modules.chainequilibrium.ChainEquilibriumGraphic.class),
            mod("colloid","No description yet",etomica.modules.colloid.ColloidGraphic.class),
            mod("crystalviewer","Permits viewing of static 3-dimensional lattices and the planes they define. Various elementary lattices can be selected, and any plane or surface can be viewed by specifying it via its Miller indices. Image can be rotated to permit viewing from any angle.",etomica.modules.crystalviewer.CrystalViewer.class),
            mod("dcvgcmd","No description yet",etomica.modules.dcvgcmd.DCVGCMDGraphic.class),
            mod("droplet","No description yet",etomica.modules.droplet.DropletGraphic.class),
            mod("droplet_atomic","No description yet",etomica.modules.droplet.DropletAtomicGraphic.class),
            mod("ensembles","No description yet",etomica.modules.ensembles.LJMCGraphic.class),
            mod("insertion","No description yet",etomica.modules.insertion.InsertionGraphic.class),
            mod("interfacial","No description yet",etomica.modules.interfacial.InterfacialSWGraphic.class),
            mod("ljmd","No description yet",etomica.modules.ljmd.LjmdGraphic.class),
            mod("materialfracture","No description yet",etomica.modules.materialfracture.MaterialFractureGraphic.class),
            mod("mu","No description yet",etomica.modules.mu.MuGraphic.class),
            mod("multiharmonic","No description yet",etomica.modules.multiharmonic.MultiharmonicGraphicMC.class),
            mod("osmosis","No description yet",etomica.modules.osmosis.Osmosis.class),
            mod("pistoncylinder","The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a chamber with a movable wall under external pressure.",etomica.modules.pistoncylinder.PistonCylinderGraphic.class),
            mod("reactionequilibrium","Simple reaction equlibrium involving two atomic species and the three dimeric molecules they can form. Atoms move about via 2-D molecular dynamics, and can \"react\" to form dimers. Dynamic equilibria is demonstrated through constant recombining of atoms, and equilibria can be quantified and analyzed with thermodynamic reaction equilibria models.",etomica.modules.reactionequilibrium.ReactionEquilibriumGraphic.class),
            mod("rheology","No description yet",etomica.modules.rheology.RheologyGraphic.class),
            mod("reverse_osmosis","No description yet",etomica.modules.rosmosis.ReverseOsmosisGraphic.class),
            mod("reverse_osmosis_water","No description yet",etomica.modules.rosmosis.ReverseOsmosisWaterGraphic.class),
            mod("sam","No description yet",etomica.modules.sam.SamGraphic.class),
            mod("swmd","No description yet",etomica.modules.swmd.SwmdGraphic.class),
            mod("vle","No description yet",etomica.modules.vle.VLE.class),
            mod("b2fit","No description yet",etomica.modules.vle.B2Fit.class),
    };
}

