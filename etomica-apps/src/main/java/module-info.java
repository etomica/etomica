module etomica.apps {
    exports etomica.overlap;
    exports etomica.virial.cluster;
    exports etomica.virial.overlap;
    exports etomica.AlkaneEH;
    exports etomica.association;
    exports etomica.dielectric;
    exports etomica.dimer;
    exports etomica.eam;
    exports etomica.eigenstuff;
    exports etomica.freeenergy.npath;
    exports etomica.gaussianwork;
    exports etomica.interfacial;
    exports etomica.kmc;
    exports etomica.liquidLJ;
    exports etomica.mappedRdf;
    exports etomica.mappedvirial;
    exports etomica.virial;
    exports etomica.virial.simulations;
    requires etomica.core;
    requires etomica.graph;
    requires java.desktop;
    requires jackson.databind;
    requires jama;
}
