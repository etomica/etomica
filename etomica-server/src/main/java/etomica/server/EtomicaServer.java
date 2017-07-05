package etomica.server;

import etomica.meta.ComponentIndex;
import etomica.simulation.Simulation;
import org.jooby.Jooby;

import java.util.Set;

public class EtomicaServer extends Jooby {

    {
        get("/", () -> "Hello World");
    }

    public static void main(String[] args) {
        ComponentIndex<Simulation> idx = new ComponentIndex<>(Simulation.class);
        Set<Class<? extends Simulation>> set = idx.getComponentSet();
        run(EtomicaServer::new, args);
    }
}
