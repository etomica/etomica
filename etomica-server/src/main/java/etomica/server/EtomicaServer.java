package etomica.server;

import org.jooby.Jooby;

public class EtomicaServer extends Jooby {

    {
        get("/", () -> "Hello World");
    }

    public static void main(String[] args) {
        run(EtomicaServer::new, args);
    }
}
