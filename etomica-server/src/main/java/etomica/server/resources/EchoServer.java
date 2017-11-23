package etomica.server.resources;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import etomica.server.dao.SimulationStore;

import javax.inject.Inject;
import javax.websocket.OnMessage;
import javax.websocket.OnOpen;
import javax.websocket.Session;
import javax.websocket.server.ServerEndpoint;

@ServerEndpoint("/ws/{test}")
@Metered
@Timed
public class EchoServer {
    private final SimulationStore store;

    @Inject
    public EchoServer(SimulationStore test) {
        store = test;
    }

    @OnOpen
    public void onOpen(final Session session) {
        String test = session.getPathParameters().get("test");
        session.getAsyncRemote().sendText(test);
    }

    @OnMessage
    public void onMessage(final Session session, String message) {
        session.getAsyncRemote().sendText(message.toUpperCase());
    }
}
