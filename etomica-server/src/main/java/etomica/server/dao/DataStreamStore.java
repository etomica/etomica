package etomica.server.dao;

import etomica.data.DataPump;

import javax.inject.Singleton;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

@Singleton
public class DataStreamStore extends ConcurrentHashMap<UUID, DataPump> {
}
