package etomica.box.storage;

import etomica.space.Space;

public class Tokens {
    public static final Token<VectorStorage> POSITION = defaultVectors();
    public static final Token<VectorStorage> VELOCITY = defaultVectors();
    public static final Token<VectorStorage> ANGULAR_VELOCITY = defaultVectors();
    public static final Token<VectorStorage> ANGULAR_MOMENTUM = defaultVectors();

    public static final Token<OrientationStorage> ORIENTATION_FULL = defaultOrientations(false);

    public static final Token<OrientationStorage> ORIENTATION = defaultOrientations(true);

    public static final Token<DoubleStorage> BOND_LENGTH = defaultDoubles();

    public interface Initializer<T extends Storage> {
        void init(int idx, T storage);
    }

    public static Token<VectorStorage> defaultVectors() {
        return new Token<VectorStorage>() {
            @Override
            public VectorStorage createStorage(Space space) {
                return new VectorStorage(space);
            }
        };
    }


    public static Token<OrientationStorage> defaultOrientations(boolean isAxisSymmetric) {
        return new Token<OrientationStorage>() {
            @Override
            public OrientationStorage createStorage(Space space) {
                return new OrientationStorage(space, isAxisSymmetric);
            }
        };
    }

    public static Token<DoubleStorage> defaultDoubles() {
        return new Token<DoubleStorage>() {
            @Override
            public DoubleStorage createStorage(Space space) {
                return new DoubleStorage();
            }
        };
    }

    public static Token<VectorStorage> vectors(Initializer<VectorStorage> init) {
        return new Token<VectorStorage>() {
            @Override
            public VectorStorage createStorage(Space space) {
                return new VectorStorage(space);
            }

            @Override
            public void init(int idx, VectorStorage storage) {
                init.init(idx, storage);
            }
        };
    }
    
    public static Token<DoubleStorage> doubles(Initializer<DoubleStorage> init) {
        return new Token<DoubleStorage>() {
            @Override
            public DoubleStorage createStorage(Space space) {
                return new DoubleStorage();
            }

            @Override
            public void init(int idx, DoubleStorage storage) {
                init.init(idx, storage);
            }
        };
    }
    
    public static <T> Token<ObjectStorage<T>> objects(ObjectStorage.Factory<T> factory) {
        return new Token<ObjectStorage<T>>() {
            @Override
            public ObjectStorage<T> createStorage(Space space) {
                return new ObjectStorage<>(factory);
            }

            @Override
            public void init(int idx, ObjectStorage<T> storage) {
                storage.create(idx);
            }
        };
    }
}
