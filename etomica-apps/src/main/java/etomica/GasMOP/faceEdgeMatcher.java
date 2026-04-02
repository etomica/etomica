package etomica.GasMOP;

import java.util.*;


public class faceEdgeMatcher {
    public static void main(String[] args) {
        List<Integer[]> faceVertices = new ArrayList<>();
        faceVertices.add(new Integer[]{1, 5, 3, 2});
        faceVertices.add(new Integer[]{2, 3, 4});
        faceVertices.add(new Integer[]{1, 2, 6, 4});
        faceVertices.add(new Integer[]{1, 2, 4, 6});
        faceVertices.add(new Integer[]{1, 2, 6, 4, 9});

        List<Integer[]> listEdges = new ArrayList<>();
        listEdges.add(new Integer[]{2, 3});
        listEdges.add(new Integer[]{6, 2});
        faceEdgeMatcher faceEdgeMatcher = new faceEdgeMatcher();
      //  Map<String, List<Integer>> edgeToFaces = faceEdgeMatcher.mapEdgesToFaces(faceVertices, listEdges);

        // Print nicely
       /* for (Map.Entry<String, List<Integer>> entry : edgeToFaces.entrySet()) {
            System.out.println("Edge " + entry.getKey() + " found in faces: " + entry.getValue());
        }*/
    }

    public Map<String, List<Integer>> mapEdgesToFaces(List<List<Integer>> faceVertices, List<Integer[]> listEdges) {
        Map<String, List<Integer>> edgeToFaceIndices = new HashMap<>();

        for (Integer[] edge : listEdges) {
            String edgeKey = getEdgeKey(edge[0], edge[1]);
            List<Integer> matchingIndices = new ArrayList<>();

            for (int i = 0; i < faceVertices.size(); i++) {
                List<Integer> face = faceVertices.get(i);
                Set<String> faceEdges = getAllEdgesFromFace(face);

                if (faceEdges.contains(edgeKey)) {
                    matchingIndices.add(i);
                }
            }

            edgeToFaceIndices.put(edgeKey, matchingIndices);
        }

        return edgeToFaceIndices;
    }

    // Extracts all edges from a face polygon
    private Set<String> getAllEdgesFromFace(List<Integer> face) {
        Set<String> edgeSet = new HashSet<>();
        int n = face.size();
        for (int i = 0; i < n; i++) {
            int a = face.get(i);
            int b = face.get((i + 1) % n); // wrap around
            edgeSet.add(getEdgeKey(a, b));
        }
        return edgeSet;
    }

    // Normalizes edge so [a, b] == [b, a]
    private String getEdgeKey(int a, int b) {
        return (a < b) ? a + "-" + b : b + "-" + a;
    }


 /*   public Map<String, List<Integer>> mapEdgesToFaces(List<List<Integer>> faceVertices, List<Integer[]> listEdges) {
        Map<String, List<Integer>> edgeToFaceIndices = new HashMap<>();

        for (Integer[] edge : listEdges) {
            String edgeKey = getEdgeKey(edge[0], edge[1]);
            List<Integer> matchingIndices = new ArrayList<>();

            for (int i = 0; i < faceVertices.size(); i++) {
                List<Integer> face = faceVertices.get(i);
                Set<String> faceEdges = getAllEdgesFromFace(face);

                if (faceEdges.contains(edgeKey)) {
                    matchingIndices.add(i);
                }
            }

            edgeToFaceIndices.put(edgeKey, matchingIndices);
        }

        return edgeToFaceIndices;
    }

    // Extracts all edges from a face polygon
    private Set<String> getAllEdgesFromFace(List<Integer> face) {
        Set<String> edgeSet = new HashSet<>();
        int n = face.size();
        for (int i = 0; i < n; i++) {
            int a = face.get(0);
            int b = face.get((i + 1) % n); // wrap around
            edgeSet.add(getEdgeKey(a, b));
        }
        return edgeSet;
    }

    // Normalizes edge so [a, b] == [b, a]
    private String getEdgeKey(int a, int b) {
        return (a < b) ? a + "-" + b : b + "-" + a;
    }*/
}