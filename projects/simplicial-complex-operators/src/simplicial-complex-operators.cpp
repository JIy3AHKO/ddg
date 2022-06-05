// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

typedef Eigen::Triplet<double> Trip;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    auto x = SparseMatrix<size_t>(mesh->nVertices(), mesh->nEdges());
    
    std::vector<Trip> triples;

    for (auto v: this->mesh->vertices()) {
        for (auto e: v.adjacentEdges()) {
            triples.push_back(Trip(v.getIndex(), e.getIndex(), 1));
        }
    }

    x.setFromTriplets(triples.begin(), triples.end());

    return x.transpose();
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    auto x = SparseMatrix<size_t>(mesh->nFaces(), mesh->nEdges());
    
    std::vector<Trip> triples;

    for (auto f: this->mesh->faces()) {
        for (auto e: f.adjacentEdges()) {
            triples.push_back(Trip(f.getIndex(), e.getIndex(), 1));
        }
    }

    x.setFromTriplets(triples.begin(), triples.end());

    return x;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    auto x = Vector<size_t>(mesh->nVertices());

    x.fill(0);


    for (auto v: subset.vertices) {
        x[v] = 1;
    }

    return x;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    auto x = Vector<size_t>(mesh->nEdges());

    x.fill(0);


    for (auto e: subset.edges) {
        x[e] = 1;
    }

    return x;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    auto x = Vector<size_t>(mesh->nFaces());

    x.fill(0);


    for (auto f: subset.faces) {
        x[f] = 1;
    }

    return x;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    auto edges = this->buildEdgeVector(subset);
    auto vertices = this->buildVertexVector(subset);


    Vector<size_t> adj_e_f = A1 * edges;

    Vector<size_t> adj_v_e = A0 * vertices;
    Vector<size_t> adj_v_f = A1 * adj_v_e;


    MeshSubset starMesh = subset.deepCopy();

    for (Vector<size_t>::InnerIterator it(adj_e_f, 0); it; ++it)
    {
        if (it.value() != 0) {
            starMesh.addFace(it.row());
        }    }
    for (Vector<size_t>::InnerIterator it(adj_v_f, 0); it; ++it)
    {
        if (it.value() != 0) {
            starMesh.addFace(it.row());
        }   
    }

    for (Vector<size_t>::InnerIterator it(adj_v_e, 0); it; ++it)
    {
        if (it.value() != 0) {
            starMesh.addEdge(it.row());
        }   
    }


    return starMesh;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    auto clMesh = subset.deepCopy();

    auto faces = this->buildFaceVector(clMesh);


    Vector<size_t> adj_f_e = A1.transpose() * faces;


    for (Vector<size_t>::InnerIterator it(adj_f_e, 0); it; ++it)
    {
        if (it.value() != 0) {
            clMesh.addEdge(it.row());
        }   
    }

    auto edges = this->buildEdgeVector(clMesh);


    Vector<size_t> adj_e_v = A0.transpose() * edges;


    for (Vector<size_t>::InnerIterator it(adj_e_v, 0); it; ++it)
    {
        if (it.value() != 0) {
            clMesh.addVertex(it.row());
        }   
    }

    return clMesh; // placeholderQ
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    auto x1 = this->closure(this->star(subset));
    auto x2 = this->star(this->closure(subset));

    x1.deleteSubset(x2); 

    return x1; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    auto cl = this->closure(subset);
    return cl.equals(subset);
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    bool vertex_has_no_edge = false;

    for (auto v_id: subset.vertices) {
        auto v = mesh->vertex(v_id);
        bool edge_found = false;

        for (auto edge: v.adjacentEdges()) {
            if (subset.edges.find(edge.getIndex()) != subset.edges.end()) {
                edge_found = true;
                break;
            }
        }
        if (!edge_found) {
            vertex_has_no_edge = true;
            break;
        }
    }

    bool edge_has_no_face = false;


    for (auto e_id: subset.edges) {
        auto e = mesh->edge(e_id);
        bool face_found = false;

        for (auto f: e.adjacentFaces()) {
            if (subset.faces.find(f.getIndex()) != subset.faces.end()) {
                face_found = true;
                break;
            }
        }
        if (!face_found) {
            edge_has_no_face = true;
            break;
        }
    }

    if (subset.vertices.size() > 0) {
        if (subset.faces.size() > 0 && !edge_has_no_face && !vertex_has_no_edge) {
            return 2;
        }

        if (subset.faces.size() == 0 && !vertex_has_no_edge) {
            return 1;
        }


        if (subset.edges.size() == 0 && subset.faces.size() == 0) {
            return 0;
        }

    }
    
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    auto degree = isPureComplex(subset);
    MeshSubset boundarySet;
    MeshSubset interiourSet;

    if (degree == 2) {
        for (auto f_id: subset.faces) {
            for (auto edge: mesh->face(f_id).adjacentEdges()) {
                auto e_id = edge.getIndex();
                if (interiourSet.edges.find(e_id) != interiourSet.edges.end()) {
                    continue;
                }

                if (boundarySet.edges.find(e_id) != boundarySet.edges.end()) {
                    boundarySet.deleteEdge(e_id);
                    interiourSet.addEdge(e_id);
                    continue;
                }
                boundarySet.addEdge(e_id);
            }

        }

        for (auto e_id: boundarySet.edges) {
            auto e = mesh->edge(e_id);
            for (auto v: e.adjacentVertices()) {
                boundarySet.addVertex(v.getIndex());
            }
        }
    }

    if (degree == 1) {
        for (auto e_id: subset.edges) {
            for (auto vertex: mesh->edge(e_id).adjacentVertices()) {
                auto v_id = vertex.getIndex();
                if (interiourSet.vertices.find(v_id) != interiourSet.vertices.end()) {
                    continue;
                }

                if (boundarySet.vertices.find(v_id) != boundarySet.vertices.end()) {
                    boundarySet.deleteVertex(v_id);
                    interiourSet.addVertex(v_id);
                    continue;
                }
                boundarySet.addVertex(v_id);
            }
        }


    }    



    return boundarySet; // placeholder
}