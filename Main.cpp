#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <stack> // For std::stack
#include <queue> // For std::queue
using namespace std;

const int MAX_NODES = 100;
const int INF = 1000000;

class Graph {
public:
    Graph(const char *filename);
    void printMatrix();
    void printAdjList();
    void dfs(int start);
    void bfs(int start);
    void prim(int start);
    void dijkstra(int start);
    void kruskal();

private:
    int adjMatrix[MAX_NODES][MAX_NODES];
    int adjList[MAX_NODES][MAX_NODES];
    int adjListSize[MAX_NODES];
    int numNodes;
    void resetVisited();
    bool visited[MAX_NODES];

    struct Edge {
        int src, dest, weight;
    };
    void sortEdges(Edge edges[], int edgeCount);
};

Graph::Graph(const char *filename) {
    for (int i = 0; i < MAX_NODES; ++i) {
        for (int j = 0; j < MAX_NODES; ++j) {
            adjMatrix[i][j] = 0;
            adjList[i][j] = -1;
        }
        visited[i] = false;
        adjListSize[i] = 0;
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        cout << "File could not be opened\n";
        exit(3);
    }

    int numEdges;
    fscanf(file, "%d %d", &numNodes, &numEdges);

    for (int i = 0; i < numEdges; i++) {
        int src, dest, weight;
        fscanf(file, "%d %d %d", &src, &dest, &weight);
        adjMatrix[src - 1][dest - 1] = weight;
        adjMatrix[dest - 1][src - 1] = weight;
        adjList[src - 1][adjListSize[src - 1]++] = dest - 1;
        adjList[dest - 1][adjListSize[dest - 1]++] = src - 1;
    }

    fclose(file);

    for (int i = 0; i < numNodes; i++) {
        std::sort(adjList[i], adjList[i] + adjListSize[i]);
    }
}

void Graph::resetVisited() {
    for (int i = 0; i < MAX_NODES; i++) {
        visited[i] = false;
    }
}

void Graph::printMatrix() {
    cout << "Adjacency Matrix:\n";
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void Graph::printAdjList() {
    cout << "Adjacency List:\n";
    for (int i = 0; i < numNodes; i++) {
        cout << i + 1 << ": ";
        for (int j = 0; j < adjListSize[i]; j++) {
            cout << adjList[i][j] + 1 << " ";
        }
        cout << endl;
    }
}

void Graph::dfs(int start) {
    resetVisited();
    stack<int> s;
    s.push(start - 1);

    cout << "DFS Node Order: ";
    while (!s.empty()) {
        int current = s.top();
        s.pop();

        if (!visited[current]) {
            visited[current] = true;
            cout << current + 1 << " ";

            for (int i = adjListSize[current] - 1; i >= 0; i--) {
                int neighbor = adjList[current][i];
                if (!visited[neighbor]) {
                    s.push(neighbor);
                }
            }
        }
    }
    cout << endl;
}

void Graph::bfs(int start) {
    resetVisited();
    queue<int> q;
    q.push(start - 1);

    cout << "BFS Node Order: ";
    while (!q.empty()) {
        int current = q.front();
        q.pop();

        if (!visited[current]) {
            visited[current] = true;
            cout << current + 1 << " ";

            for (int i = 0; i < adjListSize[current]; i++) {
                int neighbor = adjList[current][i];
                if (!visited[neighbor]) {
                    q.push(neighbor);
                }
            }
        }
    }
    cout << endl;
}

void Graph::dijkstra(int start) {
    int dist[MAX_NODES];
    bool finalized[MAX_NODES];

    // Initialize distances and finalized set
    for (int i = 0; i < numNodes; i++) {
        dist[i] = INF;  // Initially set all distances to infinity
        finalized[i] = false;
    }
    dist[start - 1] = 0; // Distance to the start node is 0

    cout << "Dijkstra's Shortest Paths from node " << start << ":\n";
    for (int count = 0; count < numNodes; count++) {
        // Find the unfinalized node with the smallest distance
        int u = -1;
        for (int i = 0; i < numNodes; i++) {
            if (!finalized[i] && (u == -1 || dist[i] < dist[u])) {
                u = i;
            }
        }

        finalized[u] = true; // Mark the node as finalized
        cout << "Node " << u + 1 << ": " << dist[u] << endl;

        // Update distances for the neighbors
        for (int v = 0; v < numNodes; v++) {
            if (adjMatrix[u][v] != 0 && !finalized[v] && dist[u] + adjMatrix[u][v] < dist[v]) {
                dist[v] = dist[u] + adjMatrix[u][v];
            }
                }
        }
}
void Graph::prim(int start) {
    bool inMST[MAX_NODES];  // Track nodes in the MST
    int parent[MAX_NODES];   // Parent of each node in MST
    int key[MAX_NODES];      // Key values for the nodes

    // Initialize all values
    for (int i = 0; i < numNodes; i++) {
        key[i] = INF;        // Initially set all keys to infinity
        inMST[i] = false;    // No node is in MST initially
    }
    key[start - 1] = 0;       // Set the key of the starting node to 0

    cout << "Prim's Minimum Spanning Tree starting from node " << start << ":\n";
    int totalCost = 0;

    // Run Prim's algorithm for all nodes
    for (int count = 0; count < numNodes; count++) {
        // Find the node with the smallest key value
        int u = -1;
        for (int i = 0; i < numNodes; i++) {
            if (!inMST[i] && (u == -1 || key[i] < key[u])) {
                u = i;
            }
        }

        inMST[u] = true;  // Include node u in MST
        totalCost += key[u];

        // Update key values of neighbors
        for (int v = 0; v < numNodes; v++) {
            if (adjMatrix[u][v] != 0 && !inMST[v] && adjMatrix[u][v] < key[v]) {
                key[v] = adjMatrix[u][v];  // Update key for v
                parent[v] = u;              // Set u as parent of v
            }
        }
    }

    // Print the MST edges and total cost
    for (int i = 0; i < numNodes; i++) {
        if (i != start - 1) {  // Avoid printing the starting node itself
            cout << "Edge: " << parent[i] + 1 << " - " << i + 1 << " (Cost: " << adjMatrix[i][parent[i]] << ")\n";
        }
    }
    cout << "Total MST Cost: " << totalCost << endl;
}


void Graph::kruskal() {
    struct DisjointSet {
        int parent[MAX_NODES];
        int rank[MAX_NODES];

        void makeSet() {
            for (int i = 0; i < MAX_NODES; i++) {
                parent[i] = i;
                rank[i] = 0;
            }
        }

        int find(int x) {
            if (parent[x] != x)
                parent[x] = find(parent[x]);
            return parent[x];
        }

        void unionSets(int x, int y) {
            int rootX = find(x);
            int rootY = find(y);

            if (rootX != rootY) {
                if (rank[rootX] > rank[rootY]) {
                    parent[rootY] = rootX;
                } else if (rank[rootX] < rank[rootY]) {
                    parent[rootX] = rootY;
                } else {
                    parent[rootY] = rootX;
                    rank[rootX]++;
                }
            }
        }
    };

    DisjointSet ds;
    ds.makeSet();

    Edge edges[MAX_NODES * MAX_NODES];
    int edgeCount = 0;

    // Collect all edges
    for (int i = 0; i < numNodes; i++) {
        for (int j = i + 1; j < numNodes; j++) {
            if (adjMatrix[i][j] != 0) {
                edges[edgeCount++] = {i, j, adjMatrix[i][j]};
            }
        }
    }

    // Sort edges based on weight using bubble sort (manual sort)
    for (int i = 0; i < edgeCount - 1; i++) {
        for (int j = i + 1; j < edgeCount; j++) {
            if (edges[i].weight > edges[j].weight) {
                Edge temp = edges[i];
                edges[i] = edges[j];
                edges[j] = temp;
            }
        }
    }

    cout << "Kruskal's Minimum Spanning Tree:\n";
    int mstCost = 0;
    for (int i = 0; i < edgeCount; i++) {
        int u = edges[i].src;
        int v = edges[i].dest;
        if (ds.find(u) != ds.find(v)) {
            ds.unionSets(u, v);
            mstCost += edges[i].weight;
            cout << "Edge: " << u + 1 << " - " << v + 1 << " (Cost: " << edges[i].weight << ")\n";
        }
    }
    cout << "Total MST Cost: " << mstCost << endl;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "File was not entered\n";
        return 2;
    }

    Graph g(argv[1]);

    cout << "Graph loaded successfully.\n";
    g.printMatrix();
    g.printAdjList();

    int startNode;

    cout << "Enter starting node for DFS: ";
    cin >> startNode;
    g.dfs(startNode);

    cout << "Enter starting node for BFS: ";
    cin >> startNode;
    g.bfs(startNode);

    // Dijkstra
    cout << "Enter starting node for Dijkstra: ";
    scanf("%d", &startNode);
    g.dijkstra(startNode);

        // Prim's Algorithm
        cout << "Enter starting node for Prim's algorithm: ";
        scanf("%d", &startNode);
        g.prim(startNode);

        // Kruskal
    g.kruskal();


    return 0;
}
