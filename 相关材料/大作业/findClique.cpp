#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

using std::vector, std::ifstream, std::string, std::cout, std::endl, std::ofstream;

// read graph from file, and store it in adjacency list
vector<vector<int>> readGraphInList(string &filename)
{
    ifstream fin(filename);
    int n, e;
    fin >> n >> e;

    vector<vector<int>> graph(n);
    int u, v;
    while (fin >> u >> v)
    {
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    return graph;
}

// find a clique
vector<int> findClique(vector<vector<int>> &graph)
{
    const int n = graph.size();

    vector<int> clique;

    return clique;
}

// check if the found clique is a clique
bool isClique(const vector<vector<int>> graph, const vector<int> &clique)
{
    for (int i = 0; i < clique.size(); i++)
    {
        for (int j = i + 1; j < clique.size(); j++)
        {
            if (std::find(graph[clique[i]].begin(), graph[clique[i]].end(), clique[j]) == graph[clique[i]].end())
            {
                return false;
            }
        }
    }
    return true;
}

int main()
{
    // read the graph and store it in adjacency list
    string file = "frb100-40.txt";
    auto graph = readGraphInList(file);

    auto clique = findClique(graph);

    if (isClique(graph, clique))
    {
        cout << "Not a clique." << endl;
    }

    // output the clique
    ofstream fout("clique.txt");
    fout << clique.size() << endl;
    for (auto v : clique)
    {
        fout << v << " ";
    }
    fout.close();

    return 0;
}
