#include <algorithm>
// #include <bmiintrin.h>
#include <x86gprintrin.h>
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <limits.h>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>

/* Size maximums for allocation. */
#define MAXVARS 100
#define MAXDOMAIN 40
#define MAXVARNEIGHBORS 35
#define MAXVIOLATIONS 1000

/* Vectorized domain structure. */
#define WORDDOMAINMEMBERS 10 /* number of domain members to pack in to each 64-bit word (need not exactly fill it) */

#define BITSPERDOMAINMEMBER 6
#define DOMAINMEMBERMASK 0x3fLL /* should reflect BITSPERDOMAINMEMBER */
#define MAXDOMAINMEMBERVIOLATIONS 0x1fLL /* should be < DOMAINMEMBERMASK, and should be all-1's in binary */

std::vector<std::vector<int>> GRAPH;
std::vector<int> CLIQUE;
/* Static instance data. */
int vars, domainwords, domain;
unsigned long long int whichconstraint[MAXVARS][MAXDOMAIN][MAXVARS];
char varneighbors[MAXVARS][MAXVARNEIGHBORS];
char whichvarneighbor[MAXVARS][MAXVARNEIGHBORS];
/* my varneighbor index in var neighbor's neighbors */
char numvarneighbors[MAXVARS];
__m256i varneighborvector[MAXVARS][MAXDOMAIN][MAXVARNEIGHBORS];
__m256i highbits, lowbits;

/* Dynamic search state data. */
char varval[MAXVARS];
/* which domain member is this var currently assigned to */
char violations[MAXVIOLATIONS][2];
/* store the two variables in the violated constraint */
int numviolations;


/* Assumes each CSP variable is represented by a clique of sequentially numbered nodes in a DIMACS format Maximum
 * Independent Set instance. */
void loadinstance(char* filename)
{
  for (int i = 0; i < vars; i++)
  {
    for (int j = 0; j < domain; j++)
    {
      for (int k = 0; k < vars; k++)
      {
        whichconstraint[i][j][k] = 0;
      }
    }
  }

  FILE* fp = fopen(filename, "r");
  char line[1000];
  int readflag = 1;
  while (readflag)
  {
    char* f = fgets(line, 1000, fp);
    if (f == NULL)
    {
      readflag = 0;
    }
    else
    {
      if (line[0] == 'e')
      {
        int v1, v2;
        sscanf(&(line[2]), "%d %d", &v1, &v2);
        whichconstraint[(v1 - 1) / domain][(v1 - 1) % domain][(v2 - 1) / domain]
          |=
          (1ULL << ((v2 - 1) % domain));
        /* Assume node numbering starts at 1. */
        whichconstraint[(v2 - 1) / domain][(v2 - 1) % domain][(v1 - 1) / domain]
          |=
          (1ULL << ((v1 - 1) % domain));
      }
    }
  }
  fclose(fp);

  int tempflag[MAXVARS][MAXVARS];
  for (int c = 0; c < vars; c++)
  {
    numvarneighbors[c] = 0;
    for (int d = 0; d < vars; d++)
    {
      tempflag[c][d] = 0;
      for (int i = 0; i < domain; i++)
      {
        for (int j = 0; j < domain; j++)
        {
          if ((c != d) && (whichconstraint[c][i][d] & (1ULL << j)) && (
            tempflag[c][d] == 0))
          {
            varneighbors[c][numvarneighbors[c]++] = d;
            tempflag[c][d] = 1;
          }
        }
      }
    }
  }
  for (int i = 0; i < vars; i++)
  {
    for (int j = 0; j < numvarneighbors[i]; j++)
    {
      int n = varneighbors[i][j];
      for (int k = 0; k < numvarneighbors[n]; k++)
      {
        if (varneighbors[n][k] == i)
          whichvarneighbor[i][j] = k;
      }
    }
  }
  for (int c = 0; c < vars; c++)
  {
    for (int i = 0; i < domain; i++)
    {
      for (int j = 0; j < numvarneighbors[c]; j++)
      {
        for (int k = 0; k < domainwords; k++)
        {
          varneighborvector[c][i][j][k] = 0LL;
        }
        for (int k = 0; k < domain; k++)
        {
          if (whichconstraint[c][i][varneighbors[c][j]] & (1ULL << k))
            varneighborvector[c][i][j][k / WORDDOMAINMEMBERS] |=
              (1LL << ((k % WORDDOMAINMEMBERS) * (BITSPERDOMAINMEMBER)));
        }
      }
    }
  }
  for (int k = 0; k < domainwords; k++)
  {
    highbits[k] = 0;
    lowbits[k] = 0;
  }
  for (int k = 0; k < domain; k++)
  {
    highbits[k / WORDDOMAINMEMBERS] |=
    (((1ULL << ((k % WORDDOMAINMEMBERS) * (BITSPERDOMAINMEMBER)))) << (
      BITSPERDOMAINMEMBER - 1));
    lowbits[k / WORDDOMAINMEMBERS] |= (
      1ULL << ((k % WORDDOMAINMEMBERS) * (BITSPERDOMAINMEMBER)));
  }
}


void printcursol()
{
  printf("sol %d: ", numviolations);
  for (int i = 0; i < vars; i++)
  {
    CLIQUE.push_back(i * MAXDOMAIN + varval[i]);
    printf("%d ", i * MAXDOMAIN + varval[i]);
  }
  printf("\n");
}


void setvar(int u, int val)
{
  varval[u] = val;
  int n = numvarneighbors[u];
  unsigned long long int* p = whichconstraint[u][val];
  for (int j = 0; j < n; j++)
  {
    int c = varneighbors[u][j];
    int r = varval[c];
    if ((r >= 0) && (p[c] & (1ULL << r)))
    {
      violations[numviolations][0] = u;
      violations[numviolations][1] = c;
      numviolations++;
    }
  }
}

void resetvar(int u, int val, int eindex, int dr, int dr2)
{
  varval[u] = val;

  numviolations--;
  violations[eindex][0] = violations[numviolations][0];
  violations[eindex][1] = violations[numviolations][1];
  int numsat = 1;
  for (int i = 0; numsat < dr; i++)
  {
    if (((violations[i][0]) == u) || ((violations[i][1]) == u))
    {
      numviolations--;
      violations[i][0] = violations[numviolations][0];
      violations[i][1] = violations[numviolations][1];
      i--;
      numsat++;
    }
  }

  int numv = 0;
  unsigned long long int* p = whichconstraint[u][val];
  for (int j = 0; numv < dr2; j++)
  {
    int c = varneighbors[u][j];
    int r = varval[c];
    if (p[c] & (1ULL << r))
    {
      violations[numviolations][0] = u;
      violations[numviolations][1] = c;
      numviolations++;
      numv++;
    }
  }
}


/* Return best replacement score. */
inline int best_replacements(int u, int* uscore, __m256i* bitset)
{
  __m256i dp256;
  __m256i ymm0 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
  for (int j = 0; j < numvarneighbors[u]; j++)
  {
    int which = whichvarneighbor[u][j];
    int nbr = varneighbors[u][j];
    int rep = varval[nbr];
    __m256i ymm1 = _mm256_load_si256(&(varneighborvector[nbr][rep][which]));
    ymm0 = _mm256_add_epi64(ymm0, ymm1);
  }
  int rmem = varval[u];
  _mm256_store_si256(&dp256, ymm0);
  int wnum = rmem / WORDDOMAINMEMBERS;
  int wshift = (rmem % WORDDOMAINMEMBERS) * BITSPERDOMAINMEMBER;
  *uscore = ((dp256[wnum]) >> wshift) & DOMAINMEMBERMASK;
  dp256[wnum] |= (MAXDOMAINMEMBERVIOLATIONS << wshift);
  /* max out its degree so we don't re-choose it */
  ymm0 = _mm256_load_si256(&dp256);

  *bitset = _mm256_sub_epi64(highbits, ymm0);
  int thisscore = *uscore;
  while (_mm256_testz_si256(*bitset, highbits))
  {
    *bitset = _mm256_add_epi64(*bitset, lowbits);
    thisscore--;
  }
  return (thisscore);
}

/* Variant of best_replacements, for use when variables are not yet set. */
int best_nonreplacements(int u, __m256i* bitset)
{
  __m256i ymm0 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
  for (int j = 0; j < numvarneighbors[u]; j++)
  {
    int which = whichvarneighbor[u][j];
    int nbr = varneighbors[u][j];
    int rep = varval[nbr];
    if (rep >= 0)
    {
      __m256i ymm1 = _mm256_load_si256(&(varneighborvector[nbr][rep][which]));
      ymm0 = _mm256_add_epi64(ymm0, ymm1);
    }
  }

  *bitset = _mm256_sub_epi64(highbits, ymm0);
  int thisscore = 0;
  while (_mm256_testz_si256(*bitset, highbits))
  {
    *bitset = _mm256_add_epi64(*bitset, lowbits);
    thisscore--;
  }
  return (thisscore);
}


inline void add_to_candidates(__m256i* bitset, int* candidates,
                              int* numcandidates)
{
  __m256i dp256;
  *bitset = _mm256_and_si256(*bitset, highbits);
  _mm256_store_si256(&dp256, *bitset);
  int v = 0;
  for (int k = 0; k < domainwords; k++)
  {
    while (dp256[k] > 0)
    {
      long long unsigned int num = _tzcnt_u64(dp256[k]);
      candidates[(*numcandidates)++] = v + (num / BITSPERDOMAINMEMBER);
      dp256[k] = _blsr_u64(dp256[k]);
    }
    v += WORDDOMAINMEMBERS;
  }
}


void setrandvar()
{
  int allunset[MAXVARS];
  int numallunset = 0;
  for (int i = 0; i < vars; i++)
  {
    if (varval[i] < 0)
      allunset[numallunset++] = i;
  }
  int unsetvar = allunset[rand() % numallunset];

  __m256i bitset;
  int bestscore = best_nonreplacements(unsetvar, &bitset);

  int bestsval[MAXDOMAIN];
  int numbests = 0;
  add_to_candidates(&bitset, bestsval, &numbests);
  setvar(unsetvar, bestsval[rand() % numbests]);
}

void init()
{
  for (int v = 0; v < vars; v++)
  {
    varval[v] = -1;
  }
  int numsetvars = 0;
  numviolations = 0;
  while (numsetvars < vars)
  {
    setrandvar();
    numsetvars++;
  }
}


/* Supports 1<=successgap<=3.  Return 1 for success (within successgap), 0 for failure. */
int checkpartialsuccess(int successgap)
{
  int d[MAXVARS];
  for (int i = 0; i < vars; i++)
  {
    d[i] = 0;
  }
  for (int i = 0; i < numviolations; i++)
  {
    d[(violations[i][1])]++;
    d[(violations[i][0])]++;
  }

  int two[MAXVARS];
  int numtwo = 0;
  for (int i = 0; i < vars; i++)
  {
    if (numviolations - d[i] + 1 <= successgap)
      return 1;
    if (d[i] >= 2)
      two[numtwo++] = i;
  }
  if (successgap >= 2)
  {
    for (int i = 0; i < numtwo; i++)
    {
      for (int j = i + 1; j < numtwo; j++)
      {
        int dtot = d[two[i]] + d[two[j]];
        int intra = 0;
        if (whichconstraint[two[i]][varval[two[i]]][two[j]] & (
          0x1LL << (varval[two[j]])))
          intra = 1;
        if (numviolations - dtot + intra + 2 <= successgap)
          return 1;
        if (successgap >= 3)
        {
          for (int k = j + 1; k < numtwo; k++)
          {
            int dtot3 = dtot + d[two[k]];
            int intra3 = intra;
            if (whichconstraint[two[i]][varval[two[i]]][two[k]] & (
              0x1LL << (varval[two[k]])))
              intra3++;
            if (whichconstraint[two[j]][varval[two[j]]][two[k]] & (
              0x1LL << (varval[two[k]])))
              intra3++;
            if (numviolations - dtot3 + intra3 + 3 <= successgap)
              return 1;
          }
        }
      }
    }
  }
  return 0;
}


/* Return #updates for success, or -1 for failure as of maxupdates. */
/* Only check successgap when numviolations<=successgapcheckviolations. */
/* Print new bests if at least as good as printviolations. */
long long int sls(long long maxupdates, int successgap,
                  int successgapcheckviolations, int printviolations)
{
  int bestviolations = INT_MAX;
  int prevvar = -1;
  long long numupdates = 1;
  long long lastupdate[MAXVARS];
  for (int j = 0; j < vars; j++)
  {
    lastupdate[j] = 0;
  }
  while (numupdates <= maxupdates)
  {
    if (numviolations == 0)
    {
      return (numupdates);
    }
    else if ((successgap > 0) && (
      numviolations <= successgapcheckviolations))
    {
      if (checkpartialsuccess(successgap))
        return (numupdates);
    }
    int vindex = rand() % numviolations;
    int var = violations[vindex][0];
    int backup = violations[vindex][1];
    if (((lastupdate[backup] < lastupdate[var]) && (backup != prevvar)) || (
      var == prevvar))
    {
      backup = var;
      var = violations[vindex][1];
    }

    int origvar = var;
    int bestsval[MAXDOMAIN * 2];
    int firstbackup = INT_MAX; /* first index of backup in bestsval */
    int numbests = 0;
    int dr, drb;
    __m256i ymm3, ymm3b;
    int bestscore = best_replacements(var, &dr, &ymm3);
    int backupscore = bestscore - 1; /* defaults to ignoring backup */
    if ((bestscore < 0) && (backup != prevvar))
    {
      /* best replacement is not at least iso; evaluate backup */
      var = backup;
      backupscore = best_replacements(var, &drb, &ymm3b);
    }
    if (bestscore >= backupscore)
      add_to_candidates(&ymm3, bestsval, &numbests);
    if (backupscore >= bestscore)
    {
      bestscore = backupscore;
      firstbackup = numbests;
      add_to_candidates(&ymm3b, bestsval, &numbests);
    }

    int whichbest = rand() % numbests;
    if (whichbest >= firstbackup)
    {
      dr = drb;
      prevvar = var;
    }
    else
    {
      prevvar = origvar;
    }
    resetvar(prevvar, bestsval[whichbest], vindex, dr, dr - bestscore);

    if (numviolations < bestviolations)
    {
      bestviolations = numviolations;
      if (bestviolations <= printviolations)
        printcursol();
      printf("bestconflicts %d iteration %lld numconflicts %d\n",
             bestviolations, numupdates, numviolations);
    }

    lastupdate[prevvar] = numupdates;
    numupdates++;
  }
  return -1; /* no success */
}

using std::vector, std::ifstream, std::string, std::cout, std::endl,
  std::ofstream;

// read graph from file, and store it in adjacency list
vector<vector<int>> readGraphInList(string& filename)
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

  GRAPH = graph;
  return graph;
}

// Function to compute the complement of a graph
vector<vector<int>> complementGraph(const vector<vector<int>>& graph)
{
  int n = graph.size();
  vector<vector<int>> complement(n);

  for (int u = 0; u < n; ++u)
  {
    vector<bool> isNeighbor(n, false);
    for (int v : graph[u])
    {
      isNeighbor[v] = true;
    }
    for (int v = 0; v < n; ++v)
    {
      if (u != v && !isNeighbor[v])
      {
        complement[u].push_back(v);
      }
    }
  }

  return complement;
}

// Function to save a graph in DIMACS format
void saveGraphInDIMACS(const vector<vector<int>>& graph, const string& filename)
{
  int n = graph.size();
  int e = 0;

  // Count the number of edges
  for (const auto& neighbors : graph)
  {
    e += neighbors.size();
  }
  // e /= 2; // Each edge is counted twice

  ofstream fout(filename);
  fout << "p edge " << n << " " << e << "\n";

  // Write each edge
  for (int u = 0; u < n; ++u)
  {
    for (int v : graph[u])
    {
      if (u < v)
      {
        // Ensure each edge is written only once
        fout << "e " << u + 1 << " " << v + 1 << "\n";
      }
    }
  }

  fout.close();
}

// check if the found clique is a clique
bool isClique(const vector<vector<int>> graph, const vector<int>& clique)
{
  for (int i = 0; i < clique.size(); i++)
  {
    for (int j = i + 1; j < clique.size(); j++)
    {
      if (std::find(graph[clique[i]].begin(), graph[clique[i]].end(), clique[j])
        == graph[clique[i]].end())
      {
        return false;
      }
    }
  }
  return true;
}

// Function to find and output the remaining nodes if they form a complete graph
vector<int> removeUseless(const vector<vector<int>>& graph, const vector<int>& clique, int n)
{
  if (n == 2)
  {
    for (int i = 0; i < clique.size() - 1; ++i)
    {
      for (int j = i + 1; j < clique.size(); ++j)
      {
        vector<int> remaining;
        for (int k = 0; k < clique.size(); ++k)
        {
          if (k != i && k != j)
          {
            remaining.push_back(clique[k]);
          }
        }
        if (isClique(graph, remaining))
        {
          cout << "Remaining nodes form a clique." << endl;
          cout << "Remaining nodes: ";
          for (int k : remaining)
          {
            cout << k << " ";
          }
          cout << endl;
          return remaining;
        }
        // else
        // {
        //   cout << "Remaining nodes don't form a clique." << endl;
        //   cout << "Remaining nodes: ";
        //   for (int k : remaining)
        //   {
        //     cout << k << " ";
        //   }
        //   cout << endl;
        // }
      }
    }
  }
  vector<int> ret;
  return ret;
}

// find a clique
vector<int> findClique(vector<vector<int>>& graph)
{
  const int n = graph.size();

  vector<int> clique;

  unsigned int seed = 123;
  srand(seed);
  vars = 100;
  domain = 40;
  domainwords = (domain + (WORDDOMAINMEMBERS - 1)) / WORDDOMAINMEMBERS;
  /* round up */
  if (domainwords > 4)
  {
    printf("FATAL ERROR: unsupported domain size");
    fflush(stdout);
    exit(1);
  }
  long long int runlength = 300000000000;
  int numruns = 1;
  int successgap = 2, successgapcheckviolations = 8, printviolations = -1;

  string tmp_file = "frb.dimacs";
  saveGraphInDIMACS(complementGraph(graph), tmp_file);
  loadinstance(tmp_file.data());

  for (int i = 0; i < numruns; i++)
  {
    init();
    long long int successval = sls(runlength, successgap,
                                   successgapcheckviolations, printviolations);
    if (successval >= 0)
    {
      printf("success at %lld\n", successval);
      printcursol();
      clique = removeUseless(graph, CLIQUE, 2);
    }
    else
    {
      printf("unsuccessful\n");
    }
  }

  return clique;
}

int main()
{
  auto start_time = std::chrono::high_resolution_clock::now();

  // read the graph and store it in adjacency list
  string file = "frb100-40.txt";
  auto graph = readGraphInList(file);

  auto clique = findClique(graph);

  if (!isClique(graph, clique))
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
  auto end_time = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
  cout << "time: " << time << "ms" << endl;

  return 0;
}
