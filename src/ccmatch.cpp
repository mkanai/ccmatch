#include <limits>
#include <vector>
#include <queue>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace std;
using namespace Rcpp;

#define INF numeric_limits<double>::infinity()
#define RET_NO_SOLUTION -1
#define RET_ABORTED -2

typedef pair<double, int> Pair;
typedef struct
{
    int to;
    int capacity;
    double cost;
    int rev;
} edge_t;
typedef vector<vector<edge_t> > graph_t;

// add edge to graph
void add_edge(graph_t &graph, int from, int to, int capacity, double cost)
{
    graph[from].push_back((edge_t){to, capacity, cost, graph[to].size()});
    graph[to].push_back((edge_t){from, 0, -cost, graph[from].size() - 1});
}


// solve minimun cost flow problem from s to t to flow f
double min_cost_flow(graph_t &graph, int s, int t, int f, bool display_progress)
{
    int V = graph.size(); // number of vertices
    vector<double> h(V, 0); // potential
    vector<double> dist(V, INF); // minimum distance from s
    vector<int> prevv(V, 0); // previous vertex
    vector<int> preve(V, 0); // previous edge

    double ret = 0;

    Progress progress(f, display_progress);
    while (f > 0)
    {
        if (Progress::check_abort())
            return RET_ABORTED;
        
        // Dijkstra algorithm
        priority_queue<Pair, vector<Pair>, greater<Pair> > que;
        fill(dist.begin(), dist.end(), INF);
        dist[s] = 0;
        que.push(Pair(0, s));      
        while (!que.empty())
        {
            Pair p = que.top();
            int v = p.second;
            double d = p.first;
            que.pop();
            if (dist[v] < d)
                continue;
            for (int i = 0; i < (int)graph[v].size(); i++)
            {
                edge_t &e = graph[v][i];
                d = dist[v] + e.cost + h[v] - h[e.to];
                if (dist[v] > d) d = dist[v]; // for rounding error
                if (e.capacity > 0 && dist[e.to] > d)
                {
                    dist[e.to] = d;
                    prevv[e.to] = v;
                    preve[e.to] = i;
                    que.push(Pair(dist[e.to], e.to));
                }
            }
        }

        // no solution
        if (dist[t] == INF)
            return RET_NO_SOLUTION;

        // update potential
        for (int v = 0; v < V; v++)
            h[v] += dist[v];

        int d = f;
        for (int v = t; v != s; v = prevv[v])
            d = min(d, graph[prevv[v]][preve[v]].capacity);
        f -= d;
        progress.increment(d);

        ret += d * h[t];
        for (int v = t; v != s; v = prevv[v])
        {
            edge_t &e = graph[prevv[v]][preve[v]];
            e.capacity -= d;
            graph[v][e.rev].capacity += d;
        }
    }
    return ret;
}

Function asDataFrame("as.data.frame");

// [[Rcpp::export(".ccmatch")]]
DataFrame ccmatch(NumericMatrix x, int N, bool display_progress) {
    int ncase = x.nrow(), ncontrol = x.ncol();

    int s = ncase + ncontrol, t = s + 1;
    int V = t + 1;
    graph_t graph(V, vector<edge_t>());

    for (int i = 0; i < ncase; i++)
    {
        add_edge(graph, s, i, N, 0);
        for (int j = 0; j < ncontrol; j++)
            add_edge(graph, i, ncase + j, 1, x(i, j));
    }
    for (int j = 0; j < ncontrol; j++)
        add_edge(graph, ncase + j, t, 1, 0);

    double min_cost = min_cost_flow(graph, s, t, ncase * N, display_progress);
    if (min_cost == RET_NO_SOLUTION)
    {
        stop("No solution found.");
    } else if (min_cost == RET_ABORTED)
    {
        stop("Aborted.");
    }

    NumericMatrix mat(ncase, N + 2);
    for (int i = 0; i < ncase; i++)
    {
        mat(i, 0) = i + 1;
        double sum = 0;
        int count = 0;
        for (int j = 0; j < ncontrol; j++)
        {
            if (graph[ncase + j][i].capacity  == 1)
            {
                mat(i, ++count) = j + 1;
                sum += x(i, j);
                if (count == N) break;
            }
        }
        mat(i, N + 1) = sum;
    }

    CharacterVector colnames = CharacterVector::create("Case");
    for (char i = '1'; i - '1' < N; i++)
    {
        string str = "Control"; str += i;
        colnames.push_back(str);
    }
    colnames.push_back("Distance");
    NumericVector rownames(ncase);
    for (int i = 0; i < ncase; i++)
        rownames[i] = i + 1;
    List dimnames = List::create(rownames, colnames);
    mat.attr("dimnames") = dimnames;

    DataFrame ret = asDataFrame(mat);
    return ret;
}