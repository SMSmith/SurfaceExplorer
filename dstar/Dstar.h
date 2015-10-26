/* Dstar.h
* Author: James Neufeld (neufeld@cs.ualberta.ca)
* Modified: Stephen Smith (smsmith@cmu.edu) - 10-20-15
*/

#ifndef DSTAR_H
#define DSTAR_H

#include <math.h>
#include <stdio.h>
#include <stack>
#include <queue>
#include <list>
#include <tr1/unordered_map>
#include <Eigen>
#include <iostream>

using namespace std;
using namespace __gnu_cxx;
using namespace Eigen;

#define OCCLUSION 3
#define OBSTACLE 1
#define BUFFER_OBS 2
#define BUFFER_SAFE 4
#define END_POINT 5
#define PATH_POINT 7

/* class state
 * ----------------------------------------------------------------------------
 * This class is the planning state, with RHS, G, and position information
 */
class state {
    public:
        Vector2d pos;
        Vector2d k;

    bool operator == (const state &s2) const {
        return ((pos(0) == s2.pos(0)) && (pos(1) == s2.pos(1)));
    }

    bool operator != (const state &s2) const {
        return ((pos(0) != s2.pos(0)) || (pos(1) != s2.pos(1)));
    }

    bool operator > (const state &s2) const {
        if (k(0)-0.00001 > s2.k(0)) return true;
        else if (k(0) < s2.k(0)-0.00001) return false;
        return k(1) > s2.k(1);
    }

    bool operator <= (const state &s2) const {
        if (k(0) < s2.k(0)) return true;
        else if (k(0) > s2.k(0)) return false;
        return k(1) < s2.k(1) + 0.00001;
    }

    bool operator < (const state &s2) const {
        if (k(0) + 0.000001 < s2.k(0)) return true;
        else if (k(0) - 0.000001 > s2.k(0)) return false;
        return k(1) < s2.k(1);
    }
};

// State variables can be printed
inline ostream &operator<<(ostream &os, state &s) {
    return os << "(" << s.pos(0) << "," << s.pos(1) << ") ";
}

/* cellInfo
* ----------------------------------------------------------------------------
* Standard values for D* planning:
* g is the total cost so far
* rhs is the cost to the parent node + the last leg cost to the current node
* cost is just the independent cost of that cell
*/
struct cellInfo {
    double g;
    double rhs;
    double cost;
};

/* class state_hash
* ----------------------------------------------------------------------------
* Used for looking up hash values (fast lookup for states)
*/
class state_hash {
    public:
        size_t operator()(const state &s) const {
            return s.pos(0) + 34245*s.pos(1);
        }
};

typedef priority_queue<state, vector<state>, greater<state> > ds_pq;
typedef tr1::unordered_map<state,cellInfo, state_hash, equal_to<state> > ds_ch;
typedef tr1::unordered_map<state, float, state_hash, equal_to<state> > ds_oh;

/* class Dstar
* ----------------------------------------------------------------------------
* The main class with public methods for updating, and replanning on the map
* and private methods for computing the plan
*/
class Dstar {

public:
    Dstar();
        bool        init(Vector2d s, Vector2d g);
        bool        updateCell(Vector2d c, double val);
        bool        updateStart(Vector2d s);
        bool        updateGoal(Vector2d g);
        bool        replan();
        bool        printPath();
        bool        setMap(MatrixXd m);
        MatrixXd    getPath(int layer);

    private:

        list<state> path;

        double C1;
        double k_m;
        state s_start, s_goal, s_last;
        int maxSteps;
        MatrixXd matrixMap;

        ds_pq openList;
        ds_ch cellHash;
        ds_oh openHash;

        bool        close(double x, double y);
        bool        makeNewCell(state u);
        double      getG(state u);
        double      getRHS(state u);
        bool        setG(state u, double g);
        double      setRHS(state u, double rhs);
        double      eightConnectedDist(state a, state b);
        int         computeShortestPath();
        bool        updateVertex(state u);
        bool        insert(state u);
        bool        remove(state u);
        double      trueDist(state a, state b);
        double      heuristic(state a, state b);
        state       calculateKey(state u);
        bool        getSucc(state u, list<state> &s);
        bool        getPred(state u, list<state> &s);
        double      cost(state a, state b); 
        bool        occupied(state u);
        bool        isValid(state u);
        float       keyHashCode(state u);
};

#endif
