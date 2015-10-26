/* Dstar.cpp
 * Author: Stephen Smith (smsmith@cmu.edu) - 10-20-15
 */

#include "Dstar.h"

/* Dstar::Dstar()
 * ----------------------------------------------------------------------------
 * Constructor sets algorithm parameters
 */
Dstar::Dstar() { 
    maxSteps = 80000;  // node expansions before we give up
    C1       = 1;      // cost of an unseen cell
}

/* keyHashCode(state u) 
 * ---------------------------------------------------------------------------- 
 * Returns the key hash code for the state u, this is used to compare
 * a state that have been updated
 */
float Dstar::keyHashCode(state u) {
    return (float)(u.k(0) + 1193*u.k(1));
}

/* isValid(state u) 
 * ----------------------------------------------------------------------------
 * Returns true if state u is on the open list or not by checking if
 * it is in the hash table.
 */
bool Dstar::isValid(state u) {
    ds_oh::iterator cur = openHash.find(u);
    if (cur == openHash.end()) return false;
    if (!close(keyHashCode(u), cur->second)) return false;
    return true;
}

/* getPath() 
 * ----------------------------------------------------------------------------
 * Returns the path created by replan() in an Nx3 matrix where N is the number
 * of points.  All points lie in a single layer (int layer).
 */
MatrixXd Dstar::getPath(int layer) {
    MatrixXd r(path.size(),3);
    int i=0;
    for(std::list<state>::const_iterator it = path.begin(), end = path.end(); it != end; ++it) {
        Vector3d v;
        v << it->pos, layer;
        r.block(i,0,1,3) = v.transpose();
        i++;
    }
    return r;
}

/* occupied(state u)
 * ----------------------------------------------------------------------------
 * Returns true if the cell is occupied (non-traversable), false
 * otherwise. Obstacles are marked with a cost < 0. 
 */
bool Dstar::occupied(state u) {
    ds_ch::iterator cur = cellHash.find(u);

    if (cur == cellHash.end()) return false;
    return (cur->second.cost < 0);
}

/* init(Vector2d s, Vector2d g)
 * ----------------------------------------------------------------------------
 * Initialize the algorithm with start and goal points. Set rhs and g to INF
 */
bool Dstar::init(Vector2d s, Vector2d g) {
    cellHash.clear();
    path.clear();
    openHash.clear();
    while(!openList.empty()) openList.pop();

    k_m = 0;

    s_start.pos = s;
    s_goal.pos  = g;

    cellInfo tmp;
    tmp.g = tmp.rhs =  0;
    tmp.cost = C1;

    cellHash[s_goal] = tmp;

    tmp.g = tmp.rhs = heuristic(s_start,s_goal);
    tmp.cost = C1;
    cellHash[s_start] = tmp;
    s_start = calculateKey(s_start);

    s_last = s_start;

    return true;
}

/* makeNewCell(state u)
 * ----------------------------------------------------------------------------
 * Checks if a cell is in the hash table, if not it adds it in.
 */
bool Dstar::makeNewCell(state u) {
    if (cellHash.find(u) != cellHash.end()) return;

    cellInfo tmp;
    tmp.g       = tmp.rhs = heuristic(u,s_goal);
    tmp.cost    = C1;
    cellHash[u] = tmp;

    return true;
}

/* getG(state u)
 * ----------------------------------------------------------------------------
 * Returns the G value for state u.
 */
double Dstar::getG(state u) {
    if (cellHash.find(u) == cellHash.end()) 
    return heuristic(u,s_goal);
    return cellHash[u].g;
}

/* getRHS(state u)
 * ----------------------------------------------------------------------------
 * Returns the rhs value for state u.
 */
double Dstar::getRHS(state u) {
    if (u == s_goal) return 0;  

    if (cellHash.find(u) == cellHash.end()) 
    return heuristic(u,s_goal);
    return cellHash[u].rhs;
}

/* setG(state u, double g)
 * ----------------------------------------------------------------------------
 * Sets the G value for state u
 */
void Dstar::setG(state u, double g) {
    makeNewCell(u);  
    cellHash[u].g = g; 
}

/* setRHS(state u, double rhs)
 * ----------------------------------------------------------------------------
 * Sets the rhs value for state u
 */
double Dstar::setRHS(state u, double rhs) {
    makeNewCell(u);
    cellHash[u].rhs = rhs;
}

/* eightConnectedDist(state a, state b) 
 * ----------------------------------------------------------------------------
 * Returns the 8-connected distance between state a and state b. (sqrt(2) on 
 * the diagonals and 1 on the cardinals)
 */
double Dstar::eightConnectedDist(state a, state b) {
    double temp;
    double min = abs(a.pos(0) - b.pos(0));
    double max = abs(a.pos(1) - b.pos(1));
    if (min > max) {
        double temp = min;
        min = max;
        max = temp;
    }
    return ((M_SQRT2-1.0)*min + max);
}

/* computeShortestPath()
 * ----------------------------------------------------------------------------
 * Stops planning after maxSteps, usually caused by obstacle containment.  
 * Returns a status integer based on the planning outcome:
 * -1 for took too long
 * 1 for success
 * 0 for
*/
int Dstar::computeShortestPath() {
    list<state> s;
    list<state>::iterator i;

    if (openList.empty()) return 1;

    int k=0;
    while ((!openList.empty()) && 
           (openList.top() < (s_start = calculateKey(s_start))) || 
           (getRHS(s_start) != getG(s_start))) {

        if (k++ > maxSteps) {
            fprintf(stderr, "At maxsteps\n");
            return -1;
        }

        state u;

        bool test = (getRHS(s_start) != getG(s_start));

        // lazy remove
        while(1) { 
            if (openList.empty()) return 1;
            u = openList.top();
            openList.pop();

            if (!isValid(u)) continue;
            if (!(u < s_start) && (!test)) return 2;
            break;
        }

        ds_oh::iterator cur = openHash.find(u);
        openHash.erase(cur);

        state k_old = u;

        if (k_old < calculateKey(u)) { // u is out of date
            insert(u);
        } else if (getG(u) > getRHS(u)) { // needs update (got better)
            setG(u,getRHS(u));
            getPred(u,s);
            for (i=s.begin();i != s.end(); i++) {
                updateVertex(*i);
            }
        } else {   // g <= rhs, state has got worse
            setG(u,INFINITY);
            getPred(u,s);
            for (i=s.begin();i != s.end(); i++) {
                updateVertex(*i);
            }
            updateVertex(u);
        }
    }
    return 0;
}

/* close(double x, double y) 
 * ----------------------------------------------------------------------------
 * Returns true if x and y are within floating point error, false otherwise.
 * Serves as equality for floating point numbers.
 */
void Dstar::close(double x, double y) {
    if (std::isinf(x) && std::isinf(y)) return true;
    return (fabs(x-y) < 0.00001);
}

/* updateVertex(state u)
 * ----------------------------------------------------------------------------
 * Handles RHS values that are not equal to g (meaning the map has been 
 * updated).
 */
void Dstar::updateVertex(state u) {
    list<state> s;
    list<state>::iterator i;

    if (u != s_goal) {
        getSucc(u,s);
        double tmp = INFINITY;
        double tmp2;

        for (i=s.begin();i != s.end(); i++) {
            tmp2 = getG(*i) + cost(u,*i);
            if (tmp2 < tmp) tmp = tmp2;
        }

        if (!close(getRHS(u),tmp)) setRHS(u,tmp);
    }

    if (!close(getG(u),getRHS(u))) insert(u);
}

/* insert(state u) 
 * ----------------------------------------------------------------------------
 * Inserts the state u into the openList and the openHash.
 */
void Dstar::insert(state u) {
    ds_oh::iterator cur;
    float csum;

    u = calculateKey(u);
    cur = openHash.find(u);
    csum = keyHashCode(u);

    openHash[u] = csum;
    openList.push(u);
} 

/* remove(state u)
 * ----------------------------------------------------------------------------
 * Removes state u from openHash. The state is removed from the
 * openList in replan to save computation.
 */
void Dstar::remove(state u) {
    ds_oh::iterator cur = openHash.find(u);
    if (cur == openHash.end()) return;
    openHash.erase(cur);
}


/* trueDist(state a, state b) 
 * ----------------------------------------------------------------------------
 * Euclidean cost between state a and state b (thanks Eigen)
 */
double Dstar::trueDist(state a, state b) {
return (a.pos-b.pos).norm();
}

/* heuristic(state a, state b)
 * ----------------------------------------------------------------------------
 * The heristic we use is the 8-way distance
 * scaled by a constant C1 (should be set to <= min cost).
 */
double Dstar::heuristic(state a, state b) {
return eightConnectedDist(a,b)*C1;
}

/* calculateKey(state u)
* ----------------------------------------------------------------------------
* As per [S. Koenig, 2002]
*/
state Dstar::calculateKey(state u) {

double val = fmin(getRHS(u),getG(u));

u.k(0)  = val + heuristic(u,s_start) + k_m;
u.k(1) = val;

return u;

}

/* double Dstar::cost(state a, state b)
* ----------------------------------------------------------------------------
* Returns the cost of moving from state a to state b. This could be
* either the cost of moving off state a or onto state b, we went with
* the former. This is also the 8-way cost.
*/
double Dstar::cost(state a, state b) {

Vector2d magDif = (a.pos-b.pos).cwiseAbs();
double scale = 1;

if (magDif.sum()>1) scale = M_SQRT2;

if (cellHash.count(a) == 0) return scale*C1;
return scale*cellHash[a].cost;

}
/* void Dstar::updateCell(int x, int y, double val)
* ----------------------------------------------------------------------------
* As per [S. Koenig, 2002]
*/
void Dstar::updateCell(Vector2d c, double val) {

state u;

u.pos = c;
if ((u == s_start) || (u == s_goal)) return;

makeNewCell(u); 
cellHash[u].cost = val;

updateVertex(u);
}

/* void Dstar::getSucc(state u,list<state> &s)
* ----------------------------------------------------------------------------
* Returns a list of successor states for state u, since this is an
* 8-way graph this list contains all of a cells neighbours. Unless
* the cell is occupied in which case it has no successors. 
*/
void Dstar::getSucc(state u,list<state> &s) {

s.clear();
u.k << -1, -1;

if (occupied(u)) return;

u.pos(0) += 1;
s.push_front(u);
u.pos(1) += 1;
s.push_front(u);
u.pos(0) -= 1;
s.push_front(u);
u.pos(0) -= 1;
s.push_front(u);
u.pos(1) -= 1;
s.push_front(u);
u.pos(1) -= 1;
s.push_front(u);
u.pos(0) += 1;
s.push_front(u);
u.pos(0) += 1;
s.push_front(u);

}

/* void Dstar::getPred(state u,list<state> &s)
* ----------------------------------------------------------------------------
* Returns a list of all the predecessor states for state u. Since
* this is for an 8-way connected graph the list contails all the
* neighbours for state u. Occupied neighbours are not added to the
* list.
*/
void Dstar::getPred(state u,list<state> &s) {

s.clear();
u.k << -1, -1;

u.pos(0) += 1;
if (!occupied(u)) s.push_front(u);
u.pos(1) += 1;
if (!occupied(u)) s.push_front(u);
u.pos(0) -= 1;
if (!occupied(u)) s.push_front(u);
u.pos(0) -= 1;
if (!occupied(u)) s.push_front(u);
u.pos(1) -= 1;
if (!occupied(u)) s.push_front(u);
u.pos(1) -= 1;
if (!occupied(u)) s.push_front(u);
u.pos(0) += 1;
if (!occupied(u)) s.push_front(u);
u.pos(0) += 1;
if (!occupied(u)) s.push_front(u);

}

/* void Dstar::updateStart(int x, int y)
* ----------------------------------------------------------------------------
* Update the position of the robot, this does not force a replan.
*/
void Dstar::updateStart(Vector2d s) {

s_start.pos = s;

k_m += heuristic(s_last,s_start);

s_start = calculateKey(s_start);
s_last  = s_start;

}

/* void Dstar::updateGoal(int x, int y)
* ----------------------------------------------------------------------------
* This is somewhat of a hack, to change the position of the goal we
* first save all of the non-empty on the map, clear the map, move the
* goal, and re-add all of non-empty cells. Since most of these cells
* are not between the start and goal this does not seem to hurt
* performance too much. Also it free's up a good deal of memory we
* likely no longer use.
*/
void Dstar::updateGoal(Vector2d g) {

list< pair<Vector2d, double> > toAdd;
pair<Vector2d, double> tp;

ds_ch::iterator i;
list< pair<Vector2d, double> >::iterator kk;

for(i=cellHash.begin(); i!=cellHash.end(); i++) {
if (!close(i->second.cost, C1)) {
tp.first = i->first.pos;
tp.second = i->second.cost;
toAdd.push_back(tp);
}
}

cellHash.clear();
openHash.clear();

while(!openList.empty())
openList.pop();

k_m = 0;

s_goal.pos = g;

cellInfo tmp;
tmp.g = tmp.rhs =  0;
tmp.cost = C1;

cellHash[s_goal] = tmp;

tmp.g = tmp.rhs = heuristic(s_start,s_goal);
tmp.cost = C1;
cellHash[s_start] = tmp;
s_start = calculateKey(s_start);

s_last = s_start;    

for (kk=toAdd.begin(); kk != toAdd.end(); kk++) {
Vector2d c = kk->first;
updateCell(c, kk->second);
}

}

/* bool Dstar::replan()
* ----------------------------------------------------------------------------
* Updates the costs for all cells and computes the shortest path to
* goal. Returns true if a path is found, false otherwise. The path is
* computed by doing a greedy search over the cost+g values in each
* cells. In order to get around the problem of the robot taking a
* path that is near a 45 degree angle to goal we break ties based on
*  the metric euclidean(state, goal) + euclidean(state,start). 
*/
bool Dstar::replan() {

path.clear();

int res = computeShortestPath();
//printf("res: %d ols: %d ohs: %d tk: [%f %f] sk: [%f %f] sgr: (%f,%f)\n",res,openList.size(),openHash.size(),openList.top().k(0),openList.top().k(1), s_start.k(0), s_start.k(1),getRHS(s_start),getG(s_start));
if (res < 0) {
fprintf(stderr, "NO PATH TO GOAL\n");
return false;
}
list<state> n;
list<state>::iterator i;

state cur = s_start; 

if (std::isinf(getG(s_start))) {
// fprintf(stderr, "NO PATH TO GOAL\n");
return false;
}

while(cur != s_goal) {

path.push_back(cur);
getSucc(cur, n);

if (n.empty()) {
// fprintf(stderr, "NO PATH TO GOAL\n");
return false;
}

double cmin = INFINITY;
double tmin;
state smin;

for (i=n.begin(); i!=n.end(); i++) {

//if (occupied(*i)) continue;
double val  = cost(cur,*i);
double val2 = trueDist(*i,s_goal) + trueDist(s_start,*i); // (Euclidean) cost to goal + cost to pred
val += getG(*i);

if (close(val,cmin)) {
if (tmin > val2) {
tmin = val2;
cmin = val;
smin = *i;
}
} else if (val < cmin) {
tmin = val2;
cmin = val;
smin = *i;
}
}
n.clear();
cur = smin;
}
path.push_back(s_goal);
return true;
}

/* bool Dstar::setMap(MatrixXd m)
* --------------------------------------------------------------------------------------------
* Sets the map for the planning algorithm.
* Obstacles are 1, 0's are open
*/
bool Dstar::setMap(MatrixXd m) {
matrixMap.resize(m.rows(),m.cols());
matrixMap = m;
for(int i=0; i<m.rows(); i++) {
for(int j=0; j<m.cols(); j++) {
if(m(i,j)==OBSTACLE || m(i,j)==OCCLUSION || m(i,j)==BUFFER_OBS) {
Vector2d v(j,i);
updateCell(v,-1);
}
}
}
return true;
}

/* bool Dstar::printPath()
* --------------------------------------------------------------------------------------------
* Prints every cell in the optimal path and
* the final map with 7's as the path, 1's 
* as the obstacles and 0's as not traveresed
*/
bool Dstar::printPath() {
list<state>::iterator iter;

for(iter=path.begin(); iter != path.end(); iter++) {
cout << *iter;
if(matrixMap.size()>0) {
matrixMap(iter->pos(1),iter->pos(0)) = 7;
}
}
cout << endl;
cout << matrixMap << endl;

return true;
}
