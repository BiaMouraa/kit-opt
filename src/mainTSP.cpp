#include "data.cpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

typedef struct Solution {
    vector<int> sequence;
    double cost;
} Solution;

typedef struct InsertionInfo{
    int nInserted;
    int edgeRemoved;
    double cost;
}InsertionInfo;

void showSolution(Solution *s){
    int n = s->sequence.size() - 1;
    for(int i = 0; i < n; i++){
        cout << s->sequence[i] << " -> ";
    }
    cout << s->sequence.back() << endl;
}

void calculateSolutionCost(Solution *s, Data& matrixAdj){
    s->cost = 0;
    int n = s->sequence.size() - 1;
    for(int i = 0; i < n; i++)
        s->cost += matrixAdj.getDistance(s->sequence[i],s->sequence[i+1]);
}

vector<InsertionInfo> calculateInsertionCost(Solution& s, vector<int> CL, Data& matrizAdj){
    vector<InsertionInfo> insertionCost ((s.sequence.size()-1)*CL.size());
    int l = 0, count = 0;
    int n = s.sequence.size() - 1;
    for(int a = 0, b = 0; count < n; a++, b++, count++){
        int i = s.sequence[a];
        int j = s.sequence[b];
        for(auto k : CL){
            insertionCost[l].cost = matrizAdj.getDistance(i, k) + matrizAdj.getDistance(k, j) - matrizAdj.getDistance(i,j);
            insertionCost[l].nInserted = k;
            insertionCost[l].edgeRemoved = a;
            l++;
        }
    }
    return insertionCost;
}

bool compare(InsertionInfo m, InsertionInfo n){
    return (m.cost < n.cost);
}

Solution Contrucao(Data& matrizAdj){
    Solution s;

    //creating CL
    size_t n = matrizAdj.getDimension();
    vector<int> CL;
    for(size_t i = 2; i <= n; i++) 
        CL.push_back(i);

    s.sequence.push_back(1);

    //choose 3 random from CL
    for(int i =0 ; i < 3; i++){
        int pos = rand() + CL.size();
        int value = CL[pos];
        s.sequence.push_back(value);
        CL.erase(CL.begin() + pos);
    }

    
    while(!CL.empty()){
        vector<InsertionInfo> insertionCost = calculateInsertionCost(s, CL, matrizAdj);
        sort(insertionCost.begin(), insertionCost.end(), compare);
        double alpha = (double) rand() / (RAND_MAX + 1.0);
        alpha += 0.000000000001;
        int selected = rand() % ((int) ceil(alpha*insertionCost.size()));

        s.sequence.insert(s.sequence.begin() + insertionCost[selected].edgeRemoved, insertionCost[selected].nInserted);

        vector<int>::iterator clPos = find(CL.begin(), CL.end(), insertionCost[selected].nInserted);
        CL.erase(clPos);
    }

    s.sequence.push_back(1);
    s.cost = 0;

    return s;
}

int main(int argc, char** argv) {

    auto data = Data(argc, argv[1]);
    data.readData();
    size_t n = data.getDimension();

    cout << "Dimension: " << n << endl;
    cout << "DistanceMatrix: " << endl;
    data.printMatrixDist();


    cout << "Exemplo de Solucao s = ";
    double cost = 0.0;
    for (size_t i = 1; i < n; i++) {
        cout << i << " -> ";
        cost += data.getDistance(i, i+1);
    }
    cost += data.getDistance(n, 1);
    cout << n << " -> " << 1 << endl;
    cout << "Custo de S: " << cost << endl;

    return 0;
}