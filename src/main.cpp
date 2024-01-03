#include "data.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono>

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

void calculateSolutionCost(Solution *s, Data& matrizAdj){
    s->cost = 0;
    int n = s->sequence.size() - 1;
    for(int i = 0; i < n; i++)
        s->cost += matrizAdj.getDistance(s->sequence[i],s->sequence[i+1]);
}

vector<InsertionInfo> calculateInsertionCost(Solution s, vector<int> CL, Data& matrizAdj){
    vector<InsertionInfo> insertionCost ((s.sequence.size()-1)*CL.size());

    int l = 0, count = 0;
    int n = s.sequence.size() - 1;
    for(int a = 0, b = 1; count < n; a++, b++){
        int i = s.sequence[a];
        int j = s.sequence[b];
        //for(int i = 0; i < CL.size() - 1; i++){
            for(auto k : CL){
            insertionCost[l].cost = matrizAdj.getDistance(i, k) + matrizAdj.getDistance(k, j) - matrizAdj.getDistance(i,j);
            insertionCost[l].nInserted = k;
            insertionCost[l].edgeRemoved = a;
            l++;
        }
        count++;
    }
    return insertionCost;
}

bool compare(InsertionInfo m, InsertionInfo n){
    return (m.cost < n.cost);
}

Solution Construction(Data& matrizAdj){
    Solution s;


    //creating CL
    size_t n = matrizAdj.getDimension();
    vector<int> CL;
    for(size_t i = 2; i <= n; i++) {
        CL.push_back(i);
        }

    s.sequence.push_back(1);
    //choose 3 random from CL and increment in s
    for(int i = 0 ; i < 3; i++){
        int pos = rand() % CL.size();
        int value = CL[pos];
        s.sequence.push_back(value);
        CL.erase(CL.begin() + pos);
    }

    s.sequence.push_back(1);

    while(!CL.empty()){
        vector<InsertionInfo> insertionCost = calculateInsertionCost(s, CL, matrizAdj);
        sort(insertionCost.begin(), insertionCost.end(), compare);

        double alpha = (double) rand() / (RAND_MAX + 1.0);
        alpha += 0.000000000001;
        if(alpha == 0)
            cout << alpha << endl;
        int selected = rand() % ((int) ceil(alpha*insertionCost.size()));
        
        s.sequence.insert(s.sequence.begin() + insertionCost[selected].edgeRemoved + 1, insertionCost[selected].nInserted);
        
        int clPos = 0;
        while(true){
            if(CL[clPos] == insertionCost[selected].nInserted){
                
                CL.erase(CL.begin() + clPos);
                break;
            }
            clPos++;
        }
        
    }
    calculateSolutionCost(&s, matrizAdj);

    return s;
}

bool bestImprovementSwap(Solution *s, Data& matrizAdj){


    double bestDelta = 0;
    int best_i, best_j;
    int n = s->sequence.size() - 1;


    for(int i = 2; i < n; i++){

        int ni = s->sequence[i];
        int ni_prev = s->sequence[i-1];
        int ni_next = s->sequence[i+1];

        for(int j = i + 1; j < n; j++){
            int nj = s->sequence[j];
            int nj_prev = s->sequence[j-1];
            int nj_next = s->sequence[j+1];

            double delta;

            if(j == i + 1)
                delta = -matrizAdj.getDistance(ni_prev, ni) - matrizAdj.getDistance(nj, nj_next) + matrizAdj.getDistance(ni_prev, nj) + matrizAdj.getDistance(ni, nj_next);
            else
                delta = -matrizAdj.getDistance(ni_prev, ni) - matrizAdj.getDistance(ni, ni_next) - matrizAdj.getDistance(nj_prev, nj) - matrizAdj.getDistance(nj, nj_next) + matrizAdj.getDistance(ni_prev, nj) + matrizAdj.getDistance(nj, ni_next) + matrizAdj.getDistance(nj_prev, ni) + matrizAdj.getDistance(ni, nj_next);

            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }


    if(bestDelta < 0){
        swap(s->sequence[best_i], s->sequence[best_j]);
        s->cost += bestDelta;

        return true;
    }

    return false;
}

bool bestImprovementeOrOpt(Solution *s, int type, Data& matrizAdj){
    double bestDelta = 0;
    int best_i, best_j;
    int n = s->sequence.size() - 1;
    double delta_i, delta_j;

    switch(type){
        case 1: //Reinsertion
            for(int i = 1; i < n; i++){

            int ni = s->sequence[i];
            int ni_prev = s->sequence[i-1];
            int ni_next = s->sequence[i+1];

            delta_i = matrizAdj.getDistance(ni_prev, ni_next) - matrizAdj.getDistance(ni_prev, ni) - matrizAdj.getDistance(ni, ni_next);

                for(int j = i + 1; j < n; j++){
                    int nj = s->sequence[j];
                    int nj_next = s->sequence[j+1];

                        delta_j = matrizAdj.getDistance(ni, nj_next) + matrizAdj.getDistance(ni, nj) - matrizAdj.getDistance(nj, nj_next);

                    if(delta_i + delta_j < bestDelta){
                        bestDelta = delta_i + delta_j;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
    
            if(bestDelta < 0){
                int reinsert = s->sequence[best_i];
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j, reinsert);

                s->cost += bestDelta;
                return true;
            }
            return false;

        case 2: //Or-Opt2
            for(int i = 1; i < n - 1; i++){

            int ni1 = s->sequence[i];
            int ni2= s->sequence[i+1];
            int ni_prev = s->sequence[i-1];
            int ni_next = s->sequence[i+2];

            delta_i = matrizAdj.getDistance(ni_prev, ni_next) - matrizAdj.getDistance(ni_prev, ni1) - matrizAdj.getDistance(ni2, ni_next);

                for(int j = i + 2; j < n; j++){
                    int nj = s->sequence[j];
                    int nj_next = s->sequence[j+1];

                        delta_j = matrizAdj.getDistance(ni2, nj_next) + matrizAdj.getDistance(ni1, nj) - matrizAdj.getDistance(nj, nj_next);

                    if(delta_i + delta_j < bestDelta){
                        bestDelta = delta_i + delta_j;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            if(bestDelta < 0){
                int reinsert1 = s->sequence[best_i];
                int reinsert2 = s->sequence[best_i + 1];

                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j - 1, reinsert2);
                s->sequence.insert(s->sequence.begin() + best_j - 1, reinsert1);

                s->cost += bestDelta;
                return true;
            }
            return false;

        case 3: //Or-Opt3
            for(int i = 1; i < n - 2; i++){

                int ni1 = s->sequence[i];
                int ni3= s->sequence[i+2];
                int ni_prev = s->sequence[i-1];
                int ni_next = s->sequence[i+3];
                

                delta_i = matrizAdj.getDistance(ni_prev, ni_next) - matrizAdj.getDistance(ni_prev, ni1) - matrizAdj.getDistance(ni3, ni_next);

                for(int j = i + 3; j < n; j++){
                    int nj = s->sequence[j];
                    int nj_next = s->sequence[j+1];

                        delta_j = matrizAdj.getDistance(ni3, nj_next) + matrizAdj.getDistance(ni1, nj) - matrizAdj.getDistance(nj, nj_next);

                    if(delta_i + delta_j < bestDelta){

                        bestDelta = delta_i + delta_j;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            
            if(bestDelta < 0){
                int reinsert1 = s->sequence[best_i];
                int reinsert2 = s->sequence[best_i + 1];
                int reinsert3 = s->sequence[best_i + 2];

                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert3);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert2);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert1);

                s->cost += bestDelta;
                return true;
            }
            return false;
    }
    return false;
}


bool bestImprovement2Opt(Solution *s, Data& matrizAdj){
    double bestDelta = 0;
    int best_i, best_j;
    int n = s->sequence.size() - 1;

    for(int i = 1; i < n; i++){

        int ni = s->sequence[i];
        int ni_next = s->sequence[i+1];

        double delta;

        for(int j = i + 2; j < n - 1; j++){
            int nj = s->sequence[j];
            int nj_next = s->sequence[j+1];

                delta = -matrizAdj.getDistance(ni, ni_next) - matrizAdj.getDistance(nj, nj_next) + matrizAdj.getDistance(ni, nj) + matrizAdj.getDistance(ni_next, nj_next);

            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(bestDelta < 0){
        for(int i = best_i + 1, j = best_j; i < j; i++, j--){
            int aux = s->sequence[i];
            s->sequence[i] = s->sequence[j];
            s->sequence[j] = aux;
        }

        s->cost += bestDelta;
        return true;
    }
    return false;
}

void LocalSearch(Solution *s, Data& matrizAdj){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;
    int a, b, c, d, e;
    a = b = c = d = e = 0;


    while(!NL.empty()){
        int n = rand() % NL.size();
        switch (NL[n])
        {
        case 1:
            improved = bestImprovementSwap(s, matrizAdj);
            break;
        case 2:
            improved = bestImprovement2Opt(s, matrizAdj);
            break;
        case 3:
            improved = bestImprovementeOrOpt(s, 1, matrizAdj);
            break;
        case 4: 
            improved = bestImprovementeOrOpt(s, 2, matrizAdj);
            break;
        case 5:
            improved = bestImprovementeOrOpt(s, 3, matrizAdj);
            break;
        }

        if(improved)
            NL = {1, 2, 3, 4, 5};
        else
            NL.erase(NL.begin() + n);
    
    }

}

void Perturbation(Solution* s){

    s->sequence.erase(s->sequence.begin());
    s->sequence.erase(s->sequence.end() - 1);


    int n = s->sequence.size();
    int tMax = n/10;

    int seg1First, seg2First = -1, insertSeg1, insertSeg2;
    vector <int> seg1 = {}, seg2 = {};

    //random definition of the segment sizes
    int tSeg1 = rand() % tMax + 2;
    int tSeg2 = rand() % tMax + 2;


    seg1First = rand() % (n - 1 - tSeg1);

     for(int i = 0; i < tSeg1; i++){
        seg1.push_back(s->sequence[seg1First + i]); //set segment1
    }
    

    while(seg2First == -1){
        seg2First = rand() % (n - 1 - tSeg2); //random initial element
        
        //check if seg1 and seg2 are overlapping
        for(int i = 0; i < seg1.size(); i++){ 
            for(int j = 0; j < tSeg2; j++){
                if(s->sequence[seg2First + j] == seg1[i]){
                    seg2First = -1;
                    break;
                } 
            }
        }

     }

    for(int i = 0; i < tSeg2; i++){
        seg2.push_back(s->sequence[seg2First + i]); //set segment2  
    }


    //switch seg2 elements
    for(int i = 0; i < seg2.size(); i++){
        if(seg2First < seg1First){ //seg2 before seg1
            vector<int>::iterator index = find(s->sequence.begin(), s->sequence.end(), seg2[i]);
            s->sequence.erase(index);
            s->sequence.insert(s->sequence.begin()+seg1First-1, seg2[i]);
        }else{ //seg2 after seg1
            vector<int>::iterator index = find(s->sequence.begin(), s->sequence.end(), seg2[tSeg2 - 1 - i]);
            s->sequence.erase(index);
            s->sequence.insert(s->sequence.begin()+seg1First, seg2[tSeg2 - 1 - i]);
        }
    }
    

    //switch seg1 elements
    for(int i = 0; i < seg1.size(); i++){
        if(seg2First < seg1First){ //seg2 before seg1
            vector<int>::iterator index = find(s->sequence.begin(), s->sequence.end(), seg1[tSeg1 - 1 - i]);
            s->sequence.erase(index);
            s->sequence.insert(s->sequence.begin()+seg2First, seg1[tSeg1 - 1 - i]);
        } else{
            vector<int>::iterator index = find(s->sequence.begin(), s->sequence.end(), seg1[i]);
            s->sequence.erase(index);
            s->sequence.insert(s->sequence.begin()+seg2First+1, seg1[i]);
        }
    }

    s->sequence.insert(s->sequence.begin(), 1);
    s->sequence.push_back(1);
    
}

Solution ILS(int maxIter, int maxIterILS, Data& matrizAdj){
    Solution bestOfAll;
    bestOfAll.cost = INFINITY;


    for(int i = 0; i < maxIter; i++){
        Solution s = Construction(matrizAdj);
        Solution best = s;

        int iterILS = 0;
        while(iterILS <= maxIterILS){
            LocalSearch(&s, matrizAdj);
            calculateSolutionCost(&s, matrizAdj);
            if(s.cost < best.cost){
                best = s;
                iterILS = 0;
            }
            Perturbation(&best);
            iterILS++;
        }
        if(best.cost < bestOfAll.cost){
            bestOfAll = best;
        }
    }
    return bestOfAll;

}

int main(int argc, char** argv) {

    srand((unsigned)time(NULL));

    auto data = Data(argc, argv[1]);
    data.readData();
    size_t n = data.getDimension();

    int maxIter, maxIterILS;
    Solution best;
    maxIter = 50;
    double sumCost = 0, sumTime = 0;

    

    if(n >= 150){
            maxIterILS = n / 2;
        } else{
            maxIterILS = n;
        }
 

    for(int i = 0; i < 10; i++){
        auto begin = chrono::high_resolution_clock::now();

        best = ILS(maxIter, maxIterILS, data);

        auto end = chrono::high_resolution_clock::now();
        auto time = chrono::duration_cast<chrono::milliseconds>(end - begin);


        sumCost += best.cost;
        sumTime += (time.count()/1000.0);

        cout << i << "- " << sumTime << " - " << best.cost << endl;
        
            
    }


    sumCost /= 10.0;
    sumTime /= 10.0;

    cout << sumTime << " " << sumCost << endl;
    

    return 0;
}