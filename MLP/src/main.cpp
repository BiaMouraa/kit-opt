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
    double accumulated;
} Solution;

typedef struct InsertionInfo{
    int nInserted;
    int edgeRemoved;
    double cost;
}InsertionInfo;

typedef struct Subsequence{
    double T, C; //duracao e custo acumulado
    int W; //custo de atraso
    int first, last; //primeiro e ultimo nos da subsequencia
    inline static Subsequence Concatenate(Subsequence& sigma1, Subsequence& sigma2, Data& data){
        Subsequence sigma;
        double temp = data.getDistance(sigma1.last,sigma2.first);
        sigma.W = sigma1.W + sigma2.W;
        sigma.T = sigma1.T + temp + sigma2.T; 
        sigma.C = sigma1.C + sigma2.W*(sigma1.T+temp) + sigma2.C;
        sigma.first = sigma1.first;
        sigma.last = sigma2.last;
        return sigma;
    }
} Subsequence;

//n = numero de nos da instancia
//s = solucao atual
//subseq_matrix[i][j] = subseq que começa no nó i e termina no nó j de s
//define todas as possibilidades de subsequencias da solucao
void UptadeAllSubseq(Solution *s, vector<vector<Subsequence>>& subseq_matrix, Data& data, int begin = -1, int end = -1){
    int n = s->sequence.size();

    //subsequencias de um unico nó
    for(int i = 0; i < n; i++){
        int v = s->sequence[i];
        if(i == 0)
            subseq_matrix[i][i].W = 0;
        else
            subseq_matrix[i][i].W = 1;
        subseq_matrix[i][i].C = 0;
        subseq_matrix[i][i].T = 0;
        subseq_matrix[i][i].first = v;
        subseq_matrix[i][i].last = v;
    }

    if(begin != -1 && end != -1){ //atualiza apenas as subsequencias que foram alteradas
        for(int i = 0; i <= end; i++){
            for(int j = begin; j <= end; j++){
                if(i >= j)
                    j = i + 1;
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j], data);
            }
        }

        for(int i = n - 1; i >= begin; i--){    
            for(int j = end; j >= 0; j--){
                if(i <= j)
                    j = i - 1;
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j], data);
            }
        }

    }else{ //define todas as subsequencias [i][j]
        for (int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j], data);
            }
        }

        //subsequencias invertidas (necessario para o 2-opt)
        for(int i = n - 1; i >= 0; i--){    
            for(int j = i - 1; j >= 0; j--){
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j], data);
            }
        }
    }
    s->accumulated = subseq_matrix[0][n].C;
}

void showSolution(Solution *s){
    int n = s->sequence.size() - 1;
    for(int i = 0; i < n; i++){
        cout << s->sequence[i] << " -> ";
    }
    cout << s->sequence.back() << endl;
}

void calculateSolutionCost(Solution *s, Data& data){
    s->cost = 0;
    int n = s->sequence.size() - 1;
    for(int i = 0; i < n; i++)
        s->cost += data.getDistance(s->sequence[i],s->sequence[i+1]);
}

vector<InsertionInfo> calculateInsertionCost(int r, vector<int> CL, Data& data){
    vector<InsertionInfo> insertionCost (CL.size());

    //calcula a distancia do ultimo no (r) a todos os elementos de CL
    int count = 0;
    for(auto k : CL){
        insertionCost[count].cost = data.getDistance(r, k);
        insertionCost[count].nInserted = k;
        count++;
    }
    return insertionCost;
}

bool compare(InsertionInfo m, InsertionInfo n){
    return (m.cost < n.cost);
}

Solution Construction(Data& data){
    Solution s = {{1,1},0};
    int r = 1;


    //creating CL
    size_t n = data.getDimension();
    vector<int> CL;
    for(size_t i = 2; i <= n; i++) {
        CL.push_back(i);
        }

    
    // //choose 3 random from CL and increment in s
    // for(int i = 0 ; i < 3; i++){
    //     int pos = rand() % CL.size();
    //     int value = CL[pos];
    //     s.sequence.push_back(value);
    //     CL.erase(CL.begin() + pos);
    // }



    while(!CL.empty()){
        vector<InsertionInfo> insertionCost = calculateInsertionCost(r, CL, data);
        sort(insertionCost.begin(), insertionCost.end(), compare);

        double alpha = (double) rand() / (RAND_MAX + 1.0);
        alpha += 0.000000000001;
        if(alpha == 0)
            cout << alpha << endl;

        //define e insere o no na sequencia
        int selected = rand() % ((int) ceil(alpha*insertionCost.size()));
        r = insertionCost[selected].nInserted; 
        s.sequence.insert(s.sequence.begin()+s.sequence.size()-1, r);

        int clPos = 0;
        while(true){
            if(CL[clPos] == insertionCost[selected].nInserted){
                
                CL.erase(CL.begin() + clPos);
                break;
            }
            clPos++;
        }
        
    }
    calculateSolutionCost(&s, data);

    return s;
}

bool bestImprovementSwap(Solution *s, vector<vector<Subsequence>>& subseq_matrix, Data& data){
    //cout << "swap";
    int best_i, best_j;
    int n = s->sequence.size() - 1;
    int totalCost = subseq_matrix[0][n].C;
    int bestCost = totalCost;

    Subsequence sigma1, sigma2, sigma3, sigma;

    for(int i = 1; i < n; i++){
        for(int j = i + 1; j < n; j++){

            if(j == i + 1){
                sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j], data);
                sigma2 = Subsequence::Concatenate(sigma1, subseq_matrix[i][i], data);
                sigma = Subsequence::Concatenate(sigma2, subseq_matrix[j+1][n], data);
            
            }else{
                sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j], data);
                sigma2 = Subsequence::Concatenate(sigma1, subseq_matrix[i+1][j-1], data);
                sigma3 = Subsequence::Concatenate(sigma2, subseq_matrix[i][i], data);
                sigma = Subsequence::Concatenate(sigma3, subseq_matrix[j+1][n], data);
            }
            if(sigma.C < bestCost){
                bestCost = sigma.C;
                best_i = i;
                best_j = j;
            }
        }
    }


    if(bestCost < totalCost){
        swap(s->sequence[best_i], s->sequence[best_j]);
        UptadeAllSubseq(s, subseq_matrix, data, best_i, best_j);
        return true;
    }

    return false;
}

bool bestImprovementOrOpt(Solution *s, vector<vector<Subsequence>>& subseq_matrix, int type, Data& data){
    int best_i, best_j;
    int n = s->sequence.size() - 1;
    double delta_i, delta_j;
    double totalCost = subseq_matrix[0][n].C;
    double bestCost = totalCost;

    Subsequence sigma1, sigma2, sigma;

    switch(type){
        case 1: //Reinsertion
            for(int i = 1; i < n; i++){
                for(int j = i + 1; j < n; j++){
                    sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+1][j], data);
                    sigma2 = Subsequence::Concatenate(sigma1, subseq_matrix[i][i], data);
                    sigma = Subsequence::Concatenate(sigma2, subseq_matrix[j+1][n], data);

                    if(sigma.C < bestCost){
                        bestCost = sigma.C;
                        best_i = i;
                        best_j = j;
                    }
                }
            }

            if(bestCost < totalCost){
                int reinsert = s->sequence[best_i];
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j, reinsert);

                UptadeAllSubseq(s, subseq_matrix, data, best_i, best_j);
                return true;
            }
            return false;

        case 2: //Or-Opt2
            for(int i = 1; i < n - 1; i++){
                for(int j = i + 2; j < n; j++){
                    
                    sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+2][j], data);
                    sigma2 = Subsequence::Concatenate(sigma1, subseq_matrix[i][i+1], data);
                    sigma = Subsequence::Concatenate(sigma2, subseq_matrix[j+1][n], data);

                    if(sigma.C < bestCost){
                        bestCost = sigma.C;
                        best_i = i;
                        best_j = j;
                    }
                }
            }

            
            if(bestCost < totalCost){
                int reinsert1 = s->sequence[best_i];
                int reinsert2 = s->sequence[best_i + 1];

                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j - 1, reinsert2);
                s->sequence.insert(s->sequence.begin() + best_j - 1, reinsert1);

                UptadeAllSubseq(s, subseq_matrix, data, best_i, best_j);
                return true;
            }
            return false;

        case 3: //Or-Opt3
            for(int i = 1; i < n - 2; i++){
                for(int j = i + 3; j < n; j++){
                    sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+3][j], data);
                    sigma2 = Subsequence::Concatenate(sigma1, subseq_matrix[i][i+2], data);
                    sigma = Subsequence::Concatenate(sigma2, subseq_matrix[j+1][n], data);

                    if(sigma.C < bestCost){

                        bestCost = sigma.C;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            
            if(bestCost < totalCost){
                int reinsert1 = s->sequence[best_i];
                int reinsert2 = s->sequence[best_i + 1];
                int reinsert3 = s->sequence[best_i + 2];

                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.erase(s->sequence.begin() + best_i);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert3);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert2);
                s->sequence.insert(s->sequence.begin() + best_j - 2, reinsert1);

                UptadeAllSubseq(s, subseq_matrix, data, best_i, best_j);
                return true;
            }
            return false;
    }
    return false;
}


bool bestImprovement2Opt(Solution *s, vector<vector<Subsequence>>& subseq_matrix, Data& data){
    int best_i, best_j;
    int n = s->sequence.size() - 1;
    double totalCost = subseq_matrix[0][n].C;
    double bestCost = totalCost;
    

    Subsequence sigma1, sigma;

    for(int i = 1; i < n; i++){
        for(int j = i + 2; j < n; j++){

            sigma1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][i], data);
            sigma = Subsequence::Concatenate(sigma1, subseq_matrix[j+1][n], data);
            
            if(sigma.C < bestCost){
                bestCost = sigma.C;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(bestCost < totalCost){
        for(int i = best_i, j = best_j; i < j; i++, j--){
            int aux = s->sequence[i];
            s->sequence[i] = s->sequence[j];
            s->sequence[j] = aux;
        }

        UptadeAllSubseq(s, subseq_matrix, data);
        return true;
    }
    return false;
}

void LocalSearch(Solution& s, vector<vector<Subsequence>>& subseq_matrix, Data& data){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;


    while(!NL.empty()){
        int n = rand() % NL.size();
        switch (NL[n])
        {
        case 1:
           improved = bestImprovementSwap(&s, subseq_matrix, data);
            break;
        case 2:
            improved = bestImprovement2Opt(&s, subseq_matrix, data);
            break;
        case 3:
            improved = bestImprovementOrOpt(&s, subseq_matrix, 1, data);
            break;
        case 4: 
            improved = bestImprovementOrOpt(&s, subseq_matrix, 2, data);
            break;
        case 5:
            improved = bestImprovementOrOpt(&s, subseq_matrix, 3, data);
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

Solution ILS(int maxIter, int maxIterILS, Data& data){
    Solution bestOfAll;
    bestOfAll.cost = INFINITY;
    double bestCostAcum;
    double bestOfAllCost;
    size_t n = data.getDimension();
    vector<vector<Subsequence>> subseq_matrix(n + 1, vector<Subsequence>(n + 1));
 
    for(int i = 0; i < maxIter; i++){
        Solution s = Construction(data);
        UptadeAllSubseq(&s, subseq_matrix, data);
        Solution best = s;
        bestCostAcum = subseq_matrix[0][n].C;
        cout << " C ";
        
        int iterILS = 0;
        while(iterILS <= maxIterILS){
            LocalSearch(s, subseq_matrix, data);
            UptadeAllSubseq(&s, subseq_matrix, data);
            calculateSolutionCost(&s, data);

            cout << " L ";

            if(subseq_matrix[0][n].C < bestCostAcum){
                best = s;
                iterILS = 0;
                bestCostAcum = subseq_matrix[0][n].C;
            } 

            Perturbation(&best);

            cout << " P ";

            UptadeAllSubseq(&best, subseq_matrix, data);
            iterILS++;
        }
        if(bestCostAcum < bestOfAllCost){
            bestOfAll = best;
            bestOfAllCost = bestCostAcum;
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
        vector<vector<Subsequence>> subseq_matrix(n + 1, vector<Subsequence>(n + 1)); 
        auto begin = chrono::high_resolution_clock::now();

        cout << maxIter << " - " << maxIterILS << endl;
        best = ILS(maxIter, maxIterILS, data);
        UptadeAllSubseq(&best, subseq_matrix, data);
        cout << "saindo" << endl;

        auto end = chrono::high_resolution_clock::now();
        auto time = chrono::duration_cast<chrono::milliseconds>(end - begin);


        sumCost += subseq_matrix[0][n].C;
        sumTime += (time.count()/1000.0);

        cout << i << "- " << sumTime << " - " << subseq_matrix[0][n].C << endl;
        
            
    }


    sumCost /= 10.0;
    sumTime /= 10.0;

    cout << sumTime << " " << sumCost << endl;
    

    return 0;
}